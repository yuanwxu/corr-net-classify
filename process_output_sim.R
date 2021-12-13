# Process results of simulated data from BAnOCC and DGCNN

library(tidyverse)
library(banocc)
library(reticulate)
library(tidytext) # for order factor within facets in ggplot
library(patchwork)
theme_set(theme_minimal())
options(bitmapType='cairo') # fix "unnable to open connection to X11 display"

# Read banocc fits for all simulated control and cases
read_bfits <- function(sample_sizes, case_number){
    out <- vector("list", length(sample_sizes))
    names(out) <- as.character(sample_sizes)
    
    for (ss in as.character(sample_sizes)) {
        out[[ss]] <- vector("list", length(case_number))
        names(out[[ss]]) <- as.character(case_number)
        for (cn in as.character(case_number)) {
            ctrl_name <- paste0(paste("fit", "ctrl", ss, cn, sep = "_"), ".rds")
            case_name <- paste0(paste("fit", "case", ss, cn, sep = "_"), ".rds")
            ctrl_bfit <- readRDS(ctrl_name)
            case_bfit <- readRDS(case_name)
            out[[ss]][[cn]] <- vector("list")
            out[[ss]][[cn]][['Ctrl']] <- ctrl_bfit
            out[[ss]][[cn]][['Case']] <- case_bfit
        }   
    }
    out
}


# Filter non-zero correlations by checking if 0 is included in the Bayesian credible interval
# conf_alpha --- 0.05 correponds to 95% credible interval
pos_nonzero_corr <- function(bfit, conf_alpha=0.05, eval_conver=TRUE){
    pos_summary <- get_banocc_output(bfit, conf_alpha = conf_alpha, eval_convergence = eval_conver)
    
    # Helper to convert a named symmetric matrix to a dataframe, using a trick 
    # to keep only upper triangular matrix values
    to_df <- function(A, val_name){
        A[lower.tri(A, diag = TRUE)] <- NA  # assign lower tri to NA, diag inclusive
        reshape2::melt(A, varnames = c("source", "target"), value.name = val_name) %>%
            na.omit() %>%
            as_tibble()
    }
    lower <- to_df(pos_summary$CI.hpd$lower, "lower")
    upper <- to_df(pos_summary$CI.hpd$upper, "upper")
    med <- to_df(pos_summary$Estimates.median, "median")
    df <- inner_join(lower, upper, by = c("source", "target")) %>%
        inner_join(med, by = c("source", "target"))
    # Filter non-zero correlations by checking if 0 is included in the credible interval
    df %>%
        filter(upper < 0 | lower > 0)
}


# Combine 'read_bfits' and 'pos_nonzero_corr'
get_nonzero_corr_all <- function(sample_sizes, case_number, conf_alpha = 0.05){
    out <- read_bfits(sample_sizes, case_number)
    for (ss in as.character(sample_sizes)) {
        for (cn in as.character(case_number)) {
            out[[ss]][[cn]][['Ctrl']] <- pos_nonzero_corr(out[[ss]][[cn]][['Ctrl']], conf_alpha)
            out[[ss]][[cn]][['Case']] <- pos_nonzero_corr(out[[ss]][[cn]][['Case']], conf_alpha)
            # Combine control and case group results
            out[[ss]][[cn]] <- bind_rows(out[[ss]][[cn]], .id = "group")
        }
        out[[ss]] <- bind_rows(out[[ss]], .id = "case_number")
    }
    out <- bind_rows(out, .id = "sample_size")
    out
}


# Read DGCNN node importance results computed using the greedy search heuristic
source_python("load_imp_res.py")
read_node_imp <- function(sample_sizes, case_number){
    out <- vector("list", length(sample_sizes))
    names(out) <- as.character(sample_sizes)
    
    for (ss in as.character(sample_sizes)) {
        out[[ss]] <- vector("list", length(case_number))
        names(out[[ss]]) <- as.character(case_number)
        for (cn in as.character(case_number)) {
            imp_res <- load_imp_res(paste0(paste("imp", ss, cn, sep="_"), ".pkl"))
            n_runs <- length(imp_res) # DGCNN runs per sample size per case number
            
            out[[ss]][[cn]] <- vector("list", n_runs)
            for (i in seq_along(imp_res)) {
                n_top_nodes <- length(imp_res[[i]]) # top N important nodes 
                
                out[[ss]][[cn]][[i]] <- vector("list", n_top_nodes)
                for (nn in seq_along(imp_res[[i]])) {
                    node_comb <- map_chr(imp_res[[i]][[nn]][[1]], ~paste(., collapse = "|"))
                    lor <- array_branch(imp_res[[i]][[nn]][[2]], 2) 
                    
                    # k-node results
                    out[[ss]][[cn]][[i]][[nn]] <- tibble(k = nn, node_comb = node_comb, lor = lor)
                }
                # Combine all k-node results for k=1...N
                out[[ss]][[cn]][[i]] <- bind_rows(out[[ss]][[cn]][[i]])
            }
            # Combinie all DGCNN runs
            out[[ss]][[cn]] <- bind_rows(out[[ss]][[cn]], .id = "run") 
        }
        # Combine all cases
        out[[ss]] <- bind_rows(out[[ss]], .id = "case_number")
    }
    # Combine all sample sizes
    out <- bind_rows(out, .id = "sample_size")
    out
}

# Show 1-node importance in a boxplot of absolute LOR
plot_1node_imp <- function(node_imp_dgcnn, sample_size, case_number, ntop = 5){
    ss <- as.character(sample_size)
    cn <- as.character(case_number)
    
    imp_pooled <- node_imp_dgcnn %>%
        filter(sample_size == ss, case_number == cn) %>%
        filter(k == 1) %>%
        group_by(node_comb) %>% 
        summarise(lor_pooled = list(flatten_dbl(lor))) %>%
        mutate(med_abs_lor = map_dbl(lor_pooled, ~ median(abs(.x)))) %>%
        slice_max(med_abs_lor, n = ntop) %>%
        mutate(node_comb = factor(node_comb, levels = node_comb)) %>%
        unnest(lor_pooled)
    
    ggplot(imp_pooled, aes(node_comb, abs(lor_pooled))) +
        geom_boxplot() +
        labs(x = "", y = "abs(LOR)", 
             title = paste("Sample size:", sample_size, " Case:", case_number))
}


# Show n-node frequency across all DGCNN runs in a barplot
plot_freq_nnode <- function(node_imp_dgcnn, sample_size, case_number, n){
    if (n > max(node_imp_dgcnn$k)) {
        stop("Increase N in N-node importance calculation.")
    }
    
    ss <- as.character(sample_size)
    cn <- as.character(case_number)
    
    df <- node_imp_dgcnn %>%
        filter(sample_size == ss, case_number == cn) %>%
        # Best n-node is contained in the first n nodes of (n+1)-node
        filter(k == n+1) %>%
        mutate(node_comb = str_split(node_comb, "\\|")) %>%
        group_by(run) %>%
        slice_head(n = 1) %>%
        summarise(n_node = map(node_comb, ~ .x[1:n])) %>%
        unnest(n_node)
    
    ggplot(df, aes(n_node)) +
        geom_bar() +
        labs(x = "", title = paste("Sample size:", sample_size, " Case:", case_number),
             subtitle = paste0(n, "-node"))
}

node_imp_dgcnn <- read_node_imp(sample_sizes = c(100, 500, 1000), case_number = 1:3)
p_1node_100_1 <- plot_1node_imp(node_imp_dgcnn, 100, 1)
p_freq_nnode_100_1 <- plot_freq_nnode(node_imp_dgcnn, 100, 1, 3)

# 1-node boxplot across all samples sizes and case numbers
sample_sizes <- c(100, 500, 1000)
case_number <- 1:3
p_1node <- vector("list", length(sample_sizes)*length(case_number))
i <- 1
for (ss in sample_sizes) {
    for (cn in case_number) {
        p_1node[[i]] <- plot_1node_imp(node_imp_dgcnn, ss, cn, ntop = 5) +
            theme(axis.text.x=element_text(angle = 90, vjust = 0.5, face = "bold"))
        i <- i + 1
    }
}
wrap_plots(p_1node, nrow = 3)
ggsave("p_1node.png", width = 10, height = 9, units = "in", dpi = 350)


# n-node frequency barplot across all samples sizes and case numbers
p_freq_nnode <- vector("list", length(sample_sizes)*length(case_number))
i <- 1
for (ss in sample_sizes) {
    # Highlight nodes used for association spike-in
    color_tbl <- tibble(node = paste0("Feature", 1:100), color = "grey30")
    if (ss == 100) {
        color_tbl$color[c(98, 49, 27, 82, 41)] <- "yellow"
    }
    if (ss == 500) {
        color_tbl$color[c(30, 36, 50, 57, 68)] <- "forestgreen"
    }
    if (ss == 1000) {
        color_tbl$color[c(61, 7, 23, 62, 97)] <- "blue"
    }
    
    for (cn in case_number) {
        if (cn == 1) {
            p_freq_nnode[[i]] <- plot_freq_nnode(node_imp_dgcnn, ss, cn, n = 3)
        }
        else if (cn == 2) {
            p_freq_nnode[[i]] <- plot_freq_nnode(node_imp_dgcnn, ss, cn, n = 4)
        }
        else {
            p_freq_nnode[[i]] <- plot_freq_nnode(node_imp_dgcnn, ss, cn, n = 5)
        }
        p_freq_nnode[[i]] <- p_freq_nnode[[i]] +
            geom_bar(aes(fill = n_node)) +
            scale_fill_manual(breaks = color_tbl$node, values = color_tbl$color) +
            theme(axis.text.x=element_text(angle = 90, vjust = 0.5, face = "bold"),
                  legend.position = "none")
        i <- i + 1
    }
}
wrap_plots(p_freq_nnode, nrow = 3)
ggsave("p_freq_nnode.png", width = 10, height = 9, units = "in", dpi = 350)



# Inspect BAnOCC infered correlations -------------------------------------

library(igraph)

# Visualize posterior median correlation networks inferred from BAnOCC
# Edge width prop to corr strength
# Edge color indicates corr sign
viz_corr_network <- function(nonzero_cor, sample_size, case_number, corr_thresh, grp = "Case"){
    ss <- as.character(sample_size)
    cn <- as.character(case_number)
    
    ig_df <- nonzero_cor %>%
        filter(sample_size == ss, case_number == cn) %>%
        filter(abs(median) > corr_thresh) 

    if (!(grp %in% ig_df$group)) {
        stop("No edges present, try lower corr threshold.")
    }
    
    ig <- ig_df %>% 
        filter(group == grp) %>%
        select(source, target, median) %>%
        graph_from_data_frame(directed = FALSE)
        
    ew <- abs(E(ig)$median*5) # edge width
    l <- layout_nicely(ig)
    # l <- layout_with_graphopt(ig, charge = 0.05, mass = 50) # layout
    
    plot(ig, edge.width = ew, edge.color = ifelse(E(ig)$median > 0, "lightblue2", "tomato"),
         vertex.label.color="black", vertex.label.dist = 2, vertex.size = 17,
         layout = l, main = paste("Sample size:", sample_size, " Case number:", case_number))
        # vertex.label.cex=1.5
}


nonzero_cor <- get_nonzero_corr_all(sample_sizes = c(100, 500, 1000), case_number = 1:3)
viz_corr_network(nonzero_cor, 100, 1, 0.3)

png("p_corr_net_sim.png", width = 10.5, height = 9, units = "in", bg = "transparent", res = 350)
graphics::layout(mat = matrix(1:9, nrow = 3, byrow = TRUE))
for (ss in sample_sizes) {
    for (cn in case_number) {
        viz_corr_network(nonzero_cor, ss, cn, corr_thresh = 0.2)
    }
}
dev.off()


# graphics::layout(mat = matrix(c(1,2,3), nrow = 1))
# viz_corr_network(nonzero_cor, runid = "r12e", corr_thresh = 0.4, title = "S1")
# viz_corr_network(nonzero_cor, runid = "r12f", corr_thresh = 0.4, title = "S2")
# viz_corr_network(nonzero_cor, runid = "r12g", corr_thresh = 0.4, title = "S3")
# dev.off()



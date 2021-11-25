# Process results of simulated data from BAnOCC and DGCNN

library(tidyverse)
library(banocc)
library(reticulate)
library(tidytext) # for order factor within facets in ggplot
library(patchwork)

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
        labs(x = "", y = "abs(LOR)")
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
        labs(x = "")
}

node_imp_dgcnn <- read_node_imp(sample_sizes = c(100, 500, 1000), case_number = 1:3)
p_1node_100_1 <- plot_1node_imp(node_imp_dgcnn, 100, 1)
p_freq_nnode_100_1 <- plot_freq_nnode(node_imp_dgcnn, 100, 1, 3)




# Inspect BAnOCC infered correlations -------------------------------------

library(igraph)

viz_corr_network <- function(nonzero_cor, sample_size, case_number, corr_thresh, 
                             grp = "Case", show_comm = FALSE, title = ""){
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
    l <- layout_with_graphopt(ig, charge = 0.05, mass = 50) # layout
    if (show_comm) { # Run community detection
        clp <- cluster_optimal(ig)
        plot(clp, ig, edge.width = ew, 
             edge.color = ifelse(E(ig)$median > 0, "blue", "orangered"),
             vertex.label.color="black", layout = l, main = title, frame = TRUE)
    }
    else {
        # Edge width prop to strength of corr
        # Edge color indicates pos or neg corr
        plot(ig, edge.width = ew,
             edge.color = ifelse(E(ig)$median > 0, "lightblue", "salmon"),
             vertex.label.color="black", layout = l, main = title, frame = TRUE)
            # vertex.label.cex=1.5
    }
}


nonzero_cor <- get_nonzero_corr_all(sample_sizes = c(100, 500, 1000), case_number = 1:3)
viz_corr_network(nonzero_cor, 100, 1, 0.3)


# graphics::layout(mat = matrix(c(1,2,3), nrow = 1))
# viz_corr_network(nonzero_cor, runid = "r12e", corr_thresh = 0.4, title = "S1")
# viz_corr_network(nonzero_cor, runid = "r12f", corr_thresh = 0.4, title = "S2")
# viz_corr_network(nonzero_cor, runid = "r12g", corr_thresh = 0.4, title = "S3")
# dev.off()



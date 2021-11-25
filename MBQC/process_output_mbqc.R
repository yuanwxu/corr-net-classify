# Process results of MBQC data from BAnOCC and DGCNN

library(tidyverse)
library(banocc)
library(reticulate)
library(tidytext) # for order factor within facets in ggplot
library(patchwork)

# Read banocc fits for control and case (CRC) 
read_bfits <- function(){
    out <- vector("list")
    out[['Ctrl']] <- readRDS("fit_mbqc_ctrl.rds")
    out[['Case']] <- readRDS("fit_mbqc_case.rds")
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
get_nonzero_corr_all <- function(conf_alpha = 0.05){
    out <- read_bfits()
    out[['Ctrl']] <- pos_nonzero_corr(out[['Ctrl']], conf_alpha)
    out[['Case']] <- pos_nonzero_corr(out[['Case']], conf_alpha)
    out <- bind_rows(out, .id = "group")
    out
}


# Read DGCNN node importance results computed using the greedy search heuristic
source_python("load_imp_res.py")
read_node_imp <- function(){
    imp_res <- load_imp_res("imp_mbqc.pkl")
    n_runs <- length(imp_res) # DGCNN runs
    out <- vector("list", n_runs)
    
    for (i in seq_along(imp_res)) {
        n_top_nodes <- length(imp_res[[i]]) # top N important nodes 
        
        out[[i]] <- vector("list", n_top_nodes)
        for (nn in seq_along(imp_res[[i]])) {
            node_comb <- map_chr(imp_res[[i]][[nn]][[1]], ~paste(., collapse = "|"))
            lor <- array_branch(imp_res[[i]][[nn]][[2]], 2) 
            
            # k-node results
            out[[i]][[nn]] <- tibble(k = nn, node_comb = node_comb, lor = lor)
        }
        # Combine all k-node results for k=1...N
        out[[i]] <- bind_rows(out[[i]])
    }
    
    # Combinie all DGCNN runs
    out <- bind_rows(out, .id = "run") 
    out
}


# Show 1-node importance in a boxplot of absolute LOR
plot_1node_imp <- function(node_imp_dgcnn, ntop = 5){
    imp_pooled <- node_imp_dgcnn %>%
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
plot_freq_nnode <- function(node_imp_dgcnn, n){
    if (n > max(node_imp_dgcnn$k)) {
        stop("Increase N in N-node importance calculation.")
    }
    
    df <- node_imp_dgcnn %>%
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

node_imp_dgcnn <- read_node_imp()
p_1node <- plot_1node_imp(node_imp_dgcnn, ntop = 5)
p_freq_3node <- plot_freq_nnode(node_imp_dgcnn, n = 3)




# Inspect BAnOCC infered correlations -------------------------------------

library(igraph)

# Visualize posterior median correlation networks involving important nodes selected by the
# greedy search heuristic. These nodes (specified via 'nodes_to_show') are highlighted (yellow)
# Node size prop to sqrt of degree.
# Edge width prop to corr strength
# Edge color indicates corr sign
viz_corr_network <- function(nonzero_cor, nodes_to_show, corr_thresh, 
                             grp = "Case", show_comm = FALSE, title = ""){
    ig_df <- nonzero_cor %>%
        filter((source %in% nodes_to_show) | (target %in% nodes_to_show)) %>%
        filter(abs(median) > corr_thresh) 

    if (!(grp %in% ig_df$group)) {
        stop("No edges present, try lower corr threshold.")
    }
    
    ig <- ig_df %>% 
        filter(group == grp) %>%
        select(source, target, median) %>%
        graph_from_data_frame(directed = FALSE)
    
    vc <- ifelse(V(ig)$name %in% nodes_to_show, "yellow", "white") # node color
    vs <- 6*sqrt(degree(ig)) # node size    
    ew <- abs(E(ig)$median*5) # edge width
    l <- layout_with_graphopt(ig, charge = 0.05, mass = 50) # layout
    if (show_comm) { # Run community detection
        clp <- cluster_optimal(ig)
        plot(clp, ig, col = vc, edge.width = ew, 
             edge.color = ifelse(E(ig)$median > 0, "blue", "orangered"), vertex.size = vs,
             vertex.label.color="black", layout = l, main = title)
    }
    else {
        plot(ig, col = vc, edge.width = ew, 
             edge.color = ifelse(E(ig)$median > 0, "blue", "orangered"), vertex.size = vs,
             vertex.label.color="black", layout = l, main = title)
            # vertex.label.cex=1.5
    }
}


nonzero_cor <- get_nonzero_corr_all()
viz_corr_network(nonzero_cor, corr_thresh = 0.4, grp = "Case", show_comm = TRUE)
viz_corr_network(nonzero_cor, corr_thresh = 0.4, grp = "Ctrl", show_comm = TRUE)



# Microbial correlation network classification
# Simulated data with known correlation structure

library(tidyverse)
library(SparseDOSSA2)

N_SAMPLE <- c(100, 500, 1000)
N_FEATURE <- 100

# Simulate control group feature abundances
control_sim <- function(){
    out <- vector("list", length(N_SAMPLE))
    for (i in seq_along(N_SAMPLE)) {
        set.seed(N_SAMPLE[i] + 1)
        sim <- SparseDOSSA2("Stool", n_sample = N_SAMPLE[i], n_feature = N_FEATURE) 
        out[[i]] <- sim
    }
    names(out) <- as.character(N_SAMPLE)
    out
}

# Simulate case group feature abundances
case_sim <- function(n_cases = 3){
    # Features to spike depend on sample size: selecting most abund features 
    # present in at least 50% samples for each sample size
    features_to_spike <- function(nsamp){
        cat("Baseliine simulation with", nsamp, "samples. \n")
        set.seed(nsamp + 1)
        null_sim <- SparseDOSSA2("Stool", n_sample = nsamp, n_feature = N_FEATURE)
        null_mat <- null_sim$simulated_matrices$a_null
        features <- null_mat %>% 
            as_tibble(rownames = NA) %>%
            tibble::rownames_to_column("Feature") %>%
            rowwise("Feature") %>%
            mutate(total_prev = sum(c_across() > 0.001), total_abund = sum(c_across())) %>%
            ungroup() %>%
            filter(total_prev > nsamp/2) %>%
            arrange(desc(total_abund)) %>%
            pull(Feature)
        features
    }
    
    # Construct spike-in configuration data frame matching SparseDOSSA API
    get_spike_df <- function(to_spike, case_number, eff_size = 1){
        n_spiked <- case_number + 1  # number of spiked associations
        feature_pairs <- str_split(paste(to_spike[1], to_spike[-1]), pattern = " ")
        eff_spec <- rep(c(1,1), n_spiked)
        eff_spec[seq(4, length(eff_spec), 4)] <- -1
        data.frame(metadata_datum = rep(1:n_spiked, each = 2, times = 2),
                   feature_spiked = rep(flatten_chr(feature_pairs), 2),
                   associated_property = rep(c("abundance", "prevalence"), each=2*n_spiked),
                   effect_size = rep(eff_spec, 2))
    }
    
    out <- vector("list", length(N_SAMPLE))
    names(out) <- as.character(N_SAMPLE)
    for (i in seq_along(N_SAMPLE)) {
        out[[i]] <- vector("list", n_cases)
        names(out[[i]]) <- paste0("Case", 1:n_cases)
        features <- features_to_spike(N_SAMPLE[i])
        for (j in seq_len(n_cases)) {
            # Pick first few features to spike: 
            # one feature will be chosen to be correlated with the rest
            to_spike <- features[1:(j+2)]
            
            cat("Simulation with spiked features", "[", to_spike, "]",
                "for case", j, "with", N_SAMPLE[i], "samples... \n")
            
            # Get metadata matrix
            # SparseDOSSA can't handle identical covariates: when (A, B) spiked with z,
            # (A, C) cannot be spiked with the same z. We use a trick by adding small Gaussian 
            # noise to z to form new covariates.
            z <- rnorm(N_SAMPLE[i])
            Z <- cbind(z, z + array(rnorm(N_SAMPLE[i]*j, sd = 0.01), dim = c(N_SAMPLE[i], j)))
            spike_df <- get_spike_df(to_spike, j)
            
            set.seed(N_SAMPLE[i] + 1)
            sim <- SparseDOSSA2("Stool", n_sample = N_SAMPLE[i], n_feature = N_FEATURE,
                                spike_metadata = spike_df,
                                metadata_matrix = Z) 
            out[[i]][[j]] <- sim
        }
    }
    out
}

# Apply filter to remove features with rel abundance less than some threshold in at least
# x% samples. Combine control and case group first then filter as want to ensure same 
# features in both groups. Allow writing results to disk on the go.
filter_feature <- function(sim0, sim1, write_results = TRUE, 
                           min_rel = 0.0001, min_perc_preval = 0.3){
    stopifnot(length(sim0) == length(sim1))
    for (i in seq_along(sim1)) { # over sample sizes
        for (j in seq_along(sim1[[i]])) { # over cases
            rel0 <- sim0[[i]]$simulated_matrices$rel %>% t()
            rel1 <- sim1[[i]][[j]]$simulated_matrices$rel %>% t()
            new_features <- data.frame(rbind(rel0, rel1)) %>%
                select(where(~ sum(.x > min_rel) > min_perc_preval*length(.x))) %>%
                names()

            # Recompute relative abundances with new features, add filtered abs and rel 
            # abundances to relevant slot in the simulation object
            a0_filtered <- sim0[[i]]$simulated_matrices$a_null[new_features, ]
            sim1[[i]][[j]]$simulated_matrices[["a_null_filtered"]] <- a0_filtered
            rel0_filtered <- t(t(a0_filtered) / colSums(a0_filtered))
            
            a1_filtered <- sim1[[i]][[j]]$simulated_matrices$a_spiked[new_features, ]
            sim1[[i]][[j]]$simulated_matrices[["a_spiked_filtered"]] <- a1_filtered
            rel1_filtered <- t(t(a1_filtered) / colSums(a1_filtered))
            sim1[[i]][[j]]$simulated_matrices[["rel_filtered"]] <- rel1_filtered
            
            if (write_results) {
                X_ctrl <- rel0_filtered %>% t() %>% data.frame()
                X_case <- rel1_filtered %>% t() %>% data.frame()
                write_csv(X_ctrl, paste0("X_ctrl_", N_SAMPLE[i], "_", j, ".csv"))
                write_csv(X_case, paste0("X_case_", N_SAMPLE[i], "_", j, ".csv"))
            }
        }
    }
    sim1
}

# Extract spike configurations
extract_spike_config <- function(sim){
    config <- vector("list", length(sim))
    names(config) <- names(sim)
    
    for (i in seq_along(sim)) {
        dfs <- map(sim[[i]], c("spike_metadata", "spike_metadata"))
        config[[i]] <- bind_rows(dfs, .id = "Case")
    }
    bind_rows(config, .id = "sample_size") %>%
        mutate(sample_size = as.integer(sample_size))
}

# Simple wrapper of corrplot
plot_corr <- function(sim, what, cor.method = "spearman"){
    stopifnot(what %in% c("a_null", "a_spiked", "a_null_filtered", 
                          "a_spiked_filtered", "rel_filtered"))
    corrplot::corrplot(cor(t(sim$simulated_matrices[[what]]), method = cor.method))
}


sim0 <- control_sim()
sim1 <- case_sim()
sim1 <- filter_feature(sim0, sim1) 

plot_corr(sim1$`100`$Case1, "a_null_filtered")
plot_corr(sim1$`100`$Case1, "a_spiked_filtered")

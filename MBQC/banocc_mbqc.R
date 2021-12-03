# Run BAnOCC for MBQC baseline data

library(tidyverse)
library(phyloseq)
library(banocc)

physeq <- readRDS("physeq_mbqc.rds")
X <- otu_table(physeq)@.Data %>% t() %>% as.data.frame()

# Use microbe names instead of OTU number. Because some OTUs unclassifed at 'Genus' were 
# still kept rather than thrown away, their names will by replaced by the next available rank.
# For example, if 'Family' is known, it will be used.
rename_taxa <- function(X, physeq){
    stopifnot(all(str_detect(names(X), "OTU")))
    stopifnot(all(names(X) == taxa_names(physeq)))
    
    new_names <- as.character(tax_table(physeq)[, "Genus"])
    ranks <- rank_names(physeq)
    i <- which(ranks == "Genus")
    unc_ind <- which(str_detect(new_names, "unclassified"))
    while (length(unc_ind) > 0) { # move to next higher rank
        i <- i - 1
        if (i < 1) {
            stop("Found unknown taxa")
        }
        new_names[unc_ind] <- as.character(tax_table(physeq)[unc_ind, ranks[i]])
        unc_ind <- which(str_detect(new_names, "unclassified"))
    }
    if (anyDuplicated(new_names)) {
        warning("Found identical taxa")
    }
    names(X) <- new_names
    X
}

X <- rename_taxa(X, physeq)

# Normalize
X <- X / rowSums(X)

# Add labels
X$status <- get_variable(physeq, "health_status")

# Split into sub-tables of Healthy and CRC
X_mbqc_ctrl <- X %>% filter(status == "Healthy") %>% select(-status)
X_mbqc_case <- X %>% filter(status == "CRC case") %>% select(-status)




# Run BAnOCC for MBQC control and cases. 
# banocc_params --- list of banocc model parameters
# hmc_params --- list Hamiltonian Monte Carlo parameters
# save_fits --- whether to save banocc fits
run_banocc_mbqc <- function(X_ctrl, X_case, banocc_params, hmc_params, save_fits = TRUE){
    compiled_banocc_model <- rstan::stan_model(model_code = banocc::banocc_model) 
    cat("Running banocc for MBQC control", "\n") 
    time0 <- Sys.time()
    fit_ctrl <- run_banocc(compiled_banocc_model = compiled_banocc_model,
                           C = X_ctrl,
                           n = rep(banocc_params$n, ncol(X_ctrl)),
                           L = banocc_params$L * diag(ncol(X_ctrl)),
                           a = banocc_params$a,
                           b = banocc_params$b,
                           iter = hmc_params$iter,
                           thin = hmc_params$thin,
                           warmup = hmc_params$warmup,
                           cores = 4,
                           control = list(adapt_delta = 0.95))
    time1 <- Sys.time()
    cat("Finished banocc for MBQC control", "[", time1-time0, "]", "\n")
    
    cat("Running banocc for MBQC case", "\n")
    time0 <- Sys.time()
    fit_case <- run_banocc(compiled_banocc_model = compiled_banocc_model,
                           C = X_case,
                           n = rep(banocc_params$n, ncol(X_case)),
                           L = banocc_params$L * diag(ncol(X_case)),
                           a = banocc_params$a,
                           b = banocc_params$b,
                           iter = hmc_params$iter,
                           thin = hmc_params$thin,
                           warmup = hmc_params$warmup,
                           cores = 4,
                           control = list(adapt_delta = 0.95))
    time1 <- Sys.time()
    cat("Finished banocc for MBQC case", "[", time1-time0, "]", "\n")
    
    if (save_fits) {
        saveRDS(fit_ctrl, "fit_mbqc_ctrl.rds")
        saveRDS(fit_case, "fit_mbqc_case.rds")
    }
            
    W_ctrl <- get_pos_W(fit_ctrl)
    W_case <- get_pos_W(fit_case)
    write_csv(W_ctrl, "W_mbqc_ctrl.csv")
    write_csv(W_case, "W_mbqc_case.csv")
}


# Get posterior log-basis correlation networks
# convert --- if TRUE, will return a data frame with 4 columns: graph index, source,
#             target, weight. If False, will return a 3-d array with the first dim 
#             sample index and the remaining dims equal to that of the corr. matrix
get_pos_W <- function(b_fit, convert = TRUE){
    pos_samples <- rstan::extract(b_fit$Fit)
    # This line from the source of banocc::get_banocc_output
    W <- aperm(array(apply(pos_samples$O, 1, function(Oi){
        cov2cor(solve(matrix(Oi, ncol=sqrt(length(Oi)))))
    }), dim=dim(pos_samples$O)[c(3, 2, 1)]),
    perm=c(3, 2, 1))
    
    dimnames(W) <- list(NULL, colnames(b_fit$Data$C), colnames(b_fit$Data$C))
    if (!convert) {
        return(W)
    }
    # Helper to convert a named symmetric matrix to a dataframe, using a trick 
    # to keep only upper triangular matrix values
    to_df <- function(A){
        A[lower.tri(A, diag = TRUE)] <- NA  # assign lower tri to NA, diag inclusive
        reshape2::melt(A, varnames = c("source", "target"), value.name = "weight") %>%
            na.omit() %>%
            as_tibble()
    }
    map_dfr(array_branch(W, 1), to_df, .id = "graph_id")
}


# Run BAnOCC
banocc_params <- list(n = 0, # will be translated to rep(0, P) where P number of features
                      L = 30, # translated to 30 * diag(P)
                      a = 0.5,
                      b = 3)
hmc_params <- list(iter = 8000,
                   warmup = 4000,
                   thin = 20) # 800 posterior corr networks for control and 800 for case
run_banocc_mbqc(X_mbqc_ctrl, X_mbqc_case, banocc_params, hmc_params)

print(warnings())

# change adapt_delta parameter in stan to 0.95 ? (should change for simulated data )

# Run BAnOCC for simulated data

library(tidyverse)
library(banocc)

# Load synthetic OTU table with relative abundances
load_synthetic_otu_table <- function(sample_sizes, case_number){
    out <- vector("list", length(sample_sizes))
    names(out) <- as.character(sample_sizes)
    
    for (ss in as.character(sample_sizes)) {
        out[[ss]] <- vector("list", length(case_number))
        names(out[[ss]]) <- as.character(case_number)
        for (cn in as.character(case_number)) {
            ctrl_name <- paste0(paste("X", "ctrl", ss, cn, sep = "_"), ".csv")
            case_name <- paste0(paste("X", "case", ss, cn, sep = "_"), ".csv")
            ctrl_df <- read_csv(ctrl_name, show_col_types = FALSE)
            case_df <- read_csv(case_name, show_col_types = FALSE)
            
            # Check dimensions match between control and case
            stopifnot(nrow(ctrl_df) == nrow(case_df), ncol(ctrl_df) == ncol(case_df))
            
            out[[ss]][[cn]] <- vector("list")
            out[[ss]][[cn]][['Ctrl']] <- ctrl_df
            out[[ss]][[cn]][['Case']] <- case_df
        }   
    }
    out
}

# Run BAnOCC for all simulated control and cases. 
# Write posterior log-basis correlation networks on the go.
# sim_data --- object returned by load_synthetic_otu_table() containing all simulated data
# banocc_params --- list of banocc model parameters
# hmc_params --- list Hamiltonian Monte Carlo parameters
# save_fits --- whether to save banocc fits
run_banocc_sim <- function(sim_data, banocc_params, hmc_params, save_fits = TRUE){
    compiled_banocc_model <- rstan::stan_model(model_code = banocc::banocc_model) 
    for (i in seq_along(sim_data)) {
        for (j in seq_along(sim_data[[i]])) {
            cat("Running banocc for control group --- Sample size:", names(sim_data)[i], 
                "Case number:", names(sim_data[[i]])[j], "\n") 
            time0 <- Sys.time()
            C0 <- sim_data[[i]][[j]]$Ctrl
            fit_ctrl <- run_banocc(compiled_banocc_model = compiled_banocc_model,
                                   C = C0,
                                   n = rep(banocc_params$n, ncol(C0)),
                                   L = banocc_params$L * diag(ncol(C0)),
                                   # This won't work: 'C' is a function name in {stats}
                                   # C = sim_data[[i]][[j]]$Ctrl,
                                   # n = rep(banocc_params$n, ncol(C)),
                                   # L = banocc_params$L * diag(ncol(C)),
                                   a = banocc_params$a,
                                   b = banocc_params$b,
                                   iter = hmc_params$iter,
                                   thin = hmc_params$thin,
                                   warmup = hmc_params$warmup,
                                   cores = 4)
            time1 <- Sys.time()
            cat("Finished banocc for control group --- Sample size:", names(sim_data)[i], 
                "Case number:", names(sim_data[[i]])[j], "[", time1-time0, "]", "\n")
            cat("Running banocc for case group --- Sample size:", names(sim_data)[i], 
                "Case number:", names(sim_data[[i]])[j], "\n")
            time0 <- Sys.time()
            C1 <- sim_data[[i]][[j]]$Case
            fit_case <- run_banocc(compiled_banocc_model = compiled_banocc_model,
                                   C = C1,
                                   n = rep(banocc_params$n, ncol(C1)),
                                   L = banocc_params$L * diag(ncol(C1)),
                                   a = banocc_params$a,
                                   b = banocc_params$b,
                                   iter = hmc_params$iter,
                                   thin = hmc_params$thin,
                                   warmup = hmc_params$warmup,
                                   cores = 4)
            time1 <- Sys.time()
            cat("Finished banocc for case group --- Sample size:", names(sim_data)[i], 
                "Case number:", names(sim_data[[i]])[j], "[", time1-time0, "]", "\n")
            
            if (save_fits) {
                saveRDS(fit_ctrl, paste0("fit_ctrl_", names(sim_data)[i], "_",
                                         names(sim_data[[i]])[j], ".rds"))
                saveRDS(fit_case, paste0("fit_case_", names(sim_data)[i], "_",
                                         names(sim_data[[i]])[j], ".rds"))
            }
            
            W_ctrl <- get_pos_W(fit_ctrl)
            W_case <- get_pos_W(fit_case)
            write_csv(W_ctrl, paste0("W_ctrl_", names(sim_data)[i], "_",
                                     names(sim_data[[i]])[j], ".csv"))
            write_csv(W_case, paste0("W_case_", names(sim_data)[i], "_",
                                     names(sim_data[[i]])[j], ".csv"))
        }
    }
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


sim_otu_tab <- load_synthetic_otu_table(sample_sizes = c(100, 500, 1000), case_number = 1:3)
banocc_params <- list(n = 0, # will be translated to rep(0, P) where P number of features
                      L = 30, # translated to 30 * diag(P)
                      a = 0.5,
                      b = 5)
hmc_params <- list(iter = 5000,
                   warmup = 2500,
                   thin = 20)
run_banocc_sim(sim_otu_tab, banocc_params, hmc_params)
                      
# change adapt_delta parameter in stan to 0.95 (divergent transitions encountered)
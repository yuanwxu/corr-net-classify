# Process integrated OTU table of MBQC baseline data

library(tidyverse)
library(data.table)
library(phyloseq)

mbqc <- fread("mbqc_integrated_otus_transpose.csv")

# Keep CRC and Healthy samples
mbqc <- mbqc %>% filter(health_status %in% c("CRC case", "Healthy"))

# Use 'Fresh' and 'Freeze-dried' specimen types
mbqc <- mbqc %>% filter(specimen_type_collapsed %in% c("Fresh", "Freeze-dried"))

# Select taxa and some relevant variables
mbqc <- mbqc %>% select(starts_with("k__"), sample, health_status, age)


get_tax_mat <- function(tax_names){
  rk_names <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
  tax_names <- tax_names %>% 
    str_split('\\|') %>% 
    tibble::enframe() %>%
    filter(lengths(value) %in% c(7, 8)) %>% # remove taxa with unknown genus
    mutate(value = map(value, ~ .x[1:7]))
                       
  tax_mat <- tax_names %>%
    mutate(name = paste0("OTU", name)) %>%
    mutate(value = map(value, ~ setNames(., rk_names))) %>%
    unnest_wider(value) %>%
    column_to_rownames("name") %>%
    as.matrix()
  tax_mat
}

get_otu_mat <- function(mbqc_df, tax_names){
  otu_mat <- mbqc_df %>% 
    select(all_of(tax_names)) %>%
    as.matrix() %>%
    t()
  colnames(otu_mat) <- mbqc_df$sample
  rownames(otu_mat) <- paste0("OTU", 1:nrow(otu_mat))
  otu_mat
}

tax_names <- setdiff(names(mbqc), c("sample", "health_status", "age"))
otu_mat <- get_otu_mat(mbqc, tax_names)
tax_mat <- get_tax_mat(tax_names)
physeq <- phyloseq(otu_table(otu_mat, taxa_are_rows = TRUE), tax_table(tax_mat))
physeq


# Prevalance and abundance filtering, first agglomerate over taxa to 'Genus' level where possible,
# if a taxon is unidentified at 'Genus' level it is still kept.
preproc_physeq <- function(physeq, min_rel_abund = 0.0001, min_perc_sam = 0.4, glom = "Genus"){
  nsam <- nsamples(physeq)
  ntax <- ntaxa(physeq)
  # Remove taxa absent from all samples
  physeq <- physeq %>%
    filter_taxa(function(x) any(x>0), prune = TRUE) 
  
  # Remove samples with zero sum of taxa abundance
  physeq <- prune_samples(sample_sums(physeq) > 0, physeq) 
  
  ## Agglomerate taxa
  physeq <- tax_glom(physeq, glom, NArm = FALSE)
  
  ## Keep features present in at least 'min_perc_sam' samples and rel abund greater than 
  # 'min_rel_abund'
  tokeep <- physeq %>%
    transform_sample_counts(function(x) x / sum(x)) %>% # rel abundance
    filter_taxa(function(x) sum(x>=min_rel_abund) >= min_perc_sam*length(x))
  physeq <- prune_taxa(tokeep, physeq)
  
  # Prune again samples with zero sum of taxa
  physeq <- prune_samples(sample_sums(physeq) > 0, physeq)
  
  if (nsamples(physeq) < nsam) {
    warning(nsam - nsamples(physeq), " samples removed for which (nearly) all taxa are absent.")
  }
  physeq
}

physeq <- preproc_physeq(physeq, min_rel_abund = 0.0001, min_perc_sam = 0.4)


# Add metadata
add_metadata_physeq <- function(physeq, mbqc_df){
  sampledata <- mbqc_df %>% 
    select(all_of(c("sample", "health_status", "age"))) %>%
    column_to_rownames("sample") %>%
    phyloseq::sample_data()
  
  merge_phyloseq(physeq, sampledata)
}

physeq <- add_metadata_physeq(physeq, mbqc)


# Save results
saveRDS(physeq, "physeq_mbqc.rds")


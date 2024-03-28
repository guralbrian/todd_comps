# Test for differential gene expression with DESeq2 
# Include composition covariates to test if they mediate DE

# Load libs
libs <- c("tidyverse", "compositions", "reshape2", "stats", "DESeq2") # list libraries here
lapply(libs, require, character.only = T)
rm(libs)

#### Loading and formatting of data ####
# Load compositions, expression, and phenotype data
decon.whole <- read.csv("data/processed/compositions/whole_samples.csv")
phenotypes <- read.csv("data/processed/bulk/phenos.csv")
bulk.all <- read.csv("data/processed/bulk/all_counts.csv", row.names = 1, check.names = F)


# Unmelt for clr transformation 
decon.wide <- decon.whole  |> 
  dplyr::select(id, CellType, Prop) |> 
  pivot_wider(names_from = "id", values_from = "Prop") |> 
  column_to_rownames("CellType") %>%
  mutate_all(as.numeric)
#decon.wide[decon.wide==0] <- 0.01

# use clr or ilr to incorp the compositions into the design matrix
comps.clr <- compositions::clr(t(decon.wide))

colnames(comps.clr) <- paste0("clr.", colnames(comps.clr))

# Prepare for DESeq2
bulk <- mutate_all(bulk.all, function(x) round(as.numeric(as.character(x)), digits = 0)) # round to integers

# Set reference factor levels for phenotypes
pheno.reorder <- phenotypes |> 
  mutate(id = as.factor(id),
         treatment = as.factor(treatment))

# Prepare sample information
sample_info <- data.frame(
  row.names = pheno.reorder$id,
  treatment = pheno.reorder$treatment
)

sample_info$treatment <- relevel(sample_info$treatment, ref = "PBS")

# Add clr transforms to sample info
sample.clr <- cbind(sample_info, comps.clr[colnames(bulk),])

#### Run DESeq2 ####
## DESeq with compositions
# Create a DESeqDataSet
dds.clr <- DESeqDataSetFromMatrix(
  countData = bulk,
  colData = sample.clr,
  design = ~ treatment  #+ clr.Macrophage
)

# Filter out lowly expressed genes
smallestGroupSize <- 3
keep <- rowSums(counts(dds.clr) >= 10) >= smallestGroupSize
dds.clr <- dds.clr[keep,]

# remove outlier
keep <- colnames(dds.clr)[!colnames(dds.clr) %in% c("PE5", "PE11", "PE8", "PE10")]
dds.clr <- dds.clr[,keep]

# Run DESeq
dds.clr <- DESeq(dds.clr)

# Save the results 
#Save outputs
if(!dir.exists("data/processed/models")){
  dir.create("data/processed/models")
}
saveRDS(dds.clr, "data/processed/models/unadjusted_de_exclude.RDS")

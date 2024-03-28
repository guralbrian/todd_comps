# Clean up txt file from todd
# make into csv with genes as row names
# check if transformed
libs <- c("tidyverse") # list libraries here
lapply(libs, function(x){library(x, character.only = T, quietly = T, verbose = F)})
rm(libs)



dat <- read.table("data/processed/bulk/2022-03-08_normalized_count_matrix.txt")

colnames(dat) <- dat[1,]
dat <- dat[-1,]

# Convert to numeric
dat[,c(1:14)] <- apply(dat[,c(1:14)], 2, as.numeric)

# Select to gene variant where duplicated
dat <- dat %>%
  #group_by(ensembl_gene_id_version) %>% 
  mutate(total = rowSums(dplyr::select(., where(is.numeric)), na.rm = TRUE)) |> 
  group_by(gene.symbol) |> 
  arrange(desc(total)) |> 
  filter(row_number() == 1) |> 
  ungroup() |> 
  subset(!is.na(gene.symbol)) |> 
  as.data.frame() |> 
  select(-c(total, ensembl.id))

# Label rows with the genes, remove the column of gene names
rownames(dat) <- dat$gene.symbol
dat <- dat |> select(-gene.symbol)

write.csv(dat, "data/processed/bulk/all_counts.csv")

# Make pheno table

pheno <- data.frame(id = colnames(dat)) |> 
  mutate(treatment = case_when(
    str_detect(id, "PBS") ~ "PBS",
    str_detect(id, "PE") ~ "PE"
  ))

# Save pheno table
write.csv(pheno, "data/processed/bulk/phenos.csv", row.names = F)

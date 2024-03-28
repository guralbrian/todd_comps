# ClusterProfiler analysis
libs <- c("tidyverse", "clusterProfiler", "DESeq2", 
          "org.Mm.eg.db", "stringr") # list libraries here
lapply(libs, require, character.only = T)
rm(libs)

# Load the results and expression matrices
res <- readRDS(paste0("data/processed/models/unadjusted_de_exclude.RDS"))  
names <- resultsNames(res)[-1]

# make a list of each result
res <- lapply(names, function(x){
  results(res, name=x) |> 
    as.data.frame()
})

names(res) <- names

# Function to perform GO enrichment with clusterProfiler
runGo <- function(data, onto){
# Get a list of significant genes
sig.genes <- data |> 
  filter(padj < 0.005 & abs(log2FoldChange) > 0.583) |> 
  row.names()

# Convert common gene names to ENSEMBLE IDs for clusterProfiler
gene.df <- bitr(sig.genes, fromType = "SYMBOL",
                toType = c("ENSEMBL"),
                OrgDb = org.Mm.eg.db)

gene.list <- bitr(row.names(data), fromType = "SYMBOL",
                toType = c("ENSEMBL"),
                OrgDb = org.Mm.eg.db)

# Check for enrichment within biological process gene clusters
ego <- enrichGO(gene          = gene.df$ENSEMBL,
                OrgDb         = org.Mm.eg.db,
                universe      = gene.list$ENSEMBL,
                keyType       = 'ENSEMBL',
                ont           = onto,
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05,
                readable      = TRUE) #|> 
        #clusterProfiler::simplify(cutoff = 0.6)
}

# Apply the function
go <- lapply(res, function(x){runGo(x, "BP")})

# Save 
if(!dir.exists("data/processed/pathway_genesets")){
  dir.create("data/processed/pathway_genesets")
}

saveRDS(go, paste0("data/processed/pathway_genesets/go_unadj_exclude_fib_0005.RDS"))

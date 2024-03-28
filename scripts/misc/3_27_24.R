# Apply the function
go <- readRDS("data/processed/pathway_genesets/go_unadj_exclude_fib_0005.RDS")

write.csv(go[[1]]@result, "data/processed/pathway_genesets/tables/unadj_exclude_fib/pe.csv")


# Load the results and expression matrices
res <- readRDS(paste0("data/processed/models/unadjusted_de_exclude.RDS"))  
names <- resultsNames(res)[-1]

# make a list of each result
res <- lapply(names, function(x){
  results(res, name=x) |> 
    as.data.frame()
})
results <- res[[1]] |>  arrange(padj)
write.csv(results, "data/processed/pathway_genesets/tables/unadj_exclude_fib/deseq.csv")

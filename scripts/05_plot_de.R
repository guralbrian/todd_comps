# Visualize DE 
libs <- c("tidyverse", "wesanderson", "ggrepel", "DESeq2","patchwork") # list libraries here
lapply(libs, require, character.only = T)
rm(libs)

# Load the results 
# Load the results and expression matrices
res <- readRDS(paste0("data/processed/models/unadjusted_de_exclude.RDS")) 
names <- resultsNames(res)[-1]

# make a list of each result
res <- lapply(names, function(x){
  results(res, name=x) |> 
    as.data.frame() |> 
    rownames_to_column(var = "gene")
})
names(res) <- names

# Set reference factor levels for phenotypes
pheno.reorder <- phenotypes |> 
      mutate(id = as.factor(id),
              treatment = as.factor(treatment))
#####
# Patch to plot from list of results contrasts 
####
plotDE <- function(x, title){
# Find top DE genes
top <- x |> 
  arrange(padj) |>
  slice_head(n = 10)

# Add conditional color formatting for significance
x <- x |> 
  mutate(significant = case_when(
    padj >  0.05 ~ FALSE,
    padj <= 0.05 & abs(log2FoldChange) >= 0.583  ~ TRUE,
    .default = FALSE
  ))

plot.fib <- ggplot(x, aes(x = log2FoldChange, y = -log10(padj), color = significant)) +
  geom_point(alpha = 1, size = 6) + 
  scale_color_manual(values = c("#999999", "#ed9209")) +
  geom_text_repel(data = top, aes(label = gene), size = 12, box.padding = unit(0.35, "lines"), 
                  force = 20, segment.linetype = 2, segment.size = 0.6, color = "black", force_pull = 0.01, min.segment.length = 0) + # label most sig genes 
  theme_minimal() +
  xlim(-5, 5) +
  labs(x = "log2(Fold Change)", y = "-log10(adjusted p-value)", 
       color = "p < 0.05 and\nFold Change > 1.5", title = title) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") + 
  geom_vline(xintercept = 0, linetype = "solid") + # add a line for p=0.05
  theme(legend.position = "bottom",
        axis.text = element_text(color = "black", size = 25),
        axis.ticks = element_blank(),
        title = element_text(size = 25),
        legend.text = element_text(size = 25),
        legend.key.size = unit(0.7, 'in'),
        legend.title = element_text(size = 28, vjust = 0.7),
        axis.title = element_text(color = "black", size = 28),
        panel.background = element_rect(color="black"),
        plot.background = element_rect(color="black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin = unit(c(1,1,1,1), "cm"))
}

p.de <- lapply(1:length(res), function(n){plotDE(res[[n]], names[[n]])})

# Save 
if(!dir.exists("results/05_plot_de")){
  dir.create("results/05_plot_de")
}
png(file = "results/05_plot_de/volcano_exclude_unadj.png",
    width = 12, 
    height = 10,
    units = "in",
    res = 600)

wrap_plots(p.de)

dev.off()



# Plot the PCA of samples
# Run PCA

# Convert to CPM
cpm <- apply(bulk,2, function(x) (x/sum(x))*1000000) |> 
  as.data.frame()

# Run PCA
pca <- prcomp(cpm)

# Get percent of variance explained by each PC
PoV <- pca$sdev^2/sum(pca$sdev^2) * 100  
PoV <- round(PoV, digits = 1)
# Merge with sample info, then plot
pca <- pca$rotation |> 
  as.data.frame() |> 
  dplyr::select(PC1, PC2) 
pca$id <- row.names(pca)

my_palette <- c("#A6CEE3", "#1F78B4", "#FDBF6F", "#FF7F00")
legend.names <- c("Sham_1","Sham_2", "MI_1", "MI_2")

pca <- pca |> left_join(pheno.reorder) |> 
  mutate(gene_treat = paste(genotype, treatment)) |> 
  filter(type == "whole") 
pca.plot <- pca |> 
  ggplot(aes(x = PC1, y = PC2, color = gene_treat)) +
  geom_point(size = 8, color = "black") +
  geom_point(size = 7) +
  geom_label_repel(aes(fill = gene_treat),
                   label = pca$id, color = "black", alpha = 0.5,
                   box.padding = 0.5, segment.curvature = -0.2,
                   segment.ncp = 3, segment.angle = 20, force = 10, 
                   max.overlaps = 8, force_pull = 0.001, size = 6,
                   nudge_x = 0.015, nudge_y = 0.02) +
  scale_color_manual(values = my_palette) +
  scale_fill_manual(values = my_palette) +
  scale_x_continuous(expand = expansion(mult = 0.3), name = paste0("PC1", " (", PoV[1], " % of total variance)")) +
  scale_y_continuous(expand = expansion(mult = 0.3), name = paste0("PC2", " (", PoV[2], " % of total variance)")) +
  theme(axis.text.x = element_text(vjust = 0.5),
        axis.ticks = element_blank(),
        legend.position = "none",
        legend.justification = c("right", "top"),
        legend.box.just = "right",
        legend.margin = margin(6, 6, 6, 6),
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color=NA),
        panel.grid.major = element_line(color = "darkgrey"),
        panel.grid.minor = element_blank(),
        legend.background = element_rect(fill='transparent'),
        legend.box.background = element_rect(fill='transparent'),
        plot.margin = unit(c(1,1,1,1), units = "cm"),
        text = element_text(size = 25))
# Save plot to results 
png(file = "results/10_plot_de/pca.png",
    width = 10, 
    height = 8,
    units = "in",
    res = 800)

pca.plot

dev.off()


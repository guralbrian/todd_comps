# Load libs
libs <- c("tidyverse", "RColorBrewer", "reshape2") # list libraries here
lapply(libs, require, character.only = T)
rm(libs)

#### Loading and formatting of data ####
# Load compositions
decon.whole <- read.csv("data/processed/compositions/whole_samples.csv")

#### Plot whole samples ####
#brewer.pal(n=8,"Paired")

my_palette <- c( "#1F78B4", "#FF7F00")

ordered_df <- decon.whole %>%
  filter(CellType == "Cardiomyocytes") %>%  # Filter rows where CellType is "Cardiomyocytes"
  arrange(Prop) %>%  # Sort these rows by Prop in descending order
  pull(id)  # Extract the new.id values in this order
# Use the ordered new.id values to reorder the original data frame
decon.whole <- decon.whole %>%
  mutate(id = factor(id, levels = ordered_df)) %>%
  arrange(id)
decon.whole$CellType_wrap = str_wrap(decon.whole$CellType, width = 12)


# Use the ordered new.id values to reorder the original data frame
decon.whole <- decon.whole %>%
  mutate(id = factor(id, levels = ordered_df)) %>%
  arrange(id) |> 
  mutate(treatment = factor(treatment, levels = c("PBS", "PE")))

cell.type.order <-  decon.whole %>%
  group_by(CellType_wrap) |> 
  mutate(mean = mean(Prop)) |> 
  arrange(desc(mean)) |>
  pull(CellType_wrap) |> 
  unique()

# Generate boxplot
comp_celltype <- decon.whole   %>%
  ggplot(aes(x = factor(CellType_wrap, levels = as.character(cell.type.order)), y = Prop, fill = treatment)) +
  #geom_boxplot(position = position_dodge(0.9), width = 0.9, color = "black") +
  geom_bar(stat = "summary",position = position_dodge(0.9),  fun = mean,width = 0.9,  color = "black", alpha = 0.8) +
  geom_jitter(inherit.aes = T, 
              position = position_jitterdodge(jitter.width = 0.005, dodge.width = 0.9),
              size = 4, alpha = 0.5) +
  theme(axis.text.x = element_text(color = "black", angle =25, vjust = 0.5),
        axis.text.y = element_text(color = "black"),
        axis.ticks = element_blank(),
        legend.position = c(0.9, 1),
        legend.justification = c("right", "top"),
        legend.box.just = "right",
        legend.margin = margin(6, 6, 6, 6),
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color=NA),
        panel.grid.minor = element_blank(),
        legend.background = element_rect(fill='transparent'),
        legend.box.background = element_rect(fill='transparent'),
        axis.title.x = element_blank(),
        plot.margin = unit(c(1,1,1,1), units = "cm"),
        text = element_text(size = 40)) +
  labs(y = "Proportion", 
       fill = "Treatment") +
  scale_fill_manual(values = my_palette)

# Save plot to results 

#Save outputs
if(!dir.exists("results/03_plot_comps")){
  dir.create("results/03_plot_comps")
}

png(file = "results/03_plot_comps/sample_comps.png",
    width = 16, 
    height = 8,
    units = "in",
    res = 400)

comp_celltype

dev.off()


# refactorize the cell type levels
decon.whole <- decon.whole |> 
  mutate(cell.type = factor(cell.type))


#### Plot fractions ####

p.frac <- decon.whole |>
  ggplot(aes(x=id, y=Prop, fill=CellType))  +
  geom_bar(stat='identity',
           position = "fill",
           width = 1,
           color = "black")+
  scale_fill_brewer(name = "Cell Type",
                    palette = "Set2")  +
  facet_wrap(~treatment,
             scales = "free_x")+
  ylab("Estimated Proportion") +
  theme_minimal() +
  theme(strip.text = element_text(size = 30),
        title = element_text(size = 20),
        axis.text.y  = element_text(size = 20),
        legend.text = element_text(size = 20),
        panel.grid = element_blank(),
        legend.title = element_blank(),
        legend.position = 'right',
        legend.justification = 'center',
        legend.margin = margin(t = -5, r = 0, b = 0, l = 0, unit = "pt"),
        axis.text.x = element_text(size = 15, angle = 30),
        plot.margin = margin(0, 5, 1, 5),
        plot.caption = element_text(size = 20, hjust = -0.3, face = "bold"),
        plot.tag = element_text(size = 20, hjust = -0.3),
        axis.title.y = element_text(size = 20),
        axis.title.x =  element_text(size = 20)) +
  xlab("Bulk RNAseq Replicates") 

png(file = "results/03_plot_comps/stacked_barplot.png",
    width = 14, 
    height = 8,
    units = "in",
    res = 600)

p.frac

dev.off()

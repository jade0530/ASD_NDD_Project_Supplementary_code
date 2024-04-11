# Library
library(ggplot2)
library(hrbrthemes)
library(readxl)
library(tidyr)
library(dplyr)
library(readr)
library(RColorBrewer)

# Here to load data
go_terms <- read_excel('/media/bml/USER_DATA/ASD_CNV_Project/HEATMAP/CNVs_significant_on_cnv_New_V3_heatmap_filled.xlsx', sheet = "Heatmap_GO_TERM")
pathway <- read_excel('/media/bml/USER_DATA/ASD_CNV_Project/HEATMAP/CNVs_significant_on_cnv_New_V3_heatmap_filled.xlsx', sheet = "Heatmap_PATHWAY")

# Prepare the dataframe
# If FDR=0, set to 1E-10; If FDR=nan, set to 1


go_terms_prep <- go_terms 
go_terms_prep$GO_term_combine = paste(go_terms_prep$GO_ID, go_terms_prep$GO_term)
go_terms_prep[is.na(go_terms_prep)] = 1
go_terms_prep[go_terms_prep<0.0009766] = 0.0006905 # Set up -log2 max value to 10
go_terms_prep <- go_terms_prep[c(3:7)]
go_terms_prep <- go_terms_prep %>%
  mutate(GO_term_combine=factor(GO_term_combine, levels=rev(sort(unique(GO_term_combine))))) %>%
   pivot_longer(
      cols = !GO_term_combine, 
      names_to = "Group", 
      values_to = "FDR"
  ) %>%
  mutate(GO_term_combine=factor(as.character(GO_term_combine), levels=rev(levels(GO_term_combine))))

# Traverse FDR to -log2
go_terms_prep$FDR_log_2 = round(-log2(go_terms_prep$FDR), digits = 0)
# Round up the FDR values to 3 decimal places, change back nan values for showing 
go_terms_prep$FDR_sci = formatC(go_terms_prep$FDR, format = "e", digits = 1)
go_terms_prep$FDR_sci = na_if(go_terms_prep$FDR_sci, "1.0e+00")
# Create Heatmap with colours:
heatmap_plot_go <- ggplot(go_terms_prep, aes(x = Group, y = GO_term_combine, fill = FDR_log_2)) + 
  geom_tile(color = "black",
            lwd = 0.2,
            linetype = 1) + 
  geom_text(aes(label = FDR_sci), color = "black", size = 1) + 
  scale_fill_gradient(high = "#0073C2FF", low = "white") +
  labs(title = "Gene Ontology Terms Significance",
       x = "Groups",
       y = "Gene Ontology Terms",
       fill = "-log2 FDR",
       ) +
  theme(axis.text.x = element_text(size = 6), 
        axis.title = element_text(size = 8), 
        axis.text.y = element_text(size = 5, face = "italic"),
        plot.title = element_text(size = 10, hjust = 0.5, face = "bold"),
        plot.title.position = "plot", 
        legend.position="right", legend.direction="vertical",
        legend.text.align = 1, 
        legend.key.height=grid::unit(0.8, "cm"),
        legend.key.width=grid::unit(0.2, "cm"),
        legend.text = element_text(size = 6),
        title = element_text(size = 7))

# Create Heatmap with colours:
heatmap_plot_go_v <- ggplot(go_terms_prep, aes(x = GO_term_combine, y = Group, fill = FDR_log_2)) + 
  geom_tile(color = "black",
            lwd = 0.2,
            linetype = 1) + 
  geom_text(aes(label = FDR_sci), color = "black", size = 1) + 
  scale_fill_gradient(high = "#0073C2FF", low = "white") +
  labs(title = "Gene Ontology Terms Significance",
       y = "Groups",
       x = "Gene Ontology Terms",
       fill = "-log2 FDR",
  ) +
  theme(axis.text.y = element_text(size = 6), 
        axis.title = element_text(size = 8), 
        axis.text.x = element_text(size = 5, face = "italic", angle = 45, vjust = 1, hjust = 1),
        plot.title = element_text(size = 10, hjust = 0.5, face = "bold"),
        plot.title.position = "plot", 
        legend.position="right", legend.direction="vertical",
        legend.text.align = 1, 
        legend.key.height=grid::unit(0.8, "cm"),
        legend.key.width=grid::unit(0.2, "cm"),
        legend.text = element_text(size = 6),
        title = element_text(size = 7))

# Same procedure for pathway heatmap 
pathway_prep <- pathway
# stupid way
pathway_prep[is.na(pathway_prep)] = 1
pathway_prep$ASD[pathway_prep$ASD<0.0000009537] = 0.0000006743 # Set up -log2 max values to 25, have to since must be greater th
pathway_prep$ASD_ID_DD[pathway_prep$ASD_ID_DD<0.0000009537] = 0.0000006743 # Set up -log2 max values to 25, have to since must be greater th
pathway_prep$ID_DD[pathway_prep$ID_DD<0.0000009537] = 0.0000006743 # Set up -log2 max values to 25, have to since must be greater th
pathway_prep$SCZ[pathway_prep$SCZ<0.0000009537] = 0.0000006743 # Set up -log2 max values to 25, have to since must be greater th
pathway_prep <- pathway_prep[c(2:6)]
# pivot the dataframe 
pathway_prep <- pathway_prep %>%
  mutate(description=factor(description, levels=rev(sort(unique(description))))) %>%
  pivot_longer(
    cols = !description, 
    names_to = "Group", 
    values_to = "FDR"
  ) %>%
  mutate(description=factor(as.character(description), levels=rev(levels(description))))
# Change to log2
pathway_prep$FDR_log2 = round(-log2(pathway_prep$FDR), digits = 0) 

# Round up the FDR values to 2 decimal places, change back nan values for showing 
pathway_prep$FDR_sci = formatC(pathway_prep$FDR, format = "e", digits = 1)
pathway_prep$FDR_sci = na_if(pathway_prep$FDR_sci, "1.0e+00")

# Create Heatmap with colours:
heatmap_plot_path <- ggplot(pathway_prep, aes(x = Group, y = description, fill = FDR_log2)) + 
  geom_tile(color = "black",
            lwd = 0.2,
            linetype = 1) + 
  geom_text(aes(label = FDR_sci), color = "black", size = 1) +
  scale_fill_gradient(high = "#CD534CFF", low = "white") +
  labs(title = "Pathway Enrichment Significance",
       x = "Groups",
       y = "Pathways",
       fill = "-log2 FDR") + 
  geom_text(aes(label = FDR), color = "black", size = 1) + 
  theme(axis.text.x = element_text(size = 6), 
        axis.title = element_text(size = 8), 
        axis.text.y = element_text(size = 5, face = "italic"),
        plot.title = element_text(size = 10, hjust = 0.5, face = "bold"),
        plot.title.position = "plot", 
        legend.position="right", legend.direction="vertical",
        legend.text.align = 1, 
        legend.key.height=grid::unit(0.8, "cm"),
        legend.key.width=grid::unit(0.2, "cm"),
        legend.text = element_text(size = 6),
        title = element_text(size = 7))

heatmap_plot_path_v <- ggplot(pathway_prep, aes(x = description, y = Group, fill = FDR_log2)) + 
  geom_tile(color = "black",
            lwd = 0.2,
            linetype = 1) + 
  geom_text(aes(label = FDR_sci), color = "black", size = 1) + 
  scale_fill_gradient(high = "#CD534CFF", low = "white") +
  labs(title = "Pathway Enrichment Significance",
       x = "Pathways",
       y = "Groups",
       fill = "-log2 FDR") + 
  theme(axis.text.y = element_text(size = 6), 
        axis.title = element_text(size = 8), 
        axis.text.x = element_text(size = 5, face = "italic", angle = 45, vjust = 1, hjust = 1),
        plot.title = element_text(size = 10, hjust = 0.5, face = "bold"),
        plot.title.position = "plot", 
        legend.position="right", legend.direction="vertical",
        legend.text.align = 1, 
        legend.key.height=grid::unit(0.8, "cm"),
        legend.key.width=grid::unit(0.2, "cm"),
        legend.text = element_text(size = 6),
        title = element_text(size = 7))

# Arrange the venn diagrams to heatmaps 
heatmap_plot_path_arrange <- ggarrange(heatmap_plot_path, labels = c("A"), ggarrange(coding_venn_title, nc_venn_title, labels = c("B","C"), nrow = 2))
heatmap_plot_go_arrange <- ggarrange(heatmap_plot_go, labels = c("A"), ggarrange(coding_venn_title, nc_venn_title, labels = c("B","C"), nrow = 2))


heatmap_plot_arrange <- ggarrange(heatmap_plot_go_v, heatmap_plot_path_v, labels = c("C", "D"), nrow = 2)
venn_arrange <- ggarrange(coding_venn_title, nc_venn_title, labels = c("A","B"), nrow = 2)
heatmap_arrange_4 <- ggarrange(venn_arrange, heatmap_plot_arrange, ncol = 2, widths = c(1,2))
heatmap_arrange_4 <- heatmap_arrange_4 + theme(
  plot.margin = margin(1,1,1,1,unit = "cm")
)
heatmap_plot_arrange_1 <- ggarrange(heatmap_plot_go_v, heatmap_plot_path_v, nrow = 2)
heatmap_plot_arrange_1 <- heatmap_plot_arrange_1 + theme(
  plot.margin = margin(1,1,1,1,unit = "cm")
)
venn_arrange_1 <- ggarrange(coding_venn_title, nc_venn_title, nrow = 2)
venn_arrange_1 <- venn_arrange_1 + theme(
  plot.margin = margin(1,1,1,1,unit = "cm")
)

# Save the diagram using ggsave 
ggsave(heatmap_plot_path, filename="/media/bml/USER_DATA/ASD_CNV_Project/heatmap_pathway.pdf", height=5.5, width=6.5, units="in", dpi=200)
ggsave(heatmap_plot_go, filename="/media/bml/USER_DATA/ASD_CNV_Project/heatmap_go.pdf", height=5.5, width=6.5, units="in", dpi=200)

ggsave(heatmap_plot_path_v, filename="/media/bml/USER_DATA/ASD_CNV_Project/heatmap_pathway_v.pdf", height=5.5, width=10, units="in", dpi=200)
ggsave(heatmap_plot_go_v, filename="/media/bml/USER_DATA/ASD_CNV_Project/heatmap_go_v.pdf", height=5.5, width=12.5, units="in", dpi=200)

ggsave(heatmap_plot_path_arrange, filename="/media/bml/USER_DATA/ASD_CNV_Project/heatmap_pathway_venn.pdf", height=5.5, width=10, units="in", dpi=200)
ggsave(heatmap_plot_go_arrange, filename="/media/bml/USER_DATA/ASD_CNV_Project/heatmap_go_venn.pdf", height=5.5, width=10, units="in", dpi=200)

ggsave(heatmap_plot_arrange_1, filename="/media/bml/USER_DATA/ASD_CNV_Project/heatmap_arrange.pdf", height=7, width=10, units="in", dpi=200)
ggsave(venn_arrange_1, filename="/media/bml/USER_DATA/ASD_CNV_Project/venn_arrange.pdf", height=7, width=5.5, units="in", dpi=200)
ggsave(heatmap_arrange_4, filename="/media/bml/USER_DATA/ASD_CNV_Project/heatmap_4_venn.pdf", height=7, width=15.3, units="in", dpi=200)


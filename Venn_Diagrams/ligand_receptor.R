# Libraries required 
library(ggvenn)
library(readxl)
library(ggplot2)
library(ggpubr)

# Load subsets of the significant genes into R program
# Each group created by filters through excel, put into a new sheet under "CNV_subsets.xlsx"
ligand_receptor <- read_excel('/Users/jadeliang/Desktop/Mona_project_CNV/cell-cell_communication/ligand_receptor_pairs.xlsx')
#ligand_receptor[is.na(ligand_receptor)] <- "N"

ligand_receptor_process <- ligand_receptor[!(is.na(ligand_receptor$ASD) & is.na(ligand_receptor$ASD_ID) 
                                          & is.na(ligand_receptor$ID_DD) & is.na(ligand_receptor$SCZ)),]

ASD <- ligand_receptor_process[ligand_receptor_process$ASD=="Y", ]
ASD <- ASD$LRsig.pathway_name[!is.na(ASD$ASD)]

ASD_ID <- ligand_receptor_process[ligand_receptor_process$ASD_ID=="Y",]
ASD_ID_network <- ASD_ID[!is.na(ASD_ID$ASD_ID),]
ASD_ID <- ASD_ID$LRsig.pathway_name[!is.na(ASD_ID$ASD_ID)]

ID_DD <- ligand_receptor_process[ligand_receptor_process$ID_DD=="Y",]

ID_DD_network <- ID_DD[!is.na(ID_DD$ID_DD),]
ID_DD <- ID_DD$LRsig.pathway_name[!is.na(ID_DD$ID_DD)]

SCZ <- ligand_receptor_process[ligand_receptor_process$SCZ=="Y", ]
SCZ_network <- SCZ[!is.na(SCZ$SCZ),]
SCZ <- SCZ$LRsig.pathway_name[!is.na(SCZ$SCZ)]


# Create lists that includes 4 groups (ASD, ASD_ID_DD, ID_DD, SCZ) according to the subsets
receptors_name <- list(
  ASD = c(ASD),
  ASD_ID = c(ASD_ID),
  ID_DD = c(ID_DD),
  SCZ = c(SCZ)
)

receptors_name_three <- list(
  ASD_ID = c(ASD_ID),
  ID_DD = c(ID_DD),
  SCZ = c(SCZ)
)
library(VennDiagram)
# Visualising coding genes duplications using ggvenn
ligand_receptor_venn <- ggvenn(
  receptors_name, 
  fill_color = c("#0073C2FF","#EFC000FF", "#868686FF", "#CD534CFF"),
  stroke_size = 0.5, show_percentage = FALSE,
)

ggsave("/Users/jadeliang/Desktop/Mona_project_CNV/cell-cell_communication/ligand_receptor_venn_pathway.svg", ligand_receptor_venn, units = "in", height = 10, width = 10)

library(ggVennDiagram)


ggVennDiagram(
  x = receptors_name,
  category.names = c("ASD" ,"ASD_ID" , "ID_DD" , "SCZ"),
  label = "count",
  set_color = c("#0073C2FF","#EFC000FF", "#868686FF", "#CD534CFF")
) + 
  scale_fill_gradient(low = "white", high = "orange")

ggVennDiagram(
  x = receptors_name_three,
  category.names = c("ASD" ,"ASD_ID" , "ID_DD" , "SCZ"),
  label = "count", 
  set_color = c("#EFC000FF", "#868686FF", "#CD534CFF")
) 


# Add graph title 


ligand_receptor_heatmap <- ggplot(ligand_receptor, aes(x = L_R, y = Group, fill = pValue_log2)) + 
  geom_tile(color = "black",
            lwd = 0.2,
            linetype = 1) + 
  geom_text(aes(label = pValue_sci), color = "black", size = 1) + 
  scale_fill_gradient(high = "#CD534CFF", low = "white") +
  labs(title = "Pathway Enrichment Significance",
       x = "Pathways",
       y = "Groups",
       fill = "-log2 pValue") + 
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



library(reshape2)

ligand_receptor_melt <- ligand_receptor[!(is.na(ligand_receptor$ASD) & is.na(ligand_receptor$ASD_ID) 
                                     & is.na(ligand_receptor$ID_DD) & is.na(ligand_receptor$SCZ)),]


ligand_receptor_melt <- ligand_receptor_melt[order(ligand_receptor_melt$LRsig.pathway_name),]
ligand_receptor_melt_1 <- ligand_receptor_melt[1:73,]
ligand_receptor_melt_2 <- ligand_receptor_melt[74:127,]
ligand_receptor_melt_3 <- ligand_receptor_melt[128:202,]

preprocess_table <- function (df) {
  ligand_receptor_melt <- melt(df, id = c('L_R', 'LRsig.pathway_name'), na.rm = FALSE)
  
  
  ligand_receptor_melt$value[ligand_receptor_melt$value=="Y"] <- 1
  ligand_receptor_melt$value[is.na(ligand_receptor_melt$value)] <- 0
  
  ligand_receptor_melt_heatmap_1 <- ggplot(ligand_receptor_melt, aes(variable, LRsig.pathway_name, alpha = as.numeric(value))) + 
    geom_tile(aes(fill = variable), color = "white") +
    scale_alpha_identity(guide = "none") +
    coord_equal(expand = 0) +
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1))
  return(ligand_receptor_melt_heatmap_1)
}

h1 <- preprocess_table(ligand_receptor_melt_1)
h2 <- preprocess_table(ligand_receptor_melt_2)
h3 <- preprocess_table(ligand_receptor_melt_3)


library(ggpubr)
heatmap_arranged <- ggarrange(h1,h2,h3,ncol = 3, nrow = 1, common.legend = TRUE, legend = "bottom")

ligand_receptor_melt_heatmap <- ggplot(ligand_receptor_melt, aes(variable, L_R, fill = LRsig.pathway_name, alpha = as.numeric(value))) + 
  geom_tile(colour = "gray50") +
  scale_alpha_identity(guide = "none") +
  coord_equal(expand = 0) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("/Users/jadeliang/Desktop/Mona_project_CNV/cell-cell_communication/ligand_receptor_melt_heatmap_class.svg", heatmap_arranged, units = "in", height = 10.5, width = 15.5, dpi = 300)

l = 0
ASD_ID_network$label = 0
for (i in 1:nrow(ASD_ID_network)) {
  if(ASD_ID_network$LRsig.pathway_name[i+1]==ASD_ID_network$LRsig.pathway_name[i]) {
    ASD_ID_network$label = l
  } else {
    l <- l + 1
    ASD_ID_network$label = l
  }
  
}
library(igraph)
relations_ASD_ID <- data.frame(
  from = c(ASD_ID_network$L),
  to = c(ASD_ID_network$R)
)

relations_ID_DD <- data.frame(
  from = c(ID_DD_network$L),
  to = c(ID_DD_network$R)
)

relations_SCZ <- data.frame(
  from = c(SCZ_network$L),
  to = c(SCZ_network$R)
)

g_ASD_ID <- graph_from_data_frame(relations_ASD_ID, directed = FALSE)
g_ID_DD <- graph_from_data_frame(relations_ID_DD, directed = FALSE)
g_SCZ <- graph_from_data_frame(relations_SCZ, directed = FALSE)

plot(g, label.cex = 0.8)
ASD_ID_graphplot <- plot.igraph(g_ASD_ID, mark.col = "lightblue")

ID_DD_graphplot <- plot.igraph(g_ID_DD, mark.col = "lightblue")

SCZ_graphplot <- plot.igraph(g_SCZ, mark.col = "lightblue")

png("/Users/jadeliang/Desktop/Mona_project_CNV/cell-cell_communication/ASD_ID_graphplot.png", 1200, 1200)
plot.igraph(g_ASD_ID,mark.expand = 15, loop.size = 1.2)
dev.off()

png("/Users/jadeliang/Desktop/Mona_project_CNV/cell-cell_communication/ID_DD_graphplot.png", 1200, 1200)
plot.igraph(g_ID_DD)
dev.off()

png("/Users/jadeliang/Desktop/Mona_project_CNV/cell-cell_communication/SCZ_graphplot.png", 1200, 1200)
plot.igraph(g_SCZ)
dev.off()

write_graph(
  g_ASD_ID,
  "/Users/jadeliang/Desktop/Mona_project_CNV/cell-cell_communication/ASD_ID_graphplot.svg",
  format = c("edgelist")
)
#ggarrange(ASD_ID_graphplot, ID_DD_graphplot, SCZ_graphplot, nrow = 2, ncol = 2, labels = c("ASD_ID", "ID_DD", "SCZ"))

# Libraries required 
library(ggvenn)
library(readxl)
library(ggplot2)
library(ggpubr)

# Generate Venn diagram for unique genes 
# Load subsets of the significant genes into R program
# Each group created by filters through excel, put into a new sheet under "CNV_subsets.xlsx"
asd_coding <- read_excel('/media/bml/USER_DATA/ASD_CNV_Project/CNV_subsets.xlsx', sheet = "ASD_coding_unique")
asd_nc <- read_excel('/media/bml/USER_DATA/ASD_CNV_Project/CNV_subsets.xlsx', sheet = "ASD_nc_unique")

asd_id_dd_coding <- read_excel('/media/bml/USER_DATA/ASD_CNV_Project/CNV_subsets.xlsx', sheet = "ASD_ID_DD_coding_unique")
asd_id_dd_nc <- read_excel('/media/bml/USER_DATA/ASD_CNV_Project/CNV_subsets.xlsx', sheet = "ASD_ID_DD_nc_unique")

id_dd_coding <- read_excel('/media/bml/USER_DATA/ASD_CNV_Project/CNV_subsets.xlsx', sheet = "ID_DD_coding_unique")
id_dd_nc <- read_excel('/media/bml/USER_DATA/ASD_CNV_Project/CNV_subsets.xlsx', sheet = "ID_DD_nc_unique")

scz_coding <- read_excel('/media/bml/USER_DATA/ASD_CNV_Project/CNV_subsets.xlsx', sheet = "SCZ_coding_unique")
scz_nc <- read_excel('/media/bml/USER_DATA/ASD_CNV_Project/CNV_subsets.xlsx', sheet = "SCZ_nc_unique")

# Create lists that includes 4 groups (ASD, ASD_ID_DD, ID_DD, SCZ) according to the subsets
coding <- list(
  ASD = asd_coding$Gene_ID, 
  ASD_ID_DD = asd_id_dd_coding$Gene_ID,
  ID_DD = id_dd_coding$Gene_ID,
  SCZ = scz_coding$Gene_ID
)

nc <- list(
  ASD = asd_nc$Gene_ID, 
  ASD_ID_DD = asd_id_dd_nc$Gene_ID,
  ID_DD = id_dd_nc$Gene_ID,
  SCZ = scz_nc$Gene_ID
)

# Visualising coding genes deletions using ggvenn
coding_venn <- ggvenn(
  coding, 
  fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
  stroke_size = 0.5, set_name_size = 4, show_percentage = FALSE
)

# Visualising noncoding genes duplications using ggvenn
nc_venn <- ggvenn(
  nc, 
  fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
  stroke_size = 0.5, set_name_size = 4, show_percentage = FALSE
)

# Arrange the 2 graphs into one 1x2 page 
ggarrange(coding_venn,nc_venn,
          labels = c("A", "B"),
          ncol = 2, nrow = 1,
          vjust = 1)


# Add graph title 
nc_venn_title <- nc_venn + labs(title = "Unique Non-coding Genes Venn Diagram for the 4 Target Groups") + 
  theme(plot.title = element_text(size = 10, hjust = 0.5, face = "bold"))


# Add graph title 
coding_venn_title <- coding_venn + labs(title = "Unique Coding Genes Venn Diagram for the 4 Target Groups") + 
  theme(plot.title = element_text(size = 10, hjust = 0.5, face = "bold"),
        text = element_text(size = 5))


# save graphs
ggsave(coding_venn_title, filename="/media/bml/USER_DATA/ASD_CNV_Project/coding_venn_title.pdf", height=3.4, width=4.7, units="in", dpi=200)
ggsave(nc_venn_title, filename="/media/bml/USER_DATA/ASD_CNV_Project/nc_venn_title.pdf", height=3.4, width=4.7, units="in", dpi=200)





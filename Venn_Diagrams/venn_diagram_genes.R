# Libraries required 
library(ggvenn)
library(readxl)
library(ggplot2)
library(ggpubr)

# Load subsets of the significant genes into R program
# Each group created by filters through excel, put into a new sheet under "CNV_subsets.xlsx"
asd_coding_dup <- read_excel('/media/bml/USER_DATA/ASD_CNV_Project/CNV_subsets.xlsx', sheet = "ASD_coding_dup")
asd_coding_del <- read_excel('/media/bml/USER_DATA/ASD_CNV_Project/CNV_subsets.xlsx', sheet = "ASD_coding_del")
asd_nc_dup <- read_excel('/media/bml/USER_DATA/ASD_CNV_Project/CNV_subsets.xlsx', sheet = "ASD_noncoding_dup")
asd_nc_del <- read_excel('/media/bml/USER_DATA/ASD_CNV_Project/CNV_subsets.xlsx', sheet = "ASD_noncoding_del")

asd_id_dd_coding_dup <- read_excel('/media/bml/USER_DATA/ASD_CNV_Project/CNV_subsets.xlsx', sheet = "ASD_ID_DD_coding_dup")
asd_id_dd_coding_del <- read_excel('/media/bml/USER_DATA/ASD_CNV_Project/CNV_subsets.xlsx', sheet = "ASD_ID_DD_coding_del")
asd_id_dd_nc_dup <- read_excel('/media/bml/USER_DATA/ASD_CNV_Project/CNV_subsets.xlsx', sheet = "ASD_ID_DD_noncoding_dup")
asd_id_dd_nc_del <- read_excel('/media/bml/USER_DATA/ASD_CNV_Project/CNV_subsets.xlsx', sheet = "ASD_ID_DD_noncoding_del")

dd_coding_dup <- read_excel('/media/bml/USER_DATA/ASD_CNV_Project/CNV_subsets.xlsx', sheet = "DD_coding_dup")
dd_coding_del <- read_excel('/media/bml/USER_DATA/ASD_CNV_Project/CNV_subsets.xlsx', sheet = "DD_coding_del")
dd_nc_dup <- read_excel('/media/bml/USER_DATA/ASD_CNV_Project/CNV_subsets.xlsx', sheet = "DD_noncoding_dup")
dd_nc_del <- read_excel('/media/bml/USER_DATA/ASD_CNV_Project/CNV_subsets.xlsx', sheet = "DD_noncoding_del")

scz_coding_dup <- read_excel('/media/bml/USER_DATA/ASD_CNV_Project/CNV_subsets.xlsx', sheet = "SCZ_coding_dup")
scz_coding_del <- read_excel('/media/bml/USER_DATA/ASD_CNV_Project/CNV_subsets.xlsx', sheet = "SCZ_coding_del")
scz_nc_dup <- read_excel('/media/bml/USER_DATA/ASD_CNV_Project/CNV_subsets.xlsx', sheet = "SCZ_noncoding_dup")
scz_nc_del <- read_excel('/media/bml/USER_DATA/ASD_CNV_Project/CNV_subsets.xlsx', sheet = "SCZ_noncoding_del")

# Create lists that includes 4 groups (ASD, ASD_ID_DD, ID_DD, SCZ) according to the subsets
coding_dup <- list(
  ASD = asd_coding_dup$Gene_ID, 
  ASD_ID_DD = asd_id_dd_coding_dup$Gene_ID,
  ID_DD = dd_coding_dup$Gene_ID,
  SCZ = scz_coding_dup$Gene_ID
)

nc_dup <- list(
  ASD = asd_nc_dup$Gene_ID, 
  ASD_ID_DD = asd_id_dd_nc_dup$Gene_ID,
  ID_DD = dd_nc_dup$Gene_ID,
  SCZ = scz_nc_dup$Gene_ID
)

coding_del <- list(
  ASD = asd_coding_del$Gene_ID, 
  ASD_ID_DD = asd_id_dd_coding_del$Gene_ID,
  ID_DD = dd_coding_del$Gene_ID,
  SCZ = scz_coding_del$Gene_ID
)

nc_del <- list(
  ASD = asd_nc_del$Gene_ID, 
  ASD_ID_DD = asd_id_dd_nc_del$Gene_ID,
  ID_DD = dd_nc_del$Gene_ID,
  SCZ = scz_nc_del$Gene_ID
)

# Visualising coding genes duplications using ggvenn
coding_dup_venn <- ggvenn(
  coding_dup, 
  fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
  stroke_size = 0.5, show_percentage = FALSE
)

# Add graph title 
coding_dup_venn <- coding_dup_venn + labs(title = "Coding Genes Venn Diagram (Duplications)")

# Visualising coding genes deletions using ggvenn
coding_del_venn <- ggvenn(
  coding_del, 
  fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
  stroke_size = 0.5, show_percentage = FALSE
)
# Add graph title 
coding_del_venn <- coding_del_venn + labs(title = "Coding Genes Venn Diagram (Deletions)")

# Visualising noncoding genes duplications using ggvenn
nc_dup_venn <- ggvenn(
  nc_dup, 
  fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
  stroke_size = 0.5, show_percentage = FALSE
)
# Add graph title 
nc_dup_venn <- nc_dup_venn + labs(title = "Non-coding Genes Venn Diagram (Duplications)")

# Visualising noncoding genes deletions using ggvenn
nc_del_venn <- ggvenn(
  nc_del, 
  fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
  stroke_size = 0.5, show_percentage = FALSE
)
# Add graph title 
nc_del_venn <- nc_del_venn + labs(title = "Non-coding Genes Venn Diagram (Deletions)")

# show graphs
coding_dup_venn
coding_del_venn
nc_dup_venn
nc_del_venn

# Arrange the 4 graphs into one 2x2 page 
ggarrange(coding_dup_venn, coding_del_venn, nc_dup_venn, nc_del_venn, 
          labels = c("A", "B", "C", "D"),
          ncol = 2, nrow = 2)

Reduce(intersect,list(dd_coding_dup$Gene_ID, asd_coding_dup$Gene_ID,asd_id_dd_coding_dup$Gene_ID))
Reduce(intersect,list(dd_coding_dup$Gene_ID, scz_coding_dup$Gene_ID))

# Load necessary libraries
library(readr)

# Load the data

# change it 
setwd("/Users/jadeliang/Desktop/Mona_project_CNV/mona_code/")


filter_gene_list <- function(data) {
  #categories <- c(colnames(data[,4:length(colnames(data))]))
  #categories <- gsub(" ", "_", tolower(categories))
  #categories <- gsub("-", "_positive", categories)
  

  #colnames(data[,4:length(colnames(data))]) <- gsub(" ", "_", tolower(colnames(data)))
  #colnames(data[4,]) <- gsub("-", "_positive", colnames(data))
  
  #print(colnames(data[,4:length(colnames(data))]))
  # Initialize a dataframe to store filtered genes
  filtered_genes_df <- data.frame(Category = character(), Gene_name = character(), stringsAsFactors = FALSE)
  
  # Filter genes by category and threshold
  for (category in colnames(data[,4:length(colnames(data))])) {
    if (category %in% colnames(data)) {
      filtered_genes <- unique(data$Gene_name[data[[category]]<0.5])
      if (length(filtered_genes)>0) {
        temp_df <- data.frame(Category = category, Gene_name = filtered_genes, stringsAsFactors = FALSE)
        filtered_genes_df <- rbind(filtered_genes_df, temp_df)
      } else {
        print("This category does not have any genes that < 0.5")
        print(category)
      }
    }
  }
  return(filtered_genes_df)
}




## Here run the code with the function
scz_data <- read_csv("scz_gene_names.csv")
filtered_genes_scz_df <- filter_gene_list(scz_data)
write.csv(filtered_genes_scz_df, "filtered_genes_SCZ.csv", row.names = FALSE)

# Do the same for the other three categories 
asd_data <- read_csv("asd_gene_names.csv")
filtered_genes_asd_df <- filter_gene_list(asd_data)
write.csv(filtered_genes_asd_df, "filtered_genes_ASD.csv", row.names = FALSE)

asd_id_dd_data <- read_csv("asd_id_dd_gene_names.csv")
filtered_genes_asd_id_dd_df <- filter_gene_list(asd_id_dd_data)
write.csv(filtered_genes_asd_id_dd_df, "filtered_genes_ASD_ID_DD.csv", row.names = FALSE)

id_dd_data <- read_csv("id_dd_gene_names.csv")
filtered_genes_id_dd_df <- filter_gene_list(id_dd_data)
write.csv(filtered_genes_id_dd_df, "filtered_genes_ID_DD.csv", row.names = FALSE)
#combine the all files
library(dplyr)


genes_ASD <- read.csv("filtered_genes_ASD.csv")
genes_ASD_ID_DD <- read.csv("filtered_genes_ASD_ID_DD.csv")
genes_ID_DD <- read.csv("filtered_genes_ID_DD.csv")
genes_SCZ <- read.csv("filtered_genes_SCZ.csv")


genes_ASD$Source <- 'filtered_genes_ASD'
genes_ASD_ID_DD$Source <- 'filtered_genes_ASD_ID_DD'
genes_ID_DD$Source <- 'filtered_genes_ID_DD'
genes_SCZ$Source <- 'filtered_genes_SCZ'

combined_genes <- rbind(genes_ASD, genes_ASD_ID_DD, genes_ID_DD, genes_SCZ)

# Remove duplicated gene names
final_genes <- combined_genes %>% distinct(Gene_name, .keep_all = TRUE)

# Optionally, save the final dataset to a new CSV file
write.csv(final_genes, "combined_genes_final.csv", row.names = FALSE)


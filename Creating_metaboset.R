if (!requireNamespace("devtools", quietly = TRUE)) {
install.packages("devtools") }

devtools::install_github("antonvsdata/notame")
# Klåvus et al. 2020
# Link to the original publication in the description!

library(notame)

hilic_neg_sample
View(pData(hilic_neg_sample))
View(fData(hilic_neg_sample))

# # Data from example samples
# write_to_excel(hilic_neg_sample,file = "Hilic_neg.xlsx")
# write_to_excel(hilic_pos_sample,file = "Hilic_pos.xlsx")
# write_to_excel(rp_neg_sample,file = "RP_pos.xlsx")
# write_to_excel(rp_pos_sample,file = "RP_neg.xlsx")




# Read in the data for HILIC_neg
HILIC_neg <- read_from_excel("Hilic_neg.xlsx", name = "Hilic_neg")

View(HILIC_neg$pheno_data)
View(HILIC_neg$feature_data)

#Construct metabosets
Hilic_neg <- construct_metabosets(exprs = HILIC_neg$exprs, pheno_data = HILIC_neg$pheno_data,
                              feature_data = HILIC_neg$feature_data)

Hilic_neg <- Hilic_neg$Hilic_neg


###############
###############
# Read in all analogical methods with looping

modes <- c("Hilic_neg", "Hilic_pos", "RP_pos", "RP_neg")
data_list <- list()

for (mode in modes) {
  filename <- paste0(mode, ".xlsx")
  data <- read_from_excel(file = filename, split_by = "Split")
  
  data_list[[mode]] <- construct_metabosets(exprs = data$exprs, pheno_data = data$pheno_data,
                                            feature_data = data$feature_data)[[1]]
}

merged <- merge_metabosets(data_list)
View(pData(merged))
View(fData(merged))

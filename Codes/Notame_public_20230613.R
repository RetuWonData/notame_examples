# 20230613
# Preprocessing workshop in metabolomics society conference 2023
# RHaikonen, VKoistinen, TMeuronen and KHanhineva

#1. Notame installation

if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
}
devtools::install_github("antonvsdata/notame", dependencies = c("Depends", "Imports", "Suggests"))

# # EXTRA
# #If some packages do not install, try this minimum installation 
# #and then install dependencies separately while skipping errors:
# 
# if (!requireNamespace("devtools", quietly = TRUE)) {
#   install.packages("devtools")
# }
# devtools::install_github("antonvsdata/notame")
# 
# install_dependencies(
#   preprocessing = TRUE,
#   extra = TRUE,
#   batch_corr = TRUE,
#   misc = TRUE,
# )

# Check instructions in browser
browseVignettes("notame")

# # EXTRA
# # Version under development
# devtools::install_github("antonvsdata/notame", ref = "dev")
# 
# # EXTRA
# Install BatchCorrMetabolomics or batchCorr (for batch correction)
# 
# if (!requireNamespace("devtools", quietly = TRUE)) {
#   install.packages("devtools")
# }
# devtools::install_github("rwehrens/BatchCorrMetabolomics")
# 
# if (!requireNamespace("devtools", quietly = TRUE)) {
#   install.packages("devtools")
# }
# devtools::install_gitlab("CarlBrunius/batchCorr")


#2. Load libraries, set up path and logging
library(notame)
library(doParallel)
library(dplyr)
library(openxlsx)

ppath <- "notame_examples/"

# log will take place in your selected path
today <- as.character(format(Sys.time(), "%Y%m%e"))
init_log(log_file = paste0(ppath, "log/log", today, ".txt"))


#3. Read data and construct the MetaboSet objects
# (note: assumption is that the data contains signals from two or more modes)
data <- read_from_excel(file = paste0(ppath, "Data/Metaboset_to_metsoc.xlsx"), split_by = "Split")


#Check how it looks
names(data)
sapply(data, class)
sapply(data, dim)

# Build the metaboset
objects <- construct_metabosets(
  exprs = data$exprs,
  pheno_data = data$pheno_data,
  feature_data = data$feature_data,
  group_col = "Tissue", time_col = "Diet")


# #3b. Alternative: Read in the separate modes and construct MetaboSet objects
# # (note: here we read in individual modefiles)
# 
# modes <- c("HILIC_neg", "HILIC_pos", "RP_neg", "RP_pos")
# objects <- list()
# for (mode in modes) {
#   # Read single mode
#   filename <- paste0(ppath, "Data/", mode, "_sample.xlsx")
#   data <- read_from_excel(file = filename,
#                           name = mode)
#   objects[[mode]] <- construct_metabosets(exprs = data$exprs,
#                                           pheno_data = data$pheno_data,
#                                           feature_data = data$feature_data,
#                                           group_col = "Tissue", time_col = "Diet")[[1]]
# }



#Check how it looks
names(objects)
sapply(objects, class)

# Also check the PCA
temp <- merge_metabosets(objects)
plot_pca(temp, color = "Injection_order")
plot_pca(temp, color = "Batch")

# # EXTRA
# Replace spaces with underscores if necessary
# 
# rownames(data$exprs) <- gsub(" ", "_", rownames(data$exprs))
# rownames(data$feature_data) <- gsub(" ", "_", rownames(data$feature_data))
# data$feature_data$Feature_ID <- gsub(" ", "_", data$feature_data$Feature_ID)
# data$feature_data$Split <- gsub(" ", "_", data$feature_data$Split)


#4. Preprocessing of all the modes (NB: visualizations disabled here for saving computational time,
# remove the hashtags to enable them. Also, create a folder called figures in the working folder if you
# create visualizations)

#Take several cores into use (default: 8)
cl <- makeCluster(8)
registerDoParallel(cl)

# Initialize empty list for processed objects
processed <- list()
for (i in seq_along(objects)) {
  name <- names(objects)[i]
  mode <- objects[[i]]
  # Set all zero abundances to NA
  mode <- mark_nas(mode, value = 0)
  mode <- flag_detection(mode, qc_limit = 0.7, group_limit = 0.8)
  
  corrected <- correct_drift(mode)
  
  corrected <- corrected %>% assess_quality() %>% flag_quality()
  processed[[i]] <- corrected
}

#Stop using several cores (releases them for other use)
stopCluster(cl)

# Merge the modes and check the PCA
merged <- merge_metabosets(processed)
plot_pca(merged, color = "Injection_order")
plot_pca(merged, color = "Batch")

#visualizations(merged, prefix = paste0(ppath, "figures/_FULL"))

# # EXTRA
# Batch correction (if data contains samples from several batches)
# View(pData(merged))
# plot_pca(merged, color = "Batch")
# plot_pca(merged, color = "Injection_order")
# plot_pca(merged, color = "Injection_order", label = "Sample_ID")
# 
# Option A: BatchCorrMetabolomics
# batch_corrected <- dobc(merged, batch = "Batch", ref = "QC", ref_label = "QC")
# 
# Option B: batchCorr (developed by Carl Brunius)
# batch_corrected_B <- normalize_batches(merged, batch = "Batch", group = "QC", ref_label = "QC")
# 
# #Check how it looks: calculate Bhattacharyya distances and save batch correction plots
# 
# pca_bhattacharyya_dist(batch_corrected, "Batch", all_features = FALSE, center = TRUE,
#                        scale = "uv", nPCs = 3)
# 
# save_batch_plots(orig = merged,
#                  corrected = batch_corrected, 
#                  file = "Batch_correction_plots.pdf", 
#                  width = 14, 
#                  height = 10, 
#                  batch = "Batch_ID",
#                  color = "Batch_ID",
#                  shape = "QC",
#                  color_scale = getOption("notame.color_scale_dis"),
#                  shape_scale = scale_shape_manual(values = c(15, 21))
#                  )


#5. Imputation (note: may not be necessary especially if gap filling by compulsion was used in MS-DIAL)

# Remove the QCs for imputation
merged_no_qc <- drop_qcs(merged)
#visualizations(merged_no_qc, prefix = paste0(ppath, "figures/FULL_NO_QC"))

#Set seed number for reproducibility
set.seed(38)
imputed <- impute_rf(merged_no_qc, all_features = FALSE)
imputed <- impute_rf(imputed, all_features = TRUE)

# Save the merged, processed and imputed data
#6a. Write to Excel (note: update the object and file name according to project)
write_to_excel(imputed, file = paste0(ppath, "preprosessed_excel_", today,".xlsx"))

#6b. Save data in RDS format
saveRDS(imputed, file = paste0(ppath, "preprosessed_example_", today,".RDS"))

# Finish up the log!
finish_log()


##
####
######
# Extra flagging

#7 Low expression flag
# continue with imputed metaboset

# Threshold for mean expression
threshold <- 5000

means <- apply(exprs(imputed), 1, finite_mean)
condition <- (means < threshold) & is.na(flag(imputed))
summary(condition)
flag(imputed)[condition] <- "Below 5K abundance"

# Threshold for smallest value
threshold <- 500
means <- apply(exprs(imputed), 1, min)
condition <- (means < threshold) & is.na(flag(imputed))
summary(condition)
flag(imputed)[condition] <- "Missing values"

# lets have a look of flagged ones
View(fData(imputed)[,c("Feature_ID","Metabolite_name","Flag")])

#8	No ms2
# Check the name of msms column
colnames(fData(imputed))

condition <- is.na(fData(imputed)$"MS_MS_spectrum") & is.na(flag(imputed))
summary(condition)
flag(imputed)[condition] <- "No MSMS"


#9 Create "goodness" calculator
# Count based on p value and/or fold change (Differently abundant ones get higher value!)

# Take specific groups
temp <- imputed[,imputed$Diet %in% c("Rye", "Wheat")]
pData(temp) <- droplevels(pData(temp)) 

# Start parallel backend
cl <- makeCluster(8)
registerDoParallel(cl)

f_change <- fold_change(temp, group = "Diet")
stat_test <- perform_t_test(temp, formula_char = "Feature ~ Diet")
colnames(stat_test)[ncol(stat_test)] <- "1_vs_2_t_test_q"
stat_test <- stat_test[,c(1,7,8)]
cohen <- cohens_d(temp, group = "Diet")
# Stop parallel backend
stopCluster(cl)

temp <- join_fData(temp, f_change)
temp <- join_fData(temp, stat_test)
temp <- join_fData(temp, cohen)


count_test <- function(object, sel_test, thresh){
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("Package \"dplyr\" needed for this function to work. Please install it.", 
         call. = FALSE)
  }
  if(is.null(fData(object)$good_flag)){
    fData(object)$good_flag <- 0
  }
  # already existing number
  original_count <- fData(object)$good_flag
  
  
  # count if fold change or cohen_D
  if (sel_test %in% c("Cohen_d")){
    thresh2 <- -thresh
    temp <- select(fData(object), contains(sel_test))
    countt <- apply(temp, 1, function(x) sum(x>thresh|x<thresh2, na.rm = T))
  }
  else if (sel_test == c("FC")){
    thresh2 <- 1/thresh
    temp <- select(fData(object), contains(sel_test))
    countt <- apply(temp, 1, function(x) sum(x>thresh|x<thresh2, na.rm = T))
  }
  
  # count if propability like t-test
  else {
    temp <- select(fData(object), contains(sel_test))
    countt <- apply(temp, 1, function(x) sum(x<thresh, na.rm = T))
  }
  # Calculate together
  comb_count <- original_count + countt
  
  fData(object)$good_flag <- comb_count
  
  object
}

temp <- count_test(temp, "t_test_p", 0.05)
temp <- count_test(temp, "t_test_q", 0.05)
temp <- count_test(temp, "FC", 2)
temp <- count_test(temp, "Cohen_d", 1)

View(fData(temp)[,c(1,6, 24,33:37)])
View(fData(temp)[is.na(fData(temp)$Flag),c("Feature_ID","Metabolite_name","Flag", "good_flag")])


# Repeat the saving the changes
#10. Write to Excel (note: update the object and file name according to project)
write_to_excel(temp, file = paste0(ppath, "Extra_flagged_example_", today,".xlsx"))
#10b. Save data in RDS format
saveRDS(temp, file = paste0(ppath, "Extra_flagged_example_", today,".RDS"))

# Plot some metabolites

notame_boxplot <- function(object, feat_name, name_col = NULL,
                           separate = group_col(object),
                           face = NULL, comp = NULL){
  metabolite <- as.data.frame(t(exprs(object[fData(object)$Feature_ID == feat_name])))
  metabolite <- cbind(metabolite, pData(object)[,separate])
  
  # Add the group to facet
  if(!is.null(face)){
    metabolite <- cbind(metabolite, pData(object)[,face])
  }
  
  
  colnames(metabolite) <- c("metabol", "sep")
  # give the name for facet group
  if(!is.null(face)){
    colnames(metabolite)[3] <- "face"
  }
  
  # add the title based on identification
  if (!is.null(name_col)){
    if (!is.na(fData(object)[feat_name, name_col])){
      otsikko <- fData(object)[feat_name, name_col]
    }
    else{
      otsikko <- feat_name
    }
  } else{
    otsikko <- feat_name
  }
  
  library(ggpubr)
  Tissue_plot1 <- ggpubr::ggboxplot(metabolite, x = "sep", y = "metabol",
                                    title = otsikko, #fill = "sep", palette = c("#009FB8", "#BB0A21", "#F6AE2D", "#9063CD", "#DC7FA1")
  ) +
    theme_bw() + 
    theme(axis.text.x = element_text(angle = 90),
          axis.title = element_blank(),
          plot.title = element_text(hjust = 0.5))
  
  if(!is.null(comp)){
    if(is.list(comp)){
      Tissue_plot1 <- Tissue_plot1 +
        geom_signif(comparisons = comp,
                    map_signif_level=sigFunc,
                    margin_top = seq(-0.1, 0.1*length(comp), length.out = length(comp))) +
        scale_y_continuous(expand = c(0.2,0))
    }
    else{
      comp <- t(combn(levels(metabolite$sep),2))
      comp <- split(comp, 1:nrow(comp))
      Tissue_plot1 <- Tissue_plot1 + 
        geom_signif(comparisons = comp,
                    map_signif_level=sigFunc,
                    margin_top = seq(-0.1, 0.1*length(comp), length.out = length(comp)),
                    vjust = 0.6)
    }
  }
  if(!is.null(face)){
    Tissue_plot1 + facet_wrap(.~face, scales = "free_y")
  }
  else{
    Tissue_plot1
  }
}

# Help function to geom_signif
sigFunc = function(x){
  # browser()
  if(x < 0.001){"***"} 
  else if(x < 0.01){"**"}
  else if(x < 0.05){"*"}
  else{NA}}

p <- notame_boxplot(imputed, feat_name = "HILIC_pos_160_1332a2_064", separate = "Diet",  comp = T)
p + ggtitle("5-AVAB")

p2 <- notame_boxplot(temp, feat_name = "HILIC_pos_138_05502_4_248", comp = T, face = "Diet")
p2 + ggtitle("Trigonelline")

# Plot otherway around
p2 <- notame_boxplot(temp, feat_name = "HILIC_pos_138_05502_4_248", separate = "Diet", comp = T, face = "Tissue")
p2 + ggtitle("Trigonelline")


###############################################################
source(file = "/labs/khatrilab/ananthg/detrimental_host_response/common_functions/perform_meta_analyses.R")
library(readr)
library(stringr)
###############################################################


###############################################################
GEO_data <- getGEOData("GSE88794")
GSE88794 <- GEO_data$originalData$GSE88794
saveRDS(GSE88794, file = "/labs/khatrilab/ananthg/detrimental_host_response/data_2024_09_27/GSE88794_raw.RDS")
###############################################################


###############################################################
#GSE88794 <- readRDS("/labs/khatrilab/ananthg/detrimental_host_response/data_2024_09_27/GSE88794_raw.RDS")
GSE88794_with_module_scores <- keep_signature_mod_genes(GSE88794)
combined_score_pheno <- cbind(GSE88794_with_module_scores$pheno, t(GSE88794_with_module_scores$expr[c("module_1_score", "module_2_score", "module_3_score", "module_4_score", "detrimental_score", "protective_score", "som_score", "tsang_ihm_score", "tsang_up_score", "tsang_down_score"), ]))
subject_info <- str_split(combined_score_pheno$source_name_ch1, ",", simplify = TRUE)
testing_info <- str_split(subject_info[, 2], " = ", simplify = TRUE)
day_info <- str_split(subject_info[, 3], " = ", simplify = TRUE)
colnames(day_info) <- c("day_of_testing", "type_of_test")
colnames(testing_info) <- c("week_of_testing", "intervention_status")
colnames(subject_info) <- c("subject_id", "week_of_testing_intervention", "day_type_testing", "time_point_testing")
combined_score_pheno <- cbind(combined_score_pheno, subject_info, testing_info, day_info)
combined_score_pheno$subject_id <- gsub("Peripheral blood mononuclear cells isolated from subject ", "", combined_score_pheno$subject_id)
write_csv(combined_score_pheno, "/labs/khatrilab/ananthg/detrimental_host_response/data_2024_09_27/GSE88794_combined_score_pheno.csv")
###############################################################

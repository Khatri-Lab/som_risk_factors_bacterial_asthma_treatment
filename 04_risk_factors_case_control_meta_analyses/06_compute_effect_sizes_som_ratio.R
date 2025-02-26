###############################################################
library(readr)
library(plyr)
library(tidyverse)
list_of_scores <- c("som_score_ratio", "detrimental_score", "protective_score", "module_1_score", "module_2_score", "module_3_score", "module_4_score")
list_of_score_names <- c("SoM", "Detrimental", "Protective", "Module 1", "Module 2", "Module 3", "Module 4")
###############################################################


###############################################################
combined_eff_size <- NULL
data_summary <- NULL
condition_name_list <- c("aging", "sex", "obesity", "diabetes")
condition_name_modified_list <- c("Age (\U2265 60 years vs. \U2264 40)", "Sex (male vs. female)", "Obesity (BMI \U2265 30 vs. \U2264 25)", "Diabetes (t1/t2 vs. healthy)")
combined_summary <- data.frame(matrix(data = NA, ncol = 3, nrow = 4), row.names = condition_name_modified_list)
colnames(combined_summary) <- c("Dataset", "Cases", "Controls")
for(i in 1:length(condition_name_list))
{
  meta_obj <- readRDS(paste0("/labs/khatrilab/ananthg/detrimental_host_response/data_2024_09_27/", condition_name_list[i], "_meta_obj.RDS"))
  condition_summary <- data.frame(Dataset = names(meta_obj$originalData), Cases = Reduce(c, lapply(meta_obj$originalData, function(x) sum(x$class == 1))), Controls = Reduce(c, lapply(meta_obj$originalData, function(x) sum(x$class == 0))))
  condition_summary$condition <- condition_name_modified_list[i]
  combined_summary[condition_name_modified_list[i], "Dataset"] <- length(names(meta_obj$originalData))
  sample_summary <- Reduce(c, lapply(meta_obj$originalData, function(x) x$class))
  combined_summary[condition_name_modified_list[i], "Cases"] <- sum(sample_summary == 1)
  combined_summary[condition_name_modified_list[i], "Controls"] <- sum(sample_summary == 0)
  score_eff_sizes <- meta_obj$metaAnalysis$pooledResults[list_of_scores, c("effectSize", "effectSizePval")]
  colnames(score_eff_sizes) <- c("eff_size", "p_val")
  score_eff_sizes$score_name <- list_of_score_names
  score_eff_sizes$condition <- condition_name_modified_list[i]
  if(is.null(combined_eff_size))
  {
    combined_eff_size <- score_eff_sizes
    data_summary <- condition_summary
  }
  else
  {
    combined_eff_size <- rbind(combined_eff_size, score_eff_sizes)
    data_summary <- rbind(data_summary, condition_summary)
  }
}
combined_eff_size$fdr <- p.adjust(combined_eff_size$p_val, method = "fdr")
###############################################################


###############################################################
row_order <- c("SoM", "Detrimental", "Protective", "Module 1", "Module 2", "Module 3", "Module 4")
clean_up_as_matrix <- function(wide_tibble, row_order, column_order)
{
  wide_matrix <- as.data.frame(wide_tibble)
  rownames_wide_matrix <- wide_tibble$score_name
  wide_matrix <- wide_matrix[, -c(1)]
  wide_matrix <- as.data.frame(sapply(wide_matrix, as.numeric))
  rownames(wide_matrix) <- rownames_wide_matrix
  wide_matrix <- t(as.matrix(wide_matrix[row_order, column_order]))
  return(wide_matrix)
}
eff_size_matrix <- clean_up_as_matrix((combined_eff_size[, c("score_name", "condition", "eff_size")] %>% pivot_wider(names_from = "condition", values_from = "eff_size")), row_order, condition_name_modified_list)
fdr_matrix <- clean_up_as_matrix((combined_eff_size[, c("score_name", "condition", "fdr")] %>% pivot_wider(names_from = "condition", values_from = "fdr")), row_order, condition_name_modified_list)
p_matrix <- clean_up_as_matrix((combined_eff_size[, c("score_name", "condition", "p_val")] %>% pivot_wider(names_from = "condition", values_from = "p_val")), row_order, condition_name_modified_list)
###############################################################


###############################################################
saveRDS(list(eff_size_matrix = eff_size_matrix, fdr_matrix = fdr_matrix, combined_summary = combined_summary), file = "/labs/khatrilab/ananthg/detrimental_host_response/data_2024_09_27/case_control_analysis_summary_som_ratio.rds")
write_csv(combined_eff_size, file = "/labs/khatrilab/ananthg/detrimental_host_response/data_2024_09_27/meta_analyses_effect_sizes_som_ratio.csv")
write_csv(combined_summary, file = "/labs/khatrilab/ananthg/detrimental_host_response/data_2024_09_27/meta_analyses_data_summary_som_ratio.csv")
write_csv(data_summary, file = "/labs/khatrilab/ananthg/detrimental_host_response/data_2024_09_27/meta_analyses_individual_dataset_summary_som_ratio.csv")
###############################################################

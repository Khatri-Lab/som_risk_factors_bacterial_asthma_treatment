###############################################################
library(readr)
library(plyr)
library(tidyverse)
source(file = "/labs/khatrilab/ananthg/detrimental_host_response/common_functions/perform_meta_analyses.R")
###############################################################


###############################################################
asthma_meta_obj <- readRDS("/labs/khatrilab/ananthg/detrimental_host_response/data_2024_09_27/asthma_GEO_datasets.RDS")
GSE134544 <- asthma_meta_obj[[2]]
GSE134544$pheno$Class <- 2
GSE134544$pheno$time_point <- sapply(GSE134544$pheno$title, function(x){str_split(x, " ")[[1]][3]})
GSE134544$pheno$Class[GSE134544$pheno$time_point == "Week" & GSE134544$pheno$`omalizumab responder status:ch1` == "Responder"] <- 0
GSE134544$pheno$Class[GSE134544$pheno$time_point == "Week" & GSE134544$pheno$`omalizumab responder status:ch1` == "Non-responder"] <- 1
GSE134544$class <- GSE134544$pheno$Class
GSE148725 <- readRDS("/labs/khatrilab/ananthg/detrimental_host_response/data_2024_09_27/asthma_GEO_GSE148725.RDS")
GSE148725$pheno$Class<- ifelse(GSE148725$pheno$`gete-responder:ch1` == "yes", 0, 1)
GSE148725$pheno$Class[GSE148725$pheno$`time point:ch1` != "at baseline"] <- 2
GSE148725$class <- GSE148725$pheno$Class
###############################################################


###############################################################
list_of_datasets <- list(GSE134544, GSE148725)
list_of_datasets_update <- list(GSE134544, GSE148725)
responder_meta_obj <- meta_analyze_datasets(list_of_datasets, list_of_datasets_update)
###############################################################


##########################################################
obtain_meta_analysis_summary <- function(meta_integrator_obj, gene_name, score_name)
{
  eff_size_summary <- data.frame(mean = meta_integrator_obj$metaAnalysis$datasetEffectSizes[gene_name, ], std_err = meta_integrator_obj$metaAnalysis$datasetEffectSizeStandardErrors[gene_name, colnames(meta_integrator_obj$metaAnalysis$datasetEffectSizes)])
  eff_size_summary$Dataset <- colnames(meta_integrator_obj$metaAnalysis$datasetEffectSizes)
  eff_size_summary <- rbind(eff_size_summary, c(meta_integrator_obj$metaAnalysis$pooledResults[gene_name, "effectSize"], meta_integrator_obj$metaAnalysis$pooledResults[gene_name, "effectSizeStandardError"], "Summary"))
  eff_size_summary$lower <- as.numeric(eff_size_summary$mean) - 1.96 * as.numeric(eff_size_summary$std_err)
  eff_size_summary$upper <- as.numeric(eff_size_summary$mean) + 1.96 * as.numeric(eff_size_summary$std_err)
  eff_size_summary <- as.data.frame(eff_size_summary)
  rownames(eff_size_summary) <- eff_size_summary$Dataset
  eff_size_summary$Controls <- 0
  eff_size_summary$Cases <- 0
  dataset_names_list <- names(meta_integrator_obj$originalData)
  for(dataset_name in dataset_names_list)
  {
    eff_size_summary[dataset_name, "Controls"] <- sum(meta_integrator_obj$originalData[[dataset_name]]$class == 0)
    eff_size_summary[dataset_name, "Cases"] <- sum(meta_integrator_obj$originalData[[dataset_name]]$class == 1)
    eff_size_summary["Summary", "Controls"] <- eff_size_summary["Summary", "Controls"] + sum(meta_integrator_obj$originalData[[dataset_name]]$class == 0)
    eff_size_summary["Summary", "Cases"] <- eff_size_summary["Summary", "Cases"] + sum(meta_integrator_obj$originalData[[dataset_name]]$class == 1)
  }
  summary_data <- data.frame(Dataset = eff_size_summary$Dataset, Controls = eff_size_summary$Controls, Cases = eff_size_summary$Cases)
  summary_data$`Effect size` <- "replace"
  summary_data$`Effect size (95% CI)` <- ifelse(is.na(eff_size_summary$mean), NA,
    sprintf("%.2f (%.2f to %.2f)",
      as.numeric(eff_size_summary$mean), as.numeric(eff_size_summary$lower), as.numeric(eff_size_summary$upper)))
  summary_data$eff_size <- eff_size_summary$mean
  summary_data$lower <- eff_size_summary$lower
  summary_data$upper <- eff_size_summary$upper
  summary_data$box_size <- ifelse(is.na(eff_size_summary$mean), NA, 1 / (as.numeric(eff_size_summary$std_err))^2)
  summary_data$box_size <- as.numeric(summary_data$box_size)
  summary_data$box_size <- 0.25 * (summary_data$box_size - min(summary_data$box_size)) / (max(summary_data$box_size) - min(summary_data$box_size)) + 0.5
  summary_data <- summary_data[, c("Dataset", "Cases", "Controls", "Effect size", "Effect size (95% CI)", "eff_size", "lower", "upper", "box_size")]
  colnames(summary_data) <- c("Dataset", "NR", "R", score_name, "Effect size (95% CI)", "eff_size", "lower", "upper", "box_size")
  return(summary_data)
}
##########################################################


###############################################################
som_summary <- obtain_meta_analysis_summary(responder_meta_obj, "som_score", "SoM")
detrimental_summary <- obtain_meta_analysis_summary(responder_meta_obj, "detrimental_score", "Detrimental")
protective_summary <- obtain_meta_analysis_summary(responder_meta_obj, "protective_score", "Protective")
module_1_summary <- obtain_meta_analysis_summary(responder_meta_obj, "module_1_score", "Module 1")
module_2_summary <- obtain_meta_analysis_summary(responder_meta_obj, "module_2_score", "Module 2")
module_3_summary <- obtain_meta_analysis_summary(responder_meta_obj, "module_3_score", "Module 3")
module_4_summary <- obtain_meta_analysis_summary(responder_meta_obj, "module_4_score", "Module 4")
###############################################################


###############################################################
combine_summaries <- function(first_module, second_module, third_module, fourth_module)
{
  colnames(first_module)[5:9] <- c("ignore", "eff_size_1", "lower_1", "upper_1", "box_size_1")
  colnames(second_module)[5:9] <- c("ignore", "eff_size_2", "lower_2", "upper_2", "box_size_2")
  combined_summary <- join(first_module, second_module, by = c("Dataset", "NR", "R"))
  combined_summary$space_1 <- "replace"
  if(!missing(third_module))
  {
    colnames(third_module)[5:9] <- c("ignore", "eff_size_3", "lower_3", "upper_3", "box_size_3")
    combined_summary <- join(combined_summary, third_module, by = c("Dataset", "NR", "R"))
    combined_summary$space_2 <- "replace"
  }
  if(!missing(fourth_module))
  {
    colnames(fourth_module)[5:9] <- c("ignore", "eff_size_4", "lower_4", "upper_4", "box_size_4")
    combined_summary <- join(combined_summary, fourth_module, by = c("Dataset", "NR", "R"))
    combined_summary$space_3 <- "replace"
  }
  return(combined_summary)
}
main_fig_summary <- combine_summaries(som_summary, detrimental_summary, protective_summary)
main_fig_summary$Dataset[main_fig_summary$Dataset == "GSE134544"] <- "Omalizumab"
main_fig_summary$Dataset[main_fig_summary$Dataset == "GSE148725"] <- "Benralizumab"
supp_fig_summary <- combine_summaries(module_1_summary, module_2_summary, module_3_summary, module_4_summary)
supp_fig_summary$Dataset[supp_fig_summary$Dataset == "GSE134544"] <- "Omalizumab"
supp_fig_summary$Dataset[supp_fig_summary$Dataset == "GSE148725"] <- "Benralizumab"
###############################################################


###############################################################
add_empty_rows <- function(summary_data)
{
  summary_data <- summary_data %>% add_row(.before = 3)
  summary_data[3, ] <- "replace"
  return(summary_data)
}
main_fig_summary <- add_empty_rows(main_fig_summary)
supp_fig_summary <- add_empty_rows(supp_fig_summary)
write_csv(main_fig_summary, "/labs/khatrilab/ananthg/detrimental_host_response/data_2024_09_27/asthma_biologics_meta_analysis_main_fig_summary.csv")
write_csv(supp_fig_summary, "/labs/khatrilab/ananthg/detrimental_host_response/data_2024_09_27/asthma_biologics_meta_analysis_supp_fig_summary.csv")
###############################################################

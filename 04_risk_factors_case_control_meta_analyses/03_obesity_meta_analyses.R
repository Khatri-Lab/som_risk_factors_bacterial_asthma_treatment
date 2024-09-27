###############################################################
source(file = "/labs/khatrilab/ananthg/detrimental_host_response/common_functions/perform_meta_analyses.R")
DHR_Obesity_metaObj <- readRDS(file = "/labs/khatrilab/jiayingt/data/DHR_Obesity/DHR_Obesity_metaAnalysis_v05_noEllison/dataset_files/05_metaAnalysisResults/DHR_Obesity_metaObj_v05_metaAnalysisResults.RDS")
###############################################################


###############################################################
list_of_datasets <- list(DHR_Obesity_metaObj$originalData$GSE18897, DHR_Obesity_metaObj$originalData$GSE41233, DHR_Obesity_metaObj$originalData$GSE53232, DHR_Obesity_metaObj$originalData$GSE110551)
list_of_datasets_update <- list(DHR_Obesity_metaObj$originalData$GSE18897, DHR_Obesity_metaObj$originalData$GSE41233, DHR_Obesity_metaObj$originalData$GSE53232, DHR_Obesity_metaObj$originalData$GSE110551)
obesity_meta_obj <- meta_analyze_datasets(list_of_datasets, list_of_datasets_update)
saveRDS(obesity_meta_obj, file = "/labs/khatrilab/ananthg/detrimental_host_response/data_2024_09_27/obesity_meta_obj.RDS")
###############################################################

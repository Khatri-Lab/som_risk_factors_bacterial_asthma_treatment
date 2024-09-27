###############################################################
source(file = "/labs/khatrilab/ananthg/detrimental_host_response/common_functions/perform_meta_analyses.R")
DHR_Aging_metaObj <- readRDS("/labs/khatrilab/jiayingt/data/DHR_Aging/DHR_Aging_metaAnalysis_v04_noEllison/dataset_files/05_metaAnalysisResults/DHR_Aging_metaObj_v04_metaAnalysisResults.RDS")
###############################################################


###############################################################
list_of_datasets <- list(DHR_Aging_metaObj$originalData$GSE59635, DHR_Aging_metaObj$originalData$GSE59654, DHR_Aging_metaObj$originalData$GSE59743, DHR_Aging_metaObj$originalData$GSE79396_Atlanta, DHR_Aging_metaObj$originalData$GSE79396_Denver, DHR_Aging_metaObj$originalData$GSE101709, DHR_Aging_metaObj$originalData$GSE101710, DHR_Aging_metaObj$originalData$GSE107990)
list_of_datasets_update <- list(DHR_Aging_metaObj$originalData$GSE59635, DHR_Aging_metaObj$originalData$GSE59654, DHR_Aging_metaObj$originalData$GSE59743, DHR_Aging_metaObj$originalData$GSE79396_Atlanta, DHR_Aging_metaObj$originalData$GSE79396_Denver, DHR_Aging_metaObj$originalData$GSE101709, DHR_Aging_metaObj$originalData$GSE101710, DHR_Aging_metaObj$originalData$GSE107990)
aging_meta_obj <- meta_analyze_datasets(list_of_datasets, list_of_datasets_update)
saveRDS(aging_meta_obj, file = "/labs/khatrilab/ananthg/detrimental_host_response/data_2024_09_27/aging_meta_obj.RDS")
###############################################################

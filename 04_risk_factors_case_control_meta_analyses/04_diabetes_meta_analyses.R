###############################################################
source(file = "/labs/khatrilab/ananthg/detrimental_host_response/common_functions/perform_meta_analyses.R")
DHR_T1DT2D_metaObj <- readRDS("/labs/khatrilab/jiayingt/data/DHR_Diabetes/DHR_combinedT1DT2D_metaAnalysis_v03/dataset_files/04_unwantedSamplesRemoved/DHR_T1DT2D_metaObj_v03_metaAnalysisResults.RDS")
###############################################################


###############################################################
sum(rownames(DHR_T1DT2D_metaObj$originalData$GSE9006_T2D$expr) != rownames(DHR_T1DT2D_metaObj$originalData$GSE9006_T1D$expr))
sum(colnames(DHR_T1DT2D_metaObj$originalData$GSE9006_T2D$pheno) != colnames(DHR_T1DT2D_metaObj$originalData$GSE9006_T1D$pheno))
GSE9006_expr <- cbind(DHR_T1DT2D_metaObj$originalData$GSE9006_T2D$expr, DHR_T1DT2D_metaObj$originalData$GSE9006_T1D$expr[, DHR_T1DT2D_metaObj$originalData$GSE9006_T1D$class == 1])
GSE9006_pheno <- rbind(DHR_T1DT2D_metaObj$originalData$GSE9006_T2D$pheno, DHR_T1DT2D_metaObj$originalData$GSE9006_T1D$pheno[DHR_T1DT2D_metaObj$originalData$GSE9006_T1D$class == 1, ])
GSE9006_class <- c(DHR_T1DT2D_metaObj$originalData$GSE9006_T2D$class, DHR_T1DT2D_metaObj$originalData$GSE9006_T1D$class[DHR_T1DT2D_metaObj$originalData$GSE9006_T1D$class == 1])
GSE9006_keys <- DHR_T1DT2D_metaObj$originalData$GSE9006_T2D$keys
GSE9006 <- list(expr = GSE9006_expr, pheno = GSE9006_pheno, class = GSE9006_class, keys = GSE9006_keys, formattedName = "GSE9006")
###############################################################


###############################################################
list_of_datasets <- list(GSE9006, DHR_T1DT2D_metaObj$originalData$GSE44314, DHR_T1DT2D_metaObj$originalData$GSE55098, DHR_T1DT2D_metaObj$originalData$GSE72376, DHR_T1DT2D_metaObj$originalData$GSE72377, DHR_T1DT2D_metaObj$originalData$GSE156035, DHR_T1DT2D_metaObj$originalData$GSE15932, DHR_T1DT2D_metaObj$originalData$GSE23561, DHR_T1DT2D_metaObj$originalData$GSE95849, DHR_T1DT2D_metaObj$originalData$GSE142153)
list_of_datasets_update <- list(GSE9006, DHR_T1DT2D_metaObj$originalData$GSE44314, DHR_T1DT2D_metaObj$originalData$GSE55098, DHR_T1DT2D_metaObj$originalData$GSE72376, DHR_T1DT2D_metaObj$originalData$GSE72377, DHR_T1DT2D_metaObj$originalData$GSE156035, DHR_T1DT2D_metaObj$originalData$GSE15932, DHR_T1DT2D_metaObj$originalData$GSE23561, DHR_T1DT2D_metaObj$originalData$GSE95849, DHR_T1DT2D_metaObj$originalData$GSE142153)
diabetes_meta_obj <- meta_analyze_datasets(list_of_datasets, list_of_datasets_update)
saveRDS(diabetes_meta_obj, file = "/labs/khatrilab/ananthg/detrimental_host_response/data_2024_09_27/diabetes_meta_obj.RDS")
###############################################################

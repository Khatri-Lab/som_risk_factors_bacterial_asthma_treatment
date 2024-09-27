###############################################################
source(file = "/labs/khatrilab/ananthg/detrimental_host_response/common_functions/perform_meta_analyses.R")
load("/labs/khatrilab/ebongen/sexDifferences2.0/isexs/sexMetaObj.RData")
###############################################################


###############################################################
flip_case_control_sex_meta_obj <- function(dataset)
{
  dataset$class <- 1 - dataset$class
  dataset$pheno$Class <- dataset$class
  return(dataset)
}
###############################################################


###############################################################
list_of_datasets <- list(sexMetaObj$originalData$gse17065, sexMetaObj$originalData$gse19151, sexMetaObj$originalData$gse21862, sexMetaObj$originalData$gse47353, sexMetaObj$originalData$gse53195, sexMetaObj$originalData$gse60491, valMetaObj$originalData$gse13485, valMetaObj$originalData$gse18323, valMetaObj$originalData$gse19442, valMetaObj$originalData$gse21311, valMetaObj$originalData$gse30453, valMetaObj$originalData$gse37069, valMetaObj$originalData$gse38484, valMetaObj$originalData$gse58137, valMetaObj$originalData$gse61821, valMetaObj$originalData$gse65219)
list_of_datasets <- lapply(list_of_datasets, flip_case_control_sex_meta_obj)
list_of_datasets_update <- list_of_datasets
sex_meta_obj <- meta_analyze_datasets(list_of_datasets, list_of_datasets_update)
saveRDS(sex_meta_obj, file = "/labs/khatrilab/ananthg/detrimental_host_response/data_2024_09_27/sex_meta_obj.RDS")
###############################################################

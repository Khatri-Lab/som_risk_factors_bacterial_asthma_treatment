###############################################################
library(MetaIntegrator)
library(psych)
library(magrittr)
library(tidyverse)
###############################################################


###############################################################
genes_mod_1 <- c("NQO2", "SLPI", "ORM1", "KLHL2", "ANXA3", "TXN", "AQP9", "BCL6", "DOK3", "PFKFB4", "TYK2")
genes_mod_2 <- c("BCL2L11", "BCAT1", "BTBD7", "CEP55", "HMMR", "PRC1", "KIF15", "CAMP", "CEACAM8", "DEFA4", "LCN2", "CTSG", "AZU1")
genes_mod_3 <- c("MAFB", "OASL", "UBE2L6", "VAMP5", "CCL2", "NAPA", "ATG3", "VRK2", "TMEM123", "CASP7")
genes_mod_4 <- c("DOK2", "HLA-DPB1", "BUB3", "SMYD2", "SIDT1", "EXOC2", "TRIB2", "KLRB1")
genes_tsang_up <- c("SLC16A10", "NDRG2", "AK5", "CD7", "RNFT2", "PHC1", "MAN1C1", "FAM102A", "RCAN3", "EIF3H", "RPS17", "RPS5", "PLXDC1", "TESPA1", "SGK223", "EIF3L", "ANAPC16", "CD27", "NPAT", "ID3", "RACK1", "APEX1", "GCNT4", "FCMR", "TCF7", "KLHL3", "AXIN2", "LY9", "RPS25", "LDHB", "PKIA", "RPL3", "N6AMT1", "GAL3ST4", "SSBP2", "CD1C", "LEF1", "RPL7", "PIK3IP1", "GPRASP1", "ABI2", "APBB1", "SPTBN1", "GPA33", "CCR9", "BCKDHB", "SCAI", "RPL4", "NOG", "TCEA3", "ETS1", "LDLRAP1", "GPR183", "ZNF548", "ZNF91", "NPM1", "MSANTD2", "KAT6B", "SLC7A6", "DCHS1", "OXNAD1", "RPS2", "RPL7A", "GRAP", "RPL23A", "RPL10A", "RPS3", "FAM175A", "RPL29", "EEF2", "EDAR", "ABLIM1", "MBLAC2", "CCR7", "ZNF573", "CAMK4", "LRRN3", "MAGI3", "RPLP2", "ZIK1", "NT5E", "FUT8", "ZNF101", "RPL34", "RPS20", "FOXP1", "ZNF550", "TSPYL2", "ATP6V0E2-AS1", "GRPEL2", "MGC57346", "SEPT1", "PRKACB", "AGMAT", "RPL11", "MYC", "ZZZ3", "RPL5")
genes_tsang_down <- c("CTRL", "NTNG2", "AP5B1", "PDCD1LG2", "DOCK4", "S100A9", "BACH1", "FAM8A1", "SECTM1", "S100A8", "TYMP", "HK3", "IL1RN", "MYD88", "REC8", "ALPK1", "SAT1", "PRKCD", "SLC26A8", "PARP9", "LMNB1", "RELT", "TAP1", "JAK2", "BRI3", "GBP2", "PLEKHO2", "ETV7", "ODF3B", "SIGLEC5", "CEACAM1", "CARD16", "ZBP1", "DDX60L", "APOL2", "CD63", "TNFAIP6", "KCNJ2", "ANKRD22", "SCARF1", "SEMA4A", "DNAJC5", "SQRDL", "HELZ2", "GADD45B", "FAS", "HSPA6", "PIK3AP1", "CLEC7A", "SERPING1", "OR52K2", "ITPRIP")
###############################################################


###############################################################
calculate_signature_score <- function(gene_list, dataset)
{
  expr_matrix <- MetaIntegrator::getSampleLevelGeneData(dataset, gene_list)
  expr_matrix[expr_matrix == 0] <- NA
  signature_score <- geometric.mean(expr_matrix)
  return(signature_score)
}
###############################################################


###############################################################
keep_signature_mod_genes <- function(dataset)
{
  new_expr <- dataset$expr
  original_expr_rows <- nrow(new_expr)
  new_expr <- rbind(new_expr, calculate_signature_score(genes_mod_1, dataset))
  new_expr <- rbind(new_expr, calculate_signature_score(genes_mod_2, dataset))
  new_expr <- rbind(new_expr, calculate_signature_score(genes_mod_3, dataset))
  new_expr <- rbind(new_expr, calculate_signature_score(genes_mod_4, dataset))
  new_expr <- rbind(new_expr, calculate_signature_score(genes_tsang_up, dataset))
  new_expr <- rbind(new_expr, calculate_signature_score(genes_tsang_down, dataset))
  rownames(new_expr)[(original_expr_rows + 1): nrow(new_expr)] <- c("module_1_score", "module_2_score", "module_3_score", "module_4_score", "tsang_up_score", "tsang_down_score")
  original_expr_rows <- nrow(new_expr)
  new_expr <- rbind(new_expr, new_expr["module_1_score", ] + new_expr["module_2_score", ])
  new_expr <- rbind(new_expr, new_expr["module_3_score", ] + new_expr["module_4_score", ])
  new_expr <- rbind(new_expr, (new_expr["module_1_score", ] + new_expr["module_2_score", ]) - (new_expr["module_3_score", ] + new_expr["module_4_score", ]))
  new_expr <- rbind(new_expr, (new_expr["module_1_score", ] + new_expr["module_2_score", ]) / (new_expr["module_3_score", ] + new_expr["module_4_score", ]))
  new_expr <- rbind(new_expr, new_expr["tsang_up_score", ] - new_expr["tsang_down_score", ])
  rownames(new_expr)[(original_expr_rows + 1): nrow(new_expr)] <- c("detrimental_score", "protective_score", "som_score", "som_score_ratio", "tsang_ihm_score")
  new_keys <- c(dataset$keys, "module_1_score", "module_2_score", "module_3_score", "module_4_score", "tsang_up_score", "tsang_down_score", "detrimental_score", "protective_score", "som_score", "som_score_ratio", "tsang_ihm_score")
  names(new_keys) <- c(names(dataset$keys), "module_1_score", "module_2_score", "module_3_score", "module_4_score", "tsang_up_score", "tsang_down_score", "detrimental_score", "protective_score", "som_score", "som_score_ratio", "tsang_ihm_score")
  dataset$keys <- new_keys
  dataset$expr <- new_expr
  return(dataset)
}
###############################################################


###############################################################
meta_analyze_datasets <- function(list_of_datasets, list_of_datasets_update)
{
  for(x in 1:length(list_of_datasets))
  {
    dataset <- list_of_datasets[[x]]
    remove_index <- which(dataset$pheno$Class == 2)
    if(length(remove_index) > 0)
    {
      dataset$pheno <- dataset$pheno[-remove_index, ]
      dataset$expr <- dataset$expr[, -remove_index]
      dataset$class <- dataset$class[-remove_index]
    }
    checkDataObject(dataset, "Dataset")
    dataset_min <- min(dataset$expr, na.rm = TRUE)
    dataset$expr[is.na(dataset$expr)] <- dataset_min
    if(dataset_min < 0)
    {
      dataset$expr <- dataset$expr - dataset_min + 1
    }
    appended_dataset <- keep_signature_mod_genes(dataset)
    list_of_datasets_update[[x]] <- appended_dataset
  }
  discovery_datasets <- list_of_datasets_update
  dataset_nms <- sapply(list_of_datasets_update, function(x) x$formattedName)
  names(discovery_datasets) <- dataset_nms
  meta_obj <- list()
  meta_obj$originalData <- discovery_datasets
  meta_analysis_obj_check <- checkDataObject(meta_obj, "Meta", "Pre-Analysis")
  if(meta_analysis_obj_check)
  {
    meta_analysis_results <- runMetaAnalysis(meta_obj, maxCores = 11)
    return(meta_analysis_results)
  }
  else
  {
    return(meta_analysis_obj_check)
  }
}
###############################################################


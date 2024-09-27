###############################################################
library(data.table)
library(psych)
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
read_process_phenotype_file <- function(file_path, subj_id_to_keep)
{
  pheno_file <- fread(file_path)
  colnames(pheno_file) <- as.character(pheno_file[1, ])
  pheno_file <- pheno_file[-1, ]
  pheno_file <- as.data.frame(sapply(pheno_file, as.numeric))
  pheno_file <- pheno_file[pheno_file$shareid %in% subj_id_to_keep, ]
  return(pheno_file)
}
###############################################################


###############################################################
compute_module_score <- function(gene_expr_matrix, gene_signature, module_title)
{
  expr_matrix <- gene_expr_matrix[gene_expr_matrix$gene_name %in% gene_signature, 2:ncol(gene_expr_matrix)]
  expr_matrix[expr_matrix == 0] <- NA
  module_score <- geometric.mean(expr_matrix)
  module_score <- as.data.frame(module_score)
  colnames(module_score) <- c(module_title)
  module_score$shareid <- rownames(module_score)
  return(module_score)
}
###############################################################


###############################################################
combine_module_scores <- function(gene_expr_matrix)
{
  module_1_score <- compute_module_score(gene_expr_matrix, genes_mod_1, "module_1_score")
  module_2_score <- compute_module_score(gene_expr_matrix, genes_mod_2, "module_2_score")
  module_3_score <- compute_module_score(gene_expr_matrix, genes_mod_3, "module_3_score")
  module_4_score <- compute_module_score(gene_expr_matrix, genes_mod_4, "module_4_score")
  tsang_up_score <- compute_module_score(gene_expr_matrix, genes_tsang_up, "tsang_up_score")
  tsang_down_score <- compute_module_score(gene_expr_matrix, genes_tsang_down, "tsang_down_score")
  score_matrix <- merge(module_1_score, module_2_score, by = "shareid")
  score_matrix <- merge(score_matrix, module_3_score, by = "shareid")
  score_matrix <- merge(score_matrix, module_4_score, by = "shareid")
  score_matrix <- merge(score_matrix, tsang_up_score, by = "shareid")
  score_matrix <- merge(score_matrix, tsang_down_score, by = "shareid")
  score_matrix$detrimental_score <- score_matrix$module_1_score + score_matrix$module_2_score
  score_matrix$protective_score <- score_matrix$module_3_score + score_matrix$module_4_score
  score_matrix$som_score <- score_matrix$detrimental_score - score_matrix$protective_score
  score_matrix$tsang_ihm_score <- score_matrix$tsang_up_score - score_matrix$tsang_down_score
  return(score_matrix)
}
###############################################################


###############################################################
obtain_module_score_matrix <- function(gene_expr_matrix)
{
  score_matrix <- combine_module_scores(gene_expr_matrix)
  score_matrix$shareid <- as.numeric(score_matrix$shareid)
  score_matrix$healthy <- 0
  score_matrix$cancer <- 0
  score_matrix$cardio <- 0
  score_matrix$diabetes <- 0
  score_matrix$fever <- 0
  score_matrix$dementia <- 0
  score_matrix$clinical_diagnosis_disease <- 0
  return(score_matrix)
}
###############################################################

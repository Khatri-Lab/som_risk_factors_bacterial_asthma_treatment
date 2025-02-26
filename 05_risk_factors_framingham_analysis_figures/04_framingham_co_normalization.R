###############################################################
library(readr)
source(file = "/labs/khatrilab/ananthg/detrimental_host_response/common_functions/perform_meta_analyses.R")
source(file = "/labs/khatrilab/ananthg/detrimental_host_response/common_functions/framingham_score_phenotype_support_functions.R")
###############################################################


###############################################################
create_meta_obj <- function(cohort_name)
{
  expr <- as.data.frame(read_csv(paste0("/local-scratch/projects/candx/guangbo/ananth_processesed_framingham_data/dhr_paper_re_processesed_2024_09_27/", cohort_name, "_expr_exon.csv")))
  rownames(expr) <- expr$gene_name
  expr <- expr[, -c(1, 2)]
  pheno_matrix <- as.data.frame(read_csv(paste0("/local-scratch/projects/candx/guangbo/ananth_processesed_framingham_data/dhr_paper_re_processesed_2024_09_27/", cohort_name, "_module_score_pheno_exon.csv")))
  rownames(pheno_matrix) <- as.character(pheno_matrix$shareid)
  pheno_meta_integrator <- pheno_matrix[match(colnames(expr), rownames(pheno_matrix)), ]
  pheno_meta_integrator$Class <- 1 - pheno_meta_integrator$healthy
  print(table(pheno_meta_integrator$Class))
  expr_keys <- rownames(expr)
  names(expr_keys) <- expr_keys
  expr_mtrx <- as.matrix(sapply(expr, as.numeric))
  rownames(expr_mtrx) <- rownames(expr)
  framingham_class <- pheno_meta_integrator$Class
  names(framingham_class) <- rownames(pheno_meta_integrator)
  framingham_obj <- list(expr = expr_mtrx, keys = expr_keys, pheno = pheno_meta_integrator, formattedName = cohort_name, class = framingham_class)
  checkDataObject(framingham_obj, "Dataset")
  return(framingham_obj)
}
###############################################################


###############################################################
off_spr_obj <- create_meta_obj("off_spr")
gen_iii_obj <- create_meta_obj("gen_iii")
framingham_meta_obj <- list(originalData = list(off_spr = off_spr_obj, gen_iii = gen_iii_obj))
framingham_conormalized_results <- coconutMetaIntegrator(framingham_meta_obj)
saveRDS(framingham_conormalized_results, file = "/local-scratch/projects/candx/guangbo/ananth_processesed_framingham_data/dhr_paper_re_processesed_2024_09_27/framingham_conormalized_by_healthy_obj.RDS")
###############################################################


###############################################################
#framingham_conormalized_results <- readRDS("/local-scratch/projects/candx/guangbo/ananth_processesed_framingham_data/dhr_paper_re_processesed_2024_09_27/framingham_conormalized_by_healthy_obj.RDS")
conorm_expr <- Reduce(cbind, list(framingham_conormalized_results$COCONUTList$off_spr$genes,
  framingham_conormalized_results$COCONUTList$gen_iii$genes,
  framingham_conormalized_results$controlList$GSEs$off_spr$genes,
  framingham_conormalized_results$controlList$GSEs$gen_iii$genes))
conorm_expr <- cbind(gene_name = rownames(conorm_expr), conorm_expr)
conorm_expr <- combine_module_scores(conorm_expr)
pheno_to_keep <- c("age", "sex", "bmi", "weight", "height", "regular_smokers", "smoking_quiters", "smoking_quitting_age", "healthy", "cancer", "cardio", "diabetes", "fever", "dementia", "clinical_diagnosis_disease", "A1C", "first_cardio_diagnosis_date", "first_cancer_diagnosis_date", "death_date", "death_cause", "c_reactive_protein_serum")
conorm_pheno <- Reduce(rbind, list(framingham_conormalized_results$COCONUTList$off_spr$pheno[, pheno_to_keep],
  framingham_conormalized_results$COCONUTList$gen_iii$pheno[, pheno_to_keep],
  framingham_conormalized_results$controlList$GSEs$off_spr$pheno[, pheno_to_keep],
  framingham_conormalized_results$controlList$GSEs$gen_iii$pheno[, pheno_to_keep]))
conorm_pheno$shareid <- rownames(conorm_pheno)
conorm_matrix <- merge(conorm_pheno, conorm_expr, by = "shareid")
write_csv(conorm_matrix, file = "/local-scratch/projects/candx/guangbo/ananth_processesed_framingham_data/dhr_paper_re_processesed_2024_09_27/framingham_conormalized_by_healthy_matrix.csv")
###############################################################

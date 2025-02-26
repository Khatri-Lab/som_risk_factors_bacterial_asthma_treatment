###############################################################
library(readr)
library(tidyverse)
source(file = "/labs/khatrilab/ananthg/detrimental_host_response/common_functions/framingham_score_phenotype_support_functions.R")
###############################################################


###############################################################
gen_iii_expr <- read_csv("/local-scratch/projects/candx/guangbo/ananth_processesed_framingham_data/dhr_paper_re_processesed_2024_09_27/gen_iii_expr_exon.csv")
gen_iii_expr <- gen_iii_expr[, -1]
gen_iii_subj <- as.numeric(colnames(gen_iii_expr)[2:ncol(gen_iii_expr)])
###############################################################


###############################################################
gen_iii_subj_pheno <- read_process_phenotype_file("/local-scratch/projects/candx/guangbo/Framingham_data_part1/PhenoGenotypeFiles/RootStudyConsentSet_phs000007.Framingham.v31.p12.c1.HMB-IRB-MDS/PhenotypeFiles/phs000007.v31.pht003099.v6.p12.c1.vr_dates_2019_a_1175s.HMB-IRB-MDS.txt", gen_iii_subj)
gen_iii_subj_pheno <- gen_iii_subj_pheno[, c("shareid", "idtype", "sex", "age2", "date2", "att3", "age3", "att1", "age1")]
gen_iii_subj_pheno$sex[gen_iii_subj_pheno$sex == 1] <- "M"
gen_iii_subj_pheno$sex[gen_iii_subj_pheno$sex == 2] <- "F"
colnames(gen_iii_subj_pheno) <- c("shareid", "idtype", "sex", "age_at_exam_2", "days_to_expr_from_exam_1", "attended_next_exam", "age_at_next_exam", "attended_previous_exam", "age_at_prervious_exam")
###############################################################


###############################################################
#cardio_events <- read_process_phenotype_file("/local-scratch/projects/candx/guangbo/Framingham_data_part1/PhenoGenotypeFiles/RootStudyConsentSet_phs000007.Framingham.v31.p12.c1.HMB-IRB-MDS/PhenotypeFiles/phs000007.v31.pht003316.v8.p12.c1.vr_survcvd_2017_a_1194s.HMB-IRB-MDS.txt", gen_iii_subj)
#cardio_events <- cardio_events[cardio_events$chd == 1 | cardio_events$chd == 1 | cardio_events$chd == 1, ]
cardio_events <- read_process_phenotype_file("/local-scratch/projects/candx/guangbo/Framingham_data_part1/PhenoGenotypeFiles/RootStudyConsentSet_phs000007.Framingham.v31.p12.c1.HMB-IRB-MDS/PhenotypeFiles/phs000007.v31.pht000309.v14.p12.c1.vr_soe_2018_a_1311s.HMB-IRB-MDS.txt", gen_iii_subj)
cardio_events <- cardio_events[cardio_events$EVENT %in% c(1, 2, 3, 4, 5, 6, 7, 10, 11, 12, 13, 16, 17, 18, 21, 22, 23, 24, 25, 26, 30, 39, 40, 41), ]
first_cardio_event <- cardio_events[, c("shareid", "DATE")] %>% group_by(shareid) %>% summarise(first_cardio = min(DATE))
gen_iii_subj_pheno <- merge(gen_iii_subj_pheno, first_cardio_event, by = "shareid", all = TRUE)
gen_iii_subj_pheno$first_cardio_diagnosis_date <- gen_iii_subj_pheno$first_cardio - gen_iii_subj_pheno$days_to_expr_from_exam_1
###############################################################


###############################################################
death_events <- read_process_phenotype_file("/local-scratch/projects/candx/guangbo/Framingham_data_part1/PhenoGenotypeFiles/RootStudyConsentSet_phs000007.Framingham.v31.p12.c1.HMB-IRB-MDS/PhenotypeFiles/phs000007.v31.pht003317.v8.p12.c1.vr_survdth_2017_a_1192s.HMB-IRB-MDS.txt", gen_iii_subj)
death_events <- death_events[!is.na(death_events$DATEDTH), ]
gen_iii_subj_pheno <- merge(gen_iii_subj_pheno, death_events[, c("shareid", "DATEDTH")], by = "shareid", all = TRUE)
gen_iii_subj_pheno$death_date <- gen_iii_subj_pheno$DATEDTH - gen_iii_subj_pheno$days_to_expr_from_exam_1
gen_iii_subj_pheno$death_cause <- NA
gen_iii_subj_pheno$death_cause[!is.na(gen_iii_subj_pheno$death_date)] <- "unknown"
gen_iii_subj_pheno$death_cause[gen_iii_subj_pheno$shareid %in% cardio_events$shareid[cardio_events$EVENT %in% c(3, 9, 21, 22, 23, 24, 25, 26)]] <- "cardio"
gen_iii_subj_pheno$death_cause[gen_iii_subj_pheno$shareid %in% cardio_events$shareid[cardio_events$EVENT %in% c(27)]] <- "cancer"
###############################################################


###############################################################
cancer_events <- read_process_phenotype_file("/local-scratch/projects/candx/guangbo/Framingham_data_part1/PhenoGenotypeFiles/RootStudyConsentSet_phs000007.Framingham.v31.p12.c1.HMB-IRB-MDS/PhenotypeFiles/phs000007.v31.pht000039.v12.p12.c1.vr_cancer_2013_a_0018s.HMB-IRB-MDS.txt", gen_iii_subj)
first_cancer_event <- cancer_events[, c("shareid", "d_date")] %>% group_by(shareid) %>% summarise(first_cancer = min(d_date))
gen_iii_subj_pheno <- merge(gen_iii_subj_pheno, first_cancer_event, by = "shareid", all = TRUE)
gen_iii_subj_pheno$first_cancer_diagnosis_date <- gen_iii_subj_pheno$first_cancer - gen_iii_subj_pheno$days_to_expr_from_exam_1
###############################################################


###############################################################
gen_iii_clinic_exam_2 <- read_process_phenotype_file("/local-scratch/projects/candx/guangbo/Framingham_data_part1/PhenoGenotypeFiles/RootStudyConsentSet_phs000007.Framingham.v31.p12.c1.HMB-IRB-MDS/PhenotypeFiles/phs000007.v31.pht003094.v6.p12.c1.e_exam_2011_m_0017s.HMB-IRB-MDS.txt", gen_iii_subj)
gen_iii_clinic_exam_2_weight_height <- gen_iii_clinic_exam_2[, c("shareid", "g3b0481", "g3b0483")]
colnames(gen_iii_clinic_exam_2_weight_height) <- c("shareid", "weight", "height")
gen_iii_clinic_exam_2_weight_height$bmi <- 703 * gen_iii_clinic_exam_2_weight_height$weight / (gen_iii_clinic_exam_2_weight_height$height ^ 2)
gen_iii_subj_pheno <- merge(gen_iii_subj_pheno, gen_iii_clinic_exam_2_weight_height, by = "shareid", all = TRUE)
gen_iii_clinic_exam_2_clinic_diagnosis <- gen_iii_clinic_exam_2[, c("shareid", "g3b0404", "g3b0405", "g3b0406", "g3b0408", "g3b0407", "g3b0410", "g3b0411", "g3b0414", "g3b0415", "g3b0416", "g3b0418", "g3b0420", "g3b0422", "g3b0424", "g3b0427", "g3b0429", "g3b0430", "g3b0431", "g3b0432", "g3b0434", "g3b0435", "g3b0436", "g3b0439", "g3b0441", "g3b0442", "g3b0444", "g3b0445", "g3b0452", "g3b0453", "g3b0454", "g3b0448", "g3b0449", "g3b0426", "g3b0446")]
gen_iii_clinic_exam_2_clinic_diagnosis[is.na(gen_iii_clinic_exam_2_clinic_diagnosis)] <- 3
gen_iii_clinic_exam_2_clinic_diagnosis$disease_status <- rowSums(gen_iii_clinic_exam_2_clinic_diagnosis[, 2:ncol(gen_iii_clinic_exam_2_clinic_diagnosis)])
gen_iii_clinic_exam_2_fever <- gen_iii_clinic_exam_2[, c("shareid", "g3b0002", "g3b0003", "g3b0004", "g3b0005", "g3b0007")]
colnames(gen_iii_clinic_exam_2_fever) <- c("shareid", "hospitalization_since_exam_7", "ER_visit_since_exam_7", "day_surgery", "major_illness_visit_to_doctor", "fever_or_infection_past_two_weeks")
gen_iii_subj_pheno <- merge(gen_iii_subj_pheno, gen_iii_clinic_exam_2_fever, by = "shareid", all = TRUE)
not_regular_smokers <- gen_iii_clinic_exam_2$shareid[gen_iii_clinic_exam_2$g3b0096 == 0 & gen_iii_clinic_exam_2$g3b0091 == 0]
gen_iii_subj_pheno$regular_smokers <- 1
gen_iii_subj_pheno$regular_smokers[gen_iii_subj_pheno$shareid %in% not_regular_smokers] <- 0
smoking_quiters <- gen_iii_clinic_exam_2$shareid[gen_iii_clinic_exam_2$g3b0097 > 0 & gen_iii_clinic_exam_2$g3b0091 == 1]
gen_iii_subj_pheno$smoking_quiters <- ifelse(gen_iii_subj_pheno$shareid %in% smoking_quiters, 1, 0)
gen_iii_subj_pheno <- merge(gen_iii_subj_pheno, gen_iii_clinic_exam_2[, c("shareid", "g3b0097")], by = "shareid")
colnames(gen_iii_subj_pheno)[colnames(gen_iii_subj_pheno) == "g3b0097"] <- "smoking_quitting_age"
###############################################################


###############################################################
gen_iii_labs <- read_process_phenotype_file("/local-scratch/projects/candx/guangbo/Framingham_data_part1/PhenoGenotypeFiles/RootStudyConsentSet_phs000007.Framingham.v31.p12.c1.HMB-IRB-MDS/PhenotypeFiles/phs000007.v31.pht002889.v4.p12.c1.l_fhslab_2011_m_0656s.HMB-IRB-MDS.txt", gen_iii_subj)
gen_iii_labs <- gen_iii_labs[, c("shareid", "HBA1C", "MO_PER", "LY_PER", "NE_PER", "MO_NUM", "LY_NUM", "NE_NUM", "CHOL", "TRIG", "GLUC", "CREAT", "HDL", "CRP", "UR_ALB", "UR_CREAT", "ALB", "AST", "ALT", "BILI", "CA", "GGT", "PHOS")]
colnames(gen_iii_labs) <- c("shareid", "HBA1C", "MO_PER", "LY_PER", "NE_PER", "MO_NUM", "LY_NUM", "NE_NUM", "total_cholestrol_plasma", "triglycerides_plasma", "glucose_plasma", "creatinine_serum", "HDL_cholestrol_plasma", "c_reactive_protein_serum", "albumin_urine", "creatinine_urine", "albumin_serum", "aspartate_aminotransferase_serum", "alanine_aminotransferase_serum", "total_bilirubin_serum", "calcium_serum", "gamma_glutamyltransferase_serum", "phosphorus_serum")
gen_iii_subj_pheno <- merge(gen_iii_subj_pheno, gen_iii_labs, by = "shareid", all = TRUE)
###############################################################


###############################################################
gen_iii_diabetes_info <- read_process_phenotype_file("/local-scratch/projects/candx/guangbo/Framingham_data_part1/PhenoGenotypeFiles/RootStudyConsentSet_phs000007.Framingham.v31.p12.c1.HMB-IRB-MDS/PhenotypeFiles/phs000007.v31.pht007775.v2.p12.c1.vr_diab_ex02_3b_0388s.HMB-IRB-MDS.txt", gen_iii_subj)
gen_iii_diabetes_info <- gen_iii_diabetes_info[, c("shareid", "CURR_DIAB2", "HX_DIAB2")]
###############################################################


###############################################################
dementia_info <- read_process_phenotype_file("/local-scratch/projects/candx/guangbo/Framingham_data_part1/PhenoGenotypeFiles/RootStudyConsentSet_phs000007.Framingham.v31.p12.c1.HMB-IRB-MDS/PhenotypeFiles/phs000007.v31.pht000691.v13.p12.c1.vr_demnp_2016_a_1088s.HMB-IRB-MDS.txt", gen_iii_subj)
dementia_info <- dementia_info[, c("shareid", "mri352")]
dementia_info <- na.omit(dementia_info)
###############################################################


###############################################################
subjects_with_cvd <- unique(cardio_events$shareid)
subjects_with_cancer <- unique(cancer_events$shareid)
subjects_with_fever <- unique(gen_iii_subj_pheno$shareid[gen_iii_subj_pheno$fever_or_infection_past_two_weeks == 1])
subjects_with_diabetes <- unique(gen_iii_diabetes_info$shareid[gen_iii_diabetes_info$CURR_DIAB2 == 1])
subjects_with_dementia <- unique(dementia_info$shareid[dementia_info$mri352 != 0])
subjects_with_clinical_diagnosis <- unique(gen_iii_clinic_exam_2_clinic_diagnosis$shareid[gen_iii_clinic_exam_2_clinic_diagnosis$disease_status != 0])
unhealthy_subjects <- unique(sort(c(subjects_with_diabetes, subjects_with_fever, subjects_with_cvd, subjects_with_dementia, subjects_with_clinical_diagnosis, subjects_with_cancer)))
healthy_subjects <- setdiff(gen_iii_subj, unhealthy_subjects)
###############################################################


###############################################################
pheno_module_score_matrix <- obtain_module_score_matrix(gen_iii_expr)
pheno_module_score_matrix$healthy[pheno_module_score_matrix$shareid %in% healthy_subjects] <- 1
pheno_module_score_matrix$cancer[pheno_module_score_matrix$shareid %in% subjects_with_cancer] <- 1
pheno_module_score_matrix$cardio[pheno_module_score_matrix$shareid %in% subjects_with_cvd] <- 1
pheno_module_score_matrix$diabetes[pheno_module_score_matrix$shareid %in% subjects_with_diabetes] <- 1
pheno_module_score_matrix$fever[pheno_module_score_matrix$shareid %in% subjects_with_fever] <- 1
pheno_module_score_matrix$dementia[pheno_module_score_matrix$shareid %in% subjects_with_dementia] <- 1
pheno_module_score_matrix$clinical_diagnosis_disease[pheno_module_score_matrix$shareid %in% subjects_with_clinical_diagnosis] <- 1
###############################################################


###############################################################
gen_iii_subj_pheno <- gen_iii_subj_pheno[, c("shareid", "sex", "age_at_exam_2", "weight", "height", "bmi", "regular_smokers", "smoking_quiters", "smoking_quitting_age", "death_date", "death_cause", "HBA1C", "MO_NUM", "LY_NUM", "NE_NUM", "MO_PER", "LY_PER", "NE_PER", "first_cardio_diagnosis_date", "first_cancer_diagnosis_date", "total_cholestrol_plasma", "triglycerides_plasma", "glucose_plasma", "creatinine_serum", "HDL_cholestrol_plasma", "c_reactive_protein_serum", "albumin_urine", "creatinine_urine", "albumin_serum", "aspartate_aminotransferase_serum", "alanine_aminotransferase_serum", "total_bilirubin_serum", "calcium_serum", "gamma_glutamyltransferase_serum", "phosphorus_serum", "c_reactive_protein_serum")]
colnames(gen_iii_subj_pheno) <- c("shareid", "sex", "age", "weight", "height", "bmi", "regular_smokers", "smoking_quiters", "smoking_quitting_age", "death_date", "death_cause", "A1C", "monocyte_count", "lymphocyte_count", "neutrophil_count", "monocyte_prop", "lymphocyte_prop", "neutrophil_prop", "first_cardio_diagnosis_date", "first_cancer_diagnosis_date", "total_cholestrol_plasma", "triglycerides_plasma", "glucose_plasma", "creatinine_serum", "HDL_cholestrol_plasma", "c_reactive_protein_serum", "albumin_urine", "creatinine_urine", "albumin_serum", "aspartate_aminotransferase_serum", "alanine_aminotransferase_serum", "total_bilirubin_serum", "calcium_serum", "gamma_glutamyltransferase_serum", "phosphorus_serum", "c_reactive_protein_serum")
pheno_module_score_matrix <- merge(pheno_module_score_matrix, gen_iii_subj_pheno, by = "shareid")
write_csv(pheno_module_score_matrix, file = "/local-scratch/projects/candx/guangbo/ananth_processesed_framingham_data/dhr_paper_re_processesed_2024_09_27/gen_iii_module_score_pheno_exon.csv")
###############################################################

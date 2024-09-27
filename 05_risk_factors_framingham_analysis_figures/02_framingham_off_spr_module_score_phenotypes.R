###############################################################
library(readr)
library(tidyverse)
source(file = "/labs/khatrilab/ananthg/detrimental_host_response/common_functions/framingham_score_phenotype_support_functions.R")
###############################################################


###############################################################
off_spr_expr <- read_csv("/local-scratch/projects/candx/guangbo/ananth_processesed_framingham_data/dhr_paper_re_processesed_2024_09_27/off_spr_expr_exon.csv")
off_spr_expr <- off_spr_expr[, -1]
off_spr_subj <- as.numeric(colnames(off_spr_expr)[2:ncol(off_spr_expr)])
###############################################################


###############################################################
off_spr_subj_pheno <- read_process_phenotype_file("/local-scratch/projects/candx/guangbo/Framingham_data_part1/PhenoGenotypeFiles/RootStudyConsentSet_phs000007.Framingham.v31.p12.c1.HMB-IRB-MDS/PhenotypeFiles/phs000007.v31.pht003099.v6.p12.c1.vr_dates_2019_a_1175s.HMB-IRB-MDS.txt", off_spr_subj)
off_spr_subj_pheno <- off_spr_subj_pheno[, c("shareid", "idtype", "sex", "age8", "date8", "att9", "age9", "att7", "age7")]
off_spr_subj_pheno$sex[off_spr_subj_pheno$sex == 1] <- "M"
off_spr_subj_pheno$sex[off_spr_subj_pheno$sex == 2] <- "F"
colnames(off_spr_subj_pheno) <- c("shareid", "idtype", "sex", "age_at_exam_8", "days_to_expr_from_exam_1", "attended_next_exam", "age_at_next_exam", "attended_previous_exam", "age_at_prervious_exam")
###############################################################


###############################################################
#cardio_events <- read_process_phenotype_file("/local-scratch/projects/candx/guangbo/Framingham_data_part1/PhenoGenotypeFiles/RootStudyConsentSet_phs000007.Framingham.v31.p12.c1.HMB-IRB-MDS/PhenotypeFiles/phs000007.v31.pht003316.v8.p12.c1.vr_survcvd_2017_a_1194s.HMB-IRB-MDS.txt", off_spr_subj)
#cardio_events <- cardio_events[cardio_events$chd == 1 | cardio_events$chd == 1 | cardio_events$chd == 1, ]
cardio_events <- read_process_phenotype_file("/local-scratch/projects/candx/guangbo/Framingham_data_part1/PhenoGenotypeFiles/RootStudyConsentSet_phs000007.Framingham.v31.p12.c1.HMB-IRB-MDS/PhenotypeFiles/phs000007.v31.pht000309.v14.p12.c1.vr_soe_2018_a_1311s.HMB-IRB-MDS.txt", off_spr_subj)
cardio_events <- cardio_events[cardio_events$EVENT %in% c(1, 2, 3, 4, 5, 6, 7, 10, 11, 12, 13, 16, 17, 18, 21, 22, 23, 24, 25, 26, 30, 39, 40, 41), ]
first_cardio_event <- cardio_events[, c("shareid", "DATE")] %>% group_by(shareid) %>% summarise(first_cardio = min(DATE))
off_spr_subj_pheno <- merge(off_spr_subj_pheno, first_cardio_event, by = "shareid", all = TRUE)
off_spr_subj_pheno$first_cardio_diagnosis_date <- off_spr_subj_pheno$first_cardio - off_spr_subj_pheno$days_to_expr_from_exam_1
###############################################################


###############################################################
death_events <- read_process_phenotype_file("/local-scratch/projects/candx/guangbo/Framingham_data_part1/PhenoGenotypeFiles/RootStudyConsentSet_phs000007.Framingham.v31.p12.c1.HMB-IRB-MDS/PhenotypeFiles/phs000007.v31.pht003317.v8.p12.c1.vr_survdth_2017_a_1192s.HMB-IRB-MDS.txt", off_spr_subj)
death_events <- death_events[!is.na(death_events$DATEDTH), ]
off_spr_subj_pheno <- merge(off_spr_subj_pheno, death_events[, c("shareid", "DATEDTH")], by = "shareid", all = TRUE)
off_spr_subj_pheno$death_date <- off_spr_subj_pheno$DATEDTH - off_spr_subj_pheno$days_to_expr_from_exam_1
off_spr_subj_pheno$death_cause <- NA
off_spr_subj_pheno$death_cause[!is.na(off_spr_subj_pheno$death_date)] <- "unknown"
off_spr_subj_pheno$death_cause[off_spr_subj_pheno$shareid %in% cardio_events$shareid[cardio_events$EVENT %in% c(3, 9, 21, 22, 23, 24, 25, 26)]] <- "cardio"
off_spr_subj_pheno$death_cause[off_spr_subj_pheno$shareid %in% cardio_events$shareid[cardio_events$EVENT %in% c(27)]] <- "cancer"
###############################################################


###############################################################
cancer_events <- read_process_phenotype_file("/local-scratch/projects/candx/guangbo/Framingham_data_part1/PhenoGenotypeFiles/RootStudyConsentSet_phs000007.Framingham.v31.p12.c1.HMB-IRB-MDS/PhenotypeFiles/phs000007.v31.pht000039.v12.p12.c1.vr_cancer_2013_a_0018s.HMB-IRB-MDS.txt", off_spr_subj)
first_cancer_event <- cancer_events[, c("shareid", "d_date")] %>% group_by(shareid) %>% summarise(first_cancer = min(d_date))
off_spr_subj_pheno <- merge(off_spr_subj_pheno, first_cancer_event, by = "shareid", all = TRUE)
off_spr_subj_pheno$first_cancer_diagnosis_date <- off_spr_subj_pheno$first_cancer - off_spr_subj_pheno$days_to_expr_from_exam_1
###############################################################


###############################################################
off_spr_clinic_exam_8 <- read_process_phenotype_file("/local-scratch/projects/candx/guangbo/Framingham_data_part1/PhenoGenotypeFiles/RootStudyConsentSet_phs000007.Framingham.v31.p12.c1.HMB-IRB-MDS/PhenotypeFiles/phs000007.v31.pht000747.v6.p12.c1.ex1_8s.HMB-IRB-MDS.txt", off_spr_subj)
off_spr_clinic_exam_8_weight_height <- off_spr_clinic_exam_8[, c("shareid", "H393", "H399")]
colnames(off_spr_clinic_exam_8_weight_height) <- c("shareid", "weight", "height")
off_spr_clinic_exam_8_weight_height$bmi <- 703 * off_spr_clinic_exam_8_weight_height$weight / (off_spr_clinic_exam_8_weight_height$height ^ 2)
off_spr_subj_pheno <- merge(off_spr_subj_pheno, off_spr_clinic_exam_8_weight_height, by = "shareid", all = TRUE)
off_spr_clinic_exam_8_clinic_diagnosis <- off_spr_clinic_exam_8[, c("shareid", "H339", "H340", "H341", "H342", "H344", "H345", "H347", "H348", "H349", "H354", "H357", "H358", "H359", "H360", "H363", "H364", "H367", "H368", "H369", "H370", "H374", "H375", "H378", "H379")]
off_spr_clinic_exam_8_clinic_diagnosis[is.na(off_spr_clinic_exam_8_clinic_diagnosis)] <- 3
off_spr_clinic_exam_8_clinic_diagnosis$disease_status <- rowSums(off_spr_clinic_exam_8_clinic_diagnosis[, 2:ncol(off_spr_clinic_exam_8_clinic_diagnosis)])
off_spr_clinic_exam_8_fever <- off_spr_clinic_exam_8[, c("shareid", "H003", "H004", "H005", "H006", "H007")]
colnames(off_spr_clinic_exam_8_fever) <- c("shareid", "hospitalization_since_exam_7", "ER_visit_since_exam_7", "day_surgery", "major_illness_visit_to_doctor", "fever_or_infection_past_two_weeks")
off_spr_subj_pheno <- merge(off_spr_subj_pheno, off_spr_clinic_exam_8_fever, by = "shareid", all = TRUE)
not_regular_smokers <- off_spr_clinic_exam_8$shareid[off_spr_clinic_exam_8$H065 == 0 & off_spr_clinic_exam_8$H060 == 0]
off_spr_subj_pheno$regular_smokers <- 1
off_spr_subj_pheno$regular_smokers[off_spr_subj_pheno$shareid %in% not_regular_smokers] <- 0
smoking_quiters <- off_spr_clinic_exam_8$shareid[off_spr_clinic_exam_8$H066 > 0 & off_spr_clinic_exam_8$H060 == 1]
off_spr_subj_pheno$smoking_quiters <- ifelse(off_spr_subj_pheno$shareid %in% smoking_quiters, 1, 0)
off_spr_subj_pheno <- merge(off_spr_subj_pheno, off_spr_clinic_exam_8[, c("shareid", "H066")], by = "shareid")
colnames(off_spr_subj_pheno)[colnames(off_spr_subj_pheno) == "H066"] <- "smoking_quitting_age"
###############################################################


###############################################################
off_spr_labs <- read_process_phenotype_file("/local-scratch/projects/candx/guangbo/Framingham_data_part1/PhenoGenotypeFiles/RootStudyConsentSet_phs000007.Framingham.v31.p12.c1.HMB-IRB-MDS/PhenotypeFiles/phs000007.v31.pht000742.v6.p12.c1.fhslab1_8s.HMB-IRB-MDS.txt", off_spr_subj)
off_spr_labs <- off_spr_labs[, c("shareid", "A1C")]
off_spr_subj_pheno <- merge(off_spr_subj_pheno, off_spr_labs, by = "shareid", all = TRUE)
###############################################################


###############################################################
off_spr_diabetes_info <- read_process_phenotype_file("/local-scratch/projects/candx/guangbo/Framingham_data_part1/PhenoGenotypeFiles/RootStudyConsentSet_phs000007.Framingham.v31.p12.c1.HMB-IRB-MDS/PhenotypeFiles/phs000007.v31.pht000041.v7.p12.c1.vr_diab_ex09_1_1002s.HMB-IRB-MDS.txt", off_spr_subj)
off_spr_diabetes_info <- off_spr_diabetes_info[, c("shareid", "CURR_DIAB8", "HX_DIAB8")]
###############################################################


###############################################################
dementia_info <- read_process_phenotype_file("/local-scratch/projects/candx/guangbo/Framingham_data_part1/PhenoGenotypeFiles/RootStudyConsentSet_phs000007.Framingham.v31.p12.c1.HMB-IRB-MDS/PhenotypeFiles/phs000007.v31.pht000691.v13.p12.c1.vr_demnp_2016_a_1088s.HMB-IRB-MDS.txt", off_spr_subj)
dementia_info <- dementia_info[, c("shareid", "mri352")]
dementia_info <- na.omit(dementia_info)
###############################################################


###############################################################
subjects_with_cvd <- unique(cardio_events$shareid)
subjects_with_cancer <- unique(cancer_events$shareid)
subjects_with_fever <- unique(off_spr_subj_pheno$shareid[off_spr_subj_pheno$fever_or_infection_past_two_weeks == 1])
subjects_with_diabetes <- unique(off_spr_diabetes_info$shareid[off_spr_diabetes_info$CURR_DIAB8 == 1])
subjects_with_dementia <- unique(dementia_info$shareid[dementia_info$mri352 != 0])
subjects_with_clinical_diagnosis <- unique(off_spr_clinic_exam_8_clinic_diagnosis$shareid[off_spr_clinic_exam_8_clinic_diagnosis$disease_status != 0])
unhealthy_subjects <- unique(sort(c(subjects_with_diabetes, subjects_with_fever, subjects_with_cvd, subjects_with_dementia, subjects_with_clinical_diagnosis, subjects_with_cancer)))
healthy_subjects <- setdiff(off_spr_subj, unhealthy_subjects)
###############################################################


###############################################################
pheno_module_score_matrix <- obtain_module_score_matrix(off_spr_expr)
pheno_module_score_matrix$healthy[pheno_module_score_matrix$shareid %in% healthy_subjects] <- 1
pheno_module_score_matrix$cancer[pheno_module_score_matrix$shareid %in% subjects_with_cancer] <- 1
pheno_module_score_matrix$cardio[pheno_module_score_matrix$shareid %in% subjects_with_cvd] <- 1
pheno_module_score_matrix$diabetes[pheno_module_score_matrix$shareid %in% subjects_with_diabetes] <- 1
pheno_module_score_matrix$fever[pheno_module_score_matrix$shareid %in% subjects_with_fever] <- 1
pheno_module_score_matrix$dementia[pheno_module_score_matrix$shareid %in% subjects_with_dementia] <- 1
pheno_module_score_matrix$clinical_diagnosis_disease[pheno_module_score_matrix$shareid %in% subjects_with_clinical_diagnosis] <- 1
###############################################################


###############################################################
off_spr_subj_pheno <- off_spr_subj_pheno[, c("shareid", "sex", "age_at_exam_8", "weight", "height", "bmi", "regular_smokers", "smoking_quiters", "smoking_quitting_age", "death_date", "death_cause", "A1C", "first_cardio_diagnosis_date", "first_cancer_diagnosis_date")]
colnames(off_spr_subj_pheno) <- c("shareid", "sex", "age", "weight", "height", "bmi", "regular_smokers", "smoking_quiters", "smoking_quitting_age", "death_date", "death_cause", "A1C", "first_cardio_diagnosis_date", "first_cancer_diagnosis_date")
pheno_module_score_matrix <- merge(pheno_module_score_matrix, off_spr_subj_pheno, by = "shareid")
write_csv(pheno_module_score_matrix, file = "/local-scratch/projects/candx/guangbo/ananth_processesed_framingham_data/dhr_paper_re_processesed_2024_09_27/off_spr_module_score_pheno_exon.csv")
###############################################################

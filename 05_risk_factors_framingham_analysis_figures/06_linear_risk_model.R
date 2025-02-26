###############################################################
library(readr)
###############################################################


###############################################################
obtain_linear_model_coeff <- function(cohort_file_name)
{
  combined_score_pheno <- read_csv(paste0("/local-scratch/projects/candx/guangbo/ananth_processesed_framingham_data/dhr_paper_re_processesed_2024_09_27/", cohort_file_name))
  combined_score_pheno$disease <- as.factor(1 - combined_score_pheno$healthy)
  mod_som_lm <- lm(formula = scale(som_score) ~ age + as.factor(sex) + bmi + disease + as.factor(regular_smokers), data = combined_score_pheno)
  mod_detrimental_lm <- lm(formula = scale(detrimental_score) ~ age + as.factor(sex) + bmi + disease + as.factor(regular_smokers), data = combined_score_pheno)
  mod_protective_lm <- lm(formula = scale(protective_score) ~ age + as.factor(sex) + bmi + disease + as.factor(regular_smokers), data = combined_score_pheno)
  mod_1_lm <- lm(formula = scale(module_1_score) ~ age + as.factor(sex) + bmi + disease + as.factor(regular_smokers), data = combined_score_pheno)
  mod_2_lm <- lm(formula = scale(module_2_score) ~ age + as.factor(sex) + bmi + disease + as.factor(regular_smokers), data = combined_score_pheno)
  mod_3_lm <- lm(formula = scale(module_3_score) ~ age + as.factor(sex) + bmi + disease + as.factor(regular_smokers), data = combined_score_pheno)
  mod_4_lm <- lm(formula = scale(module_4_score) ~ age + as.factor(sex) + bmi + disease + as.factor(regular_smokers), data = combined_score_pheno)
  scores_list <- c("SoM", "Detrimental", "Protective", "Module 1", "Module 2", "Module 3", "Module 4")
  confounder_list <- c("Intercept", "Age", "Sex", "BMI", "Disease", "Smoking")
  regression_models_list <- list(mod_som_lm, mod_detrimental_lm, mod_protective_lm, mod_1_lm, mod_2_lm, mod_3_lm, mod_4_lm)
  regression_coeff <- matrix(, nrow = length(scores_list), ncol = length(confounder_list))
  regression_p_value <- matrix(, nrow = length(scores_list), ncol = length(confounder_list))
  regression_r2 <- matrix(, nrow = length(scores_list), ncol = 1)
  for(i in 1:length(regression_models_list))
  {
    model_summary <- summary(regression_models_list[[i]])
    regression_coeff[i, ] <- model_summary$coefficients[, 1]
    regression_p_value[i, ] <- model_summary$coefficients[, 4]
    regression_r2[i, 1] <- model_summary$adj.r.squared
    p_val <- pf(model_summary$fstatistic[1], model_summary$fstatistic[2], model_summary$fstatistic[3], lower.tail = F)
    attributes(p_val) <- NULL
    print(p_val)
    if(p_val > 0.5)
    {
      print(i)
      print("Non significant model")
    }
  }
  rownames(regression_p_value) <- scores_list
  rownames(regression_coeff) <- scores_list
  rownames(regression_r2) <- scores_list
  colnames(regression_p_value) <- confounder_list
  colnames(regression_coeff) <- confounder_list
  colnames(regression_r2) <- "R2"
  return(list(regression_coeff = regression_coeff, regression_p_value = regression_p_value, regression_r2 = regression_r2))
}
###############################################################


###############################################################
combined_coconut_regression <- obtain_linear_model_coeff("framingham_conormalized_by_healthy_matrix.csv")
saveRDS(combined_coconut_regression, file = "/labs/khatrilab/ananthg/detrimental_host_response/data_2024_09_27/framingham_conormalized_linear_models_summary.rds")
###############################################################


###############################################################
obtain_linear_model_coeff_ratio <- function(cohort_file_name)
{
  combined_score_pheno <- read_csv(paste0("/local-scratch/projects/candx/guangbo/ananth_processesed_framingham_data/dhr_paper_re_processesed_2024_09_27/", cohort_file_name))
  combined_score_pheno$disease <- as.factor(1 - combined_score_pheno$healthy)
  mod_som_lm <- lm(formula = scale(som_score_ratio) ~ age + as.factor(sex) + bmi + disease + as.factor(regular_smokers), data = combined_score_pheno)
  mod_detrimental_lm <- lm(formula = scale(detrimental_score) ~ age + as.factor(sex) + bmi + disease + as.factor(regular_smokers), data = combined_score_pheno)
  mod_protective_lm <- lm(formula = scale(protective_score) ~ age + as.factor(sex) + bmi + disease + as.factor(regular_smokers), data = combined_score_pheno)
  mod_1_lm <- lm(formula = scale(module_1_score) ~ age + as.factor(sex) + bmi + disease + as.factor(regular_smokers), data = combined_score_pheno)
  mod_2_lm <- lm(formula = scale(module_2_score) ~ age + as.factor(sex) + bmi + disease + as.factor(regular_smokers), data = combined_score_pheno)
  mod_3_lm <- lm(formula = scale(module_3_score) ~ age + as.factor(sex) + bmi + disease + as.factor(regular_smokers), data = combined_score_pheno)
  mod_4_lm <- lm(formula = scale(module_4_score) ~ age + as.factor(sex) + bmi + disease + as.factor(regular_smokers), data = combined_score_pheno)
  scores_list <- c("SoM", "Detrimental", "Protective", "Module 1", "Module 2", "Module 3", "Module 4")
  confounder_list <- c("Intercept", "Age", "Sex", "BMI", "Disease", "Smoking")
  regression_models_list <- list(mod_som_lm, mod_detrimental_lm, mod_protective_lm, mod_1_lm, mod_2_lm, mod_3_lm, mod_4_lm)
  regression_coeff <- matrix(, nrow = length(scores_list), ncol = length(confounder_list))
  regression_p_value <- matrix(, nrow = length(scores_list), ncol = length(confounder_list))
  regression_r2 <- matrix(, nrow = length(scores_list), ncol = 1)
  for(i in 1:length(regression_models_list))
  {
    model_summary <- summary(regression_models_list[[i]])
    regression_coeff[i, ] <- model_summary$coefficients[, 1]
    regression_p_value[i, ] <- model_summary$coefficients[, 4]
    regression_r2[i, 1] <- model_summary$adj.r.squared
    p_val <- pf(model_summary$fstatistic[1], model_summary$fstatistic[2], model_summary$fstatistic[3], lower.tail = F)
    attributes(p_val) <- NULL
    print(p_val)
    if(p_val > 0.5)
    {
      print(i)
      print("Non significant model")
    }
  }
  rownames(regression_p_value) <- scores_list
  rownames(regression_coeff) <- scores_list
  rownames(regression_r2) <- scores_list
  colnames(regression_p_value) <- confounder_list
  colnames(regression_coeff) <- confounder_list
  colnames(regression_r2) <- "R2"
  return(list(regression_coeff = regression_coeff, regression_p_value = regression_p_value, regression_r2 = regression_r2))
}
###############################################################


###############################################################
combined_coconut_regression <- obtain_linear_model_coeff_ratio("framingham_conormalized_by_healthy_matrix.csv")
saveRDS(combined_coconut_regression, file = "/labs/khatrilab/ananthg/detrimental_host_response/data_2024_09_27/framingham_conormalized_linear_models_summary_som_ratio.rds")
###############################################################

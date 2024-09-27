###############################################################
library(readr)
library(plyr)
source("/labs/khatrilab/ananthg/detrimental_host_response/common_functions/plot_themes.R")
source("/labs/khatrilab/ananthg/detrimental_host_response/common_functions/JT_test_function.R")
###############################################################


###############################################################
tile_size <- unit(16, "pt")
ht_opt("simple_anno_size" = tile_size)
###############################################################


###############################################################
coconut_combined <- read_csv("/local-scratch/projects/candx/guangbo/ananth_processesed_framingham_data/dhr_paper_re_processesed_2024_09_27/framingham_conormalized_by_healthy_matrix.csv")
coconut_combined$disease_status <- as.factor(1 - coconut_combined$healthy)
coconut_combined$years_no_smoking <- coconut_combined$age - coconut_combined$smoking_quitting_age
###############################################################


###############################################################
case_control_summary <- readRDS("/labs/khatrilab/ananthg/detrimental_host_response/data_2024_09_27/case_control_analysis_summary.rds")
eff_size_matrix <- case_control_summary$eff_size_matrix[, c("SoM", "Detrimental", "Protective")]
fdr_matrix <- case_control_summary$fdr_matrix[, c("SoM", "Detrimental", "Protective")]
data_summary <- case_control_summary$combined_summary
eff_size_plot_title <- paste0("Case-control meta-analyses (", sum(data_summary$Dataset), " datasets, ", sum(data_summary$Cases) + sum(data_summary$Controls), " samples)")
min_eff_size <- -0.8
max_eff_size <- 0.8
eff_size_matrix_colour_pts <- c(min_eff_size, min_eff_size/2, 0, max_eff_size/2, max_eff_size)
eff_size_cols <- colorRamp2(eff_size_matrix_colour_pts, c("#762A83", "#C2A5CF", "#F7F7F7", "#ACD39E", "#1B7837"))
eff_size_matrix_legend_pts <- c(min_eff_size, max_eff_size)
eff_size_matrix_legend_labels <- c(round(min_eff_size, 2), round(max_eff_size, 2))
eff_size_legend <- Legend(col_fun = eff_size_cols,
  at = eff_size_matrix_legend_pts,
  labels = eff_size_matrix_legend_labels,
  title = "Effect size", title_position = "topcenter",
  grid_height = unit(base_text_size / 3, "pt"),
  legend_width = 4 * unit(base_text_size, "pt"),
  background = "transparent",
  title_gp = text_setting_small,
  title_gap = gap_size, gap = gap_size,
  labels_gp = text_setting_sub_small,
  direction = "horizontal")
summary_annotation <- rowAnnotation(Datasets = anno_text(data_summary$Dataset, location = unit(0.95, "npc"), just = "right", gp = text_setting_sub_small, width = tile_size),
  Cases = anno_text(data_summary$Cases, location = unit(0.95, "npc"), just = "right", gp = text_setting_sub_small, width = 1.1 * tile_size),
  Controls = anno_text(data_summary$Controls, location = unit(0.95, "npc"), just = "right", gp = text_setting_sub_small, width = 1.1 * tile_size))
ht_mp_eff_size <- Heatmap(eff_size_matrix,
  col = eff_size_cols, na_col = "#FFFFFF",
  cell_fun = function(j, i, x, y, width, height, fill)
    {
      if(fdr_matrix[i, j] <= 0.1)
      {
        grid.text("*", x, y, vjust = 0.8, gp = gpar(fontsize = base_text_size * 1.5, fontfamily = base_font_family, col = black_line_colour))
      }
    },
  rect_gp = dend_lines_gp,
  column_names_rot = 45,
  width = ncol(eff_size_matrix) * tile_size,
  name = "eff_size",
  row_names_side = "left",
  show_heatmap_legend = FALSE,
  cluster_rows = FALSE, cluster_columns = FALSE,
  gap = gap_size,
  left_annotation = summary_annotation)
###############################################################


###############################################################
clean_up_framingham_matrix <- function(wide_tibble, row_order, old_column_order, column_order)
{
  wide_tibble <- wide_tibble[, colnames(wide_tibble) != "Intercept"]
  old_column_order <- c("Age", "Sex", "BMI", "Disease", "Smoking")
  column_order <- c(paste0("Age (median (IQR) = ", round(median(coconut_combined$age)), " (", round(IQR(coconut_combined$age)), "))"), paste0("Sex (female = ", sum(coconut_combined$sex == "F"), ", male = ", sum(coconut_combined$sex == "M"), ")"), paste0("BMI (median (IQR) = ", round(median(coconut_combined$bmi, na.rm = TRUE), 2), " (", round(IQR(coconut_combined$bmi, na.rm = TRUE), 2), "))"), paste0("Disease (healthy = ", sum(coconut_combined$healthy == 1), ", not = ", sum(coconut_combined$healthy == 0), ")"), paste0("Smoking (smokers = ", sum(coconut_combined$regular_smokers == 1), ", not = ", sum(coconut_combined$regular_smokers == 0), ")"))
  expanded_condition_names <- mapvalues(colnames(wide_tibble), old_column_order, column_order)
  colnames(wide_tibble) <- expanded_condition_names
  wide_tibble <- t(as.matrix(wide_tibble[c("SoM", "Detrimental", "Protective"), column_order]))
  return(wide_tibble)
}
###############################################################


###############################################################
combined_coconut_regression <- readRDS("/labs/khatrilab/ananthg/detrimental_host_response/data_2024_09_27/framingham_conormalized_linear_models_summary.rds")
regression_coeff <- clean_up_framingham_matrix(combined_coconut_regression$regression_coeff)
regression_p_value <- clean_up_framingham_matrix(combined_coconut_regression$regression_p_value)
min_regression <- -0.35
max_regression <- 0.35
regression_matrix_colour_pts <- c(min_regression, min_regression/2, 0, max_regression/2, max_regression)
regression_cols <- colorRamp2(regression_matrix_colour_pts, c("#125A56", "#60BCE9", "#ECEADA", "#FD9A44", "#A01813"))
regression_matrix_legend_pts <- c(min_regression, max_regression)
regression_matrix_legend_labels <- c(round(min_regression, 2), round(max_regression, 2))
regression_legend <- Legend(col_fun = regression_cols,
  at = regression_matrix_legend_pts,
  labels = regression_matrix_legend_labels,
  title = "Regression\ncoefficient", title_position = "topcenter",
  grid_height = unit(base_text_size / 3, "pt"),
  legend_width = 4 * unit(base_text_size, "pt"),
  background = "transparent",
  title_gp = text_setting_small,
  title_gap = gap_size, gap = gap_size,
  labels_gp = text_setting_sub_small,
  direction = "horizontal")
ht_mp_regression <- Heatmap(regression_coeff,
  col = regression_cols, na_col = "#FFFFFF",
  cell_fun = function(j, i, x, y, width, height, fill)
    {
      if(regression_p_value[i, j] < 0.05)
      {
        grid.text("*", x, y, vjust = 0.8, gp = gpar(fontsize = base_text_size * 1.5, fontfamily = base_font_family, col = black_line_colour))
      }
    },
  rect_gp = dend_lines_gp,
  column_names_rot = 45,
  width = ncol(regression_coeff) * tile_size,
  name = "regression",
  row_names_side = "left",
  show_heatmap_legend = FALSE,
  cluster_rows = FALSE, cluster_columns = FALSE,
  gap = gap_size)
###############################################################


###############################################################
compute_risk_factors <- function(framingham_score_pheno_matrix)
{
  framingham_score_pheno_matrix$disease_risk_factors <- rowSums(framingham_score_pheno_matrix[, c("cancer", "cardio", "diabetes", "fever", "dementia", "clinical_diagnosis_disease")])
  framingham_score_pheno_matrix$age_risk <- as.integer(framingham_score_pheno_matrix$age >= 60)
  framingham_score_pheno_matrix$sex_risk <- as.integer(framingham_score_pheno_matrix$sex == "M")
  framingham_score_pheno_matrix$bmi_risk <- as.integer(framingham_score_pheno_matrix$bmi >= 30)
  framingham_score_pheno_matrix$demographic_risk_factors <- rowSums(framingham_score_pheno_matrix[, c("age_risk", "sex_risk", "bmi_risk", "regular_smokers")])
  framingham_score_pheno_matrix$risk_factors <- as.integer(rowSums(framingham_score_pheno_matrix[, c("disease_risk_factors", "demographic_risk_factors")]))
  framingham_score_pheno_matrix <- na.omit(framingham_score_pheno_matrix[, c("som_score", "detrimental_score", "protective_score", "module_1_score", "module_2_score", "module_3_score", "module_4_score", "tsang_up_score", "tsang_down_score", "tsang_ihm_score", "demographic_risk_factors", "disease_risk_factors", "risk_factors")])
  framingham_score_pheno_matrix$risk_factors_grouped <- "1 - 3"
  framingham_score_pheno_matrix$risk_factors_grouped[framingham_score_pheno_matrix$risk_factors == 0] <- "0"
  framingham_score_pheno_matrix$risk_factors_grouped[framingham_score_pheno_matrix$risk_factors > 3] <- "\u2265 4"
  framingham_score_pheno_matrix$risk_factors_grouped <- ordered(framingham_score_pheno_matrix$risk_factors_grouped, levels = c("0", "1 - 3", "\u2265 4"))
  return(framingham_score_pheno_matrix)
}
coconut_combined_risk_factors <- compute_risk_factors(coconut_combined)
###############################################################


###############################################################
plot_each_score_risk_factor <- function(score_pheno_frame, score_name, score_name_plot, print_NS = FALSE)
{
  risk_factor_subjects <- table(score_pheno_frame$risk_factors_grouped)
  score_pheno_frame[[score_name]] <- scale(score_pheno_frame[[score_name]])
  limits <- quantile(score_pheno_frame[[score_name]], c(0.025, 0.975))
  score_pheno_frame[[score_name]][score_pheno_frame[[score_name]] < limits[1]] <- NA
  score_pheno_frame[[score_name]][score_pheno_frame[[score_name]] > limits[2]] <- NA
  score_pheno_frame$scaled_score <- score_pheno_frame[[score_name]]
  test_significance <- score_pheno_frame %>%
    wilcox_test(scaled_score ~ risk_factors_grouped, comparisons = list(c("0", "1 - 3"), c("1 - 3", "\u2265 4"))) %>%
    add_xy_position(x = "risk_factors_grouped")
  test_significance$y.position <- 1.03 * test_significance$y.position
  test_significance$xmin <- 1.02 * test_significance$xmin
  test_significance$xmax <- 0.98 * test_significance$xmax
  if(!print_NS)
  {
    test_significance$p_print <- as.character(ifelse(test_significance$p < 0.01, format(test_significance$p, digits = 1, scientific = TRUE), round(test_significance$p, 2)))
  }
  else
  {
    test_significance$p_print <- as.character(ifelse(test_significance$p < 0.01, format(test_significance$p, digits = 1, scientific = TRUE), ifelse(test_significance$p >= 0.055, "NS", round(test_significance$p, 2))))
  }
  score_summary <- ggplot(score_pheno_frame, aes(x = risk_factors_grouped, y = scaled_score)) +
    ylab(score_name_plot) +
    geom_boxplot(aes(fill = risk_factors_grouped), lwd = line_size_main_mm, width = 0.6, outlier.shape = NA, colour = black_text_colour) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1)), breaks = c(-1, 0, 1, 2)) +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank(),
      legend.position = "bottom", legend.direction = "horizontal", legend.margin = margin(1, 1, 1, 1), legend.box.margin = margin(-5, 5, 3, 3)) +
    stat_pvalue_manual(test_significance, label = "p_print",
      step.increase = 0.1, tip.length = 0,
      size = text_size_labels, colour = black_text_colour, bracket.size = line_size_text_repel_mm) +
    scale_fill_manual(values = c("0" = "#ECE6F6", "1 - 3" = "#D1C4E9", "\u2265 4" = "#9474CC"), labels = c("0" = paste0("0 (n = ", risk_factor_subjects["0"], ")"), "1 - 3" = paste0("1 - 3 (n = ", risk_factor_subjects["1 - 3"], ")"), "\u2265 4" = paste0("\u2265 4 (n = ", risk_factor_subjects["\u2265 4"], ")"))) +
    labs(fill = "Risk factors") + guides(colour = guide_legend(nrow = 1, override.aes = list(size = text_size_labels * 6 / 9)))
  return(score_summary)
}
risk_factors_detrimental_plot <- plot_each_score_risk_factor(coconut_combined_risk_factors, "detrimental_score", "Detrimental (z score)")
risk_factors_protective_plot <- plot_each_score_risk_factor(coconut_combined_risk_factors, "protective_score", "Protective (z score)") + theme(legend.position = "none")
risk_factors_som_plot <- plot_each_score_risk_factor(coconut_combined_risk_factors, "som_score", "SoM (z score)") + theme(legend.position = "none")
risk_factors_combined <- plot_grid(risk_factors_som_plot, risk_factors_detrimental_plot, risk_factors_protective_plot, nrow = 1, align = "hv", axis = "tblr")
###############################################################


###############################################################
smokers_only <- coconut_combined[coconut_combined$regular_smokers == 1, c("age", "sex", "bmi", "disease_status", "years_no_smoking", "module_1_score", "module_2_score", "module_3_score", "module_4_score", "detrimental_score", "protective_score", "som_score", "tsang_up_score", "tsang_down_score", "tsang_ihm_score", "smoking_quiters")]
smokers_only$years_no_smoking_discrete <- ifelse(smokers_only$smoking_quiters == 0, "Current smokers", ifelse(smokers_only$years_no_smoking >= 5, "Past smokers", "Recent smokers"))
smokers_only <- na.omit(smokers_only[, c("age", "sex", "bmi", "disease_status", "years_no_smoking_discrete", "module_1_score", "module_2_score", "module_3_score", "module_4_score", "detrimental_score", "protective_score", "som_score", "tsang_up_score", "tsang_down_score", "tsang_ihm_score", "smoking_quiters")])
smokers_only$years_no_smoking_discrete <- ordered(smokers_only$years_no_smoking_discrete, levels = c("Current smokers", "Recent smokers", "Past smokers"))
###############################################################


###############################################################
plot_each_score_smoking <- function(score_pheno_frame, score_name, score_name_plot, print_NS = FALSE)
{
  smokers_by_subject <- table(score_pheno_frame$years_no_smoking_discrete)
  score_pheno_frame[[score_name]] <- scale(score_pheno_frame[[score_name]])
  lm_formula <- as.formula(paste0(score_name, " ~ age + as.factor(sex) + bmi + disease_status"))
  linear_model <- lm(lm_formula, data = score_pheno_frame)
  score_pheno_frame$adjusted_score <- NA
  score_pheno_frame$adjusted_score[!is.na(names(linear_model$residuals))] <- linear_model$residuals
  group_colours <- rev(c("#F4A582", "#D6604D", "#B2182B"))
  test_significance <- score_pheno_frame %>%
    wilcox_test(adjusted_score ~ years_no_smoking_discrete, comparisons = list(c("Current smokers", "Recent smokers"), c("Recent smokers", "Past smokers"))) %>%
    add_xy_position(x = "years_no_smoking_discrete")
  test_significance$y.position <- 1.03 * test_significance$y.position
  test_significance$xmin <- 1.02 * test_significance$xmin
  test_significance$xmax <- 0.98 * test_significance$xmax
  if(!print_NS)
  {
    test_significance$p_print <- as.character(ifelse(test_significance$p < 0.01, format(test_significance$p, digits = 1, scientific = TRUE), round(test_significance$p, 2)))
  }
  else
  {
    test_significance$p_print <- as.character(ifelse(test_significance$p < 0.01, format(test_significance$p, digits = 1, scientific = TRUE), ifelse(test_significance$p >= 0.055, "NS", round(test_significance$p, 2))))
  }
  score_summary <- ggplot(score_pheno_frame, aes(x = years_no_smoking_discrete, y = adjusted_score, colour = years_no_smoking_discrete)) +
    ylab(score_name_plot) +
    geom_beeswarm(cex = cex_val * 1.1, size = point_size) +
    geom_boxplot(fill = NA, lwd = line_size_main_mm, width = 0.4, outlier.shape = NA, colour = black_text_colour) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)), breaks = c(-2, 0, 2, 4)) +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank(),
      legend.position = "bottom", legend.direction = "horizontal", legend.margin = margin(1, 1, 1, 1), legend.box.margin = margin(-8, 5, 3, 3)) +
    stat_pvalue_manual(test_significance, label = "p_print",
      step.increase = 0.075, tip.length = 0,
      size = text_size_labels, colour = black_text_colour, bracket.size = line_size_text_repel_mm) +
    scale_colour_manual(values = group_colours, labels = c(`Current smokers` = paste0("Current (n = ", smokers_by_subject["Current smokers"], ")"), `Recent smokers` = paste0("Recent (n = ", smokers_by_subject["Recent smokers"], ")"), `Past smokers` = paste0("Past (n = ", smokers_by_subject["Past smokers"], ")"))) +
    labs(colour = NULL) + guides(colour = guide_legend(nrow = 1, override.aes = list(size = text_size_labels * 6 / 9)))
  return(score_summary)
}
smoking_detrimental_plot <- plot_each_score_smoking(smokers_only, "detrimental_score", "Adjusted detrimental score")
smoking_protective_plot <- plot_each_score_smoking(smokers_only, "protective_score", "Adjusted protective score")
smoking_som_plot <- plot_each_score_smoking(smokers_only, "som_score", "Adjusted SoM score")
smoking_legend <- as_ggplot(ggpubr::get_legend(smoking_som_plot))
smoking_combined <- plot_grid(smoking_som_plot + theme(legend.position = "none"), smoking_detrimental_plot + theme(legend.position = "none"), smoking_protective_plot + theme(legend.position = "none"), nrow = 1, align = "hv", axis = "tblr")
smoking_tsang_ihm_plot <- plot_each_score_smoking(smokers_only, "tsang_ihm_score", "Adjusted IHM score", print_NS = TRUE) + theme(legend.position = "right", legend.title.position = "top") + guides(colour = guide_legend(nrow = 3, override.aes = list(size = text_size_labels * 6 / 9)))
###############################################################


###############################################################
diabetes_only <- coconut_combined[coconut_combined$diabetes == 1, c("age", "sex", "bmi", "disease_status", "regular_smokers", "module_1_score", "module_2_score", "module_3_score", "module_4_score", "detrimental_score", "protective_score", "som_score", "tsang_up_score", "tsang_down_score", "tsang_ihm_score", "A1C")]
diabetes_only$a1c_discrete <- ifelse(diabetes_only$A1C < 7, "A1C < 7", ifelse(diabetes_only$A1C > 10, "A1C > 10", "7 <= A1C <= 10"))
diabetes_only <- na.omit(diabetes_only)
diabetes_only$a1c_discrete <- ordered(diabetes_only$a1c_discrete, levels = c("A1C < 7", "7 <= A1C <= 10", "A1C > 10"))
###############################################################


###############################################################
plot_each_score_a1c <- function(score_pheno_frame, score_name, score_name_plot, print_NS = FALSE)
{
  score_pheno_frame[[score_name]] <- scale(score_pheno_frame[[score_name]])
  lm_formula <- as.formula(paste0(score_name, " ~ age + as.factor(sex) + bmi + as.factor(regular_smokers)"))
  linear_model <- lm(lm_formula, data = score_pheno_frame)
  score_pheno_frame$adjusted_score <- NA
  score_pheno_frame$adjusted_score[!is.na(names(linear_model$residuals))] <- linear_model$residuals
  group_colours <- c("#ACD39E", "#5AAE61", "#1B7837")
  test_significance <- score_pheno_frame %>%
    wilcox_test(adjusted_score ~ a1c_discrete, comparisons = list(c("A1C < 7", "7 <= A1C <= 10"), c("7 <= A1C <= 10", "A1C > 10"))) %>%
    add_xy_position(x = "a1c_discrete")
  test_significance$y.position <- 1.03 * test_significance$y.position
  test_significance$xmin <- 1.02 * test_significance$xmin
  test_significance$xmax <- 0.98 * test_significance$xmax
  if(!print_NS)
  {
    test_significance$p_print <- as.character(ifelse(test_significance$p < 0.01, format(test_significance$p, digits = 1, scientific = TRUE), round(test_significance$p, 2)))
  }
  else
  {
    test_significance$p_print <- as.character(ifelse(test_significance$p < 0.01, format(test_significance$p, digits = 1, scientific = TRUE), ifelse(test_significance$p >= 0.055, "NS", round(test_significance$p, 2))))
  }
  score_summary <- ggplot(score_pheno_frame, aes(x = a1c_discrete, y = adjusted_score, colour = a1c_discrete)) +
    ylab(score_name_plot) +
    geom_beeswarm(cex = cex_val * 1.1, size = point_size) +
    geom_boxplot(fill = NA, lwd = line_size_main_mm, width = 0.4, outlier.shape = NA, colour = black_text_colour) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)), breaks = c(-2, 0, 2, 4)) +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank(),
      legend.position = "bottom", legend.direction = "horizontal", legend.margin = margin(1, 1, 1, 1), legend.box.margin = margin(-8, 5, 3, 3)) +
    stat_pvalue_manual(test_significance, label = "p_print",
      step.increase = 0.075, tip.length = 0,
      size = text_size_labels, colour = black_text_colour, bracket.size = line_size_text_repel_mm) +
    scale_colour_manual(values = group_colours, labels = c(`A1C < 7` = paste0("< 7 (n = ", sum(score_pheno_frame$a1c_discrete == "A1C < 7"), ")"), `7 <= A1C <= 10` = paste0("[7, 10] (n = ", sum(score_pheno_frame$a1c_discrete == "7 <= A1C <= 10"), ")"), `A1C > 10` = paste0("> 10 (n = ", sum(score_pheno_frame$a1c_discrete == "A1C > 10"), ")"))) +
    labs(colour = "A1C") + guides(colour = guide_legend(nrow = 1, override.aes = list(size = text_size_labels * 6 / 9)))
  return(score_summary)
}
a1c_detrimental_plot <- plot_each_score_a1c(diabetes_only, "detrimental_score", "Adjusted detrimental score")
a1c_protective_plot <- plot_each_score_a1c(diabetes_only, "protective_score", "Adjusted protective score")
a1c_som_plot <- plot_each_score_a1c(diabetes_only, "som_score", "Adjusted SoM score")
a1c_legend <- as_ggplot(ggpubr::get_legend(a1c_som_plot))
a1c_combined <- plot_grid(a1c_som_plot + theme(legend.position = "none"), a1c_detrimental_plot + theme(legend.position = "none"), a1c_protective_plot + theme(legend.position = "none"), nrow = 1, align = "hv", axis = "tblr")
a1c_tsang_ihm_plot <- plot_each_score_a1c(diabetes_only, "tsang_ihm_score", "Adjusted IHM score", print_NS = TRUE) + theme(legend.position = "right", legend.title.position = "top") + guides(colour = guide_legend(nrow = 3, override.aes = list(size = text_size_labels * 6 / 9)))
###############################################################


###############################################################
GSE88794_diet_restriction <- read_csv("/labs/khatrilab/ananthg/detrimental_host_response/data_2024_09_27/GSE88794_combined_score_pheno.csv")
fasting_ogtt_data_only <- GSE88794_diet_restriction[GSE88794_diet_restriction$type_of_test == "OGTT" & GSE88794_diet_restriction$`sampling time challenge test (h):ch1` == 0, ]
subject_ids_two_time_point <- names(table(fasting_ogtt_data_only$subject_id)[table(fasting_ogtt_data_only$subject_id) == 2])
fasting_ogtt_data_only <- fasting_ogtt_data_only[fasting_ogtt_data_only$subject_id %in% subject_ids_two_time_point, ]
intervention_group_names <- c(paste0("Control\n(n = ", length(subject_ids_two_time_point) - sum(fasting_ogtt_data_only$`intervention group assignment (er/control):ch1` == "ER") / 2, ")"), paste0("Calorie restriction\n(n = ", sum(fasting_ogtt_data_only$`intervention group assignment (er/control):ch1` == "ER") / 2, ")"))
fasting_ogtt_data_only$intervention_group_name <- intervention_group_names[1]
fasting_ogtt_data_only$intervention_group_name[fasting_ogtt_data_only$`intervention group assignment (er/control):ch1` == "ER"] <- intervention_group_names[2]
fasting_ogtt_data_only$intervention_group_name <- ordered(fasting_ogtt_data_only$intervention_group_name, levels = intervention_group_names)
fasting_ogtt_data_only$intervention_status <- ifelse(fasting_ogtt_data_only$intervention_status == "before intervention", "Baseline", "Week 12")
fasting_ogtt_data_only$intervention_status <- ordered(fasting_ogtt_data_only$intervention_status, levels = c("Baseline", "Week 12"))
###############################################################


###############################################################
plot_each_score_diet <- function(score_pheno_frame, score_name, score_name_plot, print_NS = FALSE)
{
  score_pheno_frame$scaled_score <- scale(score_pheno_frame[[score_name]])
  group_colours <- c("Baseline" = "#EE3377", "Week 12" = "#EE7733")
  test_significance <- score_pheno_frame %>% group_by(intervention_group_name) %>%
    wilcox_test(scaled_score ~ intervention_status, comparisons = list(c("Baseline", "Week 12")), paired = TRUE) %>%
    add_xy_position(x = "intervention_status")
  test_significance$y.position <- 1.03 * test_significance$y.position
  test_significance$xmin <- 1.02 * test_significance$xmin
  test_significance$xmax <- 0.98 * test_significance$xmax
  if(!print_NS)
  {
    test_significance$p_print <- as.character(ifelse(test_significance$p < 0.01, format(test_significance$p, digits = 1, scientific = TRUE), round(test_significance$p, 2)))
  }
  else
  {
    test_significance$p_print <- as.character(ifelse(test_significance$p < 0.01, format(test_significance$p, digits = 1, scientific = TRUE), ifelse(test_significance$p >= 0.055, "NS", round(test_significance$p, 2))))
  }
  score_summary <- ggplot(score_pheno_frame, aes(x = intervention_status, y = scaled_score, colour = intervention_status)) +
    ylab(score_name_plot) +
    geom_beeswarm(cex = cex_val * 1.1, size = point_size) +
    geom_line(aes(group = subject_id), position = position_beeswarm(cex = cex_val * 1.1), colour = black_line_colour, alpha = 0.2, lwd = line_size_text_repel_mm) +
    geom_boxplot(fill = NA, lwd = line_size_main_mm, width = 0.4, outlier.shape = NA, colour = black_text_colour) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank(),
      legend.position = "bottom", legend.direction = "horizontal", legend.margin = margin(1, 1, 1, 1), legend.box.margin = margin(-8, 5, 3, 3)) +
    stat_pvalue_manual(test_significance, label = "p_print",
      tip.length = 0,
      size = text_size_labels, colour = black_text_colour, bracket.size = line_size_text_repel_mm) +
    scale_colour_manual(values = group_colours) +
    labs(colour = NULL) + guides(colour = guide_legend(nrow = 1, override.aes = list(size = text_size_labels * 6 / 9))) +
    facet_wrap(vars(intervention_group_name), nrow = 1, scales = "free")
  return(score_summary)
}
GSE88794_module_1_plot <- plot_each_score_diet(fasting_ogtt_data_only, "module_1_score", "Module 1 (z score)", print_NS = TRUE)
GSE88794_module_2_plot <- plot_each_score_diet(fasting_ogtt_data_only, "module_2_score", "Module 2 (z score)")
GSE88794_som_plot <- plot_each_score_diet(fasting_ogtt_data_only, "som_score", "SoM (z score)", print_NS = TRUE)
GSE88794_legend <- as_ggplot(ggpubr::get_legend(GSE88794_som_plot))
GSE88794_combined <- plot_grid(GSE88794_module_1_plot + theme(legend.position = "none"), GSE88794_som_plot + theme(legend.position = "none"), nrow = 1, align = "hv", axis = "tblr")
GSE88794_tsang_ihm_plot <- plot_each_score_diet(fasting_ogtt_data_only, "tsang_ihm_score", "IHM (z score)", print_NS = TRUE)
###############################################################


###############################################################
cairo_pdf("/labs/khatrilab/ananthg/detrimental_host_response/figures_2024_09_27/fig_4_risk_factors.pdf", width = 7.2, height = 8, bg = "transparent")

pushViewport(viewport(layout = grid.layout(nrow = 100, ncol = 100)))

grid.rect(x = unit(0.275, "npc"), y = unit(0.98, "npc"), width = unit(0.495, "npc"), height = unit(0.0225, "npc"), gp = gpar(fill = "#EFEFEF", col = "transparent"))
grid.rect(x = unit(0.775, "npc"), y = unit(0.98, "npc"), width = unit(0.4125, "npc"), height = unit(0.0225, "npc"), gp = gpar(fill = "#EFEFEF", col = "transparent"))
grid.rect(x = unit(0.29, "npc"), y = unit(0.7, "npc"), width = unit(0.4025, "npc"), height = unit(0.0225, "npc"), gp = gpar(fill = "#EFEFEF", col = "transparent"))
grid.rect(x = unit(0.8, "npc"), y = unit(0.7, "npc"), width = unit(0.32, "npc"), height = unit(0.0225, "npc"), gp = gpar(fill = "#EFEFEF", col = "transparent"))
grid.rect(x = unit(0.24, "npc"), y = unit(0.3, "npc"), width = unit(0.37, "npc"), height = unit(0.0225, "npc"), gp = gpar(fill = "#EFEFEF", col = "transparent"))
grid.rect(x = unit(0.76, "npc"), y = unit(0.3, "npc"), width = unit(0.44, "npc"), height = unit(0.0225, "npc"), gp = gpar(fill = "#EFEFEF", col = "transparent"))

pushViewport(viewport(layout.pos.row = 1:23, layout.pos.col = 1:55))
draw(ht_mp_eff_size, height = nrow(eff_size_matrix) * tile_size, gap = gap_size, background = "transparent",
  column_title = eff_size_plot_title, column_title_gp = gpar(fontsize = base_text_size, fontfamily = base_font_family),
  newpage = FALSE)
decorate_annotation("Datasets", {
    grid.text("Datasets", rot = 45, y = unit(-0.05, "npc"), just = "right", gp = text_setting_small)
})
decorate_annotation("Cases", {
    grid.text("Cases", rot = 45, y = unit(-0.05, "npc"), just = "right", gp = text_setting_small)
})
decorate_annotation("Controls", {
    grid.text("Controls", rot = 45, y = unit(-0.05, "npc"), just = "right", gp = text_setting_small)
})
upViewport(1)

pushViewport(viewport(layout.pos.row = 21:24, layout.pos.col = 1:15))
draw(eff_size_legend)
upViewport(1)
grid.text(x = unit(0.18, "npc"), y = unit(0.78, "npc"), label = "* p \U2264 0.05", gp = gpar(fontsize = base_text_size, fontface = "plain", col = black_text_colour))

pushViewport(viewport(layout.pos.row = 1:26, layout.pos.col = 56:100))
draw(ht_mp_regression, height = nrow(regression_coeff) * tile_size, gap = gap_size, background = "transparent",
  column_title = "Framingham cohort regression (5321 subjects)", column_title_gp = gpar(fontsize = base_text_size, fontfamily = base_font_family),
  newpage = FALSE)
upViewport(1)

pushViewport(viewport(layout.pos.row = 21:24, layout.pos.col = 56:70))
draw(regression_legend)
upViewport(1)
grid.text(x = unit(0.75, "npc"), y = unit(0.78, "npc"), label = "* p < 0.05", gp = gpar(fontsize = base_text_size, fontface = "plain", col = black_text_colour))

grid.text(x = unit(0.29, "npc"), y = unit(0.70, "npc"), label = "Effect of number of risk factors (Framingham)", gp = gpar(fontsize = base_text_size, fontface = "plain", col = black_text_colour))
print(risk_factors_combined, vp = viewport(layout.pos.row = 32:67, layout.pos.col = 1:58))

grid.text(x = unit(0.8, "npc"), y = unit(0.70, "npc"), label = "12 week diet restriction (GSE88794)", gp = gpar(fontsize = base_text_size, fontface = "plain", col = black_text_colour))
print(GSE88794_module_2_plot, vp = viewport(layout.pos.row = 32:67, layout.pos.col = 61:100))

grid.text(x = unit(0.24, "npc"), y = unit(0.3, "npc"), label = "Current vs. former smokers (Framingham)", gp = gpar(fontsize = base_text_size, fontface = "plain", col = black_text_colour))
print(smoking_combined, vp = viewport(layout.pos.row = 72:97, layout.pos.col = 1:49))
print(smoking_legend, vp = viewport(layout.pos.row = 98:100, layout.pos.col = 1:49))

grid.text(x = unit(0.76, "npc"), y = unit(0.3, "npc"), label = "Controlled vs. uncontrolled diabetes (Framingham)", gp = gpar(fontsize = base_text_size, fontface = "plain", col = black_text_colour))
print(a1c_combined, vp = viewport(layout.pos.row = 72:97, layout.pos.col = 52:100))
print(a1c_legend, vp = viewport(layout.pos.row = 98:100, layout.pos.col = 52:100))

grid.text(x = unit(0.015, "npc"), y = unit(0.985, "npc"), label = "A", gp = gpar(fontsize = base_text_size + 2, fontface = "bold", col = black_text_colour))
grid.text(x = unit(0.555, "npc"), y = unit(0.985, "npc"), label = "B", gp = gpar(fontsize = base_text_size + 2, fontface = "bold", col = black_text_colour))
grid.text(x = unit(0.02, "npc"), y = unit(0.7075, "npc"), label = "C", gp = gpar(fontsize = base_text_size + 2, fontface = "bold", col = black_text_colour))
grid.text(x = unit(0.61, "npc"), y = unit(0.7075, "npc"), label = "D", gp = gpar(fontsize = base_text_size + 2, fontface = "bold", col = black_text_colour))
grid.text(x = unit(0.02, "npc"), y = unit(0.3075, "npc"), label = "E", gp = gpar(fontsize = base_text_size + 2, fontface = "bold", col = black_text_colour))
grid.text(x = unit(0.51, "npc"), y = unit(0.3075, "npc"), label = "F", gp = gpar(fontsize = base_text_size + 2, fontface = "bold", col = black_text_colour))

dev.off()
###############################################################


###############################################################
som_ihm_corr <- ggplot(coconut_combined, aes(x = som_score, y = tsang_ihm_score)) +
  xlab("SoM score") + ylab("IHM score") +
  geom_point(size = point_size) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  stat_smooth(method = "lm", formula = y ~ x, geom = "smooth", lwd = line_size_main_mm) +
  stat_cor(p.accuracy = 0.01, r.accuracy = 0.01, size = text_size_labels, label.y = Inf, label.x = 0, hjust = 1, vjust = 1.5)
###############################################################


###############################################################
proteomic_regression_coeff <- matrix(c(-0.1047, 0.0147, 0.0389, -0.0070, -0.0238, 0.0023, 0.0066, 0.0634, 0.1363, 0.0006, 0.0172, -0.0061, 0.1817, 0.0008, 0.0199, 0.0290, 0.0309, 0.0042, 0.0151, 0.0412),
  nrow = 5, byrow = TRUE)
rownames(proteomic_regression_coeff) <- c("ADM", "GRN", "CEACAM8", "AZU1", "OLR1")
colnames(proteomic_regression_coeff) <- c("Sex", "Age", "BMI", "Disease")
proteomic_regression_p_value <- matrix(c(2.42E-08, 1.56E-74, 4.90E-144, 8.20E-01, 9.51E-02, 8.65E-05, 1.06E-09, 7.34E-03, 1.70E-08, 5.80E-01, 1.12E-20, 8.78E-01, 1.81E-08, 5.33E-01, 6.97E-16, 5.87E-01, 1.69E-01, 7.26E-06, 2.23E-18, 2.70E-01),
  nrow = 5, byrow = TRUE)
proteomic_regression_coeff <- proteomic_regression_coeff[, c(2, 1, 3, 4)]
proteomic_regression_p_value <- proteomic_regression_p_value[, c(2, 1, 3, 4)]
proteomic_min_regression <- -0.19
proteomic_max_regression <- 0.19
proteomic_regression_matrix_colour_pts <- c(proteomic_min_regression, proteomic_min_regression/2, 0, proteomic_max_regression/2, proteomic_max_regression)
proteomic_regression_cols <- colorRamp2(proteomic_regression_matrix_colour_pts, c("#125A56", "#60BCE9", "#ECEADA", "#FD9A44", "#A01813"))
proteomic_regression_matrix_legend_pts <- c(proteomic_min_regression, proteomic_max_regression)
proteomic_regression_matrix_legend_labels <- c(round(proteomic_min_regression, 2), round(proteomic_max_regression, 2))
module_annotation <- rowAnnotation(Module = anno_text(c(1, 1, 2, 2, 2), location = unit(0.6, "npc"), just = "right", gp = text_setting_sub_small, width = tile_size))
proteomic_regression_legend <- Legend(col_fun = proteomic_regression_cols,
  at = proteomic_regression_matrix_legend_pts,
  labels = proteomic_regression_matrix_legend_labels,
  title = "Regression\ncoefficient", title_position = "topcenter",
  grid_height = unit(base_text_size / 3, "pt"),
  legend_width = 4 * unit(base_text_size, "pt"),
  background = "transparent",
  title_gp = text_setting_small,
  title_gap = gap_size, gap = gap_size,
  labels_gp = text_setting_sub_small,
  direction = "horizontal")
ht_mp_proteomic_regression <- Heatmap(proteomic_regression_coeff,
  col = proteomic_regression_cols, na_col = "#FFFFFF",
  cell_fun = function(j, i, x, y, width, height, fill)
    {
      if(proteomic_regression_p_value[i, j] < 0.05)
      {
        grid.text("*", x, y, vjust = 0.8, gp = gpar(fontsize = base_text_size * 1.5, fontfamily = base_font_family, col = black_line_colour))
      }
    },
  rect_gp = dend_lines_gp,
  column_names_rot = 45,
  width = ncol(proteomic_regression_coeff) * tile_size,
  name = "regression",
  row_names_side = "left",
  show_heatmap_legend = FALSE,
  cluster_rows = FALSE, cluster_columns = FALSE,
  gap = gap_size,
  left_annotation = module_annotation)
###############################################################


###############################################################
cairo_pdf("/labs/khatrilab/ananthg/detrimental_host_response/figures_2024_09_27/fig_S4_risk_factors.pdf", width = 7.2, height = 8, bg = "transparent")

pushViewport(viewport(layout = grid.layout(nrow = 100, ncol = 100)))

print(GSE88794_combined, vp = viewport(layout.pos.row = 1:28, layout.pos.col = 1:100))
print(GSE88794_legend, vp = viewport(layout.pos.row = 29:31, layout.pos.col = 1:100))

print(som_ihm_corr, vp = viewport(layout.pos.row = 34:66, layout.pos.col = 1:40))

print(GSE88794_tsang_ihm_plot, vp = viewport(layout.pos.row = 34:66, layout.pos.col = 50:100))

print(smoking_tsang_ihm_plot + theme(legend.position = "none"), vp = viewport(layout.pos.row = 69:100, layout.pos.col = 1:21))
print(as_ggplot(ggpubr::get_legend(smoking_tsang_ihm_plot)), vp = viewport(layout.pos.row = 69:83, layout.pos.col = 22:37))

print(a1c_tsang_ihm_plot + theme(legend.position = "none"), vp = viewport(layout.pos.row = 69:100, layout.pos.col = 40:60))
print(as_ggplot(ggpubr::get_legend(a1c_tsang_ihm_plot)), vp = viewport(layout.pos.row = 86:100, layout.pos.col = 24:40))

pushViewport(viewport(layout.pos.row = 70:89, layout.pos.col = 62:99))
draw(ht_mp_proteomic_regression, height = nrow(proteomic_regression_coeff) * tile_size, gap = gap_size, background = "transparent",
  column_title = "ISB cohort regression\n(proteomics, 2487 subjects)", column_title_gp = gpar(fontsize = base_text_size, fontfamily = base_font_family),
  newpage = FALSE)
decorate_annotation("Module", {
    grid.text("Module", rot = 45, y = unit(-0.05, "npc"), just = "right", gp = text_setting_small)
})
upViewport(1)

pushViewport(viewport(layout.pos.row = 93:97, layout.pos.col = 62:89))
draw(proteomic_regression_legend)
upViewport(1)
grid.text(x = unit(0.9, "npc"), y = unit(0.05, "npc"), label = "* p < 0.05", gp = gpar(fontsize = base_text_size, fontface = "plain", col = black_text_colour))

grid.text(x = unit(0.01, "npc"), y = unit(0.985, "npc"), label = "A", gp = gpar(fontsize = base_text_size + 2, fontface = "bold", col = black_text_colour))
grid.text(x = unit(0.51, "npc"), y = unit(0.985, "npc"), label = "B", gp = gpar(fontsize = base_text_size + 2, fontface = "bold", col = black_text_colour))
grid.text(x = unit(0.01, "npc"), y = unit(0.675, "npc"), label = "C", gp = gpar(fontsize = base_text_size + 2, fontface = "bold", col = black_text_colour))
grid.text(x = unit(0.51, "npc"), y = unit(0.675, "npc"), label = "D", gp = gpar(fontsize = base_text_size + 2, fontface = "bold", col = black_text_colour))
grid.text(x = unit(0.01, "npc"), y = unit(0.325, "npc"), label = "E", gp = gpar(fontsize = base_text_size + 2, fontface = "bold", col = black_text_colour))
grid.text(x = unit(0.4, "npc"), y = unit(0.325, "npc"), label = "F", gp = gpar(fontsize = base_text_size + 2, fontface = "bold", col = black_text_colour))
grid.text(x = unit(0.62, "npc"), y = unit(0.325, "npc"), label = "G", gp = gpar(fontsize = base_text_size + 2, fontface = "bold", col = black_text_colour))

dev.off()
###############################################################

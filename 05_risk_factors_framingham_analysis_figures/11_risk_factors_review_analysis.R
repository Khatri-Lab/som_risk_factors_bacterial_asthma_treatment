###############################################################
library(readr)
library(plyr)
source(file = "/labs/khatrilab/ananthg/detrimental_host_response/common_functions/perform_meta_analyses.R")
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
correlation_data <- coconut_combined
correlation_data$sex <- as.numeric(as.factor(correlation_data$sex))
correlation_data$disease_status <- as.numeric(correlation_data$disease_status)
correlation_data$regular_smokers <- as.numeric(correlation_data$regular_smokers)
risk_factors <- c("age", "sex", "bmi", "disease_status", "regular_smokers")
corr_coeff <- data.frame(matrix(, nrow = length(risk_factors), ncol = length(risk_factors)), row.names = risk_factors)
colnames(corr_coeff) <- risk_factors
corr_p_value <- data.frame(matrix(, nrow = length(risk_factors), ncol = length(risk_factors)), row.names = risk_factors)
colnames(corr_p_value) <- risk_factors
for(i in risk_factors)
{
  for(j in risk_factors)
  {
    cor_test <- cor.test(correlation_data[[i]], correlation_data[[j]])
    corr_p_value[i, j] <- cor_test$p.value
    corr_coeff[i, j] <- cor_test$estimate
  }
}
rownames(corr_coeff) <- c("age", "sex", "bmi", "disease", "smokers")
colnames(corr_coeff) <- c("age", "sex", "bmi", "disease", "smokers")
rownames(corr_p_value) <- c("age", "sex", "bmi", "disease", "smokers")
colnames(corr_p_value) <- c("age", "sex", "bmi", "disease", "smokers")
###############################################################


###############################################################
corr_coeff_plot_title <- "Framingham: Correlation analyses"
min_corr <- -1
max_corr <- 1
corr_matrix_colour_pts <- c(min_corr, min_corr/2, 0, max_corr/2, max_corr)
corr_cols <- colorRamp2(corr_matrix_colour_pts, c("#762A83", "#C2A5CF", "#F7F7F7", "#ACD39E", "#1B7837"))
corr_matrix_legend_pts <- c(min_corr, max_corr)
corr_matrix_legend_labels <- c(round(min_corr, 2), round(max_corr, 2))
corr_legend <- Legend(col_fun = corr_cols,
  at = corr_matrix_legend_pts,
  labels = corr_matrix_legend_labels,
  title = "Correlation", title_position = "topcenter",
  grid_height = unit(base_text_size / 3, "pt"),
  legend_width = 4 * unit(base_text_size, "pt"),
  background = "transparent",
  title_gp = text_setting_small,
  title_gap = gap_size, gap = gap_size,
  labels_gp = text_setting_sub_small,
  direction = "horizontal")
ht_mp_corr <- Heatmap(as.matrix(corr_coeff),
  col = corr_cols, na_col = "#FFFFFF",
  # cell_fun = function(j, i, x, y, width, height, fill)
  #   {
  #     if(corr_p_value[i, j] <= 0.05)
  #     {
  #       grid.text("*", x, y, vjust = 0.8, gp = gpar(fontsize = base_text_size * 1.5, fontfamily = base_font_family, col = black_line_colour))
  #     }
  #   },
  rect_gp = dend_lines_gp,
  column_names_rot = 45,
  width = ncol(corr_coeff) * tile_size,
  name = "corr",
  row_names_side = "left",
  show_heatmap_legend = FALSE,
  cluster_rows = FALSE, cluster_columns = FALSE,
  gap = gap_size)
###############################################################


###############################################################
colinear_data <- coconut_combined
colinear_data$age_discrete <- "[40, 60]"
colinear_data$age_discrete[colinear_data$age < 40] <- "< 40"
colinear_data$age_discrete[colinear_data$age > 60] <- "> 60"
colinear_data$age_discrete <- ordered(colinear_data$age_discrete, levels = c("< 40", "[40, 60]", "> 60"))
colinear_data$bmi_discrete <- "[25, 30]"
colinear_data$bmi_discrete[colinear_data$bmi < 25] <- "< 25"
colinear_data$bmi_discrete[colinear_data$bmi > 30] <- "> 30"
colinear_data$bmi_discrete <- ordered(colinear_data$bmi_discrete, levels = c("< 25", "[25, 30]", "> 30"))
colinear_data$disease <- ordered(ifelse(colinear_data$healthy == 1, "Healthy", "Disease"), levels = c("Healthy", "Disease"))
colinear_data$smoker <- ordered(ifelse(colinear_data$regular_smokers == 1, "Yes", "No"), levels = c("Yes", "No"))
colinear_data <- colinear_data[, c("age_discrete", "sex", "bmi_discrete", "disease", "smoker")]
colinear_test <- as.data.frame(colinear_data %>% group_by(age_discrete, sex, bmi_discrete, disease, smoker) %>% summarise(count = n()))
# There are no men < 40 years of age who are healthy, obese, and smoke. All other groups have at least one subject represented.
###############################################################


###############################################################
framingham_conormalized_results <- readRDS("/local-scratch/projects/candx/guangbo/ananth_processesed_framingham_data/dhr_paper_re_processesed_2024_09_27/framingham_conormalized_by_healthy_obj.RDS")
conorm_expr <- Reduce(cbind, list(framingham_conormalized_results$COCONUTList$off_spr$genes,
  framingham_conormalized_results$COCONUTList$gen_iii$genes,
  framingham_conormalized_results$controlList$GSEs$off_spr$genes,
  framingham_conormalized_results$controlList$GSEs$gen_iii$genes))
conorm_expr <- data.frame(t(conorm_expr[c(genes_mod_1, genes_mod_2), ]))
conorm_expr$shareid <- rownames(conorm_expr)
framingham_expr_data <- merge(coconut_combined, conorm_expr, by = 'shareid', all.x = TRUE, all.y = FALSE)
###############################################################


###############################################################
smokers_only <- framingham_expr_data[framingham_expr_data$regular_smokers == 1, c("age", "sex", "bmi", "disease_status", "years_no_smoking", "module_1_score", "module_2_score", "module_3_score", "module_4_score", "detrimental_score", "protective_score", "som_score", "tsang_up_score", "tsang_down_score", "tsang_ihm_score", "smoking_quiters", genes_mod_1, genes_mod_2)]
smokers_only$years_no_smoking_discrete <- ifelse(smokers_only$smoking_quiters == 0, "Current smokers", ifelse(smokers_only$years_no_smoking >= 5, "Past smokers", "Recent smokers"))
smokers_only <- na.omit(smokers_only[, c("age", "sex", "bmi", "disease_status", "years_no_smoking_discrete", "module_1_score", "module_2_score", "module_3_score", "module_4_score", "detrimental_score", "protective_score", "som_score", "tsang_up_score", "tsang_down_score", "tsang_ihm_score", "smoking_quiters", genes_mod_1, genes_mod_2)])
smokers_only$years_no_smoking_discrete <- ordered(smokers_only$years_no_smoking_discrete, levels = c("Current smokers", "Recent smokers", "Past smokers"))
###############################################################


###############################################################
diabetes_only <- framingham_expr_data[framingham_expr_data$diabetes == 1, c("age", "sex", "bmi", "disease_status", "regular_smokers", "module_1_score", "module_2_score", "module_3_score", "module_4_score", "detrimental_score", "protective_score", "som_score", "tsang_up_score", "tsang_down_score", "tsang_ihm_score", "A1C", genes_mod_1, genes_mod_2)]
diabetes_only$a1c_discrete <- ifelse(diabetes_only$A1C < 7, "A1C < 7", ifelse(diabetes_only$A1C > 10, "A1C > 10", "7 <= A1C <= 10"))
diabetes_only <- na.omit(diabetes_only)
diabetes_only$a1c_discrete <- ordered(diabetes_only$a1c_discrete, levels = c("A1C < 7", "7 <= A1C <= 10", "A1C > 10"))
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
#
GSE88794 <- readRDS("/labs/khatrilab/ananthg/detrimental_host_response/data_2024_09_27/GSE88794_raw.RDS")
GSE88794_module_2 <- data.frame(t(MetaIntegrator::getSampleLevelGeneData(GSE88794, c(genes_mod_1, genes_mod_2))))
GSE88794_module_2$geo_accession <- rownames(GSE88794_module_2)
fasting_ogtt_data_only_module_2 <- merge(fasting_ogtt_data_only, GSE88794_module_2, by = "geo_accession", all.x = TRUE, all.y = FALSE)
###############################################################


###############################################################
perform_diff_exp_analysis <- function(expr_data, gene_list, pheno_name, pheno_levels, summary_text)
{
  expr_data_long <- expr_data[, c(gene_list, pheno_name)] %>% pivot_longer(all_of(gene_list), names_to = "gene_name", values_to = "expr")
  expr_data_long$pheno <- as.character(expr_data_long[[pheno_name]])
  multiple_testing <- expr_data_long %>% group_by(gene_name) %>%
    wilcox_test(expr ~ pheno, comparisons = list(pheno_levels), paired = FALSE, ref.group =pheno_levels[2], p.adjust.method = "none") %>%
    adjust_pvalue(method = "fdr")
  effect_size <- expr_data_long %>% group_by(gene_name) %>%
    cohens_d(expr ~ pheno, comparisons = list(pheno_levels), paired = FALSE, ref.group =pheno_levels[1])
  multiple_testing <- merge(multiple_testing, effect_size[, c("gene_name", "effsize", "magnitude")], by = "gene_name")
  multiple_testing$summary_text <- paste0(summary_text, pheno_levels[1], " vs. ", pheno_levels[2])
  return(multiple_testing)
}
###############################################################


###############################################################
fasting_control_deg <- perform_diff_exp_analysis(fasting_ogtt_data_only_module_2[fasting_ogtt_data_only_module_2$intervention_group_name == intervention_group_names[1], ], c(genes_mod_1, genes_mod_2), "intervention_status", c("Baseline", "Week 12"), "Control (GSE88794): ")
fasting_diet_restriction_deg <- perform_diff_exp_analysis(fasting_ogtt_data_only_module_2[fasting_ogtt_data_only_module_2$intervention_group_name == intervention_group_names[2], ], c(genes_mod_1, genes_mod_2), "intervention_status", c("Baseline", "Week 12"), "Diet restriction (GSE88794): ")

smoking_current_recent_deg <- perform_diff_exp_analysis(smokers_only[smokers_only$years_no_smoking_discrete %in% c("Current smokers", "Recent smokers"), ], c(genes_mod_1, genes_mod_2), "years_no_smoking_discrete", c("Current smokers", "Recent smokers"), "Smoking (Framingham): ")
smoking_recent_past_deg <- perform_diff_exp_analysis(smokers_only[smokers_only$years_no_smoking_discrete %in% c("Recent smokers", "Past smokers"), ], c(genes_mod_1, genes_mod_2), "years_no_smoking_discrete", c("Recent smokers", "Past smokers"), "Smoking (Framingham): ")

diabetes_controlled_moderate_deg <- perform_diff_exp_analysis(diabetes_only[diabetes_only$a1c_discrete %in% c("A1C < 7", "7 <= A1C <= 10"), ], c(genes_mod_1, genes_mod_2), "a1c_discrete", c("A1C < 7", "7 <= A1C <= 10"), "Diabetes (Framingham): ")
diabetes_moderate_uncontrolled_deg <- perform_diff_exp_analysis(diabetes_only[diabetes_only$a1c_discrete %in% c("7 <= A1C <= 10", "A1C > 10"), ], c(genes_mod_1, genes_mod_2), "a1c_discrete", c("7 <= A1C <= 10", "A1C > 10"), "Diabetes (Framingham): ")

combined_data <- Reduce(rbind, list(fasting_control_deg, fasting_diet_restriction_deg, smoking_current_recent_deg, smoking_recent_past_deg, diabetes_controlled_moderate_deg, diabetes_moderate_uncontrolled_deg))
###############################################################


###############################################################
combined_data_fdr <- as.data.frame(combined_data[, c("gene_name", "p.adj", "summary_text")] %>% pivot_wider(names_from = "gene_name", values_from = "p.adj"))
rownames(combined_data_fdr) <- combined_data_fdr$summary_text
combined_data_fdr <- subset(combined_data_fdr, select = -c(summary_text))
combined_data_eff_size <- as.data.frame(combined_data[, c("gene_name", "effsize", "summary_text")] %>% pivot_wider(names_from = "gene_name", values_from = "effsize"))
rownames(combined_data_eff_size) <- combined_data_eff_size$summary_text
combined_data_eff_size <- subset(combined_data_eff_size, select = -c(summary_text))
###############################################################


###############################################################
min_eff_size <- -0.75
max_eff_size <- 0.75
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
ht_mp_eff_size <- Heatmap(as.matrix(combined_data_eff_size),
  col = eff_size_cols, na_col = "#FFFFFF",
  cell_fun = function(j, i, x, y, width, height, fill)
    {
      if(combined_data_eff_size[i, j] <= 0.05)
      {
        grid.text("*", x, y, vjust = 0.8, gp = gpar(fontsize = base_text_size * 1.5, fontfamily = base_font_family, col = black_line_colour))
      }
    },
  rect_gp = dend_lines_gp,
  column_names_rot = 45,
  width = ncol(combined_data_eff_size) * tile_size,
  name = "eff_size",
  row_names_side = "left",
  show_heatmap_legend = FALSE,
  cluster_rows = FALSE, cluster_columns = FALSE,
  gap = gap_size)
###############################################################


###############################################################
cairo_pdf("/labs/khatrilab/ananthg/detrimental_host_response/figures_2024_09_27/fig_RX_diff_exp_analysis.pdf", width = 9.25, height = 2.25, bg = "transparent")

pushViewport(viewport(layout = grid.layout(nrow = 100, ncol = 100)))

pushViewport(viewport(layout.pos.row = 1:100, layout.pos.col = 16:100))
draw(ht_mp_eff_size, height = nrow(combined_data_eff_size) * tile_size, gap = gap_size, background = "transparent",
  column_title = "Differential expression analysis", column_title_gp = gpar(fontsize = base_text_size, fontfamily = base_font_family),
  newpage = FALSE)
upViewport(1)

pushViewport(viewport(layout.pos.row = 85:97, layout.pos.col = 1:10))
draw(eff_size_legend)
upViewport(1)
grid.text(x = unit(0.15, "npc"), y = unit(0.05, "npc"), label = "* fdr < 0.05", gp = gpar(fontsize = base_text_size, fontface = "plain", col = black_text_colour))

dev.off()
###############################################################


###############################################################
cairo_pdf("/labs/khatrilab/ananthg/detrimental_host_response/figures_2024_09_27/fig_RY_feature_correlation.pdf", width = 2.25, height = 2, bg = "transparent")

pushViewport(viewport(layout = grid.layout(nrow = 100, ncol = 100)))

pushViewport(viewport(layout.pos.row = 1:100, layout.pos.col = 1:100))
draw(ht_mp_corr, height = nrow(corr_coeff) * tile_size, gap = gap_size, background = "transparent",
  column_title = "Feature correlation (Framingham)", column_title_gp = gpar(fontsize = base_text_size, fontfamily = base_font_family),
  newpage = FALSE)
upViewport(1)

pushViewport(viewport(layout.pos.row = 85:95, layout.pos.col = 10:22))
draw(corr_legend)
upViewport(1)
#grid.text(x = unit(0.15, "npc"), y = unit(0.05, "npc"), label = "* fdr < 0.05", gp = gpar(fontsize = base_text_size, fontface = "plain", col = black_text_colour))

dev.off()
###############################################################

###############################################################
library(readr)
source(file = "/labs/khatrilab/ananthg/detrimental_host_response/common_functions/plot_themes.R")
###############################################################


###############################################################
framingham_score_pheno_matrix <- read_csv("/local-scratch/projects/candx/guangbo/ananth_processesed_framingham_data/dhr_paper_re_processesed_2024_09_27/framingham_conormalized_by_healthy_matrix.csv")
framingham_score_pheno_matrix$som_score <- framingham_score_pheno_matrix$som_score_ratio
###############################################################


###############################################################
print_cox_model_summary <- function(cox_model)
{
  summary_cox_model <- summary(cox_model)
  summary_cox_model_data_frame <- as.data.frame(summary_cox_model$coefficients)[, c(2, 5)]
  colnames(summary_cox_model_data_frame) <- c("Hazard ratio", "p-value")
  summary_cox_model_data_frame$`Hazard ratio` <- paste0(round(summary_cox_model_data_frame$`Hazard ratio`, 2), " (", round(summary_cox_model$conf.int[, 3], 2), ", ", round(summary_cox_model$conf.int[, 4], 2), ")")
  summary_cox_model_data_frame$`p-value` <- as.character(ifelse(summary_cox_model_data_frame$`p-value` < 0.01, format(summary_cox_model_data_frame$`p-value`, digits = 2, scientific = TRUE), round(summary_cox_model_data_frame$`p-value`, 2)))
  return(summary_cox_model_data_frame)
}
###############################################################


###############################################################
evaluate_framingham_mortality_models <- function(n_years_cutoff, score_pheno_matrix, score_name, score_name_print, additional_variates_formula, additional_variates_name)
{
  score_pheno_matrix$status <- 1
  score_pheno_matrix$status[is.na(score_pheno_matrix$death_date) | score_pheno_matrix$death_date > 365 * n_years_cutoff] <- 0
  score_pheno_matrix$death_date[is.na(score_pheno_matrix$death_date) | score_pheno_matrix$death_date > 365 * n_years_cutoff] <- 365 * n_years_cutoff
  score_pheno_matrix$death_date <- score_pheno_matrix$death_date / 365
  score_pheno_matrix$disease_status <- as.factor(1 - score_pheno_matrix$healthy)
  print(sum(score_pheno_matrix$status == 1))
  print(sum(score_pheno_matrix$status == 0))
  #
  cox_model_all_var <- coxph(as.formula(paste0("Surv(death_date, status) ~ ", score_name, additional_variates_formula, " + age + as.factor(sex) + bmi + disease_status + as.factor(regular_smokers)")), data = score_pheno_matrix)
  cox_model_all_var_summary <- print_cox_model_summary(cox_model_all_var)
  if(length(additional_variates_name) > 0)
  {
    rownames(cox_model_all_var_summary) <- c(score_name_print, additional_variates_name, "Age", "Sex", "BMI", "Disease", "Smoking")
  }
  else
  {
    rownames(cox_model_all_var_summary) <- c(score_name_print, "Age", "Sex", "BMI", "Disease", "Smoking")
  }
  #
  lm_score_adjust <- lm(as.formula(paste0(score_name, " ~ age + as.factor(sex) + bmi + disease_status + as.factor(regular_smokers)", additional_variates_formula)), data = score_pheno_matrix)
  score_pheno_matrix$residuals <- NA
  score_pheno_matrix$residuals[as.numeric(names(lm_score_adjust$residuals))] <- lm_score_adjust$residuals
  score_pheno_matrix$score_discretized <- "High"
  score_pheno_matrix$score_discretized[score_pheno_matrix$residuals < 0] <- "Low"
  score_pheno_matrix$score_discretized[is.na(score_pheno_matrix$residuals)] <- NA
  score_pheno_matrix$score_discretized <- factor(score_pheno_matrix$score_discretized, levels = c("Low", "High"))
  #
  residual_cox_model_discretized <- coxph(Surv(death_date, status) ~ score_discretized, data = score_pheno_matrix)
  residual_cox_model_discretized_summary <- print_cox_model_summary(residual_cox_model_discretized)
  rownames(residual_cox_model_discretized_summary) <- paste0("Adjusted discretized ", score_name_print)
  #
  residual_cox_model <- coxph(Surv(death_date, status) ~ residuals, data = score_pheno_matrix)
  residual_cox_model_summary <- print_cox_model_summary(residual_cox_model)
  rownames(residual_cox_model_summary) <- paste0("Adjusted ", score_name_print)
  #
  km_fit <- survfit(Surv(death_date, status) ~ score_discretized, data = score_pheno_matrix)
  km_plot <- autoplot(km_fit, conf.int = TRUE) +
    labs(color = paste0("Adjusted ", score_name_print, "\n(discretized)"), fill = paste0("Adjusted ", score_name_print, "\n(discretized)")) +
    theme(legend.position = c(0.3, 0.25)) +
    xlab("Time (years)") + ylab("Survival rate") +
    scale_fill_manual(values = c("#1A5E1F", "#EE6C00")) +
    scale_colour_manual(values = c("#1A5E1F", "#EE6C00"))
  #
  return(list(cox_model_all_var_summary = cox_model_all_var_summary, residual_cox_model_discretized_summary = residual_cox_model_discretized_summary, residual_cox_model_summary = residual_cox_model_summary, km_plot = km_plot))
}
###############################################################


###############################################################
crp_only <- evaluate_framingham_mortality_models(7, framingham_score_pheno_matrix, "c_reactive_protein_serum", "CRP", "", c())
som_only <- evaluate_framingham_mortality_models(7, framingham_score_pheno_matrix, "som_score", "SoM", "", c())
ihm_only <- evaluate_framingham_mortality_models(7, framingham_score_pheno_matrix, "tsang_ihm_score", "IHM", "", c())
###############################################################


###############################################################
crp_vs_som <- evaluate_framingham_mortality_models(7, framingham_score_pheno_matrix, "c_reactive_protein_serum", "CRP", " + som_score", c("SoM"))
crp_vs_ihm <- evaluate_framingham_mortality_models(7, framingham_score_pheno_matrix, "c_reactive_protein_serum", "CRP", " + tsang_ihm_score", c("IHM"))
som_vs_crp <- evaluate_framingham_mortality_models(7, framingham_score_pheno_matrix, "som_score", "SoM", " + c_reactive_protein_serum", c("CRP"))
som_vs_ihm <- evaluate_framingham_mortality_models(7, framingham_score_pheno_matrix, "som_score", "SoM", " + tsang_ihm_score", c("IHM"))
ihm_vs_som <- evaluate_framingham_mortality_models(7, framingham_score_pheno_matrix, "tsang_ihm_score", "IHM", " + som_score", c("SoM"))
ihm_vs_crp <- evaluate_framingham_mortality_models(7, framingham_score_pheno_matrix, "tsang_ihm_score", "IHM", " + c_reactive_protein_serum", c("CRP"))
###############################################################


###############################################################
crp_vs_som_ihm <- evaluate_framingham_mortality_models(7, framingham_score_pheno_matrix, "c_reactive_protein_serum", "CRP", " + som_score + tsang_ihm_score", c("SoM", "IHM"))
som_vs_crp_ihm <- evaluate_framingham_mortality_models(7, framingham_score_pheno_matrix, "som_score", "SoM", " + c_reactive_protein_serum + tsang_ihm_score", c("CRP", "IHM"))
ihm_vs_som_crp <- evaluate_framingham_mortality_models(7, framingham_score_pheno_matrix, "tsang_ihm_score", "IHM", " + c_reactive_protein_serum + som_score", c("CRP", "SoM"))
###############################################################


###############################################################
som_cox_model_summaries <- evaluate_framingham_mortality_models(7, framingham_score_pheno_matrix, "som_score", "SoM", "", c())
som_cox_model_summaries_10_year <- evaluate_framingham_mortality_models(10, framingham_score_pheno_matrix, "som_score", "SoM", "", c())
som_cox_model_summaries_50_year <- evaluate_framingham_mortality_models(7, framingham_score_pheno_matrix[framingham_score_pheno_matrix$age > 50, ], "som_score", "SoM", "", c())
detrimental_cox_model_summaries <- evaluate_framingham_mortality_models(7, framingham_score_pheno_matrix, "detrimental_score", "Detrimental", "", c())
protective_cox_model_summaries <- evaluate_framingham_mortality_models(7, framingham_score_pheno_matrix, "protective_score", "Protective", "", c())
detrimental_protective_cox_model_summaries <- evaluate_framingham_mortality_models(7, framingham_score_pheno_matrix, "detrimental_score", "Detrimental", " + protective_score", c("Protective"))
ihm_vs_detrimental_protective_crp <- evaluate_framingham_mortality_models(7, framingham_score_pheno_matrix, "tsang_ihm_score", "IHM", " + c_reactive_protein_serum + detrimental_score + protective_score", c("CRP", "Detrimental", "Protective"))
###############################################################


###############################################################
table_theme <- ttheme(base_style = "blank",
  base_size = "default", base_colour = black_text_colour, padding = unit(c(5, 5), "pt"),
  colnames.style = colnames_style(size = base_text_size, face = "plain", fill = "transparent", linewidth = 0, linecolor = "transparent", parse = FALSE),
  rownames.style = rownames_style(size = base_text_size, face = "plain", fill = "transparent", linewidth = 0, linecolor = "transparent", parse = FALSE),#hjust = 0, x = 0.1
  tbody.style = tbody_style(size = base_text_size - 1, face = "plain", fill = "transparent", linewidth = 0, linecolor = "transparent", parse = FALSE))
cox_model_plot <- ggdraw() +
  draw_plot(ggtexttable(som_cox_model_summaries$cox_model_all_var_summary, theme = table_theme) %>%
    tab_add_hline(at.row = 2, row.side = "top", linewidth = line_size_text_repel_mm) %>%
    tab_add_vline(at.column = 2:3, column.side = "left", linewidth = line_size_text_repel_mm))
cox_residual_model_plot <- ggdraw() +
  draw_plot(ggtexttable(som_cox_model_summaries$residual_cox_model_summary, theme = table_theme) %>%
    tab_add_hline(at.row = 2, row.side = "top", linewidth = line_size_text_repel_mm) %>%
    tab_add_vline(at.column = 2:3, column.side = "left", linewidth = line_size_text_repel_mm))
###############################################################


###############################################################
cairo_pdf("/labs/khatrilab/ananthg/detrimental_host_response/figures_2024_09_27/fig_R5_survival.pdf", width = 7, height = 2.75, bg = "transparent")

pushViewport(viewport(layout = grid.layout(nrow = 100, ncol = 100)))

print(cox_model_plot, vp = viewport(layout.pos.row = 10:67, layout.pos.col = 3:47))
print(cox_residual_model_plot, vp = viewport(layout.pos.row = 83:93, layout.pos.col = 3:47))
print(som_cox_model_summaries$km_plot, vp = viewport(layout.pos.row = 10:100, layout.pos.col = 50:100))

grid.text(x = unit(0.02, "npc"), y = unit(0.90, "npc"), label = "A", gp = gpar(fontsize = base_text_size + 2, fontface = "bold", col = black_text_colour))
grid.text(x = unit(0.02, "npc"), y = unit(0.25, "npc"), label = "B", gp = gpar(fontsize = base_text_size + 2, fontface = "bold", col = black_text_colour))
grid.text(x = unit(0.51, "npc"), y = unit(0.90, "npc"), label = "C", gp = gpar(fontsize = base_text_size + 2, fontface = "bold", col = black_text_colour))
grid.text(x = unit(0.25, "npc"), y = unit(0.895, "npc"), label = "Cox proportional-hazards model", gp = gpar(fontsize = base_text_size, fontface = "plain", col = black_text_colour))
grid.text(x = unit(0.25, "npc"), y = unit(0.25, "npc"), label = "Cox model for adjusted SoM", gp = gpar(fontsize = base_text_size, fontface = "plain", col = black_text_colour))

grid.rect(x = unit(0.5, "npc"), y = unit(0.965, "npc"), width = unit(0.39, "npc"), height = unit(0.06, "npc"), gp = gpar(fill = "#EFEFEF", col = "transparent"))
grid.text(x = unit(0.5, "npc"), y = unit(0.965, "npc"), label = "7 year survival (Framingham)", gp = gpar(fontsize = base_text_size, fontface = "plain", col = black_text_colour))

grid.text(x = unit(0.75, "npc"), y = unit(0.18, "npc"), label = paste0("HR = ", strsplit(som_cox_model_summaries$residual_cox_model_discretized_summary[1, 1], " ")[[1]][1], "; p = ", som_cox_model_summaries$residual_cox_model_discretized_summary[1, 2]), gp = gpar(fontsize = base_text_size - 1, fontface = "plain", col = black_text_colour))

dev.off()
###############################################################


###############################################################
cox_model_plot_10_year <- ggdraw() +
  draw_plot(ggtexttable(som_cox_model_summaries_10_year$cox_model_all_var_summary, theme = table_theme) %>%
    tab_add_hline(at.row = 2, row.side = "top", linewidth = line_size_text_repel_mm) %>%
    tab_add_vline(at.column = 2:3, column.side = "left", linewidth = line_size_text_repel_mm))
cox_residual_model_plot_10_year <- ggdraw() +
  draw_plot(ggtexttable(som_cox_model_summaries_10_year$residual_cox_model_summary, theme = table_theme) %>%
    tab_add_hline(at.row = 2, row.side = "top", linewidth = line_size_text_repel_mm) %>%
    tab_add_vline(at.column = 2:3, column.side = "left", linewidth = line_size_text_repel_mm))
cox_model_plot_50_year <- ggdraw() +
  draw_plot(ggtexttable(som_cox_model_summaries_50_year$cox_model_all_var_summary, theme = table_theme) %>%
    tab_add_hline(at.row = 2, row.side = "top", linewidth = line_size_text_repel_mm) %>%
    tab_add_vline(at.column = 2:3, column.side = "left", linewidth = line_size_text_repel_mm))
cox_residual_model_plot_50_year <- ggdraw() +
  draw_plot(ggtexttable(som_cox_model_summaries_50_year$residual_cox_model_summary, theme = table_theme) %>%
    tab_add_hline(at.row = 2, row.side = "top", linewidth = line_size_text_repel_mm) %>%
    tab_add_vline(at.column = 2:3, column.side = "left", linewidth = line_size_text_repel_mm))
###############################################################


###############################################################
cox_model_plot_detrimental <- ggdraw() +
  draw_plot(ggtexttable(detrimental_cox_model_summaries$cox_model_all_var_summary, theme = table_theme) %>%
    tab_add_hline(at.row = 2, row.side = "top", linewidth = line_size_text_repel_mm) %>%
    tab_add_vline(at.column = 2:3, column.side = "left", linewidth = line_size_text_repel_mm))
cox_residual_model_plot_detrimental <- ggdraw() +
  draw_plot(ggtexttable(detrimental_cox_model_summaries$residual_cox_model_summary, theme = table_theme) %>%
    tab_add_hline(at.row = 2, row.side = "top", linewidth = line_size_text_repel_mm) %>%
    tab_add_vline(at.column = 2:3, column.side = "left", linewidth = line_size_text_repel_mm))
cox_model_plot_protective <- ggdraw() +
  draw_plot(ggtexttable(protective_cox_model_summaries$cox_model_all_var_summary, theme = table_theme) %>%
    tab_add_hline(at.row = 2, row.side = "top", linewidth = line_size_text_repel_mm) %>%
    tab_add_vline(at.column = 2:3, column.side = "left", linewidth = line_size_text_repel_mm))
cox_residual_model_plot_protective <- ggdraw() +
  draw_plot(ggtexttable(protective_cox_model_summaries$residual_cox_model_summary, theme = table_theme) %>%
    tab_add_hline(at.row = 2, row.side = "top", linewidth = line_size_text_repel_mm) %>%
    tab_add_vline(at.column = 2:3, column.side = "left", linewidth = line_size_text_repel_mm))
cox_model_plot_detrimental_protective <- ggdraw() +
  draw_plot(ggtexttable(detrimental_protective_cox_model_summaries$cox_model_all_var_summary, theme = table_theme) %>%
    tab_add_hline(at.row = 2, row.side = "top", linewidth = line_size_text_repel_mm) %>%
    tab_add_vline(at.column = 2:3, column.side = "left", linewidth = line_size_text_repel_mm))
cox_model_plot_ihm_crp_som <- ggdraw() +
  draw_plot(ggtexttable(ihm_vs_som_crp$cox_model_all_var_summary, theme = table_theme) %>%
    tab_add_hline(at.row = 2, row.side = "top", linewidth = line_size_text_repel_mm) %>%
    tab_add_vline(at.column = 2:3, column.side = "left", linewidth = line_size_text_repel_mm))
###############################################################


###############################################################
cairo_pdf("/labs/khatrilab/ananthg/detrimental_host_response/figures_2024_09_27/fig_RS5_survival.pdf", width = 7, height = 4.5, bg = "transparent")

pushViewport(viewport(layout = grid.layout(nrow = 150, ncol = 100)))

print(cox_model_plot_10_year, vp = viewport(layout.pos.row = 5:65, layout.pos.col = 2:49))
print(cox_residual_model_plot_10_year, vp = viewport(layout.pos.row = 66:77, layout.pos.col = 2:49))
print(cox_model_plot_50_year, vp = viewport(layout.pos.row = 5:65, layout.pos.col = 51:99))
print(cox_residual_model_plot_50_year, vp = viewport(layout.pos.row = 66:77, layout.pos.col = 51:99))

grid.text(x = unit(0.02, "npc"), y = unit(0.98, "npc"), label = "A", gp = gpar(fontsize = base_text_size + 2, fontface = "bold", col = black_text_colour))
grid.text(x = unit(0.51, "npc"), y = unit(0.98, "npc"), label = "B", gp = gpar(fontsize = base_text_size + 2, fontface = "bold", col = black_text_colour))
grid.text(x = unit(0.25, "npc"), y = unit(0.955, "npc"), label = "10 year survival\nCox proportional-hazards model", gp = gpar(fontsize = base_text_size, fontface = "plain", col = black_text_colour))
grid.text(x = unit(0.25, "npc"), y = unit(0.59, "npc"), label = "Cox model for adjusted SoM", gp = gpar(fontsize = base_text_size, fontface = "plain", col = black_text_colour))
grid.text(x = unit(0.75, "npc"), y = unit(0.955, "npc"), label = "7 year survival, age > 50\nCox proportional-hazards model", gp = gpar(fontsize = base_text_size, fontface = "plain", col = black_text_colour))
grid.text(x = unit(0.75, "npc"), y = unit(0.59, "npc"), label = "Cox model for adjusted SoM", gp = gpar(fontsize = base_text_size, fontface = "plain", col = black_text_colour))

print(cox_model_plot_detrimental_protective, vp = viewport(layout.pos.row = 94:150, layout.pos.col = 3:47))
print(cox_model_plot_ihm_crp_som, vp = viewport(layout.pos.row = 88:150, layout.pos.col = 53:97))

grid.text(x = unit(0.02, "npc"), y = unit(0.425, "npc"), label = "C", gp = gpar(fontsize = base_text_size + 2, fontface = "bold", col = black_text_colour))
grid.text(x = unit(0.52, "npc"), y = unit(0.425, "npc"), label = "D", gp = gpar(fontsize = base_text_size + 2, fontface = "bold", col = black_text_colour))
grid.text(x = unit(0.25, "npc"), y = unit(0.39, "npc"), label = "Cox model (detrimental & protective score)", gp = gpar(fontsize = base_text_size, fontface = "plain", col = black_text_colour))
grid.text(x = unit(0.75, "npc"), y = unit(0.415, "npc"), label = "Cox model (IHM, SoM, CRP)", gp = gpar(fontsize = base_text_size, fontface = "plain", col = black_text_colour))

dev.off()
###############################################################

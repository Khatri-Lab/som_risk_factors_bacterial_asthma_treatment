###############################################################
library(readr)
library(forestploter)
source("/labs/khatrilab/ananthg/detrimental_host_response/common_functions/perform_meta_analyses.R")
source("/labs/khatrilab/ananthg/detrimental_host_response/common_functions/plot_themes.R")
source("/labs/khatrilab/ananthg/detrimental_host_response/common_functions/JT_test_function.R")
main_fig_summary <- read_csv("/labs/khatrilab/ananthg/detrimental_host_response/data_2024_09_27/asthma_biologics_meta_analysis_main_fig_summary.csv")
main_fig_summary[main_fig_summary == "replace"] <- ""
###############################################################


###############################################################
asthma_meta_obj <- readRDS("/labs/khatrilab/ananthg/detrimental_host_response/data_2024_09_27/asthma_GEO_datasets.RDS")
GSE69683 <- asthma_meta_obj[[7]]
GSE69683_with_signatures <- keep_signature_mod_genes(GSE69683)
GSE69683_score_pheno <- as.data.frame(t(GSE69683_with_signatures$expr[c("module_1_score", "module_2_score", "module_3_score", "module_4_score", "detrimental_score", "protective_score", "som_score"), ]))
GSE69683_score_pheno$pheno <- "Severe asthma (n = 334)"
GSE69683_score_pheno$pheno[GSE69683_with_signatures$pheno$asthma_status == "Healthy, non-smokers"] <- "Healthy (n = 87)"
GSE69683_score_pheno$pheno[GSE69683_with_signatures$pheno$asthma_status == "Moderate asthma, non-smokers"] <- "Moderate asthma (n = 77)"
GSE69683_score_pheno$pheno <- ordered(GSE69683_score_pheno$pheno, levels = c("Healthy (n = 87)", "Moderate asthma (n = 77)", "Severe asthma (n = 334)"))
###############################################################


###############################################################
plot_JT_test <- function(score_pheno_frame, score_name, score_name_plot)
{
  score_pheno_frame[[score_name]] <- scale(score_pheno_frame[[score_name]])
  severity_colours <- pheno_colors <- c("Healthy (n = 87)" = "#678652", "Moderate asthma (n = 77)" = "#AFAADC", "Severe asthma (n = 334)" = "#D9A9A4")
  jt_stat <- JT.test(data = score_pheno_frame[, score_name],
    class = score_pheno_frame$pheno,
    labs = c("Healthy (n = 87)", "Moderate asthma (n = 77)", "Severe asthma (n = 334)"))
  p_value_print <- paste0("p(JT) = ", format(jt_stat$p.value, digits = 2))
  if(jt_stat$p.value == 0)
  {
    p_value_print <- "p(JT) < 2.2e-16"
  }
  score_summary <- ggplot(score_pheno_frame, aes_string(x = "pheno", y = score_name, colour = "pheno")) +
    ylab(score_name_plot) +
    geom_beeswarm(cex = cex_val * 1.1, size = point_size) +
    geom_boxplot(fill = NA, lwd = line_size_main_mm, width = 0.4, outlier.shape = NA, colour = black_text_colour) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank()) +
    annotate(geom = "text", x = Inf, y = max(score_pheno_frame[, score_name]),
      label = p_value_print,
      size = text_size_labels, hjust = 1.2, vjust = -0.5) +
    scale_colour_manual(values = severity_colours) +
    labs(colour = NULL) + guides(colour = guide_legend(ncol = 1, override.aes = list(size = text_size_labels * 6 / 9)))
  return(score_summary)
}
GSE69683_som_plot <- plot_JT_test(GSE69683_score_pheno, "som_score", "SoM")
###############################################################


###############################################################
GSE96530 <- readRDS("/labs/khatrilab/ananthg/detrimental_host_response/data_2024_09_27/asthma_GEO_GSE96530.RDS")
GSE96530_with_signatures <- keep_signature_mod_genes(GSE96530)
GSE96530_score_pheno <- as.data.frame(t(GSE96530_with_signatures$expr[c("module_1_score", "module_2_score", "module_3_score", "module_4_score", "detrimental_score", "protective_score", "som_score"), ]))
GSE96530_score_pheno$pheno <- "Asthma exacerbation (n = 19)"
GSE96530_score_pheno$pheno[GSE96530_with_signatures$pheno$visit == "Convalescent"] <- "Convalescence (n = 19)"
GSE96530_score_pheno$subject_id <- GSE96530_with_signatures$pheno$subject_id
GSE96530_score_pheno$pheno <- ordered(GSE96530_score_pheno$pheno, levels = c("Asthma exacerbation (n = 19)", "Convalescence (n = 19)"))
###############################################################


###############################################################
plot_wilcox_test <- function(score_pheno_frame, score_name, score_name_plot, direction)
{
  score_pheno_frame$score_to_plot <- scale(score_pheno_frame[[score_name]])
  severity_colours <- c("Asthma exacerbation (n = 19)" = "#44BB99", "Convalescence (n = 19)" = "#BBCC33")
  test_significance <- score_pheno_frame %>%
    wilcox_test(score_to_plot ~ pheno, paired = TRUE, comparisons = list(c("Asthma exacerbation (n = 19)", "Convalescence (n = 19)"))) %>%
    add_xy_position(x = "pheno", y = max)
  test_significance$y.position <- 1.03 * test_significance$y.position
  test_significance$xmin <- 1.02 * test_significance$xmin
  test_significance$xmax <- 0.98 * test_significance$xmax
  test_significance$p_print <- as.character(ifelse(test_significance$p < 0.01, format(test_significance$p, digits = 1, scientific = TRUE), round(test_significance$p, 2)))
  score_summary <- ggplot(score_pheno_frame, aes(x = pheno, y = score_to_plot, colour = pheno)) +
    ylab(score_name_plot) +
    geom_beeswarm(cex = cex_val * 1.1, size = point_size) +
    geom_line(aes(group = subject_id), position = position_beeswarm(cex = cex_val * 1.1), colour = black_line_colour, alpha = 0.2, lwd = line_size_text_repel_mm) +
    geom_boxplot(fill = NA, lwd = line_size_main_mm, width = 0.4, outlier.shape = NA, colour = black_text_colour) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank()) +
    stat_pvalue_manual(test_significance, label = "p_print",
      tip.length = 0,
      size = text_size_labels, colour = black_text_colour, bracket.size = line_size_text_repel_mm) +
    scale_colour_manual(values = severity_colours) +
    labs(colour = NULL) + guides(colour = guide_legend(ncol = 1, override.aes = list(size = text_size_labels * 6 / 9)))
  return(score_summary)
}
GSE96530_som_plot <- plot_wilcox_test(GSE96530_score_pheno, "som_score", "SoM", "increasing")
###############################################################


###############################################################
num_rows <- 3# GSE134544 and GSE148725, plus one empty before summary
base_gp <- gpar(fontsize = base_text_size, fontfamily = base_font_family, fontface = "plain", lwd = line_size_main, cex = 1, col = black_line_colour)
forest_theme <- forest_theme(
  base_size = base_text_size,
  base_family = base_font_family,
  ci_pch = 15,
  ci_col = black_text_colour,
  ci_alpha = 1,
  ci_fill = "#FF6718",
  ci_lty = 1,
  ci_lwd = line_size_main,
  ci_Theight = 0.2,
  legend_position = "none",
  xaxis_gp = gpar(fontsize = base_text_size - 1, fontfamily = base_font_family, fontface = "plain", lwd = line_size_main, cex = 1, col = black_text_colour),
  refline_gp = gpar(fontsize = base_text_size - 1, fontfamily = base_font_family, fontface = "plain", lwd = line_size_main, cex = 1, col = black_line_colour, lty = "dashed"),
  vertline_lwd = line_size_text_repel,
  vertline_col = black_line_colour,
  vertline_lty = "dashed",
  summary_col = black_text_colour,
  summary_fill = "#FF6718",
  summary_lwd = line_size_text_repel,
  footnote_gp = gpar(cex = 0.9, fontface = "plain", col = black_text_colour),
  title_gp = gpar(fontsize = base_text_size, fontfamily = base_font_family, fontface = "plain", lwd = line_size_main, cex = 1, col = black_text_colour),
  title_just = "center",
  arrow_gp = gpar(fontsize = base_text_size - 1, fontfamily = base_font_family, fontface = "plain", lwd = line_size_text_repel, fill = black_line_colour, col = black_line_colour),
  cols = c("Drug", "NR", "R", "SoM", "", "Detrimental", "", "Protective"),
  core = list(fg_params = list(fontsize = base_text_size - 1, fontface = "plain", hjust = 0.5, x = 0.5),
    bg_params = list(fill = "transparent")),
  colhead = list(fg_params = list(fontsize = base_text_size, fontface = "plain", hjust = 0.5, x = 0.5),
    bg_params = list(fill = "transparent")))
###############################################################


###############################################################
add_es_text_change_color <- function(forest_plot, column_no, eff_size, eff_size_lower, eff_size_upper)
{
  eff_size_z_factor <- abs(as.numeric(eff_size) * 3.92 / (as.numeric(eff_size_upper) - as.numeric(eff_size_lower)))
  eff_size_p <- exp(0 - 0.717 * eff_size_z_factor - 0.416 * eff_size_z_factor^2)
  just = "right"
  if(as.numeric(eff_size) < 0)
  {
    forest_plot <- forestploter::edit_plot(forest_plot, col = column_no, row = num_rows + 1, which = "ci", gp = gpar(fill = "#3499FF"))
    forest_plot <- forestploter::edit_plot(forest_plot, col = column_no, row = c(1 : num_rows), which = "ci", gp = gpar(fill = "#3499FF"))
    just = "left"
  }
  forest_plot <- add_text(forest_plot,
    paste0("ES=", round(as.numeric(eff_size), 2), "\np = ", round(eff_size_p, 2)),
    row = num_rows,
    col = column_no,
    part = "body",
    just = just,
    gp = gpar(fontsize = base_text_size - 1, fontfamily = base_font_family, col = black_text_colour)
  )
  return(forest_plot)
}
###############################################################


###############################################################
main_fig_summary$Dataset[main_fig_summary$Dataset == "Omalizumab"] <- "Omalizumab (GSE134544)"
main_fig_summary$Dataset[main_fig_summary$Dataset == "Benralizumab"] <- "Benralizumab (GSE148725)"
colnames(main_fig_summary)[1] <- "Monoclonal antibody"
colnames(main_fig_summary)[16] <- ""
colnames(main_fig_summary)[23] <- ""
forest_plot <- forestploter::forest(main_fig_summary[, c(1:3, 16, 4)],
  est = as.numeric(main_fig_summary$eff_size_1),
  lower = as.numeric(main_fig_summary$lower_1),
  upper = as.numeric(main_fig_summary$upper_1),
  sizes = as.numeric(main_fig_summary$box_size_1),
  is_summary = c(rep(FALSE, num_rows), TRUE),
  ci_column = 5,
  ticks_at = c(0, 0.5, 1),
  ref_line = 0,
  xlab = "Effect size",
  theme = forest_theme)
forest_plot <- forestploter::edit_plot(forest_plot, col = 1, which = "text", part = "header", hjust = 0, x = unit(3, "pt"))
forest_plot <- forestploter::edit_plot(forest_plot, col = 1, row = c(0 : num_rows + 1), which = "text", part = "body", hjust = 0, x = unit(3, "pt"))
convertHeight(forest_plot$heights, "pt", valueOnly = TRUE)
forest_plot$heights <- unit.c(unit(3, "pt"),
  #unit(base_text_size * 1.5, "pt"),
  unit(base_text_size * 1.3, "pt"),
  unit(base_text_size * 1.3, "pt"),
  unit(base_text_size * 1.3, "pt"),
  unit(base_text_size * 1.8, "pt"),
  unit(base_text_size * 1.3, "pt"),
  unit(base_text_size * 1.3, "pt"),
  unit(base_text_size * 1.6, "pt")
  )
convertWidth(forest_plot$widths, "pt", valueOnly = TRUE)
forest_plot$widths <- unit.c(unit(3, "pt"),
  unit(base_text_size * 11.5, "pt"),
  unit(base_text_size * 1.7, "pt"),
  unit(base_text_size * 1.7, "pt"),
  unit(base_text_size * 1.4, "pt"),
  unit(base_text_size * 5.5, "pt"),
  unit(3, "pt")
  )
forest_plot <- add_es_text_change_color(forest_plot, 5, main_fig_summary$eff_size_1[num_rows + 1], main_fig_summary$lower_1[num_rows + 1], main_fig_summary$upper_1[num_rows + 1])
###############################################################


###############################################################
cairo_pdf("/labs/khatrilab/ananthg/detrimental_host_response/figures_2024_09_27/fig_ats_combined_asthma.pdf", width = 4, height = 5, bg = "transparent")

pushViewport(viewport(layout = grid.layout(nrow = 105, ncol = 100)))

print(GSE69683_som_plot, vp = viewport(layout.pos.row = 4:31, layout.pos.col = 1:100))
print(GSE96530_som_plot, vp = viewport(layout.pos.row = 37:65, layout.pos.col = 1:100))
print(as_ggplot(forest_plot), vp = viewport(layout.pos.row = 73:105, layout.pos.col = 1:100))

grid.rect(x = unit(0.5, "npc"), y = unit(0.985, "npc"), width = unit(0.99, "npc"), height = unit(0.035, "npc"), gp = gpar(fill = "#EFEFEF", col = "transparent"))
grid.rect(x = unit(0.5, "npc"), y = unit(0.68, "npc"), width = unit(0.99, "npc"), height = unit(0.035, "npc"), gp = gpar(fill = "#EFEFEF", col = "transparent"))
grid.rect(x = unit(0.5, "npc"), y = unit(0.334, "npc"), width = unit(0.99, "npc"), height = unit(0.075, "npc"), gp = gpar(fill = "#EFEFEF", col = "transparent"))

grid.text(x = unit(0.5, "npc"), y = unit(0.985, "npc"), label = "Asthma (baseline, GSE69683)", gp = gpar(fontsize = base_text_size, fontface = "plain", col = black_text_colour))
grid.text(x = unit(0.5, "npc"), y = unit(0.68, "npc"), label = "Asthma (exacerbation, GSE96530)", gp = gpar(fontsize = base_text_size, fontface = "plain", col = black_text_colour))
grid.text(x = unit(0.5, "npc"), y = unit(0.334, "npc"), label = "Non-responders (NR) vs. responders (R) to\nmonoclonal antibody treatment in asthma", gp = gpar(fontsize = base_text_size, fontface = "plain", col = black_text_colour))

grid.text(x = unit(0.0175, "npc"), y = unit(0.95, "npc"), label = "A", gp = gpar(fontsize = base_text_size + 2, fontface = "bold", col = black_text_colour))
grid.text(x = unit(0.0175, "npc"), y = unit(0.645, "npc"), label = "B", gp = gpar(fontsize = base_text_size + 2, fontface = "bold", col = black_text_colour))
grid.text(x = unit(0.0175, "npc"), y = unit(0.28, "npc"), label = "C", gp = gpar(fontsize = base_text_size + 2, fontface = "bold", col = black_text_colour))

dev.off()
###############################################################

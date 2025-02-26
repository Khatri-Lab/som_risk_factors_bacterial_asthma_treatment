###############################################################
library(tidyverse)
library(data.table)
library(pROC)
library(interactions)
library(ggsurvfit)
library(metafor)
library(effsize)
library(rmeta)
library(forestploter)
library(png)
source("/labs/khatrilab/ananthg/detrimental_host_response/common_functions/plot_themes.R")
source("/labs/khatrilab/ananthg/detrimental_host_response/common_functions/JT_test_function.R")
###############################################################


###############################################################
num_rows <- 3# Burn and Vanish, plus one empty before summary
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
  core = list(fg_params = list(fontsize = base_text_size - 1, fontface = "plain", hjust = 0.5, x = 0.5),
    bg_params = list(fill = "transparent")),
  colhead = list(fg_params = list(fontsize = base_text_size, fontface = "plain", hjust = 0.5, x = 0.5),
    bg_params = list(fill = "transparent")))
###############################################################


###############################################################
severity_cols <- c("0" = "#466D2D", "1" = "#93C572", "2" = "#7B72C5", "3" = "#72BCC5", "4" = "#E1B941", "5" = "#C57B72", "6" = "#6D2D46")
severity_labels <- c("0" = "Healthy", "1" = "Asymptomatic", "2" = "Mild", "3" = "Moderate", "4" = "Serious", "5" = "Critical", "6" = "Fatal")
treatment_colours <- c(Placebo = "#CDCDCD", Hydrocortisone = "#DC3023", Anakinra = "#4B5CC4")
###############################################################


###############################################################
COCO.output.master <- readRDS(file = "/labs/khatrilab/armoore7/sepsis/conserved_host_response_figs/coco_bact_scores.RDS")
COCO.output.master$Severity_Category <- NA
COCO.output.master$Severity_Category[COCO.output.master$Severity_Scaled == 0] <- "Healthy"
COCO.output.master$Severity_Category[COCO.output.master$Severity_Scaled == 2 | COCO.output.master$Severity_Scaled == 3] <- "Non-Severe"
COCO.output.master$Severity_Category[COCO.output.master$Severity_Scaled == 4 | COCO.output.master$Severity_Scaled == 5 | COCO.output.master$Severity_Scaled == 6] <- "Severe"
df_bact <- COCO.output.master[COCO.output.master$Etiology == "Bacterial" | COCO.output.master$Etiology == "Healthy", ] %>% filter(!is.na(name))
df_bact$Severity_Scaled <- ordered(as.factor(df_bact$Severity_Scaled), levels = c("0", "1", "2", "3", "4", "5", "6"))
###############################################################


###############################################################
plot_cor_test <- function(score_pheno_frame, score_name, score_name_plot)
{
  score_pheno_frame$scaled_score <- scale(score_pheno_frame[[score_name]])
  score_summary <- ggplot(score_pheno_frame, aes(x = Severity_Scaled, y = scaled_score, colour = Severity_Scaled)) +
    ylab(score_name_plot) +
    geom_beeswarm(cex = cex_val * 0.5, size = point_size) +
    geom_boxplot(fill = NA, lwd = line_size_main_mm, width = 0.4, outlier.shape = NA, colour = black_text_colour) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank(),
          legend.position = "right", legend.direction = "horizontal",
          panel.grid.major.y = element_line(color = "grey", size = line_size_text_repel_mm)) +
    stat_cor(aes(x = as.numeric(as.character(Severity_Scaled)), y = scaled_score), method = "pearson", label.x = "0", label.y = Inf, vjust = 1.5, size = text_size_labels, family = base_font_family, colour = black_text_colour) +
    scale_colour_manual(values = severity_cols, labels = severity_labels) +
    labs(colour = NULL) + guides(colour = guide_legend(ncol = 1, override.aes = list(size = text_size_labels * 6 / 9)))
  return(score_summary)
}
bacterial_detrimental_plot <- plot_cor_test(df_bact, "detrimental_score", "Detrimental (z score)")
bacterial_module_4_plot <- plot_cor_test(df_bact, "mod4_score", "Module 4 (z score)")
bacterial_protective_plot <- plot_cor_test(df_bact, "protective_score", "Protective (z score)")
bacterial_som_plot <- plot_cor_test(df_bact, "som_score", "SoM (z score)")
bacterial_plot_legend <- as_ggplot(ggpubr::get_legend(bacterial_som_plot))
###############################################################


###############################################################
roc_df <- df_bact %>% mutate(NSvS = ifelse(Severity_Category == "Severe", 1, ifelse(Severity_Category == "Non-Severe", 0, NA)))
rocobj <- roc(roc_df$NSvS, roc_df$som_score)
auc_ci <- as.numeric(ci.auc(rocobj))
auc_plot_text <- paste0("AUC = ", format(auc_ci[2], digits = 2), " (", format(auc_ci[1], digits = 2), " - ", format(auc_ci[3], digits = 2), ")")
ci <- ci.se(rocobj)
dat.ci <- data.frame(
  x = 1 - as.numeric(rownames(ci)), 
  lower = ci[, 1], 
  upper = ci[, 3]
)
SoM_ROC <- ggroc(rocobj, color = black_line_colour, legacy.axes = TRUE) +
  ggtitle("Severe vs. non-severe") +
  xlab("1 - Specificity") + ylab("Sensitivity") +
  annotate("segment", x = 0, xend = 1, y = 0, yend = 1, color = "darkgrey", linetype = "dashed", linewidth = line_size_main_mm) +
  geom_ribbon(data = dat.ci, aes(x = x, ymin = lower, ymax = upper), fill = black_line_colour, alpha = 0.2) +
  theme(panel.grid.major.y = element_line(color = "grey", size = line_size_text_repel_mm),
        panel.grid.major.x = element_line(color = "grey", size = line_size_text_repel_mm))+
  scale_y_continuous(breaks = seq(0, 1, 0.2)) +
  scale_x_continuous(breaks = seq(0, 1, 0.2)) +
  annotate("rect", xmin = 0.11, xmax = 0.99, ymin = 0.07, ymax = 0.17, fill = "#FFFFFF", alpha = 0.5, colour = NA) +
  annotate("text", x = .55, y = .12, label = auc_plot_text, size = text_size_labels)
###############################################################


###############################################################
plot_barplot <- function(score_pheno_data, score_column_name, treatment_column_name, mortality_column_name, mortality_state_name, score_name_title, plot_title, return_odds_table = FALSE)
{
  score_pheno_data$discretized_score <- score_pheno_data[, score_column_name]
  score_pheno_data$treatment_column <- score_pheno_data[, treatment_column_name]
  score_pheno_data$mortality_column <- score_pheno_data[, mortality_column_name]
  fisher_high <- fisher.test(score_pheno_data$treatment_column[score_pheno_data$discretized_score == "\U2265 median"], score_pheno_data$mortality_column[score_pheno_data$discretized_score == "\U2265 median"])
  fisher_high_p <- paste0("p = ", round(fisher_high$p.value, 3))
  fisher_low <- fisher.test(score_pheno_data$treatment_column[score_pheno_data$discretized_score == "< median"], score_pheno_data$mortality_column[score_pheno_data$discretized_score == "< median"])
  fisher_low_p <- paste0("p = ", round(fisher_low$p.value, 3))
  prop_df <- score_pheno_data %>% 
    group_by(discretized_score, treatment_column) %>% 
    summarize(mort_prop = mean(mortality_column == mortality_state_name, na.rm = TRUE) * 100)
  p_text_space <- 0.07 * max(prop_df$mort_prop)
  protective_bar_plot <- ggplot(prop_df, aes(x = discretized_score, y = mort_prop, fill = treatment_column)) +
    geom_bar(position = "dodge", stat = "identity", width = 0.8) +
    labs(x = score_name_title, y = "Mortality") +
    scale_fill_manual(values = treatment_colours) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
    theme(legend.position = "bottom", legend.margin = margin(1, 1, 1, 1), legend.box.margin = margin(-5, 5, 3, 3)) +
    labs(fill = NULL) + guides(fill = guide_legend(nrow = 1, override.aes = list(size = text_size_labels * 6 / 9))) +
    scale_y_continuous(labels = function(x) paste0(x, '%')) +
    annotate("text", x = 2, y = p_text_space + max(prop_df$mort_prop[prop_df$discretized_score == "\U2265 median"]), label = fisher_high_p, size = text_size_labels) +
    annotate("text", x = 1, y = p_text_space + max(prop_df$mort_prop[prop_df$discretized_score == "< median"]), label = fisher_low_p, size = text_size_labels)
  if(return_odds_table)
  {
    prop_df_odds <- as.data.frame(score_pheno_data[score_pheno_data$discretized_score == "\U2265 median", ] %>% 
      group_by(treatment_column) %>% 
      summarize(alive = as.integer(sum(mortality_column != mortality_state_name)), dead = as.integer(sum(mortality_column == mortality_state_name))))
    rownames(prop_df_odds) <- prop_df_odds$treatment_column
    prop_df_odds <- subset(prop_df_odds, select = -c(treatment_column))
    return(prop_df_odds)
  }
  return(protective_bar_plot)
}
###############################################################


###############################################################
burn_trial_joined <- readRDS(file = "/labs/khatrilab/armoore7/burns/burn_shock_steroids.RDS") %>%
  mutate(SoM_score = detrimental_score - protective_score, day = ifelse(is.na(time), NA, paste0("D", time/24))) %>%
  dplyr::rename(mod1_score = clust1_score, mod2_score = clust2_score, mod3_score = clust3_score, mod4_score = clust4_score)
burn_trial_joined$treatment.ch1 <- ordered(burn_trial_joined$treatment.ch1, levels = c("Placebo", "Hydrocortisone"))

burn_trial_baseline <- burn_trial_joined %>%
  filter(sample.type.ch1 == "S1") %>%
  mutate(high_som = ifelse(SoM_score >= median(SoM_score, na.rm = TRUE), "\U2265 median", "< median"),
         high_protective = ifelse(protective_score >= median(protective_score, na.rm = TRUE), "\U2265 median", "< median"),
         high_detrimental = ifelse(detrimental_score >= median(detrimental_score, na.rm = TRUE), "\U2265 median", "< median"))
burn_trial_baseline$high_protective <- ordered(burn_trial_baseline$high_protective, levels = c("< median", "\U2265 median"))
###############################################################


###############################################################
burn_trial_baseline$survival_text <- ifelse(burn_trial_baseline$survival..d28..ch1 == "Non survivor", "Non-survivor", "Survivor")
test_significance <- burn_trial_baseline %>%
  wilcox_test(SoM_score ~ survival_text, paired = FALSE, comparisons = list(c("Non-survivor", "Survivor"))) %>%
  add_y_position(y = max)
test_significance$y.position <- 1.03 * test_significance$y.position
test_significance$xmin <- 1.02
test_significance$xmax <- 1.98
test_significance$p_print <- as.character(ifelse(test_significance$p < 0.01, format(test_significance$p, digits = 1, scientific = TRUE), round(test_significance$p, 2)))

burn_som_baseline <- ggplot(burn_trial_baseline, aes(x = survival_text, y = SoM_score))+
  geom_beeswarm(cex = cex_val, size = point_size, colour = black_text_colour) +
  geom_boxplot(lwd = line_size_main_mm, width = 0.4, outlier.shape = NA, colour = black_text_colour, fill = NA) +
  scale_y_continuous(expand = expansion(mult = c(0.025, 0.1))) +
  stat_pvalue_manual(test_significance, label = "p_print",
      tip.length = 0,
      size = text_size_labels, colour = black_text_colour, bracket.size = line_size_text_repel_mm) +
  ylab("SoM score") + ggtitle("Burn cohort (baseline)") +
  theme(legend.position = "none", axis.title.x = element_blank())
###############################################################


###############################################################
change_plot_data <- as.data.frame(burn_trial_joined %>% filter(day == "D0" | day == "D7"))
change_plot_data$day[change_plot_data$day == "D0"] <- "Day 0"
change_plot_data$day[change_plot_data$day == "D7"] <- "Day 7"
plot_change_plot <- function(change_plot_data, score_name, plot_title, round_points_0 = 2, round_points_7 = 2)
{
  change_plot_data$score_to_plot <- scale(change_plot_data[[score_name]])
  wilcox_d0 <- wilcox.test(change_plot_data$score_to_plot[change_plot_data$`treatment.ch1` == "Hydrocortisone" & change_plot_data$day == "Day 0"], change_plot_data$score_to_plot[change_plot_data$`treatment.ch1` == "Placebo" & change_plot_data$day == "Day 0"])
  wilcox_d0_p <- paste0("p = ", round(wilcox_d0$p.value, round_points_0))
  wilcox_d7 <- wilcox.test(change_plot_data$score_to_plot[change_plot_data$`treatment.ch1` == "Hydrocortisone" & change_plot_data$day == "Day 7"], change_plot_data$score_to_plot[change_plot_data$`treatment.ch1` == "Placebo" & change_plot_data$day == "Day 7"])
  wilcox_d7_p <- paste0("p = ", round(wilcox_d7$p.value, round_points_7))
  change_plot <- ggplot(change_plot_data, aes(x = day, y = score_to_plot, fill = treatment.ch1))+
    geom_boxplot(lwd = line_size_main_mm, width = 0.4, outlier.shape = NA, colour = black_text_colour) +
    scale_fill_manual(values = treatment_colours) +
    scale_y_continuous(expand = expansion(mult = c(0.025, 0.1))) +
    ylab(plot_title) +
    theme(legend.position = "none", axis.title.x = element_blank()) +
    annotate("text", x = 1, y = 1.1 * max(change_plot_data$score_to_plot[change_plot_data$day == "Day 0"]), label = wilcox_d0_p, size = text_size_labels) +
    annotate("text", x = 2, y = 1.1 * max(change_plot_data$score_to_plot[change_plot_data$day == "Day 7"]), label = wilcox_d7_p, size = text_size_labels)
  return(change_plot)
}
protective_change_plot <- plot_change_plot(change_plot_data, "protective_score", "Protective (z score)", round_points_7 = 3)
detrimental_change_plot <- plot_change_plot(change_plot_data, "detrimental_score", "Detrimental (z score)")
###############################################################


###############################################################
burns_protective_plot <- plot_barplot(burn_trial_baseline, "high_protective", "treatment.ch1", "survival..d28..ch1", "Non survivor", "Protective score", "Burn") + theme(legend.position = "none") + ggtitle("Burn")
###############################################################


###############################################################
vanish_joined <- readRDS(file = "/labs/khatrilab/armoore7/sepsis/vanish_joined.RDS") %>%
  mutate(SoM_score = detrimental_score - protective_score) %>%
  dplyr::rename(mod1_score = clust1_score, mod2_score = clust2_score, mod3_score = clust3_score, mod4_score = clust4_score) %>%
  mutate(high_som = ifelse(SoM_score >= median(SoM_score, na.rm = TRUE), "\U2265 median", "< median"),
         high_protective = ifelse(protective_score >= median(protective_score, na.rm = TRUE), "\U2265 median", "< median"),
         high_detrimental = ifelse(detrimental_score >= median(detrimental_score, na.rm = TRUE), "\U2265 median", "< median"))
vanish_joined$steroids_yn_words[vanish_joined$steroids_yn_words == "No Hydrocortisone"] <- "Placebo"
vanish_joined$high_protective <- ordered(vanish_joined$high_protective, levels = c("< median", "\U2265 median"))
vanish_joined$steroids_yn_words <- ordered(vanish_joined$steroids_yn_words, levels = c("Placebo", "Hydrocortisone"))
###############################################################


###############################################################
vanish_joined$survival_text <- ifelse(vanish_joined$Characteristics.outcome.day.28. == "Dead", "Non-survivor", "Survivor")
test_significance <- vanish_joined %>%
  wilcox_test(SoM_score ~ survival_text, paired = FALSE, comparisons = list(c("Non-survivor", "Survivor"))) %>%
  add_y_position(y = max)
test_significance$y.position <- 1.03 * test_significance$y.position
test_significance$xmin <- 1.02
test_significance$xmax <- 1.98
test_significance$p_print <- as.character(ifelse(test_significance$p < 0.01, format(test_significance$p, digits = 1, scientific = TRUE), round(test_significance$p, 2)))

vanish_som <- ggplot(vanish_joined, aes(x = survival_text, y = SoM_score))+
  geom_beeswarm(cex = cex_val, size = point_size, colour = black_text_colour) +
  geom_boxplot(lwd = line_size_main_mm, width = 0.4, outlier.shape = NA, colour = black_text_colour, fill = NA) +
  scale_y_continuous(expand = expansion(mult = c(0.025, 0.1))) +
  stat_pvalue_manual(test_significance, label = "p_print",
      tip.length = 0,
      size = text_size_labels, colour = black_text_colour, bracket.size = line_size_text_repel_mm) +
  ylab("SoM score") + ggtitle("Vanish cohort") +
  theme(legend.position = "none", axis.title.x = element_blank())
###############################################################


###############################################################
vanish_randomized <- readRDS(file = "/labs/khatrilab/armoore7/sepsis/vanish_joined.RDS") %>%
  mutate(SoM_score = detrimental_score - protective_score) %>%
  dplyr::rename(mod1_score = clust1_score, mod2_score = clust2_score, mod3_score = clust3_score, mod4_score = clust4_score) %>%
  filter(Characteristics.drug2.per.protocol. != "none")%>%
  mutate(high_som = ifelse(SoM_score >= median(SoM_score, na.rm = TRUE), "\U2265 median", "< median"),
         high_protective = ifelse(protective_score >= median(protective_score, na.rm = TRUE), "\U2265 median", "< median"),
         high_detrimental = ifelse(detrimental_score >= median(detrimental_score, na.rm = TRUE), "\U2265 median", "< median"))
vanish_randomized$steroids_yn_words[vanish_randomized$steroids_yn_words == "No Hydrocortisone"] <- "Placebo"
vanish_randomized$high_protective <- ordered(vanish_randomized$high_protective, levels = c("< median", "\U2265 median"))
vanish_randomized$steroids_yn_words <- ordered(vanish_randomized$steroids_yn_words, levels = c("Placebo", "Hydrocortisone"))
###############################################################


###############################################################
vanish_protective_plot <- plot_barplot(vanish_joined, "high_protective", "steroids_yn_words", "Characteristics.outcome.day.28.", "Dead", "Protective score", "Vanish") + ggtitle("VANISH")
vanish_detrimental_plot <- plot_barplot(vanish_joined, "high_detrimental", "steroids_yn_words", "Characteristics.outcome.day.28.", "Dead", "Detrimental score", "Vanish")
vanish_randomized_protective_plot <- plot_barplot(vanish_randomized, "high_protective", "steroids_yn_words", "Characteristics.outcome.day.28.", "Dead", "Protective score", "Vanish") + theme(legend.position = "none")
hydrocortisone_legend <- as_ggplot(ggpubr::get_legend(vanish_protective_plot))
vanish_protective_plot <- vanish_protective_plot + theme(legend.position = "none")
###############################################################


###############################################################
burns_protective_odds_frame <- plot_barplot(burn_trial_baseline, "high_protective", "treatment.ch1", "survival..d28..ch1", "Non survivor", "Protective score", "Burn", TRUE)
vanish_protective_odds_frame <- plot_barplot(vanish_joined, "high_protective", "steroids_yn_words", "Characteristics.outcome.day.28.", "Dead", "Protective score", "Vanish", TRUE)
combined_protective_odds_frame <- matrix(c(burns_protective_odds_frame[1, 1], burns_protective_odds_frame[1, 2], burns_protective_odds_frame[2, 1], burns_protective_odds_frame[2, 2], vanish_protective_odds_frame[1, 1], vanish_protective_odds_frame[1, 2], vanish_protective_odds_frame[2, 1], vanish_protective_odds_frame[2, 2]), nrow = 2, byrow = TRUE)
colnames(combined_protective_odds_frame) <- c("alive_placebo", "dead_placebo", "alive_hydrocortisone", "dead_hydrocortisone")
combined_protective_odds_ratio <- escalc(measure = "OR", ai = alive_placebo, bi = dead_placebo, ci = alive_hydrocortisone, di = dead_hydrocortisone, data = combined_protective_odds_frame)
combined_protective_odds_ratio$se <- sqrt(combined_protective_odds_ratio$vi)
combined_protective_odds_ratio$box_size <- 1 / combined_protective_odds_ratio$vi
combined_protective_odds_ratio$lower <- exp(combined_protective_odds_ratio$yi - 1.96 * combined_protective_odds_ratio$se)
combined_protective_odds_ratio$upper <- exp(combined_protective_odds_ratio$yi + 1.96 * combined_protective_odds_ratio$se)
combined_protective_odds_ratio$odds_ratio <- exp(combined_protective_odds_ratio$yi)
meta_protective_odds_ratio <- rma(yi, vi, data = combined_protective_odds_ratio)
meta_protective_odds_ratio_results <- predict(meta_protective_odds_ratio, transf = exp, digits = 5)
odds_ratio_plot_frame <- as.data.frame(combined_protective_odds_ratio[, c("odds_ratio", "lower", "upper", "box_size")])
odds_ratio_plot_frame[3, ] <- c(NA, NA, NA, NA)
odds_ratio_plot_frame[4, ] <- c(meta_protective_odds_ratio_results$pred, meta_protective_odds_ratio_results$ci.lb, meta_protective_odds_ratio_results$ci.ub, 1 / meta_protective_odds_ratio$se ^ 2)
odds_ratio_plot_frame$Cohort <- c("Burn", "VANISH", "", "Summary")
odds_ratio_plot_frame$Mortality <- c(" ", " ", " ", " ")
odds_ratio_plot_frame$box_size <- 0.25 * (odds_ratio_plot_frame$box_size - min(odds_ratio_plot_frame$box_size, na.rm = TRUE)) / (max(odds_ratio_plot_frame$box_size, na.rm = TRUE) - min(odds_ratio_plot_frame$box_size, na.rm = TRUE)) + 0.5
###############################################################


###############################################################
burn_eff_size <- cohen.d(SoM_score ~ survival_text, data = burn_trial_baseline, hedges.correction = TRUE)
vanish_eff_size <- cohen.d(SoM_score ~ survival_text, data = vanish_joined, hedges.correction = TRUE)
baseline_som_frame <- data.frame(eff_size = c(burn_eff_size$estimate, vanish_eff_size$estimate),
  lower = c(burn_eff_size$conf.int[1], vanish_eff_size$conf.int[1]),
  upper = c(burn_eff_size$conf.int[2], vanish_eff_size$conf.int[2]),
  NS = c(sum(burn_trial_baseline$survival_text == "Non-survivor"), sum(vanish_joined$survival_text == "Non-survivor")),
  S = c(sum(burn_trial_baseline$survival_text == "Survivor"), sum(vanish_joined$survival_text == "Survivor")))
colnames(baseline_som_frame) <- c("eff_size", "lower", "upper", "NS", "S")
baseline_som_frame$se <- (baseline_som_frame$upper - baseline_som_frame$lower)/3.92
baseline_som_frame$box_size <- 1 / (baseline_som_frame$se)^2
baseline_summary_eff_size <- meta.summaries(baseline_som_frame$eff_size, baseline_som_frame$se, method = "random")
baseline_som_frame <- baseline_som_frame[, c("eff_size", "lower", "upper", "NS", "S", "se", "box_size")]
baseline_som_frame[3, ] <- rep(NA, 7)
baseline_som_frame[4, ] <- c(baseline_summary_eff_size$summary, baseline_summary_eff_size$summary - 1.96 * baseline_summary_eff_size$se.summary, baseline_summary_eff_size$summary + 1.96 * baseline_summary_eff_size$se.summary, sum(baseline_som_frame$NS, na.rm = TRUE), sum(baseline_som_frame$S, na.rm = TRUE), baseline_summary_eff_size$se.summary, 1 / baseline_summary_eff_size$se.summary ^ 2)
baseline_som_frame$Cohort <- c("Burn", "VANISH", "", "Summary")
baseline_som_frame$SoM <- c(" ", " ", " ", " ")
baseline_som_frame$box_size <- 0.25 * (baseline_som_frame$box_size - min(baseline_som_frame$box_size, na.rm = TRUE)) / (max(baseline_som_frame$box_size, na.rm = TRUE) - min(baseline_som_frame$box_size, na.rm = TRUE)) + 0.5
baseline_som_frame$NS[3] <- ""
baseline_som_frame$S[3] <- ""
###############################################################


###############################################################
forest_plot <- forestploter::forest(odds_ratio_plot_frame[, c("Cohort", "Mortality")],
  est = as.numeric(log(odds_ratio_plot_frame$odds_ratio)),
  lower = as.numeric(log(odds_ratio_plot_frame$lower)),
  upper = as.numeric(log(odds_ratio_plot_frame$upper)),
  sizes = as.numeric(odds_ratio_plot_frame$box_size),
  is_summary = c(rep(FALSE, num_rows), TRUE),
  ci_column = 2,
  ref_line = 0,
  ticks_at = c(0, 2, 4),
  xlab = "Log odds ratio",
  theme = forest_theme)
forest_plot <- forestploter::edit_plot(forest_plot, col = 1, which = "text", part = "header", hjust = 0, x = unit(3, "pt"))
forest_plot <- forestploter::edit_plot(forest_plot, col = 1, row = c(0 : num_rows + 1), which = "text", part = "body", hjust = 0, x = unit(3, "pt"))
convertHeight(forest_plot$heights, "pt", valueOnly = TRUE)
forest_plot$heights <- unit.c(unit(3, "pt"),
  unit(base_text_size * 1.7, "pt"),
  unit(base_text_size * 1.4, "pt"),
  unit(base_text_size * 1.4, "pt"),
  unit(base_text_size * 2.8, "pt"),
  unit(base_text_size * 1.4, "pt"),
  unit(base_text_size * 1.4, "pt"),
  unit(3, "pt")
  )
convertWidth(forest_plot$widths, "pt", valueOnly = TRUE)
forest_plot$widths <- unit.c(unit(3, "pt"),
  unit(base_text_size * 5, "pt"),
  unit(base_text_size * 7, "pt"),
  unit(3, "pt")
  )
forest_plot <- add_text(forest_plot,
    paste0("log OR = ", format(log(odds_ratio_plot_frame$odds_ratio[4]), digits = 2), "\np = ", format(meta_protective_odds_ratio$pval, digits = 2, scientific = TRUE)),
    row = 3,
    col = 2,
    part = "body",
    just = "right",
    gp = gpar(fontsize = base_text_size - 1, fontfamily = base_font_family, col = black_text_colour)
  )
log_odds_forest_plot <- forest_plot
###############################################################


###############################################################
forest_plot <- forestploter::forest(baseline_som_frame[, c("Cohort", "NS", "S", "SoM")],
  est = as.numeric(baseline_som_frame$eff_size),
  lower = as.numeric(baseline_som_frame$lower),
  upper = as.numeric(baseline_som_frame$upper),
  sizes = as.numeric(baseline_som_frame$box_size),
  is_summary = c(rep(FALSE, num_rows), TRUE),
  ci_column = 4,
  ref_line = 0,
  ticks_at = c(0, 0.5, 1, 1.5),
  xlab = "Effect size (NS vs. S)",
  theme = forest_theme)
forest_plot <- forestploter::edit_plot(forest_plot, col = 1, which = "text", part = "header", hjust = 0, x = unit(3, "pt"))
forest_plot <- forestploter::edit_plot(forest_plot, col = 1, row = c(0 : num_rows + 1), which = "text", part = "body", hjust = 0, x = unit(3, "pt"))
convertHeight(forest_plot$heights, "pt", valueOnly = TRUE)
forest_plot$heights <- unit.c(unit(3, "pt"),
  unit(base_text_size * 1.7, "pt"),
  unit(base_text_size * 1.4, "pt"),
  unit(base_text_size * 1.4, "pt"),
  unit(base_text_size * 2.8, "pt"),
  unit(base_text_size * 1.4, "pt"),
  unit(base_text_size * 1.4, "pt"),
  unit(3, "pt")
  )
convertWidth(forest_plot$widths, "pt", valueOnly = TRUE)
forest_plot$widths <- unit.c(unit(3, "pt"),
  unit(base_text_size * 5, "pt"),
  unit(base_text_size * 2.5, "pt"),
  unit(base_text_size * 2.5, "pt"),
  unit(base_text_size * 8, "pt"),
  unit(3, "pt")
  )
baseline_som_z <- baseline_som_frame$eff_size[4] / baseline_som_frame$se[4]
baseline_som_summary_p <- exp(0 - 0.717 * baseline_som_z - 0.416 * baseline_som_z ^ 2)
forest_plot <- add_text(forest_plot,
    paste0("ES = ", format(baseline_som_frame$eff_size[4], digits = 2), "\np = ", format(baseline_som_summary_p, digits = 3)),
    row = 3,
    col = 4,
    part = "body",
    just = "right",
    gp = gpar(fontsize = base_text_size - 1, fontfamily = base_font_family, col = black_text_colour)
  )
baseline_som_forest_plot <- forest_plot
###############################################################


###############################################################
# Logistic regression for protective score
log_model_protective <- glm(d28_death_yn ~ protective_score*steroids_yn_words + Characteristics.sex. + protective_score*apache_2, data = vanish_joined, family = "binomial")
log_model_protective_summary <- summary(log_model_protective)
p_label_protective <- paste0("p = ", format(log_model_protective_summary$coefficients["protective_score:steroids_yn_words.L", "Pr(>|z|)"], digits = 1))
interaction_plot_protective <- interact_plot(model = log_model_protective, pred = protective_score, modx = steroids_yn_words, interval = TRUE, int.type = "confidence", int.width = 0.95, outcome.scale = "link", line.thickness = line_size_main_mm * 1.5) +
  ylab("Marg. effect: Mort.") + xlab("Protective score") +
  basic_gg_theme() +
  theme(legend.position = "none") +
  scale_colour_manual(values = treatment_colours) +
  scale_fill_manual(values = treatment_colours) +
  scale_linetype_manual(values = c("dashed", "solid")) +
  ggtitle("VANISH") +
  annotate("text", x = 14.5, y = 2, label = p_label_protective, size = text_size_labels)
###############################################################


###############################################################
# Logistic regression for detrimental score
log_model_detrimental <- glm(d28_death_yn ~ detrimental_score*steroids_yn_words + Characteristics.sex. + protective_score*apache_2, data = vanish_joined, family = "binomial")
log_model_detrimental_summary <- summary(log_model_detrimental)
p_label_detrimental <- paste0("p = ", format(log_model_detrimental_summary$coefficients["detrimental_score:steroids_yn_words.L", "Pr(>|z|)"], digits = 1))
interaction_plot_detrimental <- interact_plot(model = log_model_detrimental, pred = detrimental_score, modx = steroids_yn_words, interval = TRUE, int.type = "confidence", int.width = 0.95, outcome.scale = "link", line.thickness = line_size_main_mm * 1.5) +
  ylab("Marg. effect: Mort.") + xlab("Detrimental score") +
  basic_gg_theme() +
  theme(legend.position = "none") +
  scale_colour_manual(values = treatment_colours) +
  scale_fill_manual(values = treatment_colours) +
  scale_linetype_manual(values = c("dashed", "solid")) + 
  annotate("text", x = 17.5, y = 2, label = p_label_detrimental, size = text_size_labels)
###############################################################


###############################################################
aligned_main_fig <- cowplot::align_plots(bacterial_som_plot + theme(legend.position = "none"), SoM_ROC, burns_protective_plot, vanish_protective_plot, interaction_plot_protective,
  align = "hv", axis = "tblr")
cohort_description <- readPNG("/labs/khatrilab/ananthg/detrimental_host_response/figures_2024_09_27/fig_1_combined_bacterial_panel_a.png")
aligned_supp_fig <- cowplot::align_plots(bacterial_detrimental_plot + theme(legend.position = "none"), bacterial_protective_plot + theme(legend.position = "none"), protective_change_plot, detrimental_change_plot, vanish_randomized_protective_plot, vanish_detrimental_plot + theme(legend.position = "none"), interaction_plot_detrimental, bacterial_module_4_plot + theme(legend.position = "none"),
  align = "hv", axis = "tblr")
###############################################################


###############################################################
cairo_pdf("/labs/khatrilab/ananthg/detrimental_host_response/figures_2024_09_27/fig_SCCM_bacterial_infection.pdf", width = 7.5, height = 4.4, bg = "transparent")

pushViewport(viewport(layout = grid.layout(nrow = 44, ncol = 75)))
pushViewport(viewport(layout.pos.row = 2:44, layout.pos.col = 1:54))
grid.raster(cohort_description, x = 0.5, y = 0.5,)
upViewport(1)
print(as_ggplot(aligned_main_fig[[1]]), vp = viewport(layout.pos.row = 2:22, layout.pos.col = 53:75))
print(as_ggplot(aligned_main_fig[[2]]), vp = viewport(layout.pos.row = 21:42, layout.pos.col = 53:75))

grid.text(x = unit(0.5, "npc"), y = unit(0.98, "npc"), label = "Bacterial infection at presentation (17 datasets)", gp = gpar(fontsize = base_text_size, fontface = "plain", col = black_text_colour))

dev.off()
###############################################################


###############################################################
cairo_pdf("/labs/khatrilab/ananthg/detrimental_host_response/figures_2024_09_27/fig_SCCM_hydrocortisone_1.pdf", width = 7.2, height = 3, bg = "transparent")

pushViewport(viewport(layout = grid.layout(nrow = 37, ncol = 100)))

print(as_ggplot(baseline_som_forest_plot), vp = viewport(layout.pos.row = 1:36, layout.pos.col = 1:36))
print(as_ggplot(aligned_supp_fig[[3]]), vp = viewport(layout.pos.row = 4:36, layout.pos.col = 37:69))
print(as_ggplot(aligned_supp_fig[[4]]), vp = viewport(layout.pos.row = 4:36, layout.pos.col = 70:100))
print(hydrocortisone_legend, vp = viewport(layout.pos.row = 36:37, layout.pos.col = 37:100))

grid.text(x = unit(0.185, "npc"), y = unit(0.92, "npc"), label = "Non-survivor (NS) vs.\nSurvivor (S) at baseline", gp = gpar(fontsize = base_text_size, fontface = "plain", col = black_text_colour))
grid.text(x = unit(0.685, "npc"), y = unit(0.97, "npc"), label = "Burn cohort", gp = gpar(fontsize = base_text_size, fontface = "plain", col = black_text_colour))

dev.off()
###############################################################


###############################################################
cairo_pdf("/labs/khatrilab/ananthg/detrimental_host_response/figures_2024_09_27/fig_SCCM_hydrocortisone_2.pdf", width = 7.5, height = 2.6, bg = "transparent")

pushViewport(viewport(layout = grid.layout(nrow = 26, ncol = 75)))
print(as_ggplot(aligned_main_fig[[3]]), vp = viewport(layout.pos.row = 1:24, layout.pos.col = 1:18))
print(as_ggplot(aligned_main_fig[[4]]), vp = viewport(layout.pos.row = 1:24, layout.pos.col = 19:36))
print(as_ggplot(aligned_main_fig[[5]]), vp = viewport(layout.pos.row = 1:24, layout.pos.col = 38:56))
print(as_ggplot(log_odds_forest_plot), vp = viewport(layout.pos.row = 3:25, layout.pos.col = 57:74))

print(hydrocortisone_legend, vp = viewport(layout.pos.row = 25:26, layout.pos.col = 1:56))

grid.text(x = unit(0.875, "npc"), y = unit(0.9, "npc"), label = "Hydrocortisone vs. placebo\nin protective \U2265 median", gp = gpar(fontsize = base_text_size, fontface = "plain", col = black_text_colour))

dev.off()
###############################################################

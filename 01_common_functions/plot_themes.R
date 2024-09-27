###############################################################
library(ggplot2)
library(ggrepel)
library(ggpubr)
library(grid)
library(gridExtra)
library(ggbeeswarm)
library(cowplot)
library(patchwork)
library(Cairo)
library(RColorBrewer)
library(scico)
library(pals)
library(survival)
library(ggfortify)
library(gtsummary)
library(rstatix)
library(ComplexHeatmap)
library(circlize)
###############################################################


###############################################################
base_text_size <- 10
base_font_family <- "Helvetica"
text_size_labels <- (base_text_size - 1) / 2.835# in mm, ~ 9 pt
line_size_main <- 0.75
line_size_main_mm <- 0.75 / 2.835# in mm, ~0.75 in pt
line_size_text_repel <- 0.5
line_size_text_repel_mm <- 0.5 / 2.835# in mm, ~0.5 in pt
point_size <- 0.25
cex_val <- 1# For the beeswarm spread
gap_size <- unit(4, "pt")
black_text_colour <- "#252525"
white_text_colour <- "#EFEFEF"
black_line_colour <- "#505050"
grey_line_colour <- "#CCCCCC"
###############################################################


###############################################################
text_setting_sub_small_gg <- element_text(size = base_text_size - 1, family = base_font_family, colour = black_text_colour)
text_setting_small_gg <- element_text(size = base_text_size, family = base_font_family, colour = black_text_colour)
text_setting_small_title_gg <- element_text(size = base_text_size, family = base_font_family, colour = black_text_colour, hjust = 0.5, vjust = 0)
basic_gg_theme <- function()
{
  theme_bw(base_size = base_text_size, base_family = base_font_family, base_line_size = line_size_main_mm, base_rect_size = line_size_main_mm) +
  theme(text = text_setting_small_gg, axis.text = text_setting_sub_small_gg, axis.title = text_setting_small_gg,
    legend.text = text_setting_sub_small_gg, legend.title = text_setting_small_gg, legend.title.align = 0,
    plot.title = text_setting_small_title_gg,
    strip.text = text_setting_small_gg,
    line = element_line(linewidth = line_size_text_repel_mm, colour = black_line_colour),
    panel.border = element_rect(linewidth = line_size_main_mm, colour = black_line_colour, fill = "transparent"),
    axis.ticks = element_line(linewidth = line_size_text_repel_mm, colour = black_line_colour),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.background = element_blank(), plot.background = element_blank(), strip.background = element_blank(),
    legend.background = element_blank(), legend.box.background = element_blank(), legend.key = element_blank(),
    legend.key.height = unit(base_text_size, "pt"), legend.key.width = unit(base_text_size - 1, "pt"), legend.margin = margin(1, 1, 1, 1, "pt"), legend.spacing.y = unit(1.5, "pt"),
    plot.margin = margin(3, 3, 3, 3, unit = "pt")
    )
}
theme_set(basic_gg_theme())
###############################################################


###############################################################
text_setting_sub_small <- gpar(fontsize = base_text_size - 1, fontfamily = base_font_family, col = black_text_colour)
text_setting_sub_small_white <- gpar(fontsize = base_text_size - 1, fontfamily = base_font_family, col = white_text_colour)
text_setting_small <- gpar(fontsize = base_text_size, fontfamily = base_font_family, col = black_text_colour)
text_setting_small_white <- gpar(fontsize = base_text_size, fontfamily = base_font_family, col = white_text_colour)
dend_lines_gp <- gpar(col = black_line_colour, lwd = line_size_text_repel, lty = 1)
dotted_lines_gp <- gpar(col = black_line_colour, lwd = line_size_text_repel, lty = 2)
ht_opt("heatmap_row_names_gp" = text_setting_small)
ht_opt("heatmap_column_names_gp" = text_setting_small)
ht_opt("heatmap_row_title_gp" = text_setting_small)
ht_opt("heatmap_column_title_gp" = text_setting_small)
ht_opt("TITLE_PADDING" = gap_size)
ht_opt("ROW_ANNO_PADDING" = gap_size)
ht_opt("COLUMN_ANNO_PADDING" = gap_size)
ht_opt("show_parent_dend_line" = FALSE)
ht_opt("legend_border" = gpar(fill = "transparent"))
###############################################################

###############################################################
library(readr)
library(psych)
library(data.table)
library(stringr)
library(plyr)
library(tidyverse)
library(reshape2)
source(file = "/labs/khatrilab/ananthg/detrimental_host_response/common_functions/plot_themes.R")
###############################################################


###############################################################
out_folder <- "/local-scratch/projects/candx/guangbo/ananth_processesed_framingham_data/dhr_paper_re_processesed_2024_09_27/"
age_sex_gen_iii_off_spr <- read_csv("/local-scratch/projects/candx/guangbo/ananth_processesed_framingham_data/age_sex_gen_iii_off_spr.csv")
age_sex_gen_iii_off_spr$sex[age_sex_gen_iii_off_spr$sex == 2] <- "F"
age_sex_gen_iii_off_spr$sex[age_sex_gen_iii_off_spr$sex == 1] <- "M"
###############################################################


###############################################################
obtain_process_gpl <- function(gpl_file_loc)
{
  gpl_file <- fread(gpl_file_loc)
  gpl_filtered <- gpl_file[!(gpl_file$gene_assignment == "---" | is.na(gpl_file$gene_assignment) | gpl_file$gene_assignment == ""), ]
  gene_assignment_transcript_list <- strsplit(gpl_filtered$gene_assignment, split = " /// ")
  gene_assignment_list <- lapply(gene_assignment_transcript_list, function(x)
    {
      gene_name_vec <- c()
      for(i in 1:length(x))
      {
        gene_name_val <- str_split(x[i], " // ")[[1]][2]
        gene_name_vec <- append(gene_name_vec, gene_name_val)
      }
      gene_name_val <- paste(unique(gene_name_vec), sep = "", collapse = ",")
      return(gene_name_val)
    }
  )
  gene_transacript_map_frame <- data.frame(transcript_cluster_id = gpl_filtered$ID, gene_name = unlist(gene_assignment_list))
  return(gene_transacript_map_frame)
}
###############################################################


###############################################################
obtain_clean_up_gene_expr <- function(expr_file_loc, gene_transacript_map_frame, save_name)
{
  expr_subj_keys <- fread("/local-scratch/projects/candx/guangbo/Framingham_data_part1/PhenoGenotypeFiles/ChildStudyConsentSet_phs000363.Framingham.v18.p12.c1.HMB-IRB-MDS/PhenotypeFiles/phs000363.v18.pht002944.v5.p12.c1.Framingham_SABRe_Sample_Attributes.HMB-IRB-MDS.txt")
  colnames(expr_subj_keys) <- as.character(expr_subj_keys[1, ])
  expr_subj_keys <- expr_subj_keys[-1, ]
  expr_matrix <- fread(expr_file_loc)
  colnames(expr_matrix)[1] <- "transcript_cluster_id"
  expr_matrix <- merge(gene_transacript_map_frame, expr_matrix, by = "transcript_cluster_id")
  colnames(expr_matrix)[3:ncol(expr_matrix)] <- sapply(colnames(expr_matrix)[3:ncol(expr_matrix)], function(x) expr_subj_keys$shareid[expr_subj_keys$sampid == x])
  colnames_expr_matrix <- colnames(expr_matrix)
  separated_gene_names_list <- strsplit(expr_matrix$gene_name, split = ",")
  expanded_expr_matrix <- expr_matrix[rep(1:nrow(expr_matrix), sapply(separated_gene_names_list, length)), ]
  expanded_expr_matrix$gene_name <- unlist(separated_gene_names_list)
  expr_matrix <- expanded_expr_matrix[, colnames_expr_matrix[2:length(colnames_expr_matrix)]]
  gene_counts <- table(expr_matrix$gene_name)
  multiple_entry_matrix <- expr_matrix[expr_matrix$gene_name %in% names(gene_counts[gene_counts > 1]), ]
  single_entry_matrix <- expr_matrix[expr_matrix$gene_name %in% names(gene_counts[gene_counts == 1]), ]
  pruned_expr_matrix <- multiple_entry_matrix %>%
    pivot_longer(colnames(multiple_entry_matrix)[2:ncol(multiple_entry_matrix)], names_to = "share_id", values_to = "expr") %>%
    group_by(gene_name, share_id) %>%
    summarize(median_expr = median(expr)) %>%
    pivot_wider(names_from = share_id, values_from = median_expr)
  expr_matrix <- rbind(single_entry_matrix, pruned_expr_matrix[, colnames(single_entry_matrix)])
  write.csv(expr_matrix, file = paste0(out_folder, save_name))
  return(expr_matrix)
}
###############################################################


###############################################################
compute_plot_sex_score <- function(gene_expr_matrix, plot_title, save_name)
{
  male_genes = c("KDM5D", "RPS4Y1", "EIF1AY", "USP9Y", "ZBED1", "DDX3Y", "UTY", "PRKY", "ZFY", "CD99", "CYBB", "TMSB4Y")
  female_genes = c("XIST", "RPS4X", "AMMECR1", "CD40LG", "ZRSR2", "EFHC2", "CA5B", "ZFX", "EIF1AX", "MORC4", "CA5BP1", "UBA1", "SYAP1", "DDX3X", "FUNDC1", "NKRF", "ZC4H2", "PIM2", "SHROOM4", "USP9X", "SMC1A", "NUP62CL", "ERCC6L", "NAA10")
  female_gm <- geometric.mean(gene_expr_matrix[gene_expr_matrix$gene_name %in% female_genes, 2:ncol(gene_expr_matrix)])
  male_gm <- geometric.mean(gene_expr_matrix[gene_expr_matrix$gene_name %in% male_genes, 2:ncol(gene_expr_matrix)])
  sex_gm_score <- female_gm - male_gm
  sex_gm_score <- as.data.frame(sex_gm_score)
  colnames(sex_gm_score) <- c("iSEXS")
  sex_gm_score$shareid <- rownames(sex_gm_score)
  sex_gm_score_merged <- merge(sex_gm_score, age_sex_gen_iii_off_spr[, c("shareid", "sex")], by = "shareid", all.y = FALSE)
  iSEXS_plot <- ggplot(sex_gm_score_merged, aes(x = sex, y = iSEXS, fill = sex)) +
    ylab("iSEXS score") + ggtitle(plot_title) +
    geom_violin(lwd = line_size_text_repel, trim = TRUE, width = 1) +
    geom_beeswarm(cex = cex_val, size = point_size, colour = black_line_colour) +
    geom_boxplot(fill = NA, lwd = line_size_main, width = 0.2, colour = "#000000", outlier.shape = NA) +
    theme(legend.position = "none", axis.title.x = element_blank()) +
    scale_fill_manual(values = c("pink", "blue"))
  ggsave(paste0(out_folder, save_name), iSEXS_plot, width = 50, height = 50, units = "mm", bg = "transparent", device = cairo_pdf)
}
###############################################################


###############################################################
processesed_exon_gpl <- obtain_process_gpl("/local-scratch/projects/candx/guangbo/Framingham_data_part1/PhenoGenotypeFiles/ChildStudyConsentSet_phs000363.Framingham.v18.p12.c1.HMB-IRB-MDS/ExpressionFiles/phe000002.v8.FHS_SABRe_project3.marker-info.MULTI/GPL5188.txt")
off_spr_expr_exon <- obtain_clean_up_gene_expr("/local-scratch/projects/candx/guangbo/Framingham_data_part1/PhenoGenotypeFiles/ChildStudyConsentSet_phs000363.Framingham.v18.p12.c1.HMB-IRB-MDS/ExpressionFiles/phe000002.v8.FHS_SABRe_project3.expression-data-matrixfmt.c1/FinalFile_Exon_OFF_2446_Adjusted_c1.txt",
  processesed_exon_gpl, "off_spr_expr_exon.csv")
gen_iii_expr_exon <- obtain_clean_up_gene_expr("/local-scratch/projects/candx/guangbo/Framingham_data_part1/PhenoGenotypeFiles/ChildStudyConsentSet_phs000363.Framingham.v18.p12.c1.HMB-IRB-MDS/ExpressionFiles/phe000002.v8.FHS_SABRe_project3.expression-data-matrixfmt.c1/Final_Exon_GENIII_3180_c1.txt",
  processesed_exon_gpl, "gen_iii_expr_exon.csv")
compute_plot_sex_score(off_spr_expr_exon, "Off spr., exon",
  "off_spr_expr_exon_iSEXS.pdf")
compute_plot_sex_score(gen_iii_expr_exon, "Gen. III, exon",
  "gen_iii_expr_exon_iSEXS.pdf")
###############################################################

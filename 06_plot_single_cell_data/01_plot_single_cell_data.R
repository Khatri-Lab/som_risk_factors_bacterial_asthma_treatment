# /labs/khatrilab/hongzheng/singlecell/sepsis_covid/DHR.singlecell.figures.Rmd
###############################################################
library(data.table)
library(tidyverse)
library(ggplot2)
library(ggforce)
library(ggrepel)
library(ggpubr)
library(pals)
library(grid)
library(gridExtra)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(rcartocolor)
library(scales)
library(reshape2)
library(dendextend)
library(ggsci)
library(grDevices)
###############################################################


###############################################################
textsize = 10
axistextsize = 8
legendtextsize = 8
legendpointsize = 2
pointsize = 0.3
textcolor = "#252525"
axiscolor = "#252525"
linecolor = "#252525"  
pointsizen = 0.1
linesize = 0.1

mytheme = theme(
  panel.grid.major = element_line(color = "grey67", size = 0.1, linetype = 2),
  panel.grid.minor = element_blank(),
  legend.position = "right", legend.title = element_text(size = textsize, color = textcolor),
  legend.background =  element_blank(),
  legend.box.background = element_blank(),
  legend.key = element_blank(),
  legend.key.width = unit(0.2, "line"),
  legend.key.height = unit(0.2, "line"),
  panel.border = element_blank(),
  axis.line = element_line(color = linecolor, size =linesize),
  axis.ticks = element_line(color = linecolor, size =linesize),
  plot.title = element_text(face = "plain", color = textcolor, size = textsize),
  axis.text.y =  element_text(size = textsize, color = textcolor),
  axis.text.x =  element_text(size = textsize, color = textcolor, angle = 30, hjust = 1),
  axis.title = element_text(size = textsize, color = textcolor),
  strip.background = element_rect(color = NA, fill = "grey91"),
  strip.text.x = element_text(size = textsize, face = "plain", color = textcolor, margin = margin(0.1,0.04,0.1,0.04, "line")),
  strip.text.y = element_text(size = textsize, face = "plain", color = textcolor, margin = margin(0.1,0.04,0.1,0.04, "line")))
###############################################################


###############################################################
ScatterPlotRast <- function(df, featureX, featureY, groupby = NA, colorby = NA, seed =42, shuffle =T, label = NA, labeldf = data.frame(), labelfont = "plain", palette = "category", ncolors = 20, mycol = NA, colcateg = "continuous", linesize = 0.2, pointsize = 0.5, rasterdpi = 100, alphalevel = 0.4, legendpointsize = 2, labelsize = 2, legendColumn = 1, numberofrows = 1, facetdirection = "h", xlab = "", ylab = "", plottitle = "", legendposition = "right", backgroudcol = "grey97", nacol = "grey87", repel = F, addstat = F, scale = "free", facetlabel = "all", facetlabelface = "italic", showaxis = F, showgrid = F, textsize = 10, legendtextsize = 10, axistextsize = 10, contvarlower = 0.01, contvarupper = 0.99
){
  library(ggrastr)
  
  textcolor = "#252525"
  linecolor = "#252525"
  
  if(!is.na(mycol[1]) & colcateg == "discrete"){
    df[, colorby] = factor(df[, colorby], levels = names(mycol))
  }
  
  if(is.na(mycol[1]) & colcateg == "discrete"){
    mycol = GenerateCols(palette, numberofcols = ncolors)
    set.seed(seed)
    mycol =mycol[sample(1:length(mycol))]
  }
  
  if(is.na(mycol[1]) & colcateg == "continuous"){
    mycol = dichromat_pal("DarkRedtoBlue.12")(12)
  }
  
  if(colcateg == "continuous"){
    n1=quantile(df[, colorby], contvarlower)
    n2=quantile(df[, colorby], contvarupper)
    df[, colorby][which(df[, colorby]>=n2)]=n2
    df[, colorby][which(df[, colorby]<=n1)]=n1
  }
  
  if(shuffle ==T){
    set.seed(seed)
    df = df[sample(1:nrow(df)),]}
  if(is.na(groupby)){df$groupdummy =facetlabel}else{df$groupdummy = df[, groupby]}
  
  if(colcateg == "continuous"){
    p=
      ggplot(df, aes_string(x=featureX, y =featureY)) +
      theme_bw()+
      geom_point_rast(aes_string(color = colorby), size =pointsize, stroke = 0, alpha = alphalevel, shape = 16, raster.dpi =rasterdpi)+
      #geom_vline(xintercept = 250) +
      #geom_hline(yintercept = 500) +
      scale_x_continuous(breaks = pretty_breaks()) +
      scale_y_continuous(breaks = pretty_breaks()) +
      facet_wrap(.~groupdummy, nrow = numberofrows, scales = scale, dir =facetdirection) +
      labs(title =plottitle) +
      theme(text = element_text(size = textsize, color = textcolor),
            panel.grid.minor = element_blank(),
            panel.grid.major = element_line(color = "grey87", size = 0.05, linetype = 2),
            panel.border = element_blank(),
            legend.position =legendposition,
            legend.title = element_blank(),
            #legend.background =  element_rect(color =backgroudcol, fill =backgroudcol),
            legend.box.background =  element_blank(),
            legend.box.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"),
            legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"),
            legend.justification = c(0,0.5),
            legend.text = element_text(size = legendtextsize, color = textcolor),
            legend.key.width = unit(0.3, "line"), legend.key.height = unit(0.5, "line"),
            axis.line = element_line(color = linecolor, size =linesize), axis.ticks = element_line(color = linecolor, size =linesize),
            axis.text.y =  element_text(size = axistextsize, color = textcolor),
            axis.text.x =  element_text(size = axistextsize, color = textcolor, angle = 30, hjust = 1),
            axis.title = element_text(size = textsize, color = textcolor, margin = margin(0,0,0,0, "null")),
            plot.title = element_text(size = textsize, color = textcolor, face = "italic", margin = margin(0,0,0,0, "null")),
            plot.background = element_rect(color =backgroudcol, fill =backgroudcol),
            plot.margin = unit(c(0.4,0.4,0.4,0.4), "line"),
            strip.background = element_rect(color = NA, fill = "grey91"),
            strip.text.x = element_text(size = textsize, face =facetlabelface, color = textcolor, margin = margin(0.1,0.04,0.1,0.04, "line")),
            strip.text.y = element_text(size = textsize, face =facetlabelface, color = textcolor, margin = margin(0.1,0.04,0.1,0.04, "line"))
      ) +
      scale_color_gradientn(colours =mycol, name = NULL, na.value = "grey93", breaks = c(n1,n2), labels = c(round(n1, 1), round(n2, 1))) 
  }
  else{
    p=
      ggplot(df, aes_string(x=featureX, y =featureY)) +
      theme_bw()+
      geom_point_rast(aes_string(color = colorby), size =pointsize, stroke = 0, alpha = alphalevel, shape = 16, raster.dpi =rasterdpi)+
      #geom_vline(xintercept = 250) +
      #geom_hline(yintercept = 500) +
      scale_x_continuous(breaks = pretty_breaks()) +
      scale_y_continuous(breaks = pretty_breaks()) +
      facet_wrap(.~groupdummy, nrow = numberofrows, scales = scale, dir =facetdirection) +
      labs(title =plottitle) +
      theme(text = element_text(size = textsize, color = textcolor),
            panel.grid.minor = element_blank(),
            panel.grid.major = element_line(color = "grey87", size = 0.05, linetype = 2),
            panel.border = element_blank(),
            legend.position =legendposition,
            legend.title = element_blank(),
            #legend.background =  element_rect(color =, fill =backgroudcol),
            legend.box.background =  element_blank(),
            legend.box.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"),
            legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"),
            legend.key.width = unit(0.3, "line"), legend.key.height = unit(0.5, "line"),
            legend.justification = c(0,0.5),
            legend.text = element_text(size = legendtextsize, color = textcolor),
            axis.line = element_line(color = linecolor, size =linesize), axis.ticks = element_line(color = linecolor, size =linesize),
            axis.text.y =  element_text(size = axistextsize, color = textcolor),
            axis.text.x =  element_text(size = axistextsize, color = textcolor, angle = 30, hjust = 1),
            axis.title = element_text(size = textsize, color = textcolor, margin = margin(0,0,0,0, "null")),
            plot.title = element_text(size = textsize, color = textcolor, face = "plain", margin = margin(0,0,0,0, "null")),
            plot.background = element_rect(color =backgroudcol, fill =backgroudcol),
            plot.margin = unit(c(0.4,0.4,0.4,0.4), "line"),
            strip.background = element_rect(color = NA, fill = "grey91"),
            strip.text.x = element_text(size = textsize, face =facetlabelface, color = textcolor, margin = margin(0.1,0.04,0.1,0.04, "line")),
            strip.text.y = element_text(size = textsize, face =facetlabelface, color = textcolor, margin = margin(0.1,0.04,0.1,0.04, "line"))
      ) +
      scale_color_manual(values =mycol, name = NULL, na.value = "grey93") +
      guides(color = guide_legend(ncol =legendColumn, label.hjust = 0, override.aes = list(size =legendpointsize, alpha = 0.7))) 
  }
  
  
  if(!is.na(label) & repel ==T){
    p=p+
      geom_label_repel(data = df, aes(label =labeltext), size =labelsize, color = "#252525", force = 0.05, fontface =labelfont)
  }
  
  if(!is.na(label) & repel == F){
    p=p+
      geom_label(data = df, aes(label =labeltext), size =labelsize, color = "#252525", fontface =labelfont)
  }
  
  if(nrow(labeldf) > 0 & repel ==T){
    p=p+
      geom_text_repel(data =labeldf, aes(x=x, y =y, label =labeltext), size =labelsize, color = "#252525", force = 0.05, fontface =labelfont)
  }  
  
  if(nrow(labeldf) > 0  & repel == F){
    p=p+
      geom_text(data =labeldf, aes(x=x, y =y, label =labeltext), size =labelsize, color = "#252525", fontface =labelfont)
  }  
  
  if(addstat ==T){
    p=p+
      stat_cor(method = "pearson", label.x.npc= 0, label.y.npc= 0.95, size =labelsize, color = "grey67")
  }
  
  if(showaxis == F){
    p=p+theme(axis.line = element_blank(), axis.title = element_blank(), axis.ticks = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank())
  }
  
  if(showgrid == F){
    p=p+theme(panel.grid.major = element_blank())
  }
  
  return(p)
}
###############################################################


###############################################################
df4plot = fread("/labs/khatrilab/hongzheng/singlecell/sepsis_covid/rds/dtf.sepsis_covid.v2.subset1.csv")
df4plot = as.data.table(df4plot) # 724012     40

df4plot.f = df4plot[grep("lowQual|-|Erythrocyte", celltype, invert = T)][sampleid != ""][!grepl("nonSNG", sampleid)]
df4plot.f = df4plot.f[grep("other", group, invert = T)]
df4plot.f = df4plot.f[grep("Mast", celltype, invert = T)]
length(unique(df4plot.f$sampleid))

#fwrite(df4plot.f, "/labs/khatrilab/hongzheng/singlecell/sepsis_covid/rds/dtf.clean.sepsis_covid.v2.subset1.csv")

df4plot.f[,.N] # 586020
table(df4plot.f$group)
table(df4plot.f$pid)
table(df4plot.f$severity)
df4plot.f$pid = gsub("schrepping2020a|schrepping2020b", "schrepping2020", df4plot.f$pid )
df4plot.f$celltype = gsub("Pre Neutrophil", "Immat Neutrophil", df4plot.f$celltype)
###############################################################


###############################################################
### pid ----
col_pid = c("#FFDC91FF", "#4DBBD5FF", "#F39B7FFF", "#3C5488FF", "#7876B1FF"); names(col_pid) = sort(unique(df4plot.f$pid))

df4plot.f$facetvar = "study"
pointsizen = 0.001
p.umap.pid =
  ScatterPlotRast(df = as.data.frame(df4plot.f), featureX = "UMAP_1", featureY = "UMAP_2", shuffle = T, pointsize = pointsizen, alphalevel = 0.2, colorby = "pid", groupby = "facetvar", legendColumn = 4, colcateg = "discrete", numberofrows = 5, legendposition = "none", seed = 13, mycol = col_pid, backgroudcol = "white", plottitle = "", rasterdpi = 500, nacol = "grey42", facetlabelface = "plain", legendtextsize = 8) +
  labs(title = NULL)+ 
  theme(text = element_text(size = textsize),
        legend.position = c(1.01,0.5),
        legend.justification = c(0,0.5),
        legend.spacing.x = unit(0.1, "line"),
        legend.key.width = unit(0.2, "line"),
        legend.key.height = unit(0.25, "line"),
        legend.background = element_blank(),
        plot.title = element_text(size = textsize, color = textcolor, face = "plain", hjust = 0.5)) + 
  guides(color = guide_legend(ncol = 1, label.hjust = 0, override.aes = list(size = 2.7, alpha = 0.7))) 

### group ----
col_groupseverity = c(`COVID-19,severe` = "#EFC000FF",
                      `COVID-19,nonsevere` = "#B0A4E3", 
                      `COVID-19,conv` = "cyan", 
                      `bacterial sepsis,severe` = "#FF0097",
                      `bacterial sepsis,nonsevere` = "#214DC8",
                      `bacterial sepsis,conv` = "#06A5C7",
                      `other disease` = "grey67", 
                      healthy = "#1A6840"
)

df4plot.f$facetvar = "group"
pointsizen = 0.001
df4plot.f$groupseverity =factor(df4plot.f$groupseverity, levels =names(col_groupseverity))
p.umap.group=
  ScatterPlotRast(df = as.data.frame(df4plot.f), featureX = "UMAP_1", featureY = "UMAP_2", shuffle = T, pointsize = pointsizen, alphalevel = 0.2, colorby = "groupseverity", groupby = "facetvar", legendColumn = 4, colcateg = "discrete", numberofrows = 5, legendposition = "none", seed = 13, mycol = col_groupseverity, backgroudcol = "white", plottitle = "", rasterdpi = 500, nacol = "grey42", facetlabelface = "plain", legendtextsize = 8) +
  labs(title = NULL)+
  theme(text = element_text(size = textsize),
        legend.position = c(1.01,0.5),
        legend.justification = c(0,0.5),
        legend.spacing.x = unit(0.1, "line"),
        legend.key.width = unit(0.2, "line"),
        legend.key.height = unit(0.25, "line"),
        legend.background = element_blank(),
        plot.title = element_text(size = textsize, color = textcolor, face = "plain", hjust = 0.5)) +
  guides(color = guide_legend(ncol = 1, label.hjust = 0, override.aes = list(size = 2.7, alpha = 0.7)))

### celltype ----
col_celltype = c(`CD14 Monocyte` = "#8CC269", `CD16 Monocyte` = "#70F3FF", cDC = "#B0A4E3", pDC = "#CCA4E3",
                 `Neutrophil` = "#2E42B8", `Immat Neutrophil` = "#0A73DC", `Pre Neutrophil` = "#BADEFA",`Neu-Lym` = "#F39B7FFF", 
                 "Platelet" = "grey88", Erythrocyte = "#FD8CC1FF", HSPC = "#DC3023", "Mast" = "#352A87",
                 `CD4 T Naive` = "#7397AB", `CD8 T Naive` = "#D2F0F4", `CD4 T Eff/Mem` = "#33B7A0", `CD8 T Eff/Mem` = "#F86B1D",
                 Treg = "#1A6840", `Prolif T/NK` = "#F9906F", NK = "#FFA631", NKT = "#EFC000FF", 
                 B = "#F9FB0E", PB = "#FFDC91FF", "B-T" = "#003C67FF")

df.label = as.data.table(df4plot.f)[, lapply(.SD, median), by = celltype,.SDcols = c("UMAP_1", "UMAP_2")]
colnames(df.label)[2:3]= c("x", "y")
df.label$labeltext = df.label[, "celltype",with = F]
df.label = df.label[!is.na(celltype)]
df.label = as.data.table(df.label)
df.label[, x:=ifelse(celltype == "cDC", x+1, x)]
df.label[, y:=ifelse(celltype == "HSPC", y-0.5, y)]
df.label[, y:=ifelse(celltype == "CD4 T Naive", y+1.5, y)]
df.label[, y:=ifelse(celltype == "CD4 T Eff/Mem", y+0.5, y)]
df.label[, x:=ifelse(celltype == "CD8 T Eff/Mem", x+1, x)]
df.label[, x:=ifelse(celltype == "CD8 T Naive", x+1, x)]
df.label[, x:=ifelse(celltype == "CD14 Monocyte", x-1, x)]
# df.label[, y:=ifelse(celltype == "NK", y-0.5, y)]
df4plot.f$facetvar = "cell type"
pointsizen = 0.001
p.umap.celltype =
  ScatterPlotRast(df = as.data.frame(df4plot.f), featureX = "UMAP_1", featureY = "UMAP_2", shuffle = T, pointsize = 0.01, alphalevel = 0.3, colorby = "celltype", groupby = "facetvar", legendColumn = 4, colcateg = "discrete", numberofrows = 5, legendposition = "none", seed = 13, mycol = col_celltype, backgroudcol = "white", plottitle = "", rasterdpi = 300, nacol = "grey42", labeldf = df.label, facetlabelface = "plain", repel = F, legendtextsize = 8) +
  labs(title = NULL)+ 
  theme(text = element_text(size = textsize),
        legend.position = c(1.01,0.55),
        legend.justification = c(0,0.5),
        legend.spacing.x = unit(0.1, "line"),
        legend.key.width = unit(0.2, "line"),
        legend.key.height = unit(0.2, "line"),
        legend.background = element_blank(),
        plot.title = element_text(size = textsize, color = textcolor, face = "plain", hjust = 0.5)) + 
  guides(color = guide_legend(ncol = 1, label.hjust = 0, override.aes = list(size = 2.7, alpha = 0.7))) 
###############################################################


###############################################################
# UMAP ----
#df4plot.f = fread("/labs/khatrilab/hongzheng/singlecell/sepsis_covid/rds/dtf.clean.sepsis_covid.v2.subset1.csv")
#df4plot.f$celltype = gsub("Pre Neutrophil", "Immat Neutrophil", df4plot.f$celltype)
scores2plot = c(paste0("SoMmodule",1:4), "detrimental", "protective", "SoMScore")
names(scores2plot) = c(paste0("module ",1:4), "detrimental", "protective", "SoM score")
pointsizen = 0.001
df4plot.f$SoMScore = df4plot.f$SoMmodule1 + df4plot.f$SoMmodule2 - (df4plot.f$SoMmodule3 + df4plot.f$SoMmodule4)

p.umap.score = list()
for(i in names(scores2plot)){
  df4plot.f$facetvar = i
  p.umap.score[[scores2plot[i]]]=
    ScatterPlotRast(df = as.data.frame(df4plot.f), featureX = "UMAP_1", featureY = "UMAP_2", shuffle = T, pointsize = pointsizen, alphalevel = 0.2, colorby = scores2plot[i], groupby = "facetvar", legendColumn = 4, colcateg = "continuous", numberofrows = 5, legendposition = "right", seed = 13, mycol = dichromat_pal("DarkRedtoBlue.12")(12), backgroudcol = "white", plottitle = "", rasterdpi = 300, nacol = "grey42", facetlabelface = "plain") +
    theme(text = element_text(size = textsize),
          #legend.position = c(0.5,-0.06),
          legend.position = "none",
          legend.justification = c(0.5,1),
          legend.text = element_text(size = textsize),
          legend.title = element_text(size = textsize),
          legend.key.width = unit(0.5, "line"),
          legend.key.height = unit(0.2, "line"),
          legend.box.background =  element_blank(),
          legend.title.align = 0.5,
          legend.direction = "horizontal",
          legend.spacing.x = unit(.2, "line"),
          legend.spacing.y = unit(.2, "line")
    ) + labs(title = NULL)+
    scale_color_gradientn(colours = c("grey95", dichromat_pal("DarkRedtoBlue.12")(12)[c(7:12)]), name = NULL, breaks = c(1,1.5), labels = c("low", "high")) +
    guides(color = guide_colourbar(direction = "horizontal", title.position = "top", title.hjust = 0.5, title.vjust = 1, ticks = TRUE, label = TRUE, label.hjust = 0.5))
}

# module heatmap ----
matexpr = df4plot.f[, lapply(.SD, mean), by = "celltype",.SDcols = c("SoMScore", "detrimental", "protective", paste0("SoMmodule",1:4))]
mat.rownames = matexpr$celltype
mat.expr = as.matrix(matexpr[, 2:ncol(matexpr)])
rownames(mat.expr) = mat.rownames
  
mat.rownames = matexpr$celltype
mat.expr = as.matrix(matexpr[, 2:ncol(matexpr)])
rownames(mat.expr) = mat.rownames

mat.expr = scale(mat.expr, center = T)
maxvalue = 2.5
mat.expr[mat.expr> maxvalue] = maxvalue
mat.expr[mat.expr< -maxvalue]= -maxvalue
mat.expr = t(mat.expr)

mat.expr = mat.expr[, intersect(names(col_celltype), colnames(mat.expr))]

scale_values <- function(x){(x-min(x))/(max(x)-min(x))}
values = scale_values(table(df4plot.f$celltype)) + 0.001
mat.N = matrix(rep(values,7), nrow = 7, byrow = T)
colnames(mat.N) = names(values)

rownames(mat.expr) =  c("SoM score", "detrimental", "protective", paste0("module ",1:4))
rownames(mat.N) = rownames(mat.expr) 

mat.N = mat.N[, colnames(mat.expr)]
mat.N = mat.N^0.2

col_fun = col_fun = circlize::colorRamp2(seq(-maxvalue, maxvalue, length = 13), 
                                         c(dichromat_pal("DarkRedtoBlue.12")(12)[1:6], "white", dichromat_pal("DarkRedtoBlue.12")(12)[7:12]))
p.SoMmodules = Heatmap(mat.expr,
               rect_gp = gpar(type = "none"),
               col = col_fun,
               na_col = "grey93",
               cluster_columns = F, show_row_names = T,
               cluster_rows = F, show_row_dend = F, show_column_dend = F, row_dend_side = "left", column_dend_side = "top",
               clustering_method_rows = "ward.D2", clustering_method_columns  = "ward.D2",
               clustering_distance_columns  = "euclidean", clustering_distance_rows  = "euclidean",
               row_names_side = "left",
               column_names_side = "top", column_names_rot = 55, column_names_gp = gpar(fontsize = 10, fontface = "plain"),
               row_names_gp = gpar(fontsize = 8, fontface = "plain"),
               column_title = NULL, column_title_side = "top", column_title_gp = gpar(fontsize = 10, fontface = "plain"),
               show_heatmap_legend = F,
               cell_fun =  function(j, i, x, y, width, height, fill){
                 grid.rect(x = x, y = y, width = width, height =  height, gp = gpar(lwd = 0.2, col = "grey80", fill = "grey67", alpha = 0.2))
                 #grid.circle(x = x, y = y, r = 1.2 * mat.N[i,j]  * min(unit.c(width, height)), gp = gpar(fill = col_fun(mat.expr[i, j]), col = NA, alpha = 0.9))},
                 grid.circle(x = x, y = y, r = 1.1 * min(unit.c(width, height)), gp = gpar(fill = col_fun(mat.expr[i, j]), col = NA, alpha = 0.9))},
               heatmap_legend_param=list(title = "expression", title_gp = gpar( fontsize = 10), labels_gp = gpar(fontsize = 8), legend_width = unit(1.8, "line"), legend_height = unit(.2, "line"), grid_width = unit(1.5, "line"), grid_height = unit(.1, "line"), at = c(1, max(mat.expr)), labels = c("low", "high"), by_row = F, direction = "horizontal", title_position = "lefttop"))

###############################################################


###############################################################
cairo_pdf(file = "/labs/khatrilab/ananthg/detrimental_host_response/figures_2024_09_27/DHR.singlecell.figure.pdf",width = 5.8, height = 6.5)
pushViewport(viewport(layout =grid.layout(nrow= 110, ncol =48)))

print(p.umap.score$detrimental, vp = viewport(layout.pos.row = 42:60, layout.pos.col = 13:21))
print(p.umap.score$protective, vp = viewport(layout.pos.row = 59:77, layout.pos.col = 13:21))

print(p.umap.score$SoMmodule1, vp = viewport(layout.pos.row = 42:60, layout.pos.col = 24:32))
print(p.umap.score$SoMmodule2, vp = viewport(layout.pos.row = 42:60, layout.pos.col = 34:42))
print(p.umap.score$SoMmodule3, vp = viewport(layout.pos.row = 59:77, layout.pos.col = 24:32))
print(p.umap.score$SoMmodule4, vp = viewport(layout.pos.row = 59:77, layout.pos.col = 34:42))

print(p.umap.score$SoMScore, vp = viewport(layout.pos.row = 50:68, layout.pos.col = 2:10))
print(p.umap.celltype, vp = viewport(layout.pos.row = 2:41, layout.pos.col = 19:39))
print(p.umap.pid, vp = viewport(layout.pos.row = 2:21, layout.pos.col = 2:11))
print(p.umap.group, vp = viewport(layout.pos.row = 22:41, layout.pos.col = 2:11))

pushViewport(viewport(layout.pos.row = 76:110, layout.pos.col = 4:40))
draw(p.SoMmodules, merge_legend =T, heatmap_legend_side = "bottom", newpage = FALSE)
upViewport(1)

col_fun.label1 = circlize::colorRamp2(seq(1,3, length =7), c("grey93", dichromat_pal("DarkRedtoBlue.12")(12)[c(7:12)]))
lgd = Legend(at = c(1.1,2.9), labels = c("low", "high"),
             col_fun = col_fun.label1, title = "score",
             grid_width = unit(0.1, "line"), grid_height = unit(0.1, "line"),
             legend_height = unit(3, "line"),
             direction = "vertical", title_position = "leftcenter-rot", title_gp = gpar(fontsize = 10, fontface = "plain"), labels_gp = gpar(fontsize = 10))
draw(lgd, x = unit(0.92, "npc"), y = unit(0.47, "npc"))

col_fun.label2 = circlize::colorRamp2(seq(-1,1, length = 13), c(dichromat_pal("DarkRedtoBlue.12")(12)[1:6], "white", dichromat_pal("DarkRedtoBlue.12")(12)[7:12]))
lgd = Legend(at = c(-0.7,0.7), labels = c("low", "high"),
             col_fun = col_fun.label2, title = "scaled score",
             grid_width = unit(0.1, "line"), grid_height = unit(0.1, "line"),
             legend_height = unit(3, "line"),
             direction = "vertical", title_position = "leftcenter-rot", title_gp = gpar(fontsize = 10, fontface = "plain"), labels_gp = gpar(fontsize = 10))
#draw(lgd, x = unit(0.1, "npc"), y = unit(0.23, "npc"))
draw(lgd, x = unit(0.88, "npc"), y = unit(0.08, "npc"))

grid.text(x= unit(0.02, "npc"), y = unit(0.98, "npc"), label = "A", gp=gpar(fontsize = textsize, fontface = "bold"))
grid.text(x= unit(0.02, "npc"), y = unit(0.8, "npc"), label = "B", gp=gpar(fontsize = textsize, fontface = "bold"))
grid.text(x= unit(0.37, "npc"), y = unit(0.98, "npc"), label = "C", gp=gpar(fontsize = textsize, fontface = "bold"))
grid.text(x= unit(0.02, "npc"), y = unit(0.61, "npc"), label = "D", gp=gpar(fontsize = textsize, fontface = "bold"))
grid.text(x= unit(0.02, "npc"), y = unit(0.28, "npc"), label = "E", gp=gpar(fontsize = textsize, fontface = "bold"))

grid.lines(x = unit(c(0.345, 0.355), "npc"), y = unit(c(0.47, 0.47), "npc"), gp = gpar(fill = "grey42", lwd = 1.5, col = "grey17"))
grid.text(x= unit(0.23, "npc"), y = unit(c(0.47), "npc"), label = "= ", gp=gpar(fontsize = 14, fontface = "plain", col = "grey17"))
grid.text(x= unit(0.46, "npc"), y = unit(0.39, "npc"), label = "= ", gp=gpar(fontsize = 14, fontface = "plain", col = "grey17"))
grid.text(x= unit(0.46, "npc"), y = unit(0.54, "npc"), label = "= ", gp=gpar(fontsize = 14, fontface = "plain", col = "grey17"))
grid.text(x= unit(0.68, "npc"), y = unit(0.39, "npc"), label = "+", gp=gpar(fontsize = 14, fontface = "plain", col = "grey17"))
grid.text(x= unit(0.68, "npc"), y = unit(0.54, "npc"), label = "+", gp=gpar(fontsize = 14, fontface = "plain", col = "grey17"))

dev.off()
###############################################################


###############################################################
library(Matrix)
mat = readMM("/labs/khatrilab/hongzheng/singlecell/sepsis_covid/rds/mat.sepsis_covid.v2.sub1.mtx")
matrownames = fread("/labs/khatrilab/hongzheng/singlecell/sepsis_covid/rds/mat.sepsis_covid.v2.sub1.mtx.rownames.txt", header = F)
rownames(mat) = matrownames$V1

SoMmodule1 = c("NQO2","SLPI","ORM1","KLHL2","ANXA3","TXN","AQP9","BCL6","DOK3","PFKFB4","TYK2")
SoMmodule2 = c("BCL2L11","BCAT1","BTBD7","CEP55","HMMR","PRC1","KIF15","CAMP","CEACAM8","DEFA4","LCN2","CTSG","AZU1")
SoMmodule3 = c("MAFB","OASL","UBE2L6","VAMP5","CCL2","NAPA","ATG3","VRK2","TMEM123","CASP7")
SoMmodule4 = c("DOK2","HLA-DPB1","BUB3","SMYD2","SIDT1","EXOC2","TRIB2","KLRB1")

selectedgenes = intersect(c(SoMmodule1,SoMmodule2,SoMmodule3,SoMmodule4), rownames(mat))
mat.selectedgenes =  mat[selectedgenes,] #  42 586020

dtfmeta = fread("/labs/khatrilab/hongzheng/singlecell/sepsis_covid/rds/dtf.clean.sepsis_covid.v2.subset1.csv")
dtfmetagenes = cbind(dtfmeta, as.matrix(t(mat.selectedgenes)))

dtfmetagenes = dtfmetagenes[grep("Mast", celltype, invert = T)]
dtfmetagenes$celltype = gsub("Pre Neutrophil", "Immat Neutrophil", dtfmetagenes$celltype)

genes = intersect(c(SoMmodule1,SoMmodule2,SoMmodule3,SoMmodule4), rownames(mat.selectedgenes))
matexpr = dtfmetagenes[, lapply(.SD, mean), by = "celltype",.SDcols =genes]
#matexpr = dtfmetagenes[, lapply(.SD, geomMean), by = "celltype",.SDcols =genes]
mat.rownames = matexpr$celltype
mat.expr = as.matrix(matexpr[,2:ncol(matexpr)])
rownames(mat.expr) = mat.rownames
  
mat.rownames = matexpr$celltype
mat.expr = as.matrix(matexpr[,2:ncol(matexpr)])
rownames(mat.expr) = mat.rownames

mat.expr = scale(mat.expr, center = T)
maxvalue = 2.5
mat.expr[mat.expr> maxvalue] = maxvalue
mat.expr[mat.expr< -maxvalue]= -maxvalue
mat.expr = t(mat.expr)

mat.expr = mat.expr[, intersect(names(col_celltype), colnames(mat.expr))]

scale_values <- function(x){(x-min(x))/(max(x)-min(x))}
values = scale_values(table(dtfmetagenes$celltype)) + 0.001
mat.N = matrix(rep(values,42), nrow = 42, byrow = T)
colnames(mat.N) = names(values)
rownames(mat.N) = rownames(mat.expr) 

mat.N = mat.N[rownames(mat.expr), colnames(mat.expr)]
mat.N = mat.N^0.2

annodf = data.table(gene = c(SoMmodule1,SoMmodule2,SoMmodule3,SoMmodule4),
                    module = c(rep("module1", length(SoMmodule1)), rep("module2", length(SoMmodule2)), rep("module3", length(SoMmodule3)), rep("module4", length(SoMmodule4)))
                    )

annodf = annodf[match(rownames(mat.expr), annodf$gene),]
modulecol = c("module1" = "#EFC002", "module2" = "#F39B7F", "module3" = "#00A087", "module4" = "#3C5488")

row_ha = rowAnnotation(
  module = factor(annodf$module),
  simple_anno_size = unit(.3, "line"),
  annotation_name_gp = gpar(fontsize = 10),
  col =list(module = modulecol),
  show_annotation_name = F,
   show_legend = F,
  annotation_legend_param = list(
    module = list(
      nrow = 4,
      title_gp = gpar(fontsize = 10),
      legend_width = unit(0.2, "line"), legend_height = unit(0.2, "line"),
      grid_width = unit(0.2, "line"), grid_height = unit(0.5, "line"),
      labels_gp = gpar(fontsize = 10, legend_width = unit(0.5, "line")))))


col_fun = col_fun = circlize::colorRamp2(seq(-maxvalue, maxvalue, length = 13), c(dichromat_pal("DarkRedtoBlue.12")(12)[1:6], "white", dichromat_pal("DarkRedtoBlue.12")(12)[7:12]))
p.42 = Heatmap(mat.expr,
               rect_gp = gpar(type = "none"),
               col = col_fun,
               na_col = "grey93",
               right_annotation = row_ha,
               cluster_columns = F, show_row_names = T,
               cluster_rows = F, show_row_dend = F, show_column_dend = F, row_dend_side = "left", column_dend_side = "top",
               clustering_method_rows = "ward.D2", clustering_method_columns  = "ward.D2",
               clustering_distance_columns  = "euclidean", clustering_distance_rows  = "euclidean",
               row_names_side = "left",
               column_names_side = "top", column_names_rot = 90, column_names_gp = gpar(fontsize = 10, fontface = "plain"),
               row_names_gp = gpar(fontsize = 8, fontface = "italic"),
               column_title = NULL, column_title_side = "top", column_title_gp = gpar(fontsize = 10, fontface = "plain"),
               show_heatmap_legend = T,
               cell_fun =  function(j, i, x, y, width, height, fill){
                 grid.rect(x = x, y = y, width = width, height =  height, gp = gpar(lwd = 0.2, col = "grey90", fill = NA, alpha = 0.2))
                 grid.circle(x = x, y = y, r = 0.8 * mat.N[i,j]  * min(unit.c(width, height)), gp = gpar(fill = col_fun(mat.expr[i, j]), col = NA, alpha = 0.9))},
                 #grid.circle(x = x, y = y, r = 0.8   * min(unit.c(width, height)), gp = gpar(fill = col_fun(mat.expr[i, j]), col = NA, alpha = 0.9))},
               heatmap_legend_param=list(title = "expression", title_gp = gpar( fontsize = 10), labels_gp = gpar(fontsize = 8), legend_width = unit(1.8, "line"), legend_height = unit(.2, "line"), grid_width = unit(1.5, "line"), grid_height = unit(.1, "line"), at = c(-2,2), labels = c("low", "high"), by_row = F, direction = "horizontal", title_position = "lefttop"),
               width = ncol(mat.expr) * unit(16, "pt"))
height_p42 <- nrow(mat.expr) * unit(16, "pt")

cairo_pdf(file =paste0("/labs/khatrilab/ananthg/detrimental_host_response/figures_2024_09_27/DHR.singlecell.figure.sup.pdf"),width = 3.4, height = 6)
pushViewport(viewport(layout =grid.layout(nrow= 100, ncol = 100)))
pushViewport(viewport(layout.pos.row = 1:100, layout.pos.col = 1:95))
draw(p.42, merge_legend =T, heatmap_legend_side = "bottom", newpage = FALSE)

grid.text(x= unit(0.995, "npc"), y = unit(0.7, "npc"), label = "module 1", rot =90, gp=gpar(fontsize = textsize, fontface = "plain", col = "#EFC002"))
grid.text(x= unit(0.995, "npc"), y = unit(0.49, "npc"), label = "module 2", rot =90, gp=gpar(fontsize = textsize, fontface = "plain", col = "#F39B7F"))
grid.text(x= unit(0.995, "npc"), y = unit(0.27, "npc"), label = "module 3", rot =90, gp=gpar(fontsize = textsize, fontface = "plain", col = "#00A087"))
grid.text(x= unit(0.995, "npc"), y = unit(0.12, "npc"), label = "module 4", rot =90, gp=gpar(fontsize = textsize, fontface = "plain", col = "#3C5488"))

upViewport(1)
dev.off()
###############################################################

# revision

###############################################################
df4plot = fread("/labs/khatrilab/hongzheng/singlecell/sepsis_covid/rds/dtf.sepsis_covid.v2.subset1.csv")
df4plot = as.data.table(df4plot) # 724012     40

df4plot.f = df4plot[grep("lowQual|-|Erythrocyte", celltype, invert = T)][sampleid != ""][!grepl("nonSNG", sampleid)]
df4plot.f = df4plot.f[grep("other", group, invert = T)]
df4plot.f = df4plot.f[grep("Mast", celltype, invert = T)]
length(unique(df4plot.f$sampleid))

#fwrite(df4plot.f, "/labs/khatrilab/hongzheng/singlecell/sepsis_covid/rds/dtf.clean.sepsis_covid.v2.subset1.csv")

df4plot.f[,.N] # 586020
table(df4plot.f$group)
table(df4plot.f$pid)
table(df4plot.f$severity)
df4plot.f$pid = gsub("schrepping2020a|schrepping2020b", "schrepping2020", df4plot.f$pid )
df4plot.f$celltype = gsub("Pre Neutrophil", "Immat Neutrophil", df4plot.f$celltype)

# update group labels
df4plot.f$groupseverity = gsub("bacterial sepsis,nonsevere", "other infection,nonsevere", df4plot.f$groupseverity)
df4plot.f$groupseverity = gsub("bacterial sepsis,severe", "bacterial sepsis", df4plot.f$groupseverity)
table(df4plot.f$groupseverity)

###############################################################


###############################################################
### pid ----
col_pid = c("#FFDC91FF", "#4DBBD5FF", "#F39B7FFF", "#3C5488FF", "#7876B1FF"); names(col_pid) = sort(unique(df4plot.f$pid))

df4plot.f$facetvar = "study"
pointsizen = 0.001
p.umap.pid =
  ScatterPlotRast(df = as.data.frame(df4plot.f), featureX = "UMAP_1", featureY = "UMAP_2", shuffle = T, pointsize = pointsizen, alphalevel = 0.2, colorby = "pid", groupby = "facetvar", legendColumn = 4, colcateg = "discrete", numberofrows = 5, legendposition = "none", seed = 13, mycol = col_pid, backgroudcol = "white", plottitle = "", rasterdpi = 500, nacol = "grey42", facetlabelface = "plain", legendtextsize = 8) +
  labs(title = NULL)+ 
  theme(text = element_text(size = textsize),
        legend.position = c(1.01,0.5),
        legend.justification = c(0,0.5),
        legend.spacing.x = unit(0.1, "line"),
        legend.key.width = unit(0.2, "line"),
        legend.key.height = unit(0.25, "line"),
        legend.background = element_blank(),
        plot.title = element_text(size = textsize, color = textcolor, face = "plain", hjust = 0.5)) + 
  guides(color = guide_legend(ncol = 1, label.hjust = 0, override.aes = list(size = 2.7, alpha = 0.7))) 

### group ----
col_groupseverity = c(`COVID-19,severe` = "#EFC000FF",
                      `COVID-19,nonsevere` = "#B0A4E3", 
                      `COVID-19,conv` = "cyan", 
                      `bacterial sepsis` = "#FF0097",
                      `other infection,nonsevere` = "#214DC8",
                      healthy = "#1A6840"
)

df4plot.f$facetvar = "group"
pointsizen = 0.001
df4plot.f$groupseverity =factor(df4plot.f$groupseverity, levels =names(col_groupseverity))
p.umap.group=
  ScatterPlotRast(df = as.data.frame(df4plot.f), featureX = "UMAP_1", featureY = "UMAP_2", shuffle = T, pointsize = pointsizen, alphalevel = 0.2, colorby = "groupseverity", groupby = "facetvar", legendColumn = 4, colcateg = "discrete", numberofrows = 5, legendposition = "none", seed = 13, mycol = col_groupseverity, backgroudcol = "white", plottitle = "", rasterdpi = 500, nacol = "grey42", facetlabelface = "plain", legendtextsize = 8) +
  labs(title = NULL)+
  theme(text = element_text(size = textsize),
        legend.position = c(1.01,0.5),
        legend.justification = c(0,0.5),
        legend.spacing.x = unit(0.1, "line"),
        legend.key.width = unit(0.2, "line"),
        legend.key.height = unit(0.25, "line"),
        legend.background = element_blank(),
        plot.title = element_text(size = textsize, color = textcolor, face = "plain", hjust = 0.5)) +
  guides(color = guide_legend(ncol = 1, label.hjust = 0, override.aes = list(size = 2.7, alpha = 0.7)))

### celltype ----
col_celltype = c(`CD14 Monocyte` = "#8CC269", `CD16 Monocyte` = "#70F3FF", cDC = "#B0A4E3", pDC = "#CCA4E3",
                 `Neutrophil` = "#2E42B8", `Immat Neutrophil` = "#0A73DC", `Pre Neutrophil` = "#BADEFA",`Neu-Lym` = "#F39B7FFF", 
                 "Platelet" = "grey88", Erythrocyte = "#FD8CC1FF", HSPC = "#DC3023", "Mast" = "#352A87",
                 `CD4 T Naive` = "#7397AB", `CD8 T Naive` = "#D2F0F4", `CD4 T Eff/Mem` = "#33B7A0", `CD8 T Eff/Mem` = "#F86B1D",
                 Treg = "#1A6840", `Prolif T/NK` = "#F9906F", NK = "#FFA631", NKT = "#EFC000FF", 
                 B = "#F9FB0E", PB = "#FFDC91FF", "B-T" = "#003C67FF")

df.label = as.data.table(df4plot.f)[, lapply(.SD, median), by = celltype,.SDcols = c("UMAP_1", "UMAP_2")]
colnames(df.label)[2:3]= c("x", "y")
df.label$labeltext = df.label[, "celltype",with = F]
df.label = df.label[!is.na(celltype)]
df.label = as.data.table(df.label)
df.label[, x:=ifelse(celltype == "cDC", x+1, x)]
df.label[, y:=ifelse(celltype == "HSPC", y-0.5, y)]
df.label[, y:=ifelse(celltype == "CD4 T Naive", y+1.5, y)]
df.label[, y:=ifelse(celltype == "CD4 T Eff/Mem", y+0.5, y)]
df.label[, x:=ifelse(celltype == "CD8 T Eff/Mem", x+1, x)]
df.label[, x:=ifelse(celltype == "CD8 T Naive", x+1, x)]
df.label[, x:=ifelse(celltype == "CD14 Monocyte", x-1, x)]
# df.label[, y:=ifelse(celltype == "NK", y-0.5, y)]
df4plot.f$facetvar = "cell type"
pointsizen = 0.001
p.umap.celltype =
  ScatterPlotRast(df = as.data.frame(df4plot.f), featureX = "UMAP_1", featureY = "UMAP_2", shuffle = T, pointsize = 0.01, alphalevel = 0.3, colorby = "celltype", groupby = "facetvar", legendColumn = 4, colcateg = "discrete", numberofrows = 5, legendposition = "none", seed = 13, mycol = col_celltype, backgroudcol = "white", plottitle = "", rasterdpi = 300, nacol = "grey42", labeldf = df.label, facetlabelface = "plain", repel = F, legendtextsize = 8) +
  labs(title = NULL)+ 
  theme(text = element_text(size = textsize),
        legend.position = c(1.01,0.55),
        legend.justification = c(0,0.5),
        legend.spacing.x = unit(0.1, "line"),
        legend.key.width = unit(0.2, "line"),
        legend.key.height = unit(0.2, "line"),
        legend.background = element_blank(),
        plot.title = element_text(size = textsize, color = textcolor, face = "plain", hjust = 0.5)) + 
  guides(color = guide_legend(ncol = 1, label.hjust = 0, override.aes = list(size = 2.7, alpha = 0.7))) 
###############################################################


###############################################################
# UMAP ----
#df4plot.f = fread("/labs/khatrilab/hongzheng/singlecell/sepsis_covid/rds/dtf.clean.sepsis_covid.v2.subset1.csv")
#df4plot.f$celltype = gsub("Pre Neutrophil", "Immat Neutrophil", df4plot.f$celltype)
scores2plot = c(paste0("SoMmodule",1:4), "detrimental", "protective", "SoMScore")
names(scores2plot) = c(paste0("module ",1:4), "detrimental", "protective", "SoM score")
pointsizen = 0.001
df4plot.f$SoMScore = df4plot.f$SoMmodule1 + df4plot.f$SoMmodule2 - (df4plot.f$SoMmodule3 + df4plot.f$SoMmodule4)

p.umap.score = list()
for(i in names(scores2plot)){
  df4plot.f$facetvar = i
  p.umap.score[[scores2plot[i]]]=
    ScatterPlotRast(df = as.data.frame(df4plot.f), featureX = "UMAP_1", featureY = "UMAP_2", shuffle = T, pointsize = pointsizen, alphalevel = 0.2, colorby = scores2plot[i], groupby = "facetvar", legendColumn = 4, colcateg = "continuous", numberofrows = 5, legendposition = "right", seed = 13, mycol = dichromat_pal("DarkRedtoBlue.12")(12), backgroudcol = "white", plottitle = "", rasterdpi = 300, nacol = "grey42", facetlabelface = "plain") +
    theme(text = element_text(size = textsize),
          #legend.position = c(0.5,-0.06),
          legend.position = "none",
          legend.justification = c(0.5,1),
          legend.text = element_text(size = textsize),
          legend.title = element_text(size = textsize),
          legend.key.width = unit(0.5, "line"),
          legend.key.height = unit(0.2, "line"),
          legend.box.background =  element_blank(),
          legend.title.align = 0.5,
          legend.direction = "horizontal",
          legend.spacing.x = unit(.2, "line"),
          legend.spacing.y = unit(.2, "line")
    ) + labs(title = NULL)+
    scale_color_gradientn(colours = c("grey95", dichromat_pal("DarkRedtoBlue.12")(12)[c(7:12)]), name = NULL, breaks = c(1,1.5), labels = c("low", "high")) +
    guides(color = guide_colourbar(direction = "horizontal", title.position = "top", title.hjust = 0.5, title.vjust = 1, ticks = TRUE, label = TRUE, label.hjust = 0.5))
}

# module heatmap ----
matexpr = df4plot.f[, lapply(.SD, mean), by = "celltype",.SDcols = c("SoMScore", "detrimental", "protective", paste0("SoMmodule",1:4))]
mat.rownames = matexpr$celltype
mat.expr = as.matrix(matexpr[, 2:ncol(matexpr)])
rownames(mat.expr) = mat.rownames
  
mat.rownames = matexpr$celltype
mat.expr = as.matrix(matexpr[, 2:ncol(matexpr)])
rownames(mat.expr) = mat.rownames

mat.expr = scale(mat.expr, center = T)
maxvalue = 2.5
mat.expr[mat.expr> maxvalue] = maxvalue
mat.expr[mat.expr< -maxvalue]= -maxvalue
mat.expr = t(mat.expr)

mat.expr = mat.expr[, intersect(names(col_celltype), colnames(mat.expr))]

scale_values <- function(x){(x-min(x))/(max(x)-min(x))}
values = scale_values(table(df4plot.f$celltype)) + 0.001
mat.N = matrix(rep(values,7), nrow = 7, byrow = T)
colnames(mat.N) = names(values)

rownames(mat.expr) =  c("SoM score", "detrimental", "protective", paste0("module ",1:4))
rownames(mat.N) = rownames(mat.expr) 

mat.N = mat.N[, colnames(mat.expr)]
mat.N = mat.N^0.2

col_fun = col_fun = circlize::colorRamp2(seq(-maxvalue, maxvalue, length = 13), 
                                         c(dichromat_pal("DarkRedtoBlue.12")(12)[1:6], "white", dichromat_pal("DarkRedtoBlue.12")(12)[7:12]))
p.SoMmodules = Heatmap(mat.expr,
               rect_gp = gpar(type = "none"),
               col = col_fun,
               na_col = "grey93",
               cluster_columns = F, show_row_names = T,
               cluster_rows = F, show_row_dend = F, show_column_dend = F, row_dend_side = "left", column_dend_side = "top",
               clustering_method_rows = "ward.D2", clustering_method_columns  = "ward.D2",
               clustering_distance_columns  = "euclidean", clustering_distance_rows  = "euclidean",
               row_names_side = "left",
               column_names_side = "top", column_names_rot = 55, column_names_gp = gpar(fontsize = 10, fontface = "plain"),
               row_names_gp = gpar(fontsize = 8, fontface = "plain"),
               column_title = NULL, column_title_side = "top", column_title_gp = gpar(fontsize = 10, fontface = "plain"),
               show_heatmap_legend = F,
               cell_fun =  function(j, i, x, y, width, height, fill){
                 grid.rect(x = x, y = y, width = width, height =  height, gp = gpar(lwd = 0.2, col = "grey80", fill = "grey67", alpha = 0.2))
                 #grid.circle(x = x, y = y, r = 1.2 * mat.N[i,j]  * min(unit.c(width, height)), gp = gpar(fill = col_fun(mat.expr[i, j]), col = NA, alpha = 0.9))},
                 grid.circle(x = x, y = y, r = 1.1 * min(unit.c(width, height)), gp = gpar(fill = col_fun(mat.expr[i, j]), col = NA, alpha = 0.9))},
               heatmap_legend_param=list(title = "expression", title_gp = gpar( fontsize = 10), labels_gp = gpar(fontsize = 8), legend_width = unit(1.8, "line"), legend_height = unit(.2, "line"), grid_width = unit(1.5, "line"), grid_height = unit(.1, "line"), at = c(1, max(mat.expr)), labels = c("low", "high"), by_row = F, direction = "horizontal", title_position = "lefttop"))

###############################################################


###############################################################
cairo_pdf(file = "/labs/khatrilab/ananthg/detrimental_host_response/figures_2024_09_27/DHR.singlecell.figure.pdf",width = 5.8, height = 6.5)
pushViewport(viewport(layout =grid.layout(nrow= 110, ncol =48)))

print(p.umap.score$detrimental, vp = viewport(layout.pos.row = 42:60, layout.pos.col = 13:21))
print(p.umap.score$protective, vp = viewport(layout.pos.row = 59:77, layout.pos.col = 13:21))

print(p.umap.score$SoMmodule1, vp = viewport(layout.pos.row = 42:60, layout.pos.col = 24:32))
print(p.umap.score$SoMmodule2, vp = viewport(layout.pos.row = 42:60, layout.pos.col = 34:42))
print(p.umap.score$SoMmodule3, vp = viewport(layout.pos.row = 59:77, layout.pos.col = 24:32))
print(p.umap.score$SoMmodule4, vp = viewport(layout.pos.row = 59:77, layout.pos.col = 34:42))

print(p.umap.score$SoMScore, vp = viewport(layout.pos.row = 50:68, layout.pos.col = 2:10))
print(p.umap.celltype, vp = viewport(layout.pos.row = 2:41, layout.pos.col = 19:39))
print(p.umap.pid, vp = viewport(layout.pos.row = 2:21, layout.pos.col = 2:11))
print(p.umap.group, vp = viewport(layout.pos.row = 22:41, layout.pos.col = 2:11))

pushViewport(viewport(layout.pos.row = 76:110, layout.pos.col = 4:40))
draw(p.SoMmodules, merge_legend =T, heatmap_legend_side = "bottom", newpage = FALSE)
upViewport(1)

col_fun.label1 = circlize::colorRamp2(seq(1,3, length =7), c("grey93", dichromat_pal("DarkRedtoBlue.12")(12)[c(7:12)]))
lgd = Legend(at = c(1.1,2.9), labels = c("low", "high"),
             col_fun = col_fun.label1, title = "score",
             grid_width = unit(0.1, "line"), grid_height = unit(0.1, "line"),
             legend_height = unit(3, "line"),
             direction = "vertical", title_position = "leftcenter-rot", title_gp = gpar(fontsize = 10, fontface = "plain"), labels_gp = gpar(fontsize = 10))
draw(lgd, x = unit(0.92, "npc"), y = unit(0.47, "npc"))

col_fun.label2 = circlize::colorRamp2(seq(-1,1, length = 13), c(dichromat_pal("DarkRedtoBlue.12")(12)[1:6], "white", dichromat_pal("DarkRedtoBlue.12")(12)[7:12]))
lgd = Legend(at = c(-0.7,0.7), labels = c("low", "high"),
             col_fun = col_fun.label2, title = "scaled score",
             grid_width = unit(0.1, "line"), grid_height = unit(0.1, "line"),
             legend_height = unit(3, "line"),
             direction = "vertical", title_position = "leftcenter-rot", title_gp = gpar(fontsize = 10, fontface = "plain"), labels_gp = gpar(fontsize = 10))
#draw(lgd, x = unit(0.1, "npc"), y = unit(0.23, "npc"))
draw(lgd, x = unit(0.88, "npc"), y = unit(0.08, "npc"))

grid.text(x= unit(0.02, "npc"), y = unit(0.98, "npc"), label = "A", gp=gpar(fontsize = textsize, fontface = "bold"))
grid.text(x= unit(0.02, "npc"), y = unit(0.8, "npc"), label = "B", gp=gpar(fontsize = textsize, fontface = "bold"))
grid.text(x= unit(0.37, "npc"), y = unit(0.98, "npc"), label = "C", gp=gpar(fontsize = textsize, fontface = "bold"))
grid.text(x= unit(0.02, "npc"), y = unit(0.61, "npc"), label = "D", gp=gpar(fontsize = textsize, fontface = "bold"))
grid.text(x= unit(0.02, "npc"), y = unit(0.28, "npc"), label = "E", gp=gpar(fontsize = textsize, fontface = "bold"))

grid.lines(x = unit(c(0.345, 0.355), "npc"), y = unit(c(0.47, 0.47), "npc"), gp = gpar(fill = "grey42", lwd = 1.5, col = "grey17"))
grid.text(x= unit(0.23, "npc"), y = unit(c(0.47), "npc"), label = "= ", gp=gpar(fontsize = 14, fontface = "plain", col = "grey17"))
grid.text(x= unit(0.46, "npc"), y = unit(0.39, "npc"), label = "= ", gp=gpar(fontsize = 14, fontface = "plain", col = "grey17"))
grid.text(x= unit(0.46, "npc"), y = unit(0.54, "npc"), label = "= ", gp=gpar(fontsize = 14, fontface = "plain", col = "grey17"))
grid.text(x= unit(0.68, "npc"), y = unit(0.39, "npc"), label = "+", gp=gpar(fontsize = 14, fontface = "plain", col = "grey17"))
grid.text(x= unit(0.68, "npc"), y = unit(0.54, "npc"), label = "+", gp=gpar(fontsize = 14, fontface = "plain", col = "grey17"))

dev.off()
###############################################################


###############################################################
c(seq(-2.3, 0, length.out = 7)[1:6], seq(0, 4.4, length.out = 7))

c("#2A0BD9", "#264EFF", "#40A1FF", "#73DAFF", "#ABF8FF", "#E0FFFF", 
"white", "#FFFFBF", "#FFE099", "#FFAD73", "#F76E5E", "#D92632", 
"#A60021")

summary(df4plot.f$SoMScore)

summary(df4plot.f$SoMmodule1)
summary(df4plot.f$SoMmodule2)
summary(df4plot.f$SoMmodule3)
summary(df4plot.f$SoMmodule4)

summary(df4plot.f$detrimental)
summary(df4plot.f$protective)

p.umap.score.bygroup = list()
for(i in names(scores2plot)){
  for(j in c("COVID-19,severe", "COVID-19,nonsevere", "COVID-19,conv", "bacterial sepsis", 
             "other infection,nonsevere", "healthy")){
    df4plot.f$facetvar = i
    p.umap.score.bygroup[[j]][[scores2plot[i]]]=
      ScatterPlotRast(df = as.data.frame(df4plot.f[groupseverity %in% j]), featureX = "UMAP_1", featureY = "UMAP_2", shuffle = T, pointsize = pointsizen, alphalevel = 0.2, colorby = scores2plot[i], groupby = "facetvar", legendColumn = 4, colcateg = "continuous", numberofrows = 5, legendposition = "right", seed = 13, mycol = c("grey95", dichromat_pal("DarkRedtoBlue.12")(12)[c(7:12)]), backgroudcol = "white", plottitle = "", rasterdpi = 300, nacol = "grey42", facetlabelface = "plain") +
      theme(text = element_text(size = textsize),
            #legend.position = c(0.5,-0.06),
            legend.position = "bottom",
            legend.justification = c(0.5,0.5),
            legend.text = element_text(size = textsize / 2),
            legend.title = element_blank(),
            legend.key.width = unit(0.5, "line"),
            legend.key.height = unit(0.2, "line"),
            legend.box.background =  element_blank(),
            #legend.title.align = 0.5,
            legend.margin = margin(-2, 1, 1, 1),
            #plot.margin = margin(1, 1, 1, 1),
            legend.direction = "horizontal",
            legend.spacing.x = unit(.2, "line"),
            legend.spacing.y = unit(.2, "line"))
  }
}
###############################################################


###############################################################
p.umap.score.bytissue = list()
schrepping_pbmc <- df4plot.f[df4plot.f$pid == "schrepping2020" & grepl("PBMC", df4plot.f$sampleid), ]
schrepping_wb <- df4plot.f[df4plot.f$pid == "schrepping2020" & !grepl("PBMC", df4plot.f$sampleid), ]
for(i in names(scores2plot)){
    schrepping_pbmc$facetvar = i
    p.umap.score.bytissue[["PBMC"]][[scores2plot[i]]]=
      ScatterPlotRast(df = as.data.frame(schrepping_pbmc), featureX = "UMAP_1", featureY = "UMAP_2", shuffle = T, pointsize = pointsizen, alphalevel = 0.2, colorby = scores2plot[i], groupby = "facetvar", legendColumn = 4, colcateg = "continuous", numberofrows = 5, legendposition = "right", seed = 13, mycol = c("grey95", dichromat_pal("DarkRedtoBlue.12")(12)[c(7:12)]), backgroudcol = "white", plottitle = "", rasterdpi = 300, nacol = "grey42", facetlabelface = "plain") +
      theme(text = element_text(size = textsize),
            #legend.position = c(0.5,-0.06),
            legend.position = "bottom",
            legend.justification = c(0.5,0.5),
            legend.text = element_text(size = textsize / 2),
            legend.title = element_blank(),
            legend.key.width = unit(0.5, "line"),
            legend.key.height = unit(0.2, "line"),
            legend.box.background =  element_blank(),
            #legend.title.align = 0.5,
            legend.margin = margin(-2, 1, 1, 1),
            #plot.margin = margin(1, 1, 2, 1),
            legend.direction = "horizontal",
            legend.spacing.x = unit(.2, "line"),
            legend.spacing.y = unit(.2, "line"))
}
for(i in names(scores2plot)){
    schrepping_wb$facetvar = i
    p.umap.score.bytissue[["WB"]][[scores2plot[i]]]=
      ScatterPlotRast(df = as.data.frame(schrepping_wb), featureX = "UMAP_1", featureY = "UMAP_2", shuffle = T, pointsize = pointsizen, alphalevel = 0.2, colorby = scores2plot[i], groupby = "facetvar", legendColumn = 4, colcateg = "continuous", numberofrows = 5, legendposition = "right", seed = 13, mycol = c("grey95", dichromat_pal("DarkRedtoBlue.12")(12)[c(7:12)]), backgroudcol = "white", plottitle = "", rasterdpi = 300, nacol = "grey42", facetlabelface = "plain") +
      theme(text = element_text(size = textsize),
            #legend.position = c(0.5,-0.06),
            legend.position = "bottom",
            legend.justification = c(0.5,0.5),
            legend.text = element_text(size = textsize / 2),
            legend.title = element_blank(),
            legend.key.width = unit(0.5, "line"),
            legend.key.height = unit(0.2, "line"),
            legend.box.background =  element_blank(),
            #legend.title.align = 0.5,
            legend.margin = margin(-2, 1, 1, 1),
            #plot.margin = margin(1, 1, 2, 1),
            legend.direction = "horizontal",
            legend.spacing.x = unit(.2, "line"),
            legend.spacing.y = unit(.2, "line"))
}
###############################################################


###############################################################
#source("/labs/khatrilab/ananthg/detrimental_host_response/common_functions/JT_test_function.R")
library(ggbeeswarm)

group_wise_prop <- df4plot.f %>% group_by(sampleid, celltype) %>% summarise(count = n(), severity = first(groupseverity), median_som = median(SoMScore)) %>% mutate(prop = count / sum(count))
group_wise_prop_neutrophils <- group_wise_prop[group_wise_prop$celltype %in% c("Immat Neutrophil", "Neutrophil"), ]
group_wise_prop_neutrophils <- group_wise_prop_neutrophils[!(group_wise_prop_neutrophils$celltype == "Neutrophil" & grepl("PBMC", group_wise_prop_neutrophils$sampleid)), ]
group_wise_prop_neutrophils <- group_wise_prop_neutrophils[!(group_wise_prop_neutrophils$celltype == "Neutrophil" & group_wise_prop_neutrophils$prop < 0.1), ]
group_wise_prop_neutrophils$severity <- as.character(group_wise_prop_neutrophils$severity)
group_wise_prop_neutrophils$severity <- ordered(group_wise_prop_neutrophils$severity, levels = c("healthy", "COVID-19,conv", "other infection,nonsevere", "COVID-19,nonsevere", "bacterial sepsis", "COVID-19,severe"))
group_wise_prop_neutrophils$numeric_severity <- as.numeric(group_wise_prop_neutrophils$severity)

plot_cor_test <- function(score_pheno_frame, score_name, score_name_plot)
{
  #score_pheno_frame$scaled_score <- scale(score_pheno_frame[[score_name]])
  score_summary <- ggplot(score_pheno_frame, aes_string(x = "severity", y = score_name, colour = "severity")) +
    ylab(score_name_plot) +
    geom_beeswarm(cex = 1.25, size = pointsizen) +
    geom_boxplot(fill = NA, lwd = linesize, width = 0.4, outlier.shape = NA, colour = "#303030") +
    #scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank(),
          legend.position = "bottom", legend.direction = "horizontal",
          panel.grid.minor = element_blank(), plot.background = element_blank()) +
    stat_cor(aes_string(x = "numeric_severity", y = score_name), method = "pearson", label.x = "healthy", label.y = Inf, vjust = 1.5, size = textsize / 5, colour = "#505050") +
    scale_colour_manual(values = col_groupseverity) +#, labels = severity_labels) +
    labs(colour = NULL) + guides(colour = guide_legend(nrow = 2, override.aes = list(size = (textsize / 6) * 6 / 9))) +
    facet_wrap(vars(celltype), nrow = 1)
  return(score_summary)
}
som_group_plot <- plot_cor_test(group_wise_prop_neutrophils, "median_som", "SoM score")
prop_group_plot <- plot_cor_test(group_wise_prop_neutrophils, "prop", "Proportion") + theme(legend.position = "none")

som_prop_plot <- ggplot(group_wise_prop_neutrophils, aes(x = prop, y = median_som, colour = severity)) +
  ylab("SoM score") + xlab("Proportion") +
  geom_point(size = pointsizen) +
  #scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
  theme(legend.position = "none",
        panel.grid.minor = element_blank(), plot.background = element_blank()) +
  stat_cor(aes(x = prop, y = median_som), method = "pearson", label.x = 0, label.y = Inf, vjust = 1.5, size = textsize / 5, colour = "#505050") +
  scale_colour_manual(values = col_groupseverity) +
  facet_wrap(vars(celltype), nrow = 1)
###############################################################


###############################################################
cairo_pdf(file = "/labs/khatrilab/ananthg/detrimental_host_response/figures_2024_09_27/DHR.singlecell.figure.sup3.pdf",width = 12, height = 22)
pushViewport(viewport(layout =grid.layout(nrow= 300, ncol = 86)))

# COVID-19, severe
print(p.umap.score.bygroup$`COVID-19,severe`$detrimental, vp = viewport(layout.pos.row = 2:20, layout.pos.col = 13:21))
print(p.umap.score.bygroup$`COVID-19,severe`$protective, vp = viewport(layout.pos.row = 20:38, layout.pos.col = 13:21))

print(p.umap.score.bygroup$`COVID-19,severe`$SoMmodule1, vp = viewport(layout.pos.row = 2:20, layout.pos.col = 24:32))
print(p.umap.score.bygroup$`COVID-19,severe`$SoMmodule2, vp = viewport(layout.pos.row = 2:20, layout.pos.col = 34:42))
print(p.umap.score.bygroup$`COVID-19,severe`$SoMmodule3, vp = viewport(layout.pos.row = 20:38, layout.pos.col = 24:32))
print(p.umap.score.bygroup$`COVID-19,severe`$SoMmodule4, vp = viewport(layout.pos.row = 20:38, layout.pos.col = 34:42))

print(p.umap.score.bygroup$`COVID-19,severe`$SoMScore, vp = viewport(layout.pos.row = 10:28, layout.pos.col = 2:10))

pushViewport(viewport(layout.pos.row = 2:38, layout.pos.col = 2:42))
grid.rect(gp = gpar(col = "#505050", fill = "transparent"), width = 0.99, height = 0.98)
grid.text(x= unit(0.1, "npc"), y = unit(0.82, "npc"), label = "COVID-19,\nsevere", gp=gpar(fontsize = textsize, fontface = "plain"))
upViewport(1)


# COVID-19, nonsevere
print(p.umap.score.bygroup$`COVID-19,nonsevere`$detrimental, vp = viewport(layout.pos.row = 39:57, layout.pos.col = 13:21))
print(p.umap.score.bygroup$`COVID-19,nonsevere`$protective, vp = viewport(layout.pos.row = 57:75, layout.pos.col = 13:21))

print(p.umap.score.bygroup$`COVID-19,nonsevere`$SoMmodule1, vp = viewport(layout.pos.row = 39:57, layout.pos.col = 24:32))
print(p.umap.score.bygroup$`COVID-19,nonsevere`$SoMmodule2, vp = viewport(layout.pos.row = 39:57, layout.pos.col = 34:42))
print(p.umap.score.bygroup$`COVID-19,nonsevere`$SoMmodule3, vp = viewport(layout.pos.row = 57:75, layout.pos.col = 24:32))
print(p.umap.score.bygroup$`COVID-19,nonsevere`$SoMmodule4, vp = viewport(layout.pos.row = 57:75, layout.pos.col = 34:42))

print(p.umap.score.bygroup$`COVID-19,nonsevere`$SoMScore, vp = viewport(layout.pos.row = 47:65, layout.pos.col = 2:10))

pushViewport(viewport(layout.pos.row = 39:75, layout.pos.col = 2:42))
grid.rect(gp = gpar(col = "#505050", fill = "transparent"), width = 0.99, height = 0.98)
grid.text(x= unit(0.1, "npc"), y = unit(0.82, "npc"), label = "COVID-19,\nnonsevere", gp=gpar(fontsize = textsize, fontface = "plain"))
upViewport(1)


#COVID-19, conv
print(p.umap.score.bygroup$`COVID-19,conv`$detrimental, vp = viewport(layout.pos.row = 76:94, layout.pos.col = 13:21))
print(p.umap.score.bygroup$`COVID-19,conv`$protective, vp = viewport(layout.pos.row = 94:112, layout.pos.col = 13:21))

print(p.umap.score.bygroup$`COVID-19,conv`$SoMmodule1, vp = viewport(layout.pos.row = 76:94, layout.pos.col = 24:32))
print(p.umap.score.bygroup$`COVID-19,conv`$SoMmodule2, vp = viewport(layout.pos.row = 76:94, layout.pos.col = 34:42))
print(p.umap.score.bygroup$`COVID-19,conv`$SoMmodule3, vp = viewport(layout.pos.row = 94:112, layout.pos.col = 24:32))
print(p.umap.score.bygroup$`COVID-19,conv`$SoMmodule4, vp = viewport(layout.pos.row = 94:112, layout.pos.col = 34:42))

print(p.umap.score.bygroup$`COVID-19,conv`$SoMScore, vp = viewport(layout.pos.row = 84:102, layout.pos.col = 2:10))

pushViewport(viewport(layout.pos.row = 76:112, layout.pos.col = 2:42))
grid.rect(gp = gpar(col = "#505050", fill = "transparent"), width = 0.99, height = 0.98)
grid.text(x= unit(0.1, "npc"), y = unit(0.82, "npc"), label = "COVID-19,\nconv", gp=gpar(fontsize = textsize, fontface = "plain"))
upViewport(1)


# bacterial sepsis
print(p.umap.score.bygroup$`bacterial sepsis`$detrimental, vp = viewport(layout.pos.row = 2:20, layout.pos.col = 56:64))
print(p.umap.score.bygroup$`bacterial sepsis`$protective, vp = viewport(layout.pos.row = 20:38, layout.pos.col = 56:64))

print(p.umap.score.bygroup$`bacterial sepsis`$SoMmodule1, vp = viewport(layout.pos.row = 2:20, layout.pos.col = 67:75))
print(p.umap.score.bygroup$`bacterial sepsis`$SoMmodule2, vp = viewport(layout.pos.row = 2:20, layout.pos.col = 77:85))
print(p.umap.score.bygroup$`bacterial sepsis`$SoMmodule3, vp = viewport(layout.pos.row = 20:38, layout.pos.col = 67:75))
print(p.umap.score.bygroup$`bacterial sepsis`$SoMmodule4, vp = viewport(layout.pos.row = 20:38, layout.pos.col = 77:85))

print(p.umap.score.bygroup$`bacterial sepsis`$SoMScore, vp = viewport(layout.pos.row = 10:28, layout.pos.col = 45:53))

pushViewport(viewport(layout.pos.row = 2:38, layout.pos.col = 45:85))
grid.rect(gp = gpar(col = "#505050", fill = "transparent"), width = 0.99, height = 0.98)
grid.text(x= unit(0.1, "npc"), y = unit(0.82, "npc"), label = "bacterial sepsis", gp=gpar(fontsize = textsize, fontface = "plain"))
upViewport(1)


# other infection, nonsevere
print(p.umap.score.bygroup$`other infection,nonsevere`$detrimental, vp = viewport(layout.pos.row = 39:57, layout.pos.col = 56:64))
print(p.umap.score.bygroup$`other infection,nonsevere`$protective, vp = viewport(layout.pos.row = 57:75, layout.pos.col = 56:64))

print(p.umap.score.bygroup$`other infection,nonsevere`$SoMmodule1, vp = viewport(layout.pos.row = 39:57, layout.pos.col = 67:75))
print(p.umap.score.bygroup$`other infection,nonsevere`$SoMmodule2, vp = viewport(layout.pos.row = 39:57, layout.pos.col = 77:85))
print(p.umap.score.bygroup$`other infection,nonsevere`$SoMmodule3, vp = viewport(layout.pos.row = 57:75, layout.pos.col = 67:75))
print(p.umap.score.bygroup$`other infection,nonsevere`$SoMmodule4, vp = viewport(layout.pos.row = 57:75, layout.pos.col = 77:85))

print(p.umap.score.bygroup$`other infection,nonsevere`$SoMScore, vp = viewport(layout.pos.row = 47:65, layout.pos.col = 45:53))

pushViewport(viewport(layout.pos.row = 39:75, layout.pos.col = 45:85))
grid.rect(gp = gpar(col = "#505050", fill = "transparent"), width = 0.99, height = 0.98)
grid.text(x= unit(0.1, "npc"), y = unit(0.82, "npc"), label = "other infection,\nnonsevere", gp=gpar(fontsize = textsize, fontface = "plain"))
upViewport(1)


#healthy
print(p.umap.score.bygroup$`healthy`$detrimental, vp = viewport(layout.pos.row = 76:94, layout.pos.col = 56:64))
print(p.umap.score.bygroup$`healthy`$protective, vp = viewport(layout.pos.row = 94:112, layout.pos.col = 56:64))

print(p.umap.score.bygroup$`healthy`$SoMmodule1, vp = viewport(layout.pos.row = 76:94, layout.pos.col = 67:75))
print(p.umap.score.bygroup$`healthy`$SoMmodule2, vp = viewport(layout.pos.row = 76:94, layout.pos.col = 77:85))
print(p.umap.score.bygroup$`healthy`$SoMmodule3, vp = viewport(layout.pos.row = 94:112, layout.pos.col = 67:75))
print(p.umap.score.bygroup$`healthy`$SoMmodule4, vp = viewport(layout.pos.row = 94:112, layout.pos.col = 77:85))

print(p.umap.score.bygroup$`healthy`$SoMScore, vp = viewport(layout.pos.row = 84:102, layout.pos.col = 45:53))

pushViewport(viewport(layout.pos.row = 76:112, layout.pos.col = 45:85))
grid.rect(gp = gpar(col = "#505050", fill = "transparent"), width = 0.99, height = 0.98)
grid.text(x= unit(0.1, "npc"), y = unit(0.82, "npc"), label = "healthy", gp=gpar(fontsize = textsize, fontface = "plain"))
upViewport(1)



# PBMC, schrepping
print(p.umap.score.bytissue$PBMC$detrimental, vp = viewport(layout.pos.row = 113:131, layout.pos.col = 13:21))
print(p.umap.score.bytissue$PBMC$protective, vp = viewport(layout.pos.row = 131:149, layout.pos.col = 13:21))

print(p.umap.score.bytissue$PBMC$SoMmodule1, vp = viewport(layout.pos.row = 113:131, layout.pos.col = 24:32))
print(p.umap.score.bytissue$PBMC$SoMmodule2, vp = viewport(layout.pos.row = 113:131, layout.pos.col = 34:42))
print(p.umap.score.bytissue$PBMC$SoMmodule3, vp = viewport(layout.pos.row = 131:149, layout.pos.col = 24:32))
print(p.umap.score.bytissue$PBMC$SoMmodule4, vp = viewport(layout.pos.row = 131:149, layout.pos.col = 34:42))

print(p.umap.score.bytissue$PBMC$SoMScore, vp = viewport(layout.pos.row = 121:139, layout.pos.col = 2:10))

pushViewport(viewport(layout.pos.row = 113:149, layout.pos.col = 2:42))
grid.text(x= unit(0.1, "npc"), y = unit(0.82, "npc"), label = "PBMC\n(schrepping2020)", gp=gpar(fontsize = textsize, fontface = "plain"))
#grid.rect(gp = gpar(col = "#505050", fill = "transparent"), width = 0.99, height = 0.98)
upViewport(1)


# WB, schrepping
print(p.umap.score.bytissue$WB$detrimental, vp = viewport(layout.pos.row = 113:131, layout.pos.col = 56:64))
print(p.umap.score.bytissue$WB$protective, vp = viewport(layout.pos.row = 131:149, layout.pos.col = 56:64))

print(p.umap.score.bytissue$WB$SoMmodule1, vp = viewport(layout.pos.row = 113:131, layout.pos.col = 67:75))
print(p.umap.score.bytissue$WB$SoMmodule2, vp = viewport(layout.pos.row = 113:131, layout.pos.col = 77:85))
print(p.umap.score.bytissue$WB$SoMmodule3, vp = viewport(layout.pos.row = 131:149, layout.pos.col = 67:75))
print(p.umap.score.bytissue$WB$SoMmodule4, vp = viewport(layout.pos.row = 131:149, layout.pos.col = 77:85))

print(p.umap.score.bytissue$WB$SoMScore, vp = viewport(layout.pos.row = 121:139, layout.pos.col = 45:53))

pushViewport(viewport(layout.pos.row = 113:149, layout.pos.col = 45:85))
grid.text(x= unit(0.1, "npc"), y = unit(0.82, "npc"), label = "WB\n(schrepping2020)", gp=gpar(fontsize = textsize, fontface = "plain"))
#grid.rect(gp = gpar(col = "#505050", fill = "transparent"), width = 0.99, height = 0.98)
upViewport(1)


pushViewport(viewport(layout.pos.row = 150:300, layout.pos.col = 47:85))
draw(p.42, merge_legend =T, heatmap_legend_side = "bottom", height = height_p42, newpage = FALSE)
grid.text(x= unit(0.95, "npc"), y = unit(0.78, "npc"), label = "module 1", rot =90, gp=gpar(fontsize = textsize, fontface = "plain", col = "#EFC002"))
grid.text(x= unit(0.95, "npc"), y = unit(0.52, "npc"), label = "module 2", rot =90, gp=gpar(fontsize = textsize, fontface = "plain", col = "#F39B7F"))
grid.text(x= unit(0.95, "npc"), y = unit(0.29, "npc"), label = "module 3", rot =90, gp=gpar(fontsize = textsize, fontface = "plain", col = "#00A087"))
grid.text(x= unit(0.95, "npc"), y = unit(0.12, "npc"), label = "module 4", rot =90, gp=gpar(fontsize = textsize, fontface = "plain", col = "#3C5488"))
upViewport(1)


print(som_group_plot, vp = viewport(layout.pos.row = 155:210, layout.pos.col = 2:44))
print(prop_group_plot, vp = viewport(layout.pos.row = 215:255, layout.pos.col = 2:44))
print(som_prop_plot, vp = viewport(layout.pos.row = 260:300, layout.pos.col = 2:44))

grid.text(x= unit(0.0075, "npc"), y = unit(0.9925, "npc"), label = "A", gp=gpar(fontsize = textsize, fontface = "bold"))
grid.text(x= unit(0.0075, "npc"), y = unit(0.6225, "npc"), label = "B", gp=gpar(fontsize = textsize, fontface = "bold"))
grid.text(x= unit(0.0075, "npc"), y = unit(0.485, "npc"), label = "C", gp=gpar(fontsize = textsize, fontface = "bold"))
grid.text(x= unit(0.0075, "npc"), y = unit(0.2875, "npc"), label = "D", gp=gpar(fontsize = textsize, fontface = "bold"))
grid.text(x= unit(0.0075, "npc"), y = unit(0.1375, "npc"), label = "E", gp=gpar(fontsize = textsize, fontface = "bold"))
grid.text(x= unit(0.55, "npc"), y = unit(0.49, "npc"), label = "F", gp=gpar(fontsize = textsize, fontface = "bold"))


dev.off()
###############################################################

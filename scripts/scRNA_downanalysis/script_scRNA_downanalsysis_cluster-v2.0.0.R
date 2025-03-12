#!/bin/Rscript
## 参数整理
library(optparse)
option_list <- list(
  make_option(c("-i", "--seurat_rds"), type = "character", help  = "Single Cell Seurat Data (RDS FILE)"),
  make_option(c("-o", "--out"), type = "character", help  = "Output dir"),
  make_option(c("-t", "--n_cpu_source"), type = "integer", default = 8, help  = "Thread count (integer;  default = 8)"),
  make_option(c("-s", "--species"), type = "character", help  = "species of human or mouse"),
  make_option(c("--memory"), type = "integer", default = 10, help  = "Memory consumption per thread (integer; default = 10)"),
  make_option(c("--num_feature_low_threshold"), type = "integer", default = 200, help  = "The number of genes detected is low threshold (integer; default = 200)"),
  make_option(c("--num_feature_high_threshold"), type = "integer", default = 10000, help  = "High threshold for the number of genes detected (integer; default = 10000)"),
  make_option(c("--mt_gene_threshold"), type = "integer", default = 30, help  = "Mitochondrial percentage threshold (integer %; default = 30)"),
  make_option(c("--hb_gene_threshold"), type = "integer", default = 5, help  = "Mitochondrial percentage threshold (integer %; default = 5)"),
  make_option(c("--marker"), type = "character", help  = "An xlsx file containing marker information (the first line is the cell type name, followed by the corresponding marker)")
)
args <-
  parse_args(
    OptionParser(option_list = option_list, usage = "Usage: %prog [options] \nDescription: Single cell analysis - broad class analysis!")
  )
# version 2:0:0
## 质控出图，添加每个样本的质控图片结果
## 增加双细胞去除DoubletFinder

# 记录开始时间
start_time <- Sys.time()

# 获取当前项目src配置文件(不能动)
source("~/project/scRNA_downanalysis/src/src.R")


using(tidyverse, ggplot2, scales, gridExtra, Matrix, Seurat, harmony, cowplot, future, DoubletFinder, openxlsx, clustree)
set.seed(2024)

species <- args$species
seurat_raw <- readRDS(args$seurat_rds)

## create dir
out_dir <- args$out
dir.create(out_dir, recursive = T)
dir_list <- c('01_figs', "02_tabs", '03_files')
lapply(dir_list, function(x){dir.create(paste(sep = "/", out_dir, x), recursive = T)})
out <- list()
out[['dir']] <- out_dir
out[['figs']] <- paste(sep = "/", out_dir, '01_figs')
out[['tabs']] <- paste(sep = "/", out_dir, '02_tabs')
out[['files']] <- paste(sep = "/", out_dir, '03_files')

## function
human_QC_analysis <- function(seurat_raw, dr_gene, hb_gene, S_gene, G2M_gene){
  message('**We are counting QC index...')
  seurat_raw[["percent.mt"]] <- PercentageFeatureSet(seurat_raw, pattern = "^MT-")
  seurat_raw[["percent.HB"]] <- PercentageFeatureSet(seurat_raw, features = hb_gene)
  seurat_raw <- CellCycleScoring(object = seurat_raw, g2m.features = G2M_gene, s.features = S_gene)
  seurat_raw <- AddModuleScore(seurat_raw, features = list(dr_gene), name = 'Dissociation_reaction_marker')
  return(seurat_raw)
}

## marker
human_common_marker_list <- lapply(as.data.frame(read.xlsx(args$marker)), function(x) x[complete.cases(x)])

## 物种判断及相关数据匹配读取
if (species == 'human') {
  message('**The defined species is human...')
  dr_gene <- read.table(config$`Others File`$Dissociation_reaction_marker)[,1]
  hb_gene <- read.table(config$`Others File`$Erythrocyte_marker)[,1]
  S_gene <- read.table(config$`Others File`$Cell_cycle_S_marker, header = T)[,1]
  G2M_gene <- read.table(config$`Others File`$Cell_cycle_G2M_marker, header = T)[,1]
  
  hb_gene <- hb_gene[hb_gene %in% rownames(seurat_raw@assays$RNA)]
  dr_gene <- dr_gene[dr_gene %in% rownames(seurat_raw@assays$RNA)]
  seurat_deal <- human_QC_analysis(seurat_raw, dr_gene, hb_gene, S_gene, G2M_gene)
}else if (species == 'mouse') {
  # seurat_raw[["percent.mt"]] <- PercentageFeatureSet(seurat_raw, pattern = "^mt-")
  stop('The mouse pipeline needs to be perfected')
}else{
  stop('species is not a parameter of either human or mouse')
}



# 解离密度图及5%过滤点
P <- ggplot(seurat_deal@meta.data, aes(x = Dissociation_reaction_marker1)) + 
  geom_density(fill = "skyblue", alpha = 0.5, color = "darkblue", size = 1) +  # 更改密度曲线和填充颜色
  # geom_vline(xintercept = 2, color = "blue", linetype = "dashed", size = 1) +  # 在 x = 2 位置添加垂直线
  labs(title = "Density Plot with Custom Vertical Lines",
       x = "Dissociation Reaction Marker 1",
       y = "Density") +
  theme_minimal(base_size = 15) +  # 使用简约主题并调整基准字号
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 20),  # 标题居中，粗体，大小调整
    axis.title.x = element_text(face = "bold", size = 15),  # x 轴标题加粗并调整大小
    axis.title.y = element_text(face = "bold", size = 15),  # y 轴标题加粗并调整大小
    axis.text = element_text(size = 12)  # 坐标轴刻度文字大小调整
  )


# 计算右百分之5的位置
x_95 <- quantile(seurat_deal@meta.data$Dissociation_reaction_marker1, 0.95)

# 在右百分之5的位置添加垂直线并标注 x 数值
P <- P + 
  geom_vline(xintercept = x_95, color = "red", linetype = "dashed", size = 1) +
  annotate("text", x = x_95, y = Inf, label = paste("x =", round(x_95, 2)), vjust = 1.5, color = "red", size = 5, fontface = "bold")
save_picture(P, out = out$figs, "01-Dissociation_reaction_score", width_len = 900, height_len = 600)



# 过滤前
P <- VlnPlot(seurat_deal, features = c('nCount_RNA', 'nFeature_RNA', 'percent.mt', 'percent.HB', 'Dissociation_reaction_marker1'), pt.size = 0)
save_picture(P, out = out$figs, "02-QC_before_vlnplot", width_len = 900, height_len = 600)


# 添加按照每个样本信息质控图 v.2.0.0
lapply(unique(seurat_deal$orig.ident), function(x){
  subseurat_temp <- subset(seurat_deal, orig.ident == x)
  P <- VlnPlot(subseurat_temp, c('nCount_RNA', 'nFeature_RNA', 'percent.mt', 'percent.HB', 'Dissociation_reaction_marker1'), pt.size = 0)
  save_picture(P, out = out$figs, paste0("02-QC-sub_before_vlnplot-", x), width_len = 900, height_len = 600)
})


# QC 过滤
message('**Filtering in progress...')
num_feature_low_threshold <- args$num_feature_low_threshold
num_feature_high_threshold <- args$num_feature_high_threshold
mt_gene_threshold <- args$mt_gene_threshold
hb_gene_threshold <- args$hb_gene_threshold

message(paste0('Filtering info: 
feature low : ', num_feature_low_threshold, '
feature high : ', num_feature_high_threshold, '
mt threshold : ', mt_gene_threshold, '
hb threshold : ', hb_gene_threshold, "
Dissociation reaction threshold : ", signif(x_95, 4)))

seurat_deal_filter <-
  subset(seurat_deal,
         subset = nFeature_RNA > num_feature_low_threshold &
           nFeature_RNA < num_feature_high_threshold & percent.mt < mt_gene_threshold & Dissociation_reaction_marker1 < x_95 & percent.HB < hb_gene_threshold)



# 过滤后
P <- VlnPlot(seurat_deal_filter, features = c('nCount_RNA', 'nFeature_RNA', 'percent.mt', 'percent.HB', 'Dissociation_reaction_marker1'), pt.size = 0)
save_picture(P, out = out$figs, "03-QC_after_vlnplot", width_len = 900, height_len = 600)

# 添加按照每个样本信息质控图 v.2.0.0
lapply(unique(seurat_deal_filter$orig.ident), function(x){
  subseurat_temp <- subset(seurat_deal_filter, orig.ident == x)
  P <- VlnPlot(subseurat_temp, c('nCount_RNA', 'nFeature_RNA', 'percent.mt', 'percent.HB', 'Dissociation_reaction_marker1'), pt.size = 0)
  save_picture(P, out = out$figs, paste0("03-QC-sub_after_vlnplot-", x), width_len = 900, height_len = 600)
})


message('**QC scatter...') 
index <- 1
features <- c('nCount_RNA', 'nFeature_RNA', 'percent.mt', 'percent.HB', 'Dissociation_reaction_marker1')
for (i in seq_along(features)) {
  for (k in (i+1):length(features)) {
    # 获取实际的特征名称
    feature1 <- features[i]
    feature2 <- features[k]
    if (is.na(feature2)) {
      break 
    }
    
    # 生成图形
    P <- (FeatureScatter(seurat_deal, feature1 = feature1, feature2 = feature2, pt.size = 0.1) + 
            ylab(paste0(feature2, " - before QC")) + xlab(paste0(feature1, " - before QC"))) + 
      (FeatureScatter(seurat_deal_filter, feature1 = feature1, feature2 = feature2, pt.size = 0.1) + 
         ylab(paste0(feature2, " - after QC")) + xlab(paste0(feature1, " - after QC")))
    
    # 保存图形
    save_picture(P, out = out$figs, Pname = paste0("04-", index, "-QC_Scatter"), width_len = 900, height_len = 400)
    
    index <- index + 1
  }
}

# 标准化数据
message('**Normalization is under way...')
seurat_deal_filter <-
  NormalizeData(seurat_deal_filter,
                normalization.method = "LogNormalize",
                scale.factor = 10000)

t <- args$n_cpu_source
mem <- args$memory
if (t == 1) {
  plan("sequential")
} else{
  plan("multiprocess", workers = t)
  options(future.globals.maxSize = as.numeric(mem) * 1000 * 1024 ^ 2)
}

message('**High variant gene screening...') 
seurat_deal_filter <-
  FindVariableFeatures(seurat_deal_filter, selection.method = "vst", nfeatures = 2000)
P <- LabelPoints(VariableFeaturePlot(seurat_deal_filter, pt.size = 0.8), head(VariableFeatures(seurat_deal_filter), 10), repel = T)
save_picture(P, out = out$figs, Pname = "05-VariableFeatures", width_len = 600, height_len = 400)

message('**It is being standardized...') 
seurat_deal_filter <- ScaleData(seurat_deal_filter, vars.to.regress = c("percent.mt", "S.Score", "G2M.Score", "percent.HB", "Dissociation_reaction_marker1"))
plan("sequential")

message('**Pca dimension reduction is in progress')
seurat_deal_filter <-
  RunPCA(seurat_deal_filter, features = VariableFeatures(object = seurat_deal_filter), verbose = FALSE)

message('**Pca dimension reduction is in progress')
P <- VizDimLoadings(seurat_deal_filter, dims = 1:2)
save_picture(P, out = out$figs, Pname = "06-PC_coef", width_len = 600, height_len = 400)
P <- DimPlot(seurat_deal_filter, reduction = "pca", pt.size = 0.8)
save_picture(P, out = out$figs, Pname = "07-PC_dim", width_len = 600, height_len = 400)
P <- ElbowPlot(seurat_deal_filter, ndims = 50)
save_picture(P, out = out$figs, Pname = "08-PC_Elbow", width_len = 600, height_len = 400)

message('**Harmony correction in progress...')
seurat_deal_filter_harmony <- RunHarmony(seurat_deal_filter, "orig.ident", plot_convergence = TRUE)



# ElbowPlot(seurat_deal_filter_harmony, ndims = 50)
# stdevs  <- Stdev(object = seurat_deal_filter_harmony, reduction = "pca")
# # 计算第一和第二差数
# first_derivative <- diff(stdevs)
# second_derivative <- diff(first_derivative)
# # 寻找第二导数的最小值的位置
# dim.use <- which.min(second_derivative) + 3
dim.use <- 30

message('**clustering analysis in progress...')
seurat_deal_filter_harmony_cluster <- seurat_deal_filter_harmony %>%
  RunUMAP(reduction = "harmony", dims = 1:dim.use) %>%
  FindNeighbors(reduction = "harmony", dims = 1:dim.use) 


message('**Resolution debugging in progress...')
seurat_temp  <- seurat_deal_filter_harmony_cluster
for (res in seq(0.1, 2, by = 0.1)) {
  if (res == 0.1) {
    rm(res_aim)
  }
  seurat_temp <- FindClusters(seurat_temp, resolution = res)
  if (length(unique(seurat_temp@active.ident)) >= 20 & length(unique(seurat_temp@active.ident)) < 30 & !(exists('res_aim'))) {
    res_aim <- res
  }
}


P <- clustree(seurat_temp, prefix = 'RNA_snn_res.') + coord_flip()
save_picture(
  P,
  out =  out$figs,
  Pname = "09-clusttree",
  width_len = 1600,
  height_len = 900
)

if (length(unique(seurat_temp@active.ident)) < 20) {
  stop("Resolution 2 does not have the right number of clusters")
}

message(paste0("res_aim : ", res_aim))

message('**TSNE added in progress...')
seurat_deal_filter_harmony_cluster <- FindClusters(seurat_deal_filter_harmony_cluster, resolution = res_aim)
seurat_deal_filter_harmony_cluster <- seurat_deal_filter_harmony_cluster %>%
  RunTSNE(reduction = "harmony", dims = 1:dim.use)

# 双细胞去除DoubletFinder v.2.0.0
seurat_deal_filter_harmony_cluster <- doubleCell_filter(seurat_deal_filter_harmony_cluster, fig_out = out$figs, tab_out = out$tabs)



## marker
message('**marker information supplement in progress...')
dir.create(paste(sep = "/", out$files, "01-common_marker"),recursive = T)
if (species == 'human') {
  for (i in 1:length(human_common_marker_list)) {
    gene<-human_common_marker_list[[i]][human_common_marker_list[[i]] %in% rownames(seurat_deal_filter_harmony_cluster@assays$RNA)]
    if (length(gene) > 0){
      P<-VlnPlot(seurat_deal_filter_harmony_cluster, features = gene, pt.size = 0)
      save_picture(P, out = paste(sep = "/",out$files,"01-common_marker"), Pname = names(human_common_marker_list)[i], width_len = 1100, height_len = 900)
    }
  }
}else if (species == 'mouse') {
  library(homologene)
  for (i in 1:length(human_common_marker_list)) {
    gene<-human2mouse(human_common_marker_list[[i]])$mouseGene[human2mouse(human_common_marker_list[[i]])$mouseGene %in% rownames(seurat_deal_filter_harmony_cluster@assays$RNA)]
    if (length(gene) > 0){
      P<-VlnPlot(seurat_deal_filter_harmony_cluster, features = gene, pt.size = 0)
      save_picture(P, out = paste(sep = "/",out$files,"01-common_marker"), Pname = names(human_common_marker_list)[i], width_len = 1100, height_len = 900)
    }
  }
}else{
  stop('species nonexits')
}

P <- DimPlot(seurat_deal_filter_harmony_cluster, reduction = "umap", raster = FALSE, pt.size = 0.8)
save_picture(
  P,
  out =  out$figs,
  Pname = "11-unAnn cluster",
  width_len = 900,
  height_len = 600
)
P <- DimPlot(seurat_deal_filter_harmony_cluster, reduction = "umap", group.by = "orig.ident", raster = FALSE, pt.size = 0.8)
save_picture(
  P,
  out =  out$figs,
  Pname = "12-unAnn cluster by sample",
  width_len = 900,
  height_len = 600
)





message('**Working environment preservation...')
save.image(paste(sep = "/", out$dir, "env.Rdata"))

# rna_edit_out <- paste0(dirname(out$dir), '/03.RNAedit')
# dir.create(rna_edit_out, recursive = T)
# step2_after_ann(paste(sep = "/", out$dir, "env.Rdata"), rna_edit_out)

# 记录结束时间
end_time <- Sys.time()

# 计算运行时间
execution_time <- end_time - start_time

# 输出运行时间
print(paste("脚本运行时间:", signif(execution_time, 4)))

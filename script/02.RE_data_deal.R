library(optparse)
option_list <- list(
  make_option(c("-s", "--src"), type = "character", default = "/mnt/sda/project/yjc1/pipeline/02_scRNA_edit/script/src.R", help  = "configuration files, src.R. default: /mnt/sda/project/yjc1/pipeline/02_scRNA_edit/script/src.R"),
  make_option(c("-g", "--group"), type = "character", help  = "group file"),
  make_option(c("-o", "--out"), type = "character", help  = "Output dir"),
  make_option(c("--annovar"), type = "character", default = "/mnt/sda/project/yjc1/software/annovar/annotate_variation.pl", help  = "annotate_variation.pl. default: /mnt/sda/project/yjc1/software/annovar/annotate_variation.pl"),
  make_option(c("--ref"), type = "character", default = "/mnt/sda/project/yjc1/software/annovar/hg38", help  = "annotate ref. default: /mnt/sda/project/yjc1/software/annovar/hg38")
)
args <-
  parse_args(
    OptionParser(option_list = option_list, usage = "Usage: %prog [options] \nDescription: Single cell analysis - broad class analysis!")
  )
library(tidyverse)
library(Seurat)
library(dplyr)
library(tidyr)
library(stringr)
library(data.table)
library(harmony)
set.seed(2023)
calculate_davies_bouldin <- function(embeddings, clusters) {
  # 计算质心
  centroids <- aggregate(embeddings, by=list(cluster=clusters), FUN=mean)
  centroids <- centroids[, -1]
  
  # 计算每个点到质心的距离
  distances <- matrix(0, nrow=length(clusters), ncol=max(clusters))
  for (i in 1:max(clusters)) {
    cluster_points <- embeddings[clusters == i, ]
    centroid <- centroids[i, ]
    distances[clusters == i, i] <- sqrt(rowSums((cluster_points - centroid)^2))
  }
  
  # 计算每个聚类内的平均距离
  avg_distances <- apply(distances, 2, mean)
  
  # 计算质心之间的距离
  centroid_distances <- as.matrix(dist(centroids))
  
  # 计算 Davies-Bouldin 指数
  DB_index <- 0
  for (i in 1:max(clusters)) {
    max_ratio <- 0
    for (j in 1:max(clusters)) {
      if (i != j) {
        ratio <- (avg_distances[i] + avg_distances[j]) / centroid_distances[i, j]
        if (ratio > max_ratio) {
          max_ratio <- ratio
        }
      }
    }
    DB_index <- DB_index + max_ratio
  }
  DB_index <- DB_index / max(clusters)
  
  return(DB_index)
}
#source("~/project/scRNA_downanalysis/src/src.R")
source(args$src)

#> list.files('/home/data/lpb1/project/28_WD/03.RNAedit')
# [1] "J23020119" "J23020120" "J23020121" "J23020122" "J23020123" "J23020126"
# [7] "J23030438" "J23080320" "J23080321" "J23080483"
#name <-  list.files('/home/data/lpb1/project/28_WD/03.RNAedit')
#name_dir <-  list.files('/home/data/lpb1/project/28_WD/03.RNAedit')
out <- paste(sep = "/", args$out, "QC")
out_result <- args$out
if (file.exists(out)) {
  unlink(out, recursive = TRUE)
}
name <- list.files(args$out)
name_dir <-  list.files(args$out)

#out <- '/home/data/lpb1/project/28_WD/03.RNAedit/QC'
#out_result <- "/home/data/lpb1/project/28_WD/03.RNAedit"
# 定义质控文件路径
#file_path <- "/home/data/lpb1/project/28_WD/03.RNAedit/QC/QC_data-test.csv"
file_path <- paste(sep = "/", args$out, "QC/QC_data-test.csv")

for (name_index in 1:length(name)) {
  if (length(name) != length(name_dir)) {
    break
  }
  out_sub <- paste(sep = "/", out, name[name_index])
  dir.create(out_sub, recursive = T)
  print(paste(sep = "/", args$out, name_dir[name_index], "result"))
  setwd(paste(sep = "/", args$out, name_dir[name_index], "result"))
  
  dir_list <- list.files(".")
  
  RNA_edit_Region_list <- list()
  RNA_edit_function_list <- list()
  for (dir_name in dir_list) {
    message(dir_name)
    print(dir_name)
    
    RNA_edit <- fread(paste(sep = "/", dir_name, 'REDItools_serial_table.txt'), sep = "\t", header = T, check.names = F)
    
    # Step 1: 拆分AllSubs列并创建新的行
    df_split <- RNA_edit %>%
      separate_rows(AllSubs, sep = " ")
    
    # Step 2: 将BaseCount[A,C,G,T]列拆分成A、C、G、T四列
    df_split <- df_split %>%
      mutate(`BaseCount[A,C,G,T]` = str_remove_all(`BaseCount[A,C,G,T]`, "[\\[\\]]")) %>%
      separate(`BaseCount[A,C,G,T]`, into = c("A", "C", "G", "T"), sep = ",", convert = TRUE)
    
    # Step 3: 计算A、C、G、T四列的和
    df_split$BaseDepth <- rowSums(subset(df_split, select = c('A', 'C', 'G', 'T')))
    
    # 30X
    df_split <- subset(df_split, BaseDepth >= 30)
    
    # Step 4: 计算AllSubs后一个字符对应列的数量并新增一列保存
    df_split <- df_split %>%
      mutate(AllSubs_LastChar = str_sub(AllSubs, -1),
             Alter_Reads = case_when(
               AllSubs_LastChar == "A" ~ A,
               AllSubs_LastChar == "C" ~ C,
               AllSubs_LastChar == "G" ~ G,
               AllSubs_LastChar == "T" ~ T,
               TRUE ~ 0
             ))
    
    
    # Alter Reads>
    df_split <- subset(df_split, Alter_Reads > 3)
    
    # type
    df_split$type <- paste(sep = "->", df_split$Reference, df_split$AllSubs_LastChar)
    # Group
    df_split$Group <- dir_name
    
    # 行整理
    df_split <- as.data.frame(df_split[, c('Region', 'Position', 'Reference', 'AllSubs_LastChar', "BaseDepth", "Alter_Reads", "type", "Group")])
    colnames(df_split) <- c('Chr', 'Position', 'Reference', 'Alter', "Depth", "Alter_Reads", "type", "Group")
    
    
    # 准备Annovar注释文件
    
    # 在Position列左边添加End列，改名Position为Start
    df_split_annovar <- df_split %>%
      mutate(End = Position) %>%
      rename(Start = Position) %>%
      select(Chr, End, Start, everything())
    fwrite(df_split_annovar, file = paste(sep = "/", dir_name, 'REDItools_annovar_input.txt'), sep = "\t", col.names = F, row.names = F, quote = F)
    
    # 定义命令
    command <-
      paste0("perl ", args$annovar, " -out ./", dir_name, "/Annovar -build hg38 ./", dir_name,"/REDItools_annovar_input.txt ", args$ref)
    
    # 使用system函数执行命令并等待其完成
    system(command, wait = TRUE)
    
    annovar_region <-
      fread(
        paste(sep = "/", dir_name, 'Annovar.variant_function'),
        sep = "\t",
        header = F,
        check.names = F
      )
    colnames(annovar_region) <-
      c(
        'Region_Ann',
        'RelatedtGene',
        'Chr',
        'State',
        'End',
        'Reference',
        'Alter',
        "Depth",
        "Alter_Reads",
        "type",
        "Group"
      )
    # annovar_function <- fread(
    #   paste(sep = "/", dir_name, 'Annovar.exonic_variant_function'),
    #   sep = "\t",
    #   header = F,
    #   check.names = F
    # )
    # colnames(annovar_function) <-
    #   c(
    #     "Line",
    #     'Function',
    #     'Gene_AlterType',
    #     'Chr',
    #     'State',
    #     'End',
    #     'Reference',
    #     'Alter',
    #     "Depth",
    #     "Alter_Reads",
    #     "type",
    #     "Group"
    #   )
    RNA_edit_Region_list[[dir_name]] <- annovar_region
    # RNA_edit_function_list[[dir_name]] <- annovar_function
  }
  
  RNA_edit_Region_merge <- do.call(rbind, RNA_edit_Region_list)
  # RNA_edit_function_merge <- do.call(rbind, RNA_edit_function_list)
  write.table(RNA_edit_Region_merge, file = paste(sep = "/", out_sub, "S01-RNA编辑位点.csv"), quote = F, sep = "\t", col.names = T,row.names = T)
  # write.table(RNA_edit_function_merge, file = paste(sep = "/", out_sub, "S02-注释RNA编辑位点.csv"), quote = F, sep = "\t", col.names = T,row.names = T)
  
  # RNA编辑位点检出数量
  P <- ggplot(RNA_edit_Region_merge, aes(x = Group)) +
    geom_bar( width = 0.7, color = "steelblue", fill = "steelblue") +
    labs(title = "RNA Edit Number by Spot", x = "Spot", y = "RNA Edit Num") +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
      axis.title.x = element_text(size = 14, face = "bold"),
      axis.title.y = element_text(size = 14, face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    ) 
  save_picture(P, out = out_sub, Pname = 'P01-编辑位点检出数量柱状图', width_len = 600, height_len = 400)
  
  # RNA编辑位点注释检出数量
  # P <- ggplot(RNA_edit_function_merge, aes(x = Group)) +
  #   geom_bar( width = 0.7, color = "steelblue", fill = "steelblue") +
  #   labs(title = "RNA Edit annotation Number by Spot", x = "Spot", y = "RNA Edit Num") +
  #   theme_minimal() +
  #   theme(
  #     plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
  #     axis.title.x = element_text(size = 14, face = "bold"),
  #     axis.title.y = element_text(size = 14, face = "bold"),
  #     axis.text.x = element_text(angle = 45, hjust = 1),
  #     legend.position = "none"
  #   ) 
  # save_picture(P, out = out_sub, Pname = 'P02-注释编辑位点检出数量柱状图', width_len = 600, height_len = 400)
}



QC_file <- list.files(out)
for (QC_index in QC_file) {
  message(QC_index)
  out_sub <- paste(sep = "/", out, QC_index)
  if (QC_index == QC_file[1]) {
    RNA_edit_Region_merge <- read.table(paste(sep = "/", out_sub, "S01-RNA编辑位点.csv"))
    RNA_edit_Region_merge$source <- QC_index
  }else{
    temp <- read.table(paste(sep = "/", out_sub, "S01-RNA编辑位点.csv"), sep = "\t", header = T, row.names = 1)
    temp$source <- QC_index
    RNA_edit_Region_merge <- rbind(RNA_edit_Region_merge, temp)
  }
}
RNA_edit_Region_merge$Group_raw <- RNA_edit_Region_merge$Group
RNA_edit_Region_merge$Group <- unlist(lapply(RNA_edit_Region_merge$Group_raw, function(x){return(str_split(x, "_spot")[[1]][1])}))
RNA_edit_Region_merge$barcode <- paste(sep = "_", RNA_edit_Region_merge$source, RNA_edit_Region_merge$Group_raw)
RNA_edit_Region_merge <- as.data.frame(RNA_edit_Region_merge)


RNA_edit_Region_merge$Cor_gene <- unlist(lapply(RNA_edit_Region_merge$RelatedtGene, function(x) {return(str_split(str_split(x, "[(]")[[1]][1], ',')[[1]][1])}))

#Group_info <- read.table("/home/data/lpb1/project/28_WD/sample.group", sep = "\t", header = T)
Group_info <- read.table(args$group, sep = "\t", header = T)
library(plyr)
RNA_edit_Region_merge$Group2 <- mapvalues(RNA_edit_Region_merge$source,
                                          from = Group_info$sample_id, 
                                          to = Group_info$Group)
fwrite(RNA_edit_Region_merge, paste(sep = "/", out_result, "RNA_edit_raw.table"), sep = "\t", col.names = T, row.names = F, quote = F)
# 构建矩阵，Group为行坐标标签，合并的五列为列坐标标签
result_matrix <- RNA_edit_Region_merge %>%
  unite("Combined", Cor_gene, Chr, State, Reference, Alter, sep = "-") %>%
  dplyr::select(barcode, Combined, Alter_Reads) %>%
  distinct() %>%
  pivot_wider(names_from = Combined, values_from = Alter_Reads, values_fn = list(Combined = length), values_fill = list(Combined = 0))
result_matrix[is.na(result_matrix)] <- 0
matrix <- as.data.frame(t(result_matrix %>% column_to_rownames('barcode')))
fwrite(matrix, paste(sep = "/", out_result, "RNA_edit.matrix"), sep = "\t", col.names = T, row.names = T, quote = F)

# title: "SingleCell_function"
# author: "三口先生"
# date: "2024/12/23"
# version: "3.0.0"
library(ggplot2)
color1<-unique(c("#9FCCE6","#DD5856","#F38D28","#F69B9A","#B5982F","#FEBF80","#DD5856","#F38D28","#F69B9A","#B5982F","#7D6D6A","#8DCA7D","#9FCCE6","#4E7CA6","#D57096","#AD7AA5","#FDBFD2","#DBA3C8","#CC5C15","#D42F7E",'#E5C06A','#E28D8F',"#80B880","#EDB1C4","#5378A1","#847C7D","#4E928D","#E7863C","#B9AAAD","#D45760","#F1AD7C","#66639E","#92B646","#157BB7","#EDBD41","#229873","#CC5C15","#D42F7E","#DBA513","#638CC9","#6F6AA7","#A1C7DC","#359837","#ED7A1C","#8A298C","#286EA1","#DEEDD0","#F49897",'#61B5BA','#D96951','#40507D','#9D5B34','#5D3C8A','#997635','#729D44','#666666','#C2362B','#EABB7B','#C1AFCE'))
col_temp<-factor(color1,levels = color1)
ggplot(data = data.frame(col=col_temp,y=rep(1,length(col_temp))),
       aes(x=col_temp,y = y))+
  geom_bar(stat='identity',fill=col_temp)+
  scale_fill_manual(values = col_temp)+
  theme(axis.text.x = element_text(angle = 45,hjust = 1,size = 5,vjust = 1))

suppressMessages(library(jsonlite))
# 获取当前项目json配置文件(不能动)
config <- fromJSON("/mnt/sda/project/yjc1/pipeline/02_scRNA_edit/script/config.json")

# library
using <- function(...) {
  packages <- as.character(match.call(expand.dots = FALSE)[[2]])
  
  if (length(packages) == 0) {
    return(invisible())
  }
  # Attempt to load packages making note of which don't load
  loaded <- sapply(packages, function(x) {
    # Try to load package
    if (suppressPackageStartupMessages(require(x, character.only = TRUE, quietly = TRUE))) {
      return(TRUE)
    }
    # Couldn't load
    return(FALSE)
  })
  
  # Give a warning if some packags couldn't be loaded
  if (!all(loaded)) {
    failed <- packages[!loaded]
    warning("\n Failed to load: ", paste(failed, collapse = ", "))
  }
  return(invisible(loaded))
}

# message('save_picture')
save_picture<-function(picture, out, Pname, width_len = 1100, height_len = 900){
  library(ggplot2)
  while (!is.null(dev.list()))  dev.off()
  png(filename = paste(sep = "/", out, paste0(Pname,".png")),width = width_len, height = height_len)
  print(picture)
  while (!is.null(dev.list()))  dev.off()
  while (!is.null(dev.list()))  dev.off()
  pdf(file = paste(sep = "/", out, paste0(Pname,".pdf")),width = as.numeric(width_len/100), height = as.numeric(height_len/100))
  print(picture)
  while (!is.null(dev.list()))  dev.off()
}

# message('ann_and_plot')
ann_and_plot<-function(mouse_deal, ann_file, meta_names){
  set.seed(2023)
  library(Seurat)
  library(future)
  library(ggplot2)
  library(tidyverse)
  if (!(file.exists(ann_file))) {
    stop("Error: ann file isn't existsed")
  }
  #读取注释文件
  ann<-read.table(file = ann_file,
                  header = F,check.names = F)
  ##注释
  all.equal(as.character(levels(mouse_deal)),as.character(ann$V1))
  new_id<-ann$V2
  names(new_id)<-levels(mouse_deal)
  new_id
  mouse_deal <- RenameIdents(mouse_deal,new_id)
  all.equal(rownames(mouse_deal@meta.data),names(mouse_deal@active.ident))
  mouse_deal@meta.data[,meta_names]<-as.factor(mouse_deal@active.ident)
  
  
  return(mouse_deal)
}

step2_after_ann <- function(Rdata_file, out_edit){
  use_threads <- as.integer(as.numeric(system("grep -c ^processor /proc/cpuinfo", intern = TRUE))*2*2/3)
  Rscript <- paste0('load("', Rdata_file, '")
#==========================
seurat_ann <- ann_and_plot(
  seurat_deal_filter_harmony_cluster,
  ann_file = paste0(out$files, "/ann.txt"),
  meta_names = "cell_type"
)
saveRDS(seurat_ann, paste(sep = "/", out$files, "seurat_ann.rds"))
# seurat_ann <- readRDS(paste(sep = "/", out$files, "seurat_ann.rds"))

#==========================
sampleRDS_dir <- paste(sep = "/", out$files, "02-sampleRDS")
dir.create(sampleRDS_dir, recursive = T)
script_dir <- paste(sep = "/", out$files, "03-RNAedit_sh")
dir.create(script_dir, recursive = T)
log_dir <- paste(sep = "/", out$files, "04-RNAedit_log")
dir.create(log_dir, recursive = T)
for (seurat_index in unique(seurat_ann$orig.ident)) {
  message(seurat_index)
  temp_seurat <- subset(seurat_ann, orig.ident == seurat_index)
  # 计算每个群体的数量
  group_counts <- table(temp_seurat@active.ident)
  # 找到数量大于等于20的群体
  valid_groups <- names(group_counts[group_counts >= 20])
  # 更新 active.ident，只保留有效的群体
  temp_seurat <- subset(temp_seurat, idents = valid_groups)
  rownames(temp_seurat@meta.data) <- unlist(lapply(rownames(temp_seurat@meta.data), function(x){return(str_split(x, "_")[[1]][2])}))
  saveRDS(temp_seurat, paste(sep = "/", out$files, "02-sampleRDS", paste0(seurat_index, ".rds")))
  
  
  script <- paste0("#!/bin/bash
  
  # 记录开始时间
  start_time=$(date +%s)
  
  # 运行 R 脚本并记录资源使用情况
  /usr/bin/time -v nohup /mnt/sda/project/yjc1/miniconda3/envs/scRE/bin/Rscript \
  #/mnt/sda/project/yjc1/pipeline/02_scRNA_edit/script/sc_script-single_Cell-group.R change by yijiacheng
  /mnt/sda/project/yjc1/pipeline/02_scRNA_edit/script/sc_script-single_Cell-group.new.R \
  -i ", paste(sep = "/", out$files, "02-sampleRDS", paste0(seurat_index, ".rds")), " \
  -o ',out_edit,'/", seurat_index, " \
  -t ',use_threads,' \
  -b ',out$work,'/01_Cellranger/", seurat_index, "/outs/possorted_genome_bam.bam \
  -g 20 \
  > output.log 2> error.log &
    
    # 等待脚本运行结束
    wait $!
    
    # 记录结束时间
    end_time=$(date +%s)
    
    # 计算总时长
    total_time=$((end_time - start_time))
    echo \'总时长: $total_time 秒\' >> output.log")
  # Split the script into lines
  lines <- unlist(strsplit(script, "\n"))
  # Combine the lines with a space, except for the first line
  modified_script <- paste(paste(lines[c(1:6)], collapse = "\n"), paste(lines[c(7:14)], collapse = " "), paste(lines[c(15:24)], collapse = "\n"), sep = "\n")
  
  file_path <- paste0(script_dir, "/", seurat_index, "_script.sh")
  writeLines(modified_script, file_path)
  # 追加写入第二个脚本内容
  
  # 构造命令
  log_file <- paste0(log_dir, "/", seurat_index, ".log")
  cmd <- paste("nohup bash", file_path, ">", log_file, "2>&1")
  
  # 执行命令
  system(cmd, wait = TRUE)
  
}')
file_path <- paste0(dirname(out$dir), "/step2_afterAnn.R")
writeLines(Rscript, file_path)
}
step2_after_ann_gatk <- function(Rdata_file, out_edit){
  use_threads <- as.integer(as.numeric(system("grep -c ^processor /proc/cpuinfo", intern = TRUE))*2*2/3)
  Rscript <- paste0('load("', Rdata_file, '")
#==========================
seurat_ann <- ann_and_plot(
  seurat_deal_filter_harmony_cluster,
  ann_file = paste0(out$files, "/ann.txt"),
  meta_names = "cell_type"
)
saveRDS(seurat_ann, paste(sep = "/", out$files, "seurat_ann.rds"))
# seurat_ann <- readRDS(paste(sep = "/", out$files, "seurat_ann.rds"))

#==========================
sampleRDS_dir <- paste(sep = "/", out$files, "02-sampleRDS")
dir.create(sampleRDS_dir, recursive = T)
script_dir <- paste(sep = "/", out$files, "03-RNAedit_sh")
dir.create(script_dir, recursive = T)
log_dir <- paste(sep = "/", out$files, "04-RNAedit_log")
dir.create(log_dir, recursive = T)
for (seurat_index in unique(seurat_ann$orig.ident)) {
  message(seurat_index)
  temp_seurat <- subset(seurat_ann, orig.ident == seurat_index)
  # 计算每个群体的数量
  group_counts <- table(temp_seurat@active.ident)
  # 找到数量大于等于20的群体
  valid_groups <- names(group_counts[group_counts >= 20])
  # 更新 active.ident，只保留有效的群体
  temp_seurat <- subset(temp_seurat, idents = valid_groups)
  rownames(temp_seurat@meta.data) <- unlist(lapply(rownames(temp_seurat@meta.data), function(x){return(str_split(x, "_")[[1]][2])}))
  saveRDS(temp_seurat, paste(sep = "/", out$files, "02-sampleRDS", paste0(seurat_index, ".rds")))
  
  
  script <- paste0("#!/bin/bash
  
  # 记录开始时间
  start_time=$(date +%s)
  
  # 运行 R 脚本并记录资源使用情况
  /usr/bin/time -v nohup /mnt/sda/project/yjc1/miniconda3/envs/scRE/bin/Rscript \
  /mnt/sda/project/yjc1/pipeline/02_scRNA_edit/script/sc_script-single_Cell-group-gatk.R \
  -i ", paste(sep = "/", out$files, "02-sampleRDS", paste0(seurat_index, ".rds")), " \
  -o ',out_edit,'/", seurat_index, " \
  -t ',use_threads,' \
  -b ',out$work,'/01_Cellranger/", seurat_index, "/outs/possorted_genome_bam.bam \
  -g 20 \
  > output.log 2> error.log &
    
    # 等待脚本运行结束
    wait $!
    
    # 记录结束时间
    end_time=$(date +%s)
    
    # 计算总时长
    total_time=$((end_time - start_time))
    echo \'总时长: $total_time 秒\' >> output.log")
  # Split the script into lines
  lines <- unlist(strsplit(script, "\n"))
  # Combine the lines with a space, except for the first line
  modified_script <- paste(paste(lines[c(1:6)], collapse = "\n"), paste(lines[c(7:14)], collapse = " "), paste(lines[c(15:24)], collapse = "\n"), sep = "\n")
  
  file_path <- paste0(script_dir, "/", seurat_index, "_script.sh")
  writeLines(modified_script, file_path)
  # 追加写入第二个脚本内容
  
  # 构造命令
  log_file <- paste0(log_dir, "/", seurat_index, ".log")
  cmd <- paste("nohup bash", file_path, ">", log_file, "2>&1")
  
  # 执行命令
  system(cmd, wait = TRUE)
  
}')
file_path <- paste0(dirname(out$dir), "/step2_afterAnn.R")
writeLines(Rscript, file_path)
}

step2_after_ann_seeksoul <- function(Rdata_file, out_edit){
  use_threads <- as.integer(as.numeric(system("grep -c ^processor /proc/cpuinfo", intern = TRUE))*2*2/3)
  Rscript <- paste0('load("', Rdata_file, '")
#==========================
seurat_ann <- ann_and_plot(
  seurat_deal_filter_harmony_cluster,
  ann_file = paste0(out$files, "/ann.txt"),
  meta_names = "cell_type"
)
saveRDS(seurat_ann, paste(sep = "/", out$files, "seurat_ann.rds"))
# seurat_ann <- readRDS(paste(sep = "/", out$files, "seurat_ann.rds"))

#==========================
sampleRDS_dir <- paste(sep = "/", out$files, "02-sampleRDS")
dir.create(sampleRDS_dir, recursive = T)
script_dir <- paste(sep = "/", out$files, "03-RNAedit_sh")
dir.create(script_dir, recursive = T)
log_dir <- paste(sep = "/", out$files, "04-RNAedit_log")
dir.create(log_dir, recursive = T)
for (seurat_index in unique(seurat_ann$orig.ident)) {
  message(seurat_index)
  temp_seurat <- subset(seurat_ann, orig.ident == seurat_index)
  # 计算每个群体的数量
  group_counts <- table(temp_seurat@active.ident)
  # 找到数量大于等于20的群体
  valid_groups <- names(group_counts[group_counts >= 20])
  # 更新 active.ident，只保留有效的群体
  temp_seurat <- subset(temp_seurat, idents = valid_groups)
  rownames(temp_seurat@meta.data) <- unlist(lapply(rownames(temp_seurat@meta.data), function(x){return(str_split(x, "_")[[1]][2])}))
  saveRDS(temp_seurat, paste(sep = "/", out$files, "02-sampleRDS", paste0(seurat_index, ".rds")))
  
  
  script <- paste0("#!/bin/bash
source ~/.bashrc
conda activate scRNA-edit

  
  # 记录开始时间
  start_time=$(date +%s)
  
  python /home/data/lpb1/project/seeksoul_bam_split_dup/add_CB_UR.py \
  -i ',dirname(out$dir),'/01.data/02.out/", seurat_index, "/step2/STAR/", seurat_index, "_SortedByCoordinate.bam \
  -o ',dirname(out$dir),'/01.data/02.out/", seurat_index, "/step2/STAR/", seurat_index, "_SortedByCoordinate_addCB.bam \
  -t ',use_threads,'

  # 运行 R 脚本并记录资源使用情况
  /usr/bin/time -v nohup /mnt/sda/project/yjc1/miniconda3/envs/scRE/bin/Rscript \
  /mnt/sda/project/yjc1/pipeline/02_scRNA_edit/script/sc_script-single_Cell-group-seeksoul.R \
  -i ", paste(sep = "/", out$files, "02-sampleRDS", paste0(seurat_index, ".rds")), " \
  -o ',out_edit,'/", seurat_index, " \
  -t ',use_threads,' \
  -b ',out$work,'/01_Cellranger/", seurat_index, "/step2/STAR/", seurat_index, "_SortedByCoordinate_addCB.bam \
  -g 20 \
  > output.log 2> error.log &
    
    # 等待脚本运行结束
    wait $!
    
    # 记录结束时间
    end_time=$(date +%s)
    
    # 计算总时长
    total_time=$((end_time - start_time))
    echo \'总时长: $total_time 秒\' >> output.log")
  # Split the script into lines
  lines <- unlist(strsplit(script, "\n"))
  # Combine the lines with a space, except for the first line
  modified_script <- paste(paste(lines[c(1:8)], collapse = "\n"), paste(lines[c(9:14)], collapse = " "), paste(lines[c(15:22)], collapse = " "), paste(lines[c(23:32)], collapse = "\n"), sep = "\n")
  file_path <- paste0(script_dir, "/", seurat_index, "_script.sh")
  writeLines(modified_script, file_path)
  # 追加写入第二个脚本内容
  
  # 构造命令
  log_file <- paste0(log_dir, "/", seurat_index, ".log")
  cmd <- paste("nohup bash", file_path, ">", log_file, "2>&1")
  
  # 执行命令
  system(cmd, wait = TRUE)
  
}')
file_path <- paste0(dirname(out$dir), "/step2_afterAnn.R")
writeLines(Rscript, file_path)
}

doubleCell_filter <- function(scRNA_harmony, fig_out, tab_out){
  paramSweep_fix <- function (seu, PCs = 1:10, sct = FALSE, num.cores = 1) {
    require(Seurat)
    require(fields)
    require(parallel)
    pK <- c(5e-04, 0.001, 0.005, seq(0.01, 0.3, by = 0.01))
    pN <- seq(0.05, 0.3, by = 0.05)
    min.cells <- round(nrow(seu@meta.data)/(1 - 0.05) - nrow(seu@meta.data))
    pK.test <- round(pK * min.cells)
    pK <- pK[which(pK.test >= 1)]
    orig.commands <- seu@commands
    if (nrow(seu@meta.data) > 10000) {
      real.cells <- rownames(seu@meta.data)[sample(1:nrow(seu@meta.data), 
                                                   10000, replace = FALSE)]
      data <- seu@assays[["RNA"]]@counts[, real.cells]
      n.real.cells <- ncol(data)
    }
    if (nrow(seu@meta.data) <= 10000) {
      real.cells <- rownames(seu@meta.data)
      data <- seu@assays[["RNA"]]@counts
      n.real.cells <- ncol(data)
    }
    if (num.cores > 1) {
      require(parallel)
      cl <- makeCluster(num.cores)
      output2 <- mclapply(as.list(1:length(pN)), FUN = parallel_paramSweep, 
                          n.real.cells, real.cells, pK, pN, data, orig.commands, 
                          PCs, sct, mc.cores = num.cores)
      stopCluster(cl)
    }
    else {
      output2 <- lapply(as.list(1:length(pN)), FUN = parallel_paramSweep, 
                        n.real.cells, real.cells, pK, pN, data, orig.commands, 
                        PCs, sct)
    }
    sweep.res.list <- list()
    list.ind <- 0
    for (i in 1:length(output2)) {
      for (j in 1:length(output2[[i]])) {
        list.ind <- list.ind + 1
        sweep.res.list[[list.ind]] <- output2[[i]][[j]]
      }
    }
    name.vec <- NULL
    for (j in 1:length(pN)) {
      name.vec <- c(name.vec, paste("pN", pN[j], "pK", pK, 
                                    sep = "_"))
    }
    names(sweep.res.list) <- name.vec
    return(sweep.res.list)
  }
  doubletFinder_fix <- function (seu, PCs, pN = 0.25, pK, nExp, reuse.pANN = FALSE, 
                                 sct = FALSE, annotations = NULL) {
    require(Seurat)
    require(fields)
    require(KernSmooth)
    if (reuse.pANN != FALSE) {
      pANN.old <- seu@meta.data[, reuse.pANN]
      classifications <- rep("Singlet", length(pANN.old))
      classifications[order(pANN.old, decreasing = TRUE)[1:nExp]] <- "Doublet"
      seu@meta.data[, paste("DF.classifications", pN, pK, nExp, 
                            sep = "_")] <- classifications
      return(seu)
    }
    if (reuse.pANN == FALSE) {
      real.cells <- rownames(seu@meta.data)
      data <- seu@assays[["RNA"]]@counts[, real.cells]
      n_real.cells <- length(real.cells)
      n_doublets <- round(n_real.cells/(1 - pN) - n_real.cells)
      print(paste("Creating", n_doublets, "artificial doublets...", 
                  sep = " "))
      real.cells1 <- sample(real.cells, n_doublets, replace = TRUE)
      real.cells2 <- sample(real.cells, n_doublets, replace = TRUE)
      doublets <- (data[, real.cells1] + data[, real.cells2])/2
      colnames(doublets) <- paste("X", 1:n_doublets, sep = "")
      data_wdoublets <- cbind(data, doublets)
      if (!is.null(annotations)) {
        stopifnot(typeof(annotations) == "character")
        stopifnot(length(annotations) == length(Cells(seu)))
        stopifnot(!any(is.na(annotations)))
        annotations <- factor(annotations)
        names(annotations) <- Cells(seu)
        doublet_types1 <- annotations[real.cells1]
        doublet_types2 <- annotations[real.cells2]
      }
      orig.commands <- seu@commands
      if (sct == FALSE) {
        print("Creating Seurat object...")
        seu_wdoublets <- CreateSeuratObject(counts = data_wdoublets)
        print("Normalizing Seurat object...")
        seu_wdoublets <- NormalizeData(seu_wdoublets, normalization.method = orig.commands$NormalizeData.RNA@params$normalization.method, 
                                       scale.factor = orig.commands$NormalizeData.RNA@params$scale.factor, 
                                       margin = orig.commands$NormalizeData.RNA@params$margin)
        print("Finding variable genes...")
        seu_wdoublets <- FindVariableFeatures(seu_wdoublets, 
                                              selection.method = orig.commands$FindVariableFeatures.RNA$selection.method, 
                                              loess.span = orig.commands$FindVariableFeatures.RNA$loess.span, 
                                              clip.max = orig.commands$FindVariableFeatures.RNA$clip.max, 
                                              mean.function = orig.commands$FindVariableFeatures.RNA$mean.function, 
                                              dispersion.function = orig.commands$FindVariableFeatures.RNA$dispersion.function, 
                                              num.bin = orig.commands$FindVariableFeatures.RNA$num.bin, 
                                              binning.method = orig.commands$FindVariableFeatures.RNA$binning.method, 
                                              nfeatures = orig.commands$FindVariableFeatures.RNA$nfeatures, 
                                              mean.cutoff = orig.commands$FindVariableFeatures.RNA$mean.cutoff, 
                                              dispersion.cutoff = orig.commands$FindVariableFeatures.RNA$dispersion.cutoff)
        print("Scaling data...")
        seu_wdoublets <- ScaleData(seu_wdoublets, features = orig.commands$ScaleData.RNA$features, 
                                   model.use = orig.commands$ScaleData.RNA$model.use, 
                                   do.scale = orig.commands$ScaleData.RNA$do.scale, 
                                   do.center = orig.commands$ScaleData.RNA$do.center, 
                                   scale.max = orig.commands$ScaleData.RNA$scale.max, 
                                   block.size = orig.commands$ScaleData.RNA$block.size, 
                                   min.cells.to.block = orig.commands$ScaleData.RNA$min.cells.to.block)
        print("Running PCA...")
        seu_wdoublets <- RunPCA(seu_wdoublets, features = orig.commands$ScaleData.RNA$features, 
                                npcs = length(PCs), rev.pca = orig.commands$RunPCA.RNA$rev.pca, 
                                weight.by.var = orig.commands$RunPCA.RNA$weight.by.var, 
                                verbose = FALSE)
        pca.coord <- seu_wdoublets@reductions$pca@cell.embeddings[, 
                                                                  PCs]
        cell.names <- rownames(seu_wdoublets@meta.data)
        nCells <- length(cell.names)
        rm(seu_wdoublets)
        gc()
      }
      if (sct == TRUE) {
        require(sctransform)
        print("Creating Seurat object...")
        seu_wdoublets <- CreateSeuratObject(counts = data_wdoublets)
        print("Running SCTransform...")
        seu_wdoublets <- SCTransform(seu_wdoublets)
        print("Running PCA...")
        seu_wdoublets <- RunPCA(seu_wdoublets, npcs = length(PCs))
        pca.coord <- seu_wdoublets@reductions$pca@cell.embeddings[, 
                                                                  PCs]
        cell.names <- rownames(seu_wdoublets@meta.data)
        nCells <- length(cell.names)
        rm(seu_wdoublets)
        gc()
      }
      print("Calculating PC distance matrix...")
      dist.mat <- fields::rdist(pca.coord)
      print("Computing pANN...")
      pANN <- as.data.frame(matrix(0L, nrow = n_real.cells, 
                                   ncol = 1))
      if (!is.null(annotations)) {
        neighbor_types <- as.data.frame(matrix(0L, nrow = n_real.cells, 
                                               ncol = length(levels(doublet_types1))))
      }
      rownames(pANN) <- real.cells
      colnames(pANN) <- "pANN"
      k <- round(nCells * pK)
      for (i in 1:n_real.cells) {
        neighbors <- order(dist.mat[, i])
        neighbors <- neighbors[2:(k + 1)]
        pANN$pANN[i] <- length(which(neighbors > n_real.cells))/k
        if (!is.null(annotations)) {
          for (ct in unique(annotations)) {
            neighbors_that_are_doublets = neighbors[neighbors > 
                                                      n_real.cells]
            if (length(neighbors_that_are_doublets) > 0) {
              neighbor_types[i, ] <- table(doublet_types1[neighbors_that_are_doublets - 
                                                            n_real.cells]) + table(doublet_types2[neighbors_that_are_doublets - 
                                                                                                    n_real.cells])
              neighbor_types[i, ] <- neighbor_types[i, 
              ]/sum(neighbor_types[i, ])
            }
            else {
              neighbor_types[i, ] <- NA
            }
          }
        }
      }
      print("Classifying doublets..")
      classifications <- rep("Singlet", n_real.cells)
      classifications[order(pANN$pANN[1:n_real.cells], decreasing = TRUE)[1:nExp]] <- "Doublet"
      seu@meta.data[, paste("pANN", pN, pK, nExp, sep = "_")] <- pANN[rownames(seu@meta.data), 
                                                                      1]
      seu@meta.data[, paste("DF.classifications", pN, pK, nExp, 
                            sep = "_")] <- classifications
      if (!is.null(annotations)) {
        colnames(neighbor_types) = levels(doublet_types1)
        for (ct in levels(doublet_types1)) {
          seu@meta.data[, paste("DF.doublet.contributors", 
                                pN, pK, nExp, ct, sep = "_")] <- neighbor_types[, 
                                                                                ct]
        }
      }
      return(seu)
    }
  }
  ## 对"scRNA_harmony"这个单细胞对象进行pN-pK参数扫描，以生成人工双细胞并计算每个细胞的pANN值
  sweep.res.list <- paramSweep_fix(scRNA_harmony, PCs = 1:20, sct = T, num.cores = 4)
  ## 对参数扫描的结果进行汇总，计算每种pN和pK组合下的双细胞检测指标。这些指标可以用来评估不同参数下的双细胞检测效果，并选择最优的参数。参数GT表示是否提供了真实的双细胞标签，此处没有提供
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE) 
  ## 对汇总结果按照真实双胞胎比例（BCreal）进行升序排序，并显示排序后的数据框。这可以帮助我们找到双胞胎检测效果最好的参数组合
  sweep.stats[order(sweep.stats$BCreal),]
  ## 根据汇总结果找到最优的pK参数。
  bcmvn <- find.pK(sweep.stats) 
  ## 提取出全局最优的pK值，储存于"pK_bcmvn"
  pK_bcmvn <- as.numeric(bcmvn$pK[which.max(bcmvn$BCmetric)]) 
  ## 估计的同源双细胞（由两个或多个相同类型的细胞融合而成的假阳性细胞，它们通常比异源双细胞更难以检测和去除）的比例     
  homotypic.prop <- modelHomotypic(scRNA_harmony$seurat_clusters) 
  ## 计算总的双细胞数量（假设双细胞形成率为 7.5%）
  nExp_poi <- round(0.075 *nrow(scRNA_harmony@meta.data)) 
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop)) # 计算异源双细胞数量
  ## 使用确定好的参数鉴定doublets
  scRNA_harmony <- doubletFinder_fix(scRNA_harmony, PCs = 1:20, pN = 0.25, pK = pK_bcmvn,
                                     nExp = nExp_poi.adj, reuse.pANN = F, sct = T)
  scRNA_harmony$DF.classifications <- scRNA_harmony@meta.data[,colnames(scRNA_harmony@meta.data)[length(colnames(scRNA_harmony@meta.data))]]
  ## 可视化
  P <- DimPlot(scRNA_harmony, reduction = "umap", group.by = "DF.classifications", raster = FALSE, pt.size = 0.8)
  save_picture(
    P,
    out =  fig_out,
    Pname = "10-Double Cell",
    width_len = 900,
    height_len = 600
  )
  write.csv(scRNA_harmony@meta.data, file = paste(sep = "/", tab_out, "01-Double_Cell_info"))
  ## 将双细胞剔除后生成新的对象scRNA_harmony.singlet，便于我们后续分析
  scRNA_harmony.singlet <- subset(scRNA_harmony, subset = DF.classifications == "Singlet")
  return(scRNA_harmony.singlet)
}
#.libPaths(c("/home/data/lpb/R/x86_64-pc-linux-gnu-library/4.2", .libPaths()))

set.seed(2023)
suppressMessages(library(tidyverse))
suppressMessages(library(Seurat))
suppressMessages(library(data.table))
suppressMessages(library(openxlsx))
suppressMessages(library(dplyr))
suppressMessages(library(parallel))
suppressMessages(library(here))
suppressMessages(library(jsonlite))

# 获取当前项目json配置文件(不能动)
config <- fromJSON("/mnt/sda/project/yjc1/pipeline/02_scRNA_edit/script/config.json")

# 定义主函数
group_barcodes_by_ident <- function(seurat_obj, group_prefix = "spot", group_num) {
  # 获取所有的条形码和它们的分类
  barcodes <- seurat_obj@meta.data %>% rownames_to_column("barcode") %>% select(barcode, cell_type)
  
  # 获取唯一分类
  levels_ident <- levels(seurat_obj@active.ident)
  
  # 初始化一个新的列，用于存储新的分组名称
  barcodes$group <- NA
  
  # 按照 levels 进行分组
  result_barcodes <- do.call(rbind, lapply(levels_ident, function(ident_level) {
    subset_barcodes <- barcodes %>% filter(cell_type == ident_level)
    split_barcodes(subset_barcodes, group_prefix, group_num)
  }))
  
  # 去掉多余的行
  result_barcodes <- result_barcodes %>% filter(!is.na(group))
  
  return(result_barcodes)
}
# 定义分组函数
split_barcodes <- function(barcodes, group_name, group_num) {
  group_num <- as.numeric(group_num)
  num_barcodes <- nrow(barcodes)
  if (num_barcodes >= group_num) {
    # 随机抽样条形码
    sampled_barcodes <- barcodes[sample(nrow(barcodes), num_barcodes), ]
    num_groups <- floor(num_barcodes / group_num)
    sampled_barcodes <- sampled_barcodes[1:(num_groups * group_num), ]
    sampled_barcodes$group <- rep(paste0(group_name, 1:num_groups), each = group_num)
  }
  return(sampled_barcodes)
}
# 定义一个函数来检查文件是否存在
check_file_exists <- function(file_path) {
  if (!file.exists(file_path) && is.null(Sys.readlink(file_path))) {
    # 如果文件不存在，输出信息并退出脚本
    message("文件不存在：", file_path, " \n")
    quit(status = 1)
  } else {
    # 如果文件存在或是有效的软链接，输出信息
    message("文件存在或是有效的软链接：", file_path, " \n")
  }
}

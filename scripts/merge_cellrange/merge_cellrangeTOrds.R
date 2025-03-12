#!/bin/Rscript
## 参数整理
library(optparse)
option_list <- list(
  make_option(c("-s", "--sample_info"), type = "character", help  = "sample info"),
  make_option(c("-o", "--out"), type = "character", help  = "Output dir"),
  make_option(c("-n", "--name"), type = "character", help  = "RDS name")
)
args <-
  parse_args(
    OptionParser(option_list = option_list, usage = "Usage: %prog [options] \nDescription: Merge data into seurat objects!")
  )
# version 1:0:0

set.seed(2023)
suppressMessages(library(tidyverse))
suppressMessages(library(Seurat))
suppressMessages(library(data.table))
suppressMessages(library(openxlsx))


# out<-"src"
# name<-"raw_seurat"
# sample_info <- read.table('./src/sample_deal.group', header = T)
out <- args$out
name <- args$name

sample_info <- read.table(args$sample_info, header = T)

message(paste0('Check whether the file exists'))
file.exists(sample_info$Path)

for (index in 1:nrow(sample_info)) {
  temp_matrix<-Read10X(sample_info$Path[index], gene.column = 2)
  colnames(temp_matrix)<-paste(sample_info$sample_id[index],'_',colnames(temp_matrix),sep = "")
  seurat_temp <- CreateSeuratObject(counts =  temp_matrix, min.cells = 3, min.features = 200)
  if (index == 1) {
    seurat <- seurat_temp
  }else{
    seurat <- merge(seurat, seurat_temp)
  }
}

library(plyr)
seurat@meta.data$Group <- mapvalues(seurat@meta.data$orig.ident,
                                    from = sample_info$sample_id, 
                                    to = sample_info$Group)

saveRDS(seurat,file = paste(out,
                            paste0(name,".rds"),sep = "/"))
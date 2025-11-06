#!/bin/Rscript
library(optparse)
option_list<-list(
  make_option(c("-i","--seurat_rds"),type = "character", help  = "Single Cell Seurat Data (RDS FILE) (\t)"),
  make_option(c("-o","--out"),type = "character", help  = "Output dir"),
  make_option(c("-t","--n_cpu_source"),type = "integer", default = 8, help  = "Use CPU resources (integer-please test, not exactly consistent with threads and cores)"),
  make_option(c("-b","--bam"),type = "character", help  = "BAM FILE"),
  make_option(c("-g","--group_num"),type = "integer", default = 20, help  = "Each Spot makes up the number of cells (integer)")
)
args <- parse_args(OptionParser(option_list=option_list,usage = "Usage: %prog [options] \nDescription: RNA-Edit Calling by Spot!"))
# version 6:0:0


# 获取当前项目src配置文件(不能动)
source("/home/data/lpb1/project/RNA_edit_pipeline/src/src.R", encoding = "UTF-8")


check_file_exists(args$seurat_rds)
check_file_exists(args$bam)
seurat_ann <- readRDS(args$seurat_rds) # 加载activate.idents 注释好的Seurat对象




# 输出目录
dir.create(args$out, recursive = T)
barcode_dir <- paste(sep ="/", normalizePath(args$out) , "barcode")
script_dir <- paste(sep ="/", normalizePath(args$out) , "script")
result_dir <- paste(sep ="/", normalizePath(args$out) , 'result')
log_dir <- paste(sep ="/", normalizePath(args$out) , 'log')
message(paste("\n", barcode_dir, script_dir, result_dir, log_dir))

dir.create(barcode_dir, recursive = T)
dir.create(script_dir, recursive = T)
dir.create(result_dir, recursive = T)
dir.create(log_dir, recursive = T)

# CPU资源使用
nthread <- args$n_cpu_source
# 检查n_cpu_source是否大于线程数
nthread_ref <- as.numeric(detectCores()*2)
if (nthread > nthread_ref) {
  # 如果n_cpu_source大于线程数，输出信息并退出脚本
  message("n_cpu_source (", nthread, ") 大于系统线程数 (", nthread_ref, ")。")
  quit(status = 1)
} else {
  # 如果n_cpu_source不大于线程数，继续执行其他代码
  message("n_cpu_source (", nthread, ") 不大于系统线程数 (", nthread_ref, ")。")
  # 这里可以添加后续的处理代码
}

result_barcodes <- group_barcodes_by_ident(seurat_ann, group_num = args$group_num)
result_barcodes$group2 <- paste(sep = "_", result_barcodes$cell_type, result_barcodes$group)

task_init <- as.numeric(system("pgrep -u $(whoami) | wc -l", intern = TRUE))

for (celltype in unique(result_barcodes$group2)) {
  
  aim_barcode <- result_barcodes$barcode[result_barcodes$group2 == celltype]
  write.table(aim_barcode, file = paste0(barcode_dir, "/", celltype, "_barcode.txt"), col.names = F, row.names = F,quote = F,sep = "\t")
  
  script <- paste0('source ~/.bashrc
conda activate scRNA-edit
', config$`Others File`$Cell2editing, ' \
        ', normalizePath(args$bam), ' \
        ', paste0(barcode_dir, "/", celltype, "_barcode.txt"), ' \
        CB:Z: \
        ', celltype, ' \
        ', result_dir, ' \
        ', config$`Others File`$ref_genome, ' \
        ', config$`Others File`$ref_dbsnp, ' \
        ', config$`Others File`$ref_simpleRepeat, ' \
        ', config$`Others File`$ref_alu, '  ', config$`Others File`$Split_bacorde_v2, '  ', config$`Others File`$duplicate.v2, '  ', config$`Others File`$red_ML)
  script2 <- paste0('conda activate gatk_env
samtools index ', result_dir, '/', celltype, '/', celltype, '.marked_duplicates.bam
bash ', config$`Others File`$gatk_script, ' \
        -i ', result_dir, '/', celltype, '/', celltype, '.marked_duplicates.bam \
        -n ', celltype, ' \
        -o ', result_dir, '/', celltype,' \
        -r ', config$`Others File`$gatk_script_ref_genome)
  # Split the script into lines
  lines <- unlist(strsplit(script, "\n"))
  lines2 <- unlist(strsplit(script2, "\n"))
  
  # Combine the lines with a space, except for the first line
  modified_script <- paste(lines[1], lines[2], paste(lines[-c(1:2)], collapse = " "), sep = "\n")
  modified_script2 <- paste(lines2[1], lines2[2], paste(lines2[-c(1:2)], collapse = " "), sep = "\n")
  
  file_path <- paste0(script_dir, "/", celltype, "_script.sh")
  writeLines(modified_script, file_path)
  # 追加写入第二个脚本内容
  cat(modified_script2, file = file_path, append = TRUE, sep = "\n")
  
  cat(readLines(file_path), sep = "\n")
  # 监控进程数量
  while (TRUE) {
    # 获取当前用户运行的后台进程数量
    cmd_count <- system("pgrep -u $(whoami) | wc -l", intern = TRUE)
    current_count <- as.numeric(cmd_count)
    
    if ((current_count-task_init) < nthread) {
      # 构造命令
      log_file <- paste0(log_dir, "/", celltype, ".log")
      cmd <- paste("nohup bash", file_path, ">", log_file, "2>&1")
      
      # 执行命令
      system(cmd, wait = FALSE)
      break
    } else {
      # 等待1分钟
      Sys.sleep(60)
    }
  }
}


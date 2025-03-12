version 1.0

workflow RNA_Edit_Workflow {
    input {
	String base_path
        String sample_list
        String sample_group
        String ref_genome
        String marker_file
        String up_out
        String status_file
    }

    # Step 0: create status file
    call create_status_file {
        input:
	    status_file=base_path + "/" + status_file
    }

    # Step 1: up analysis
    call up_analysis {
        input:
            sample_list=base_path + "/" + sample_list,
            up_out=base_path + "/" + up_out,
            ref_genome=ref_genome,
            status_file=base_path + "/" + status_file,
	    dependency_task_over=create_status_file.task_over
    }

    # Step 2: Add sample data information
    call add_sample_data_info {
        input:
            sample_group=base_path + "/" + sample_group,
            status_file=base_path + "/" + status_file,
	    base_path=base_path,
            dependency_task_over=up_analysis.task_over
    }

    # Step 3: Merge data into Seurat object
    call merge_seurat {
        input:
            status_file=base_path + "/" + status_file,
	    base_path=base_path,
            dependency_task_over=add_sample_data_info.task_over
    }

    # Step 4: scRNA analysis
    call scrna_analysis {
        input:
            marker_file=marker_file,
            status_file=base_path + "/" + status_file,
	    base_path=base_path,
            dependency_task_over=merge_seurat.task_over
    }

    output {
        File output_log = "output.log"
        File error_log = "error.log"
	File task_over = "task_over"
    }
}

# Task to create the status file in the current directory
task create_status_file {
    input {
	String status_file
    }

    command {
        # 判断 status_file 是否存在
        if [ -f ~{status_file} ]; then
            echo "Warning: ~{status_file} already exists."
        else
            touch ~{status_file}
        fi
        echo "任务通道 create_status_file 完成" > task_over
    }

    output {
	File task_over = "task_over"
    }
}

task up_analysis {
    input {
        String sample_list
        String up_out
        String ref_genome
        String status_file
	File dependency_task_over
    }

    command {
        set -e
        set -o pipefail

        if ! grep -q "up_analysis" ~{status_file}; then
            echo "*****up analysis working directory: $(pwd)" | tee -a output.log
            mkdir -p ~{up_out}
            /home/data/lpb1/project/00_pipeline_scRNA_edit/workflow/scripts/cellranger/script_cellranger_batch.sh \
                -s ~{sample_list} \
                -o ~{up_out} \
                -r ~{ref_genome} \
                -t 10 \
                -p 4 \
                >> output.log 2>> error.log
            # 如果前面的命令成功执行，则更新状态文件
            if [ $? -eq 0 ]; then
                echo "up_analysis" >> ~{status_file}
                echo "任务通道 up_analysis 完成" > task_over
            else
                echo "Task failed, not updating status file" | tee -a error.log
                exit 1
            fi
        fi
    }

    output {
        File output_log = "output.log"
        File error_log = "error.log"
	File task_over = "task_over"
    }
}

task add_sample_data_info {
    input {
        String sample_group
        String status_file
	String base_path
	File dependency_task_over
    }

    command <<<
        set -e
        set -o pipefail

        if ! grep -q "add_sample_data_info" ~{status_file}; then
            echo "*****Add sample data information working directory: $(pwd)" | tee -a output.log
            mkdir -p ~{base_path}/src
            awk -v abs_path="~{base_path}" 'NR==1 {print $0 "\tPath"}
                NR>1 {print $0 "\t" abs_path "/01.data/02.out/" $1 "/outs/filtered_feature_bc_matrix"}' \
                ~{sample_group} > ~{base_path}/src/sample_deal.group
            # 如果前面的命令成功执行，则更新状态文件
            if [ $? -eq 0 ]; then
                echo "add_sample_data_info" >> ~{status_file}
                echo "任务通道 add_sample_data_info 完成" > task_over
            else
                echo "Task failed, not updating status file" | tee -a error.log
                exit 1
            fi
        fi
    >>>
    
    output {
        File output_log = "output.log"
        File? error_log = "error.log"
	File task_over = "task_over"
    }
}

task merge_seurat {
    input {
        String status_file
	String base_path
	File dependency_task_over
    }

    command {
        set -e
        set -o pipefail

        if ! grep -q "merge_seurat" ~{status_file}; then
            echo "*****merge Seurat working directory: $(pwd)" | tee -a output.log
            /usr/local/lib/R/bin/Rscript \
                /home/data/lpb1/project/00_pipeline_scRNA_edit/workflow/scripts/merge_cellrange/merge_cellrangeTOrds.R \
                -s ~{base_path}/src/sample_deal.group \
                -o ~{base_path}/src \
                -n raw_seurat \
                >> output.log 2>> error.log
            # 如果前面的命令成功执行，则更新状态文件
            if [ $? -eq 0 ]; then
                echo "merge_seurat" >> ~{status_file}
                echo "任务通道 merge_seurat 完成" > task_over
            else
                echo "Task failed, not updating status file" | tee -a error.log
                exit 1
            fi
        fi
    }

    output {
        File output_log = "output.log"
        File error_log = "error.log"
	File task_over = "task_over"
    }
}

task scrna_analysis {
    input {
        String marker_file
        String status_file
	String base_path
	File dependency_task_over
    }

    command {
        set -e
        set -o pipefail

        if ! grep -q "scrna_analysis" ~{status_file}; then
            echo "*****scRNA analysis working directory: $(pwd)" | tee -a output.log
            mkdir -p ~{base_path}/02.scdata
            /usr/local/lib/R/bin/Rscript \
                /home/data/lpb1/project/00_pipeline_scRNA_edit/workflow/scripts/scRNA_downanalysis/script_scRNA_downanalsysis_cluster.R \
                -i ~{base_path}/src/raw_seurat.rds \
                -o ~{base_path}/02.scdata \
                -t 1 \
                -s human \
                --memory 10 \
                --num_feature_low_threshold 200 \
                --num_feature_high_threshold 10000 \
                --mt_gene_threshold 30 \
                --marker ~{marker_file} \
                >> output.log 2>> error.log
           # 如果前面的命令成功执行，则更新状态文件
            if [ $? -eq 0 ]; then
                echo "scrna_analysis" >> ~{status_file}
                echo "任务通道 scrna_analysis 完成" > task_over
            else
                echo "Task failed, not updating status file" | tee -a error.log
                exit 1
            fi
        fi
    }

    output {
        File output_log = "output.log"
        File error_log = "error.log"
	File task_over = "task_over"
    }
}

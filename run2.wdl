version 1.0

workflow RNA_Edit_Workflow2 {
    input {
        String base_path
    }
    # step 1: RNA-edit after ann.txt
    call rna_edit {
        input:
            step2_afterAnn = base_path + "/step2_afterAnn.R",
            ann_txt = base_path + "/02.scdata/03_files/ann.txt",
            base_path = base_path
    }
}

# run RNA-edit after ann.txt
task rna_edit {
    input {
        String step2_afterAnn
        String ann_txt
        String base_path
    }

    command {
        set -e
        set -o pipefail

        # 判断 step2_afterAnn 文件是否存在
        if [ ! -f ~{ann_txt} ] || [ ! -f ~{step2_afterAnn} ]; then
            echo "Error: ~{step2_afterAnn} does not exist. Exiting." | tee -a error.log
            echo "Advise: Try running run.wdl and then run2.wdl." | tee -a error.log
            exit 1  # 退出命令并返回错误状态
        else
            echo "Message: ~{step2_afterAnn} exists. Proceeding with the workflow." | tee -a output.log
            # 确保03.RNAedit目录存在
            mkdir -p ~{base_path}/03.RNAedit
            cd ~{base_path}
            /usr/local/lib/R/bin/Rscript \
            ~{step2_afterAnn} \
            >> output.log 2>> error.log
        fi
        }

    output {
        File output_log = "output.log"
        File error_log = "error.log"
    }
}

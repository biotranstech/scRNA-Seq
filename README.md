scRNA-Seq work flow
======

## scRNA-Seq overview

### What is RNA editing
RNA editing is a dynamic post-transcriptional modification with significant implications for gene regulation and disease mechanisms.

### Introduction
This is a workflow for processing single-cell RNA editing analysis, developed based on the Python package DAGflow. The workflow is divided into two main stages:
1. **up_analysis.py**: Complete the whole process from raw data processing to scRNA analysis
2. **after_anno.py**: RNA editing was performed after processing

## Workflow file directory structure

```
├── up_analysis.py             # Main workflow definition file
├── after_anno.py              # Second workflow definition file
├── common.py                  # Some workflow file processing modules
├── config.py                  # Pipeline software configuration file
├── dagflow/                   # Python DAGflow package
├── sample.list                # Sample list file
├── scripts/                   # Script directory
└── human_common_marker.xlsx   # Human marker information
```

## Installation
* REQUIREMENT
   * [Cromwell](https://github.com/broadinstitute/cromwell/releases)
   * [Cromwell](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/what-is-cell-ranger)
   * [R](https://www.r-project.org)
### Create and activate R environment
For scRNA-Seq, the R version need is over 4.2. You can create a YAML file called environment.yml with the following content:
```
name: r_environment
dependencies:
  - r-base=4.2.3
  - bioconda::bioconductor-annotationdbi
  - r-corrplot
  - r-corrplot
  - r-optparse
  - r-homologene
  - r-openxlsx
  - r-Seurat
  - r-tidyverse
  - r-data.table
  - r-plyr
  - r-ggplot
```

## Usage

### First：From raw data to scRNA analysis (up_analysis.py)

```bash
python ./up_analysis.py -i sample.xls --work_dir ./01_work --out_dir ./02_result 
```

### Second：RNA editing processing (after_anno.py)

After the first phase is complete, the second phase of processing is performed：

```bash
python ./after_anno.py --work_dir ./01_work --out_dir ./02_result -s sample.xls
```

### sample.xls
Sample information file. The first column contains sample IDs, the second and third columns contain data paths, and the fourth column contains sample grouping information: 

```json
sample1 /../demo/sample1/R22071372-J22120174-J22120174_combined_R1.fastq.gz    /../demo/sample1/R22071372-J22120174-J22120174_combined_R2.fastq.gz    CC1
sample2 /../demo/sample2/R22071373-J22120175-J22120175_combined_R1.fastq.gz    /../demo/sample2/R22071373-J22120175-J22120175_combined_R2.fastq.gz    CC1
sample5 /../demo/sample5/R22071383-J22120191-J22120191_combined_R1.fastq.gz    /../demo/sample5/R22071383-J22120191-J22120191_combined_R2.fastq.gz    CC2
sample6 /../demo/sample6/R22071384-J22120192-J22120192_combined_R1.fastq.gz    /./demo/sample6/R22071384-J22120192-J22120192_combined_R2.fastq.gz    CC2
```


## Status tracking and restart

When submitting a job using nohup, you can monitor the workflow’s progress in real time by checking the nohup.out file. When submitting the job directly in the command line, you can monitor the workflow’s progress through the messages printed on the screen. If the workflow interruption needs to be restarted, the system will automatically check the file, skip the completed task, and continue to execute from the breakpoint, improving work efficiency.


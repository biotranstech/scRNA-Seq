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
Sample information file: 

```json
{
  "RNA_Edit_Workflow.base_path": "/home/data/lpb1/project/00_pipeline_scRNA_edit",
  "RNA_Edit_Workflow.sample_list": "workflow/sample.list",
  "RNA_Edit_Workflow.sample_group": "workflow/sample.group",
  "RNA_Edit_Workflow.ref_genome": "/path/to/reference_genome",
  "RNA_Edit_Workflow.marker_file": "/path/to/marker_file",
  "RNA_Edit_Workflow.up_out": "01.data",
  "RNA_Edit_Workflow.status_file": "workflow/task_status.log"
}
```

### run_workflow_inputs2.json
Input parameter file for Phase 2 workflow (run2.wdl) :

```json
{
  "RNA_Edit_Workflow2.base_path": "/home/data/lpb1/project/00_pipeline_scRNA_edit"
}
```

## Status tracking and restart

The system records completed tasks in the task_status.log file. If the workflow interruption needs to be restarted, the system will automatically check the file, skip the completed task, and continue to execute from the breakpoint, improving work efficiency.


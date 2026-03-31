scRNA-Seq RNA Editing Workflow
======

## scRNA-Seq overview
This project provides an automated workflow for single-cell RNA editing analysis, built on the Python DAGflow framework.
It enables a complete pipeline from raw FASTQ data to RNA editing results, including:
1. Single-cell RNA-seq processing
2. Expression matrix generation
3. RNA editing detection
4. Statistical analysis
5. Result output and visualization

## Why Use This Workflow
Single-cell RNA editing analysis is typically complex and fragmented:
```
| Challenge            | Traditional Approach                |
| -------------------- | ----------------------------------- |
| Multiple tools       | cellranger + custom scripts + R     |
| Complex workflow     | Multi-step manual execution         |
| Error-prone          | Path / format / parameter issues    |
| Poor reproducibility | Different users → different results |
```
✅ Advantages of This Workflow

This pipeline provides:

1. Fully automated workflow
2. Standardized processing
3. Integrated RNA editing analysis
4. Structured outputs
5. Restartable execution (checkpoint support)

## Workflow Overview
```
FASTQ
 ↓
Cell Ranger (scRNA processing)
 ↓
Expression Matrix
 ↓
RNA Editing Detection
 ↓
Statistical Analysis
 ↓
Results
```

## Project Structure
```
├── up_analysis.py        # Stage 1: scRNA processing
├── after_anno.py        # Stage 2: RNA editing analysis
├── config.py            # Configuration file (MUST EDIT)
├── common.py            # Utility modules
├── dagflow/             # DAGflow framework
├── scripts/             # Analysis scripts
├── sample.list          # Sample input file
└── human_common_marker.xlsx
```

## Installation

### Install Conda (if not installed)
Recommended:
```
Miniconda
Anaconda
```

### Create Python Environment
```
conda create -n scrna python=3.9 -y
conda activate scrna
```

### Install R Environment (Required)
Create environment.yml:
```
name: scrna_r
channels:
  - conda-forge
  - bioconda
dependencies:
  - r-base=4.2.3
  - r-seurat
  - r-tidyverse
  - r-data.table
  - r-openxlsx
  - r-optparse
  - bioconductor-annotationdbi
```
Run:
```
conda env create -f environment.yml
conda activate scrna_r
```

### Install Required Software
🔹 Single-cell processing
Install:
```
CellRanger（https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/what-is-cell-ranger）
```
Download from official 10x Genomics website.

## Configuration (CRITICAL)
Edit:
```
config.py
```
Update the following paths:
1. cellranger path
2. Rscript path
3. reference genome path
Example:
```
cellranger_path = "/your/path/cellranger"
rscript_path = "/your/path/Rscript"
reference = "/your/path/reference"
```
⚠️ If not configured → pipeline will fail

## Input Preparation
sample.xls (Required)
```
sample1   /path/sample1_R1.fastq.gz   /path/sample1_R2.fastq.gz   group1
sample2   /path/sample2_R1.fastq.gz   /path/sample2_R2.fastq.gz   group1
sample3   /path/sample3_R1.fastq.gz   /path/sample3_R2.fastq.gz   group2
```
📌 Format Description
```
| Column | Description   |
| ------ | ------------- |
| 1      | Sample ID     |
| 2      | R1 FASTQ path |
| 3      | R2 FASTQ path |
| 4      | Group         |
```

## Usage
### Step 1: Run scRNA Processing
```
Step 1: Run scRNA Processing
```
This step performs:
1. Cell Ranger processing
2. Expression matrix generation

### Step 2: Run RNA Editing Analysis
```
python after_anno.py \
  --work_dir ./01_work \
  --out_dir ./02_result \
  -s sample.xls
```
This step performs:
1. RNA editing detection
2. Statistical analysis

## Output
```
02_result/
├── 01_Cellranger
├── 02_merge
├── 03_scdata
├── 04_RNAedit
├── Figure
└── table
```

## Monitoring & Resume
### Monitor Progress
```
nohup python up_analysis.py ... &
tail -f nohup.out
```

### Resume After Interruption
Simply rerun the same command:
```
python up_analysis.py ...
```
✔ Completed steps will be skipped
✔ Pipeline resumes automatically

## Common Issues
### Missing Software
Check installation:
```
Cell Ranger
R
```

### Incorrect Paths
Check:
```
config.py
```

### Permission Issues
```
chmod -R 755 project_directory
```

### Missing R Packages
```
conda install r-seurat
```

## Minimal Working Example
Prepare:
1. sample.xls
2. FASTQ files
3. config.py properly set

Run:
```
python up_analysis.py -i sample.xls
python after_anno.py -s sample.xls
```

## Summary
This workflow provides a complete solution for scRNA RNA editing analysis.
```
| Feature         | Description           |
| --------------- | --------------------- |
| Automation      | One-command execution |
| Standardization | Reduced human errors  |
| Reproducibility | DAGflow-controlled    |
| Scalability     | Modular design        |
```

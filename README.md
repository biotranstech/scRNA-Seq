# scRNA-Seq

This is a workflow for processing single-cell RNA editing analysis, written in the WDL language and executed based on the Cromwell workflow engine. The workflow is divided into two main stages:
1. **run.wdl**: Complete the whole process from raw data processing to scRNA analysis
2. **run2.wdl**: RNA editing was performed after processing

## Workflow file directory structure

```
├── run.wdl                    # Main workflow definition file
├── run2.wdl                   # Second workflow definition file
├── run_workflow_inputs.json   # Main workflow input parameter file
├── run_workflow_inputs2.json  # Second workflow input parameter file
├── sample.group               # Sample group file
├── sample.list                # Sample list file
├── scripts/                   # Script directory
│   ├── cellranger/            # Cellranger processing script
│   ├── merge_cellrange/       # Data merge script
│   └── scRNA_downanalysis/    # scRNA downstream analysis script
├── src/                       # Intermediate file directory
├── step2_afterAnn.R           # RNA editing post-processing R script
└── task_status.log            # Task status record file
```

## 工作流程详解

### 第一阶段：从原始数据到scRNA分析 (run.wdl)

该阶段包含五个主要步骤：

1. **创建状态文件**：初始化`task_status.log`文件，用于记录已完成的任务
   ```
   call create_status_file
   ```

2. **上游分析**：使用Cellranger处理原始数据
   ```
   call up_analysis
   ```
   - 输入：样本列表文件`sample.list`
   - 输出：位于`01.data/02.out/`目录下的Cellranger处理结果

3. **样本数据信息添加**：为每个样本添加路径信息
   ```
   call add_sample_data_info
   ```
   - 输入：样本分组文件`sample.group`
   - 输出：`src/sample_deal.group`文件，包含样本ID、分组和路径信息

4. **Seurat对象合并**：将所有样本数据合并到一个Seurat对象中
   ```
   call merge_seurat
   ```
   - 输入：`src/sample_deal.group`文件
   - 输出：`src/raw_seurat.rds`文件，包含所有样本的Seurat对象

5. **scRNA分析**：执行单细胞RNA分析
   ```
   call scrna_analysis
   ```
   - 输入：`src/raw_seurat.rds`文件和标记基因文件
   - 输出：位于`02.scdata/`目录下的分析结果，包括`03_files/ann.txt`文件

### 第二阶段：RNA编辑处理 (run2.wdl)

在第一阶段完成后，执行第二阶段处理：

1. **RNA编辑处理**：运行`step2_afterAnn.R`脚本进行RNA编辑分析
   ```
   call rna_edit
   ```
   - 输入：`step2_afterAnn.R`脚本和`02.scdata/03_files/ann.txt`文件
   - 输出：位于`03.RNAedit/`目录下的RNA编辑分析结果

## 输入文件说明

### sample.list
样本列表文件，每行包含一个样本ID，用于Cellranger处理。

示例内容：
```
sample1
sample2
sample3
```

### sample.group
样本分组文件，包含样本ID和分组信息。

示例内容：
```
Sample  Group
sample1 Control
sample2 Treatment
sample3 Control
```

### run_workflow_inputs.json
主工作流(run.wdl)的输入参数文件，包含以下关键参数：

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
第二阶段工作流(run2.wdl)的输入参数文件：

```json
{
  "RNA_Edit_Workflow2.base_path": "/home/data/lpb1/project/00_pipeline_scRNA_edit"
}
```

## 运行工作流

### 运行第一阶段工作流

```bash
cd /home/data/lpb1/project/00_pipeline_scRNA_edit
java -jar /path/to/cromwell.jar run workflow/run.wdl -i workflow/run_workflow_inputs.json
```

### 运行第二阶段工作流

确保第一阶段已完成后再运行：

```bash
cd /home/data/lpb1/project/00_pipeline_scRNA_edit
java -jar /path/to/cromwell.jar run workflow/run2.wdl -i workflow/run_workflow_inputs2.json
```

## 工作流依赖工具

1. **Cromwell** 工作流引擎
2. **Cellranger** 用于单细胞RNA-seq数据处理
3. **R** 及以下R包:
   - Seurat: 单细胞分析
   - ggplot2: 可视化
   - dplyr: 数据处理
   - 其他RNA编辑分析相关的R包

## 输出目录结构

工作流执行后将生成以下主要输出目录：

1. **01.data/**: Cellranger处理结果
   - 02.out/: 各样本的Cellranger输出
   
2. **02.scdata/**: scRNA分析结果
   - 01_files/: 质控后的文件
   - 02_files/: 聚类和降维结果 
   - 03_files/: 注释结果，包含关键的ann.txt文件
   
3. **03.RNAedit/**: RNA编辑分析结果
   - 各类RNA编辑位点信息
   - 统计图表
   
4. **src/**: 中间文件
   - sample_deal.group: 处理后的样本信息
   - raw_seurat.rds: 合并后的Seurat对象

## 状态追踪与重启

系统通过`task_status.log`文件记录已完成的任务。如果工作流中断需要重启，系统会自动检查该文件并跳过已完成的任务，从中断点继续执行，提高工作效率。

## 注意事项

1. 第一阶段工作流(run.wdl)必须先于第二阶段(run2.wdl)运行
2. 确保有足够的磁盘空间，RNA-seq分析需要大量存储
3. 对于大型数据集，可能需要调整内存参数
4. 所有路径配置应使用绝对路径以避免错误

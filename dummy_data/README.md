# RNA编辑工作流管理系统

## 概述

这是一个用于管理和执行RNA编辑分析工作流的综合系统。该系统包含了一个基于Web的用户界面和完整的后端工作流程，使研究人员能够方便地进行RNA编辑分析。

## 目录结构

```
/00_pipeline_scRNA_edit/
├── app.py                     # Flask应用主文件
├── static/                   # 静态资源目录
├── templates/                # HTML模板目录
├── workflow/                 # 工作流文件目录
    ├── run.wdl               # 主工作流定义文件
    ├── run2.wdl              # 第二工作流定义文件
    ├── run_workflow_inputs.json  # 主工作流输入文件
    ├── run_workflow_inputs2.json # 第二工作流输入文件
    ├── sample.group          # 样本分组文件
    ├── sample.list           # 样本列表文件
    ├── scripts/              # 各类脚本目录
        ├── cellranger/       # Cellranger相关脚本
        ├── merge_cellrange/  # 合并Cellranger数据的脚本
        └── scRNA_downanalysis/ # scRNA下游分析脚本
    ├── src/                  # 源代码和中间文件目录
    ├── step2_afterAnn.R      # RNA编辑后处理R脚本
    └── task_status.log       # 任务状态记录
```

## 系统依赖

本系统依赖以下软件和环境：

1. Python 3.6+
2. Flask 框架
3. R 4.0+
4. Cromwell 工作流引擎
5. Cellranger
6. Seurat R包
7. 其他RNA编辑分析所需的生物信息学工具

## 安装和配置

1. 确保系统中已安装所有依赖项
2. 克隆或复制整个目录到目标位置
3. 配置相关路径：
   - 更新`app.py`中的路径配置（如果需要）
   - 检查`workflow`目录中的WDL文件和JSON输入文件，确保路径正确

## 使用方法

### 启动Web界面

```bash
cd 00_pipeline_scRNA_edit
python app.py
```

访问 http://localhost:5000 即可打开Web界面。

### 运行工作流

#### 通过Web界面运行

1. 打开Web界面
2. 上传或选择样本文件
3. 配置工作流参数
4. 点击"运行工作流"按钮
5. 监控任务进度

#### 通过命令行运行

可以直接使用Cromwell运行WDL工作流：

```bash
java -jar /path/to/cromwell.jar run workflow/run.wdl -i workflow/run_workflow_inputs.json
```

## 工作流程说明

本系统包含两个主要工作流：

### 第一阶段 (run.wdl)

1. **准备阶段**：创建状态文件
2. **上游分析**：运行Cellranger处理原始数据
3. **样本信息处理**：添加样本数据信息
4. **Seurat对象合并**：将数据合并到Seurat对象中
5. **scRNA分析**：进行单细胞RNA分析

### 第二阶段 (run2.wdl)

1. **RNA编辑处理**：运行`step2_afterAnn.R`进行RNA编辑处理

## 输入文件

### sample.list

样本列表文件，每行包含一个样本ID。

### sample.group

样本分组文件，包含样本ID和分组信息。

### run_workflow_inputs.json

工作流输入参数配置，包括：

```json
{
  "RNA_Edit_Workflow.base_path": "工作目录路径",
  "RNA_Edit_Workflow.sample_list": "样本列表文件名",
  "RNA_Edit_Workflow.sample_group": "样本分组文件名",
  "RNA_Edit_Workflow.ref_genome": "参考基因组路径",
  "RNA_Edit_Workflow.marker_file": "标记基因文件路径",
  "RNA_Edit_Workflow.up_out": "上游分析输出目录",
  "RNA_Edit_Workflow.status_file": "状态文件名称"
}
```

## 输出目录和文件

系统将生成以下主要输出：

1. **01.data**：存放原始数据处理结果
2. **02.scdata**：存放scRNA分析结果
3. **03.RNAedit**：存放RNA编辑分析结果
4. **src**：存放中间文件和Seurat对象

## 状态追踪

系统使用`task_status.log`文件记录已完成的任务步骤，以便在重启时能够跳过已完成的步骤，提高效率。

## 注意事项

1. 确保数据目录有足够的磁盘空间，RNA编辑分析可能需要大量存储空间
2. 大型数据集分析时，考虑增加内存配置
3. 第一阶段工作流（run.wdl）必须在第二阶段（run2.wdl）之前完成

## 常见问题解答

**Q: 如何跳过已完成的任务？**
A: 系统会自动检查`task_status.log`文件，已完成的任务会被自动跳过。

**Q: 如何调整分析参数？**
A: 可以通过Web界面调整参数，或直接修改JSON输入文件。

**Q: 如何处理分析错误？**
A: 检查日志文件（output.log和error.log）以获取错误详情，解决问题后重新运行工作流。

## 维护和支持

如需进一步帮助，请联系系统管理员或研发团队。
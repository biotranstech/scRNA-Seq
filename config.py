import os.path
from collections import OrderedDict

ROOT = os.path.dirname(os.path.abspath(__file__))
BIN = os.path.join(ROOT, "bin")
SCRIPTS = os.path.join(ROOT, "script")

CELLRANGER = "/mnt/sda/project/yjc1/pipeline/02_scRNA_edit/script/cellranger"
SAMTOOLS = "/mnt/sda/project/yjc1/miniconda3/bin/samtools"
PYSAM_PYTHON = "/mnt/sda/project/yjc1/miniconda3/envs/htseq/bin/python"
SCRE_PYTHON = "/mnt/sda/project/yjc1/miniconda3/envs/scRE/bin/python"
JAVA_HOME = "/mnt/sda/project/yjc1/pipeline/02_scRNA_edit/script/jdk-17.0.17"
PICARD = "/mnt/sda/project/yjc1/miniconda3/envs/scRE/bin/picard"

RSCRIPT = "/mnt/sda/project/yjc1/miniconda3/envs/scRE/bin/Rscript"
python = "/mnt/sda/project/yjc1/miniconda3/bin/python"
PERL = "/mnt/sda/project/yjc1/miniconda3/bin/perl"

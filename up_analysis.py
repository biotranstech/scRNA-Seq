#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import json
import logging
import argparse

from collections import OrderedDict
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), "../"))
from config import *
from dagflow import DAG, Task, ParallelTask, do_dag
from common import check_path, mkdir, read_tsv, read_files, read_sample_file


LOG = logging.getLogger(__name__)

__author__ = ("Liu pingbo",)
__email__ = "..."
__version__ = "v1.0.0"

def get_soft_link(input,ref,work_dir,out_dir,job_type,thread):

    prefix, read1, read2 = read_sample_file(input)
    raw_data_dir = mkdir(os.path.join(work_dir, "00_data"))

    tasks = ParallelTask(
        id="link_file",
        work_dir=work_dir,
        type=job_type,
        option="-pe smp 4",
        script="""
ln -s {{read1}} {raw_data_dir}/{{prefix}}_S1_L001_R1_001.fastq.gz
ln -s {{read2}} {raw_data_dir}/{{prefix}}_S1_L001_R2_001.fastq.gz
{CELLRANGER} count --id={{prefix}} --transcriptome={ref} --fastqs={raw_data_dir} \\
    --sample={{prefix}} --nosecondary --localcores={thread} \\

""".format(raw_data_dir=raw_data_dir,
            scripts=SCRIPTS,
            CELLRANGER=CELLRANGER,
            thread=thread,
            ref=ref,
            out_dir=out_dir
        ),
        prefix=prefix,
        read1=read1,
        read2=read2
    )

    file1 = work_dir + "/sample_deal.group"
    fw = open(file1, 'w')
    line1 = "sample_id\tGroup\tPath\n"
    fw.write(line1)
    with open(input,'r') as f:
        for line in f:
            line = line.strip()
            ll = line.split("\t")
            out_file = work_dir + "/" + ll[0] + "/outs/filtered_feature_bc_matrix"
            line2 = ll[0]+"\t"+ll[3]+"\t"+out_file+"\n"
            fw.write(line2)
    return tasks, file1

def merge_file(group, work_dir, job_type, out_dir):

    task = Task(
        id="work_merge_gene_count",
        work_dir=work_dir,
        type=job_type,
        option="-pe smp 4",
        script="""
{Rscript} {scripts}/merge_cellrangeTOrds.R -s {group} \\
    -o {out_dir} \\
    -n raw_seurat

""".format(Rscript=RSCRIPT,
            scripts=SCRIPTS,
            group=group,
            work_dir=work_dir,
            out_dir=out_dir
        )
    )
    seurat_file = os.path.join(out_dir, 'raw_seurat.rds')
    return task, seurat_file

def scRNA_downanalsysis_cluster(seurat, thread, memory, species, mt_gene_threshold,
                low_threshold, high_threshold, marker, job_type, work_dir, out_dir):

    work_dir2 = os.path.join(work_dir, '../')
    task = Task(
        id="work_scRNA_downanalsysis_cluster",
        work_dir=work_dir,
        type=job_type,
        option="-pe smp 4",
        script="""
{Rscript} {scripts}/script_scRNA_downanalsysis_cluster.R -i {seurat} \\
    -o {out_dir} -t {thread} -s {species} --memory {memory} \\
    -w {work_dir} \\
    --num_feature_low_threshold {low_threshold} --num_feature_high_threshold {high_threshold} \\
    --mt_gene_threshold {mt_gene_threshold} --marker {marker_file}

""".format(Rscript=RSCRIPT,
            scripts=SCRIPTS,
            seurat=seurat,
            thread=thread,
            species=species,
            memory=memory,
            low_threshold=low_threshold,
            high_threshold=high_threshold,
            mt_gene_threshold=mt_gene_threshold,
            marker_file=marker,
            work_dir=work_dir2,
            out_dir=out_dir
        )
    )
    return task

def run_up_analysis(input, ref, thread, job_type, memory, species,
                concurrent, refresh, work_dir, out_dir, marker,
                mt_gene_threshold, low_threshold, high_threshold):

    work_dir = mkdir(work_dir)
    out_dir = mkdir(out_dir)
    input = os.path.abspath(input)

    dag = DAG("work_cellranger")

    cellranger_tasks, group_file  = get_soft_link(
        input = input,
        ref = ref,
        work_dir = mkdir(os.path.join(work_dir, "01_Cellranger")),
        out_dir = mkdir(os.path.join(out_dir, "01_Cellranger")),
        job_type = job_type,
        thread = thread
    )

    dag.add_task(*cellranger_tasks)

    merge_task, seurat_file = merge_file(
        group = group_file, 
        work_dir = mkdir(os.path.join(work_dir, "02_merge")), 
        out_dir = mkdir(os.path.join(out_dir, "02_merge")),
        job_type = job_type
    )

    dag.add_task(merge_task)
    merge_task.set_upstream(*cellranger_tasks)

    scRNA_cluster_task = scRNA_downanalsysis_cluster(
        seurat = seurat_file, 
        thread = thread, 
        memory = memory, 
        species = species, 
        mt_gene_threshold = mt_gene_threshold,
        low_threshold = low_threshold, 
        high_threshold = high_threshold, 
        marker = marker, 
        job_type = job_type, 
        work_dir = mkdir(os.path.join(work_dir, "03_scdata")), 
        out_dir = mkdir(os.path.join(out_dir, "03_scdata"))
    )

    dag.add_task(scRNA_cluster_task)
    scRNA_cluster_task.set_upstream(merge_task)

    do_dag(dag, concurrent, refresh)

    return None

def up_analysis(args):

    run_up_analysis(
        input=args.input,
        ref=args.reference,
        thread=args.thread,
        memory=args.memory,
        species=args.species,
        mt_gene_threshold=args.mt_gene_threshold,
        low_threshold=args.low_threshold,
        high_threshold=args.high_threshold,
        marker=args.marker,
        job_type=args.job_type,
        concurrent=args.concurrent,
        refresh=args.refresh,
        work_dir=args.work_dir,
        out_dir=args.out_dir
    )

    return 0

def add_up_analysis_args(parser):

    parser.add_argument("-i", "--input", metavar="STR", type=str, required=True,
        help="Input the name of the sample and Second generation sequencing path.")
    parser.add_argument("-ref", "--reference", metavar="STR", type=str, 
        default="/mnt/sda/project/yjc1/pipeline/02_scRNA_edit/ref/refdata-gex-GRCh38-2020-A",
        help="""Input the host's reference database dir.\
        Default: /mnt/sda/project/yjc1/pipeline/02_scRNA_edit/ref/refdata-gex-GRCh38-2020-A; \
        mmu: /mnt/sda/project/yjc1/pipeline/02_scRNA_edit/ref/refdata-gex-mm10-2020-A""")
    parser.add_argument("-s", "--species", metavar="STR", type=str, 
        default="human",
        help="""species of human or mouse.\
        Default: human; """)
    parser.add_argument("--marker", metavar="STR", type=str, 
        default="/mnt/sda/project/yjc1/pipeline/02_scRNA_edit/ref/human_common_marker.xlsx",
        help="""An xlsx file containing marker information \
        (the first line is the cell type name, followed by the corresponding marker)\
        Default: /mnt/sda/project/yjc1/pipeline/02_scRNA_edit/ref/human_common_marker.xlsx; \
        mmu: /mnt/sda/project/yjc1/pipeline/02_scRNA_edit/ref/mouse_common_marker.xlsx""")
    parser.add_argument("--memory", metavar="INT", type=int, default=10,
        help="Memory consumption per thread; default: 10")
    parser.add_argument("--low_threshold", metavar="INT", type=int, default=200,
        help="The number of genes detected is low threshold; default: 200")
    parser.add_argument("--high_threshold", metavar="INT", type=int, default=10000,
        help="High threshold for the number of genes detected; default: 10000")
    parser.add_argument("--mt_gene_threshold", metavar="INT", type=int, default=30,
        help="Mitochondrial percentage threshold (integer; default = 30)")
    parser.add_argument("--thread", metavar="INT", type=int, default=4,
        help="Analysis of the number of threads used, default: 4")
    parser.add_argument("--concurrent", metavar="INT", type=int, default=10,
        help="Maximum number of jobs concurrent  (default: 10)")
    parser.add_argument("--refresh", metavar="INT", type=int, default=30,
        help="Refresh time of log in seconds  (default: 30)")
    parser.add_argument("--job_type", choices=["sge", "local"], default="local",
        help="Jobs run on [sge, local]  (default: local)")
    parser.add_argument("--work_dir", metavar="DIR", default=".",
        help="Work directory (default: current directory)")
    parser.add_argument("--out_dir", metavar="DIR", default=".",
        help="Output directory (default: current directory)")


    return parser


def main():

    logging.basicConfig(
        stream=sys.stderr,
        level=logging.INFO,
        format="[%(levelname)s] %(message)s"
    )

    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description="""

version: %s
contact:  %s <%s>\
    """ % (__version__, " ".join(__author__), __email__))

    parser = add_up_analysis_args(parser)
    args = parser.parse_args()
    up_analysis(args)


if __name__ == "__main__":
    main()

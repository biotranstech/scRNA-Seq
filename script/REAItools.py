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


def work_REDItools(group, bam_file, barcode, ref, work_dir, job_type, out_dir):

    group_list = []
    with open(group,'r') as fr:
        for line in fr:
            line = line.strip()
            if line not in group_list:
                group_list.append(line)
    tasks = ParallelTask(
        id="work_REDItools",
        work_dir=work_dir,
        type=job_type,
        option="-pe smp 4",
        script="""
mkdir -p {out_dir}/{{name}}
perl {scripts}/Split_bacorde_v2.pl -f1 {barcode}/{{name}}_barcode.txt \\
    -f2 {bam_file} \\
    -f3 CB:Z: -out1 {out_dir}/{{name}}/{{name}}.bam
{samtools} view -H {bam_file} > {out_dir}/{{name}}/header.txt
cat {out_dir}/{{name}}/header.txt > {out_dir}/{{name}}/{{name}}.bam
cat {out_dir}/{{name}}/{{name}}.bam.1 >> {out_dir}/{{name}}/{{name}}.bam
rm -rf {out_dir}/{{name}}/{{name}}.bam.1
{samtools} sort {out_dir}/{{name}}/{{name}}.bam -o {out_dir}/{{name}}/{{name}}.sort.bam
{samtools} index {out_dir}/{{name}}/{{name}}.sort.bam

{samtools} idxstats {out_dir}/{{name}}/{{name}}.sort.bam > {out_dir}/{{name}}/{{name}}.sort.stat
sh {scripts}/stat_chr.sh {out_dir}/{{name}}/{{name}}.sort.stat {out_dir}/{{name}}/chromosome.txt

#rm -rf {out_dir}/{{name}}/run.sh {out_dir}/{{name}}/merge.path
{perl} {scripts}/get_dup.pl {out_dir}/{{name}}/chromosome.txt \\
    {out_dir}/{{name}}/run.sh \\
    {out_dir}/{{name}}/merge.path \\
    {PYSAM_PYTHON} \\
    {scripts}/duplicate.v2.py \\
    {out_dir}/{{name}} \\
    {{name}}

sh {out_dir}/{{name}}/run.sh
{samtools} merge -b {out_dir}/{{name}}/merge.path {out_dir}/{{name}}/{{name}}.marked_duplicates.bam
export JAVA_HOME={JAVA_HOME}
export CLASSPATH=.:$JAVA_HOME/lib/
export PATH=.:$JAVA_HOME/bin:$PATH
{PICARD} SortSam \\
    -I {out_dir}/{{name}}/{{name}}.marked_duplicates.bam \\
    -O {out_dir}/{{name}}/{{name}}.marked_duplicates.sort.bam \\
    -SORT_ORDER coordinate

{samtools} index {out_dir}/{{name}}/{{name}}.marked_duplicates.bam
{SCRE_PYTHON} {scripts}/reditools.py -S -q 25 -mbp 6 -Mbp 6 \\
    -f {out_dir}/{{name}}/{{name}}.marked_duplicates.bam \\
    -r {ref} \\
    -o {out_dir}/{{name}}/REDItools_serial_table.txt
""".format(samtools=SAMTOOLS,
            scripts=SCRIPTS,
            perl=PERL,
            PYSAM_PYTHON=PYSAM_PYTHON,
            SCRE_PYTHON=SCRE_PYTHON,
            bam_file=bam_file,
            JAVA_HOME=JAVA_HOME,
            PICARD=PICARD,
            barcode=barcode,
            ref=ref,
            out_dir=out_dir
        ),
        name=group_list
    )
    REDItools_list = []
    for i in group_list:
        REDItools = out_dir+"/"+i+"/"+"REDItools_serial_table.txt"
        REDItools_list.append(REDItools)
    return tasks, REDItools_list

def run_REDItools(job_type, concurrent, refresh,  
    ref, dir, sample, dir2):

    sample_list = []
    with open(sample,'r') as f:
        for line in f:
            l = line.strip().split('\t')
            sample_list.append(l[0])

    for line in sample_list:
        work_dir = mkdir(os.path.join(dir, "04_RNAedit", line, "script"))
        out_dir = mkdir(os.path.join(dir, "04_RNAedit", line, "result"))
        barcode_dir = mkdir(os.path.join(dir, "04_RNAedit", line, "barcode"))
        bam_file = os.path.join(dir2, "01_Cellranger", line, "outs/possorted_genome_bam.bam")
        group = os.path.join(dir, "04_RNAedit", line, "script/group.xls")

        job_name = "work_"+line+"_REAItools"
        dag = DAG(job_name)

        REAItools_tasks , list1 = work_REDItools(
            group = group, 
            bam_file = bam_file, 
            barcode = barcode_dir, 
            ref = ref, 
            work_dir = work_dir, 
            job_type = job_type, 
            out_dir = out_dir
        )
        dag.add_task(*REAItools_tasks)
        do_dag(dag, concurrent, refresh)


    return None

def REAItools(args):

    run_REDItools(
        job_type=args.job_type,
        concurrent=args.concurrent,
        refresh=args.refresh, 
        ref=args.reference, 
        dir=args.out_dir, 
        sample=args.sample, 
        dir2=args.cell_dir)

    return 0

def add_REAItools_args(parser):

    parser.add_argument("-ref", "--reference", metavar="STR", type=str, 
        default="/mnt/sda/project/yjc1/pipeline/02_scRNA_edit/ref/hg38.fa",
        help="""Input the host's reference database dir.\
        Default: /mnt/sda/project/yjc1/pipeline/02_scRNA_edit/ref/hg38.fa """)
    parser.add_argument("-s", "--sample", metavar="STR", type=str, 
        help="""sample file.""")
    parser.add_argument("--concurrent", metavar="INT", type=int, default=10,
        help="Maximum number of jobs concurrent  (default: 10)")
    parser.add_argument("--refresh", metavar="INT", type=int, default=30,
        help="Refresh time of log in seconds  (default: 30)")
    parser.add_argument("--job_type", choices=["sge", "local"], default="local",
        help="Jobs run on [sge, local]  (default: local)")
    parser.add_argument("--cell_dir", metavar="DIR", default=".",
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

    parser = add_REAItools_args(parser)
    args = parser.parse_args()
    REAItools(args)


if __name__ == "__main__":
    main()

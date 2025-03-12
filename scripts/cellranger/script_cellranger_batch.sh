#!/bin/bash

# Usage and help documentation
usage() {
    echo "Usage: $0 -s sample_list -o output_dir -r ref_dir -t threads -p parallel_tasks"
    echo "  -s  Path to the sample information file"
    echo "  -o  Output directory"
    echo "  -r  Reference directory"
    echo "  -t  Number of threads per task"
    echo "  -p  Number of parallel tasks"
    exit 1
}

# Parsing command line arguments
while getopts "s:o:r:t:p:" opt; do
    case ${opt} in
        s ) sample_list=$OPTARG ;;
        o ) output_dir=$OPTARG ;;
        r ) ref_dir=$OPTARG ;;
        t ) threads=$OPTARG ;;
        p ) parallel_tasks=$OPTARG ;;
        * ) usage ;;
    esac
done

# Checking if all required arguments are provided
if [ -z "$sample_list" ] || [ -z "$output_dir" ] || [ -z "$ref_dir" ] || [ -z "$threads" ] || [ -z "$parallel_tasks" ]; then
    usage
fi

# Creating necessary directories
mkdir -p ${output_dir}/01.data
mkdir -p ${output_dir}/03.logs
mkdir -p ${output_dir}/02.out

# Function to process each sample
process_sample() {
    local sample_id=$1
    local fastq_R1=$2
    local fastq_R2=$3

    # Creating soft links
    ln -s $fastq_R1 ${output_dir}/01.data/${sample_id}_S1_L001_R1_001.fastq.gz
    ln -s $fastq_R2 ${output_dir}/01.data/${sample_id}_S1_L001_R2_001.fastq.gz

    # Change to output directory
    cd ${output_dir}/02.out

    # Running cellranger count
    /home/data/lpb1/project/scRNA_cellranger/script/cellranger count \
	             --id=$sample_id \
                     --transcriptome=$ref_dir \
                     --fastqs=${output_dir}/01.data \
                     --sample=$sample_id \
                     --nosecondary \
                     --localcores=$threads \
                     &> ../03.logs/${sample_id}.log
}

# Exporting the function and variables for parallel execution
export -f process_sample
export output_dir ref_dir threads

# Reading the sample information file and processing in parallel
tail -n +2 $sample_list | \
    awk -F'\t' '{print $1, $2, $3}' | \
    xargs -n 3 -P $parallel_tasks -I {} bash -c 'process_sample {}'

# Logging resource usage information
echo "Resource usage information:" >> ${output_dir}/03.logs/resource_usage.log

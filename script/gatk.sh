#!/bin/bash
 
# Usage and help documentation function
usage() {
    echo "Usage: $0 -i sample.marked_duplicates.bam -n Sample_name -o output_dir -r ref_genome"
    echo "  -i  Path to the marked_duplicates.bam"
    echo "  -n  Sample name"
    echo "  -o  Output directory"
    echo "  -r  Path to the reference genome file"
    exit 1
}
 
# Initialize variables
input_bam=""
sample_name=""
output_dir=""
ref_genome=""

# __version__ = "1.0.0"
# __author__ = "刘品博"
# __email__ = "liupinbo1999@163.com"
set -e

# Parsing command line arguments
while getopts "i:n:o:r:" opt; do
    case ${opt} in
        i ) input_bam=$OPTARG ;;
	n ) sample_name=$OPTARG ;;
        o ) output_dir=$OPTARG ;;
        r ) ref_genome=$OPTARG ;;
        * ) usage ;;
    esac
done
 
# Check if all required arguments are provided
if [ -z "$input_bam" ] || [ -z "$output_dir" ] || [ -z "$ref_genome" ] || [ -z "$sample_name" ]; then
    usage
    exit 1  # 假设usage函数中没有包含exit语句
fi

# Ensure the sample list file exists
if [[ ! -f "$input_bam" ]]; then
    echo "Error: input marked_duplicates.bam '$sample_list_file' does not exist."
    exit 1
fi
 
# Function to process a single sample
process_sample() {
    local dup_bam=$1
    local sample_name=$2
    local out=${output_dir}
 
    # Create temporary directory for intermediate files
    local tmp_dir="${out}/tmp"
    mkdir -p "$tmp_dir"
 
    # Mark duplicates (assuming marked_duplicates.bam is already available)
    # This step should ideally be done beforehand and not included here
    # If necessary, uncomment the following line and ensure the input BAM is correct
    # gatk MarkDuplicates -I "$dup_bam" -O "${tmp_dir}/${sample_name}.deduped.bam" --CREATE_INDEX true --REMOVE_DUPLICATES true --TMP_DIR "$tmp_dir"
 

    # Add or replace read groups
    gatk AddOrReplaceReadGroups -I "$dup_bam" -O "${tmp_dir}/${sample_name}.picard.bam" -LB "$sample_name" -PL illumina -PU "$sample_name" -SM "$sample_name"
    samtools index "${tmp_dir}/${sample_name}.picard.bam"
 
    # Split NCigar reads
    gatk SplitNCigarReads -R "$ref_genome" -I "${tmp_dir}/${sample_name}.picard.bam" -O "${tmp_dir}/${sample_name}.dedup_split.bam"

    # HaplotypeCaller to generate GVCF
    gatk HaplotypeCaller -R "$ref_genome" -ERC GVCF -I "${tmp_dir}/${sample_name}.dedup_split.bam" -O "${tmp_dir}/${sample_name}.g.vcf" --native-pair-hmm-threads 1
 
    # GenotypeGVCFs to generate VCF
    gatk GenotypeGVCFs -R "$ref_genome" -V "${tmp_dir}/${sample_name}.g.vcf" -O "${tmp_dir}/${sample_name}.vcf"
 
    # Index the VCF file
    gatk IndexFeatureFile -F "${tmp_dir}/${sample_name}.vcf"
 
    # Select SNPs and indels separately
    gatk SelectVariants -V "${tmp_dir}/${sample_name}.vcf" -select-type SNP -O "${tmp_dir}/${sample_name}.snp.vcf"
    gatk SelectVariants -V "${tmp_dir}/${sample_name}.vcf" -select-type INDEL -O "${tmp_dir}/${sample_name}.indel.vcf"
 
    # Index the selected VCF files
    gatk IndexFeatureFile -F "${tmp_dir}/${sample_name}.snp.vcf"
    gatk IndexFeatureFile -F "${tmp_dir}/${sample_name}.indel.vcf"
 
    # Base recalibration (assuming we want to do this step)
    # Note: This typically requires a panel of normals or other known sites
    # For simplicity, we'll use the sample's SNPs and indels as known sites here
    gatk BaseRecalibrator -R "$ref_genome" -I "${tmp_dir}/${sample_name}.dedup_split.bam" --known-sites "${tmp_dir}/${sample_name}.snp.vcf" --known-sites "${tmp_dir}/${sample_name}.indel.vcf" -O "${tmp_dir}/${sample_name}.recal_data.table"
    # Apply BQSR
    gatk ApplyBQSR -R "$ref_genome" -I "${tmp_dir}/${sample_name}.dedup_split.bam" --bqsr-recal-file "${tmp_dir}/${sample_name}.recal_data.table" -O "${tmp_dir}/${sample_name}.bqsr.bam"
    samtools index "${tmp_dir}/${sample_name}.bqsr.bam"
 
    # Split NCigar reads again after BQSR (optional, depending on your workflow)
    gatk SplitNCigarReads -R "$ref_genome" -I "${tmp_dir}/${sample_name}.bqsr.bam" -O "${tmp_dir}/${sample_name}.bqsr.bed.bam"
    gatk HaplotypeCaller -R "$ref_genome" -ERC GVCF -I "${tmp_dir}/${sample_name}.bqsr.bed.bam" -O "${tmp_dir}/${sample_name}.bed.g.vcf" --native-pair-hmm-threads 1

    gatk GenotypeGVCFs -R "$ref_genome" -V "${tmp_dir}/${sample_name}.bed.g.vcf" -O "${tmp_dir}/${sample_name}.bed.vcf"

    # Compress and index the final VCF
    bgzip "${tmp_dir}/${sample_name}.bed.vcf"
    tabix "${tmp_dir}/${sample_name}.bed.vcf.gz"
 
    # Select and compress SNPs and indels separately
    gatk SelectVariants -V "${tmp_dir}/${sample_name}.bed.vcf.gz" -select-type SNP -O "${tmp_dir}/${sample_name}.snp.vcf.gz"
    gatk SelectVariants -V "${tmp_dir}/${sample_name}.bed.vcf.gz" -select-type INDEL -O "${tmp_dir}/${sample_name}.indel.vcf.gz"
 
    # Variant filtration (optional, depending on your workflow)
    gatk VariantFiltration -V "${tmp_dir}/${sample_name}.snp.vcf.gz" -filter "QD < 2.0" --filter-name "QD2" -filter "QUAL < 30.0" --filter-name "QUAL30" -filter "FS > 60.0" --filter-name "FS60" -filter "MQ < 40.0" --filter-name "MQ40" -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" -O "${out}/${sample_name}.vcf"
    # Clean up temporary files
    rm -rf "$tmp_dir"
}
 
# Export function for GNU Parallel (if used)
export -f process_sample
 
# Read sample information from the list file and process each sample in parallel
process_sample "$input_bam" "$sample_name"


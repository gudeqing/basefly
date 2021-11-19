#!/bin/sh

# Copyright (c) 2016-2020 Sentieon Inc. All rights reserved

# *******************************************
# Script to perform TNhaplotyper2 variant calling
# using a matched paired Tumor+normal sample with fastq
# files named normal_1.fastq.gz, normal_2.fastq.gz
# tumor_1.fastq.gz, tumor_2.fastq.gz
# *******************************************

set -eu

# Update with the fullpath location of your sample fastq
# usage: TN.sh tumor_name tumor_r1 tumor_r2 normal_name normal_r1 noraml_r2
TUMOR_SM=$1
TUMOR_RGID="rg_$TUMOR_SM" #read group ID
NORMAL_SM=$4
NORMAL_RGID="rg_$NORMAL_SM" #read group ID
PL="ILLUMINA" #or other sequencing platform
TUMOR_FASTQ_1="$2"
TUMOR_FASTQ_2="$3" #If using Illumina paired data
NORMAL_FASTQ_1="$5"
NORMAL_FASTQ_2="$6"

# Update with the location of the reference data files
FASTA_DIR="/enigma/datasets/broad-genome-references/Homo_sapiens_assembly19_1000genomes_decoy"
FASTA="$FASTA_DIR/Homo_sapiens_assembly19_1000genomes_decoy.fasta"
KNOWN_DBSNP="$FASTA_DIR/Homo_sapiens_assembly19_1000genomes_decoy.dbsnp138.vcf"
KNOWN_INDELS="$FASTA_DIR/Homo_sapiens_assembly19_1000genomes_decoy.known_indels_20120518.vcf"
KNOWN_MILLS="$FASTA_DIR/Mills_and_1000G_gold_standard.indels.b37.sites.vcf"
#CONTAMINATION_VCF="$FASTA_DIR/germline_vcf-af-only-gnomad.raw.sites.vcf"  # A VCF of germline sites to use for contamination detection
# Optional files
#PON=  # the Mutect2 panel-of-normals VCF file
#GERMLINE_VCF=  # A VCF of known germline sites

# Update with the location of the Sentieon software package and license file
SENTIEON_INSTALL_DIR=/opt/sentieon-genomics-202010.01
#export SENTIEON_LICENSE=/home/Licenses/Sentieon.lic #or using licsrvr: c1n11.sentieon.com:5443

# Other settings
NT=$(nproc) #number of threads to use in computation, set to number of cores in the server
START_DIR="$PWD/$TUMOR_SM" #Determine where the output files will be stored





# You do not need to modify any of the lines below unless you want to tweak the pipeline

# ************************************************************************************************************************************************************************



# ******************************************
# 0. Setup
# ******************************************
WORKDIR="$START_DIR"
mkdir -p $WORKDIR
LOGFILE=$WORKDIR/run.log
exec >$LOGFILE 2>&1
cd $WORKDIR

# ******************************************
# 1a. Mapping reads with BWA-MEM, sorting for tumor sample
# ******************************************
#The results of this call are dependent on the number of threads used. To have number of threads independent results, add chunk size option -K 10000000
( $SENTIEON_INSTALL_DIR/bin/sentieon bwa mem -M -R "@RG\tID:$TUMOR_RGID\tSM:$TUMOR_SM\tPL:$PL" \
    -t $NT -K 10000000 $FASTA $TUMOR_FASTQ_1 $TUMOR_FASTQ_2 || \
    { echo -n 'BWA error'; exit 1; } ) | \
    $SENTIEON_INSTALL_DIR/bin/sentieon util sort -o tumor_sorted.bam -t $NT --sam2bam -i - || \
    { echo "Alignment1 failed"; exit 1; }

# ******************************************
# 1b. Mapping reads with BWA-MEM, sorting for normal sample
# ******************************************
#The results of this call are dependent on the number of threads used. To have number of threads independent results, add chunk size option -K 10000000
( $SENTIEON_INSTALL_DIR/bin/sentieon bwa mem -M -R "@RG\tID:$NORMAL_RGID\tSM:$NORMAL_SM\tPL:$PL" \
    -t $NT -K 10000000 $FASTA $NORMAL_FASTQ_1 $NORMAL_FASTQ_2 || \
    { echo -n 'BWA error'; exit 1; } ) | \
    $SENTIEON_INSTALL_DIR/bin/sentieon util sort -o normal_sorted.bam -t $NT --sam2bam -i - || \
    { echo "Alignment2 failed"; exit 1; }

# ******************************************
# 2a. Metrics for tumor sample
# ******************************************
$SENTIEON_INSTALL_DIR/bin/sentieon driver -r $FASTA -t $NT -i tumor_sorted.bam \
    --algo MeanQualityByCycle tumor_mq_metrics.txt \
    --algo QualDistribution tumor_qd_metrics.txt --algo GCBias \
    --summary tumor_gc_summary.txt tumor_gc_metrics.txt --algo AlignmentStat \
    --adapter_seq '' tumor_aln_metrics.txt --algo InsertSizeMetricAlgo tumor_is_metrics.txt || \
    { echo "Metrics1 failed"; exit 1; }

$SENTIEON_INSTALL_DIR/bin/sentieon plot GCBias -o tumor_gc-report.pdf tumor_gc_metrics.txt
$SENTIEON_INSTALL_DIR/bin/sentieon plot QualDistribution \
    -o tumor_qd-report.pdf tumor_qd_metrics.txt
$SENTIEON_INSTALL_DIR/bin/sentieon plot MeanQualityByCycle \
    -o tumor_mq-report.pdf tumor_mq_metrics.txt
$SENTIEON_INSTALL_DIR/bin/sentieon plot InsertSizeMetricAlgo \
    -o tumor_is-report.pdf tumor_is_metrics.txt

# ******************************************
# 2b. Metrics for normal sample
# ******************************************
$SENTIEON_INSTALL_DIR/bin/sentieon driver -r $FASTA -t $NT -i normal_sorted.bam \
    --algo MeanQualityByCycle normal_mq_metrics.txt \
    --algo QualDistribution normal_qd_metrics.txt --algo GCBias \
    --summary normal_gc_summary.txt normal_gc_metrics.txt --algo AlignmentStat \
    --adapter_seq '' normal_aln_metrics.txt --algo InsertSizeMetricAlgo normal_is_metrics.txt || \
    { echo "Metrics2 failed"; exit 1; }

$SENTIEON_INSTALL_DIR/bin/sentieon plot GCBias -o normal_gc-report.pdf normal_gc_metrics.txt
$SENTIEON_INSTALL_DIR/bin/sentieon plot QualDistribution \
    -o normal_qd-report.pdf normal_qd_metrics.txt
$SENTIEON_INSTALL_DIR/bin/sentieon plot MeanQualityByCycle \
    -o normal_mq-report.pdf normal_mq_metrics.txt
$SENTIEON_INSTALL_DIR/bin/sentieon plot InsertSizeMetricAlgo \
    -o normal_is-report.pdf normal_is_metrics.txt

# ******************************************
# 3a. Remove Duplicate Reads for tumor
# sample. It is possible
# to remove instead of mark duplicates
# by adding the --rmdup option in Dedup
# ******************************************
$SENTIEON_INSTALL_DIR/bin/sentieon driver -t $NT -i tumor_sorted.bam --algo LocusCollector \
    --fun score_info tumor_score.txt || { echo "LocusCollector1 failed"; exit 1; }

$SENTIEON_INSTALL_DIR/bin/sentieon driver -t $NT -i tumor_sorted.bam --algo Dedup \
    --score_info tumor_score.txt --metrics tumor_dedup_metrics.txt tumor_deduped.bam || \
    {  echo "Dedup1 failed"; exit 1; }

# ******************************************
# 3b. Remove Duplicate Reads for normal
# sample. It is possible
# to remove instead of mark duplicates
# by adding the --rmdup option in Dedup
# ******************************************
$SENTIEON_INSTALL_DIR/bin/sentieon driver -t $NT -i normal_sorted.bam --algo LocusCollector \
    --fun score_info normal_score.txt || { echo "LocusCollector2 failed"; exit 1; }

$SENTIEON_INSTALL_DIR/bin/sentieon driver -t $NT -i normal_sorted.bam --algo Dedup \
    --score_info normal_score.txt --metrics normal_dedup_metrics.txt normal_deduped.bam || \
    { echo "Dedup2 failed"; exit 1; }

# ******************************************
# 2a. Coverage metrics for normal sample
# ******************************************
$SENTIEON_INSTALL_DIR/bin/sentieon driver -r $FASTA -t $NT -i normal_deduped.bam \
    --algo CoverageMetrics normal_coverage_metrics || \
    { echo "CoverageMetrics1 failed"; exit 1; }

# ******************************************
# 2b. Coverage metrics for tumor sample
# ******************************************
$SENTIEON_INSTALL_DIR/bin/sentieon driver -r $FASTA -t $NT -i tumor_deduped.bam \
    --algo CoverageMetrics tumor_coverage_metrics || { echo "CoverageMetrics2 failed"; exit 1; }

# ******************************************
# 4a. Base recalibration for tumor sample
# ******************************************
$SENTIEON_INSTALL_DIR/bin/sentieon driver -r $FASTA -t $NT -i tumor_deduped.bam \
    --algo QualCal -k $KNOWN_DBSNP -k $KNOWN_MILLS -k $KNOWN_INDELS tumor_recal_data.table
#$SENTIEON_INSTALL_DIR/bin/sentieon driver -r $FASTA -t $NT -i tumor_deduped.bam \
#    -q tumor_recal_data.table --algo QualCal -k $KNOWN_DBSNP -k $KNOWN_MILLS \
#    -k $KNOWN_INDELS tumor_recal_data.table.post
#$SENTIEON_INSTALL_DIR/bin/sentieon driver -t $NT --algo QualCal --plot \
#    --before tumor_recal_data.table --after tumor_recal_data.table.post tumor_recal.csv
#$SENTIEON_INSTALL_DIR/bin/sentieon plot QualCal -o tumor_recal_plots.pdf tumor_recal.csv
# ******************************************
# 4b. Base recalibration for normal sample
# ******************************************
$SENTIEON_INSTALL_DIR/bin/sentieon driver -r $FASTA -t $NT -i normal_deduped.bam \
    --algo QualCal -k $KNOWN_DBSNP -k $KNOWN_MILLS -k $KNOWN_INDELS normal_recal_data.table
#$SENTIEON_INSTALL_DIR/bin/sentieon driver -r $FASTA -t $NT -i normal_deduped.bam \
#    -q normal_recal_data.table --algo QualCal -k $KNOWN_DBSNP -k $KNOWN_MILLS \
#    -k $KNOWN_INDELS normal_recal_data.table.post
#$SENTIEON_INSTALL_DIR/bin/sentieon driver -t $NT --algo QualCal --plot \
#    --before normal_recal_data.table --after normal_recal_data.table.post normal_recal.csv
#$SENTIEON_INSTALL_DIR/bin/sentieon plot QualCal -o normal_recal_plots.pdf normal_recal.csv

# ******************************************
# 5. Somatic Variant Calling
# ******************************************
$SENTIEON_INSTALL_DIR/bin/sentieon driver -r $FASTA -t $NT -i tumor_deduped.bam \
    -i normal_deduped.bam -q tumor_recal_data.table -q normal_recal_data.table \
    --algo TNhaplotyper2 --tumor_sample $TUMOR_SM --normal_sample $NORMAL_SM \
    ${PON:+--pon $PON} ${GERMLINE_VCF:+--germline_vcf $GERMLINE_VCF} output-tnhap2-tmp.vcf.gz \
    --algo OrientationBias --tumor_sample $TUMOR_SM output-orientation \
#    --algo ContaminationModel --tumor_sample $TUMOR_SM --normal_sample $NORMAL_SM \
#    --vcf $CONTAMINATION_VCF \
#    --tumor_segments output-contamination-segments output-contamination || \
    { echo "TNhaplotyper2 failed"; exit 1; }

$SENTIEON_INSTALL_DIR/bin/sentieon driver -r $FASTA --algo TNfilter \
    -v output-tnhap2-tmp.vcf.gz --tumor_sample $TUMOR_SM --normal_sample $NORMAL_SM \
#    --contamination output-contamination --tumor_segments output-contamination-segments \
    --orientation_priors output-orientation output-tnhap2.vcf.gz || \
    { echo "TNfilter failed"; exit 1; }

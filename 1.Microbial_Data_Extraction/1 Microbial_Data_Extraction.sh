#!/bin/bash

if [ $# -ne 4 ]
then
	echo "You should provide five parameters"
	echo "1.sample_name=sample name"
	echo "2.mem=Maximum memory usage"
	echo "3.fqPath=The path of sample_fq"
	echo "4.analysis_folder=The path of analysis"
	echo "eg: sh scripts_nature.sh SRR11770300 80g /data/hour/sahmi/GSE150290-gc/SRR11770300/ /data/hour/sahmi/GSE150290-gc/invadeseq/"
	exit
fi
sample_name=$1
mem=$2
fqPath=$3
analysis_folder=$4
mkdir $analysis_folder
cd $analysis_folder
##################
/data/hour/sahmi/cellranger-6.1.2/cellranger count \
--id=$sample_name \
--transcriptome=/data/hour/refdata-gex-GRCh38-2020-A/ \
--fastqs=$fqPath >> $sample_name.cr.log
##################
/data/hour/sahmi/gatk-4.1.3.0/gatk --java-options "-Xmx${mem} -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" \
PathSeqPipelineSpark --input $analysis_folder/$sample_name/outs/possorted_genome_bam.bam \
--filter-bwa-image /data/hour/sahmi/pathseq_microbe_db/pathseq_host.fa.img \
--kmer-file /data/hour/sahmi/pathseq_microbe_db/pathseq_host.bfi \
--min-clipped-read-length 60 \
--microbe-fasta /data/hour/sahmi/pathseq_microbe_db/pathseq_microbe.fa \
--microbe-bwa-image /data/hour/sahmi/pathseq_microbe_db/pathseq_microbe.fa.img \
--taxonomy-file /data/hour/sahmi/pathseq_microbe_db/pathseq_taxonomy.db \
--output $analysis_folder/$sample_name/outs/pathseq.complete.bam \
--scores-output $analysis_folder/$sample_name/outs/pathseq.complete.csv \
--is-host-aligned false \
--filter-duplicates false \
--min-score-identity 0.7  >> $sample_name.ps.log
##################
python3 /data/hour/sahmi/invadeseq/UMI_annotator.py \
$analysis_folder/$sample_name/outs/possorted_genome_bam.bam \
$sample_name \
$analysis_folder/$sample_name/outs/filtered_feature_bc_matrix/barcodes.tsv.gz \
$analysis_folder/$sample_name/outs/pathseq.complete.bam \
$analysis_folder/$sample_name/outs/pathseq.complete.csv \
$analysis_folder/$sample_name/outs/nova.filtered_matrix.readname \
$analysis_folder/$sample_name/outs/nova.filtered_matrix.unmap_cbub.bam \
$analysis_folder/$sample_name/outs/nova.filtered_matrix.unmap_cbub.fasta \
$analysis_folder/$sample_name/outs/nova.filtered_matrix.list \
$analysis_folder/$sample_name/outs/nova.raw.filtered_matrix.readnamepath \
$analysis_folder/$sample_name/outs/nova.filtered_matrix.genus.cell \
$analysis_folder/$sample_name/outs/nova.filtered_matrix.genus.csv \
$analysis_folder/$sample_name/outs/nova.filtered_matrix.validate.csv

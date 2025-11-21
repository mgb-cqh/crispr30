#!/bin/bash
#SBATCH --job-name=star_align
#SBATCH --partition=patralab,largemem,batch
#SBATCH -N 1		
#SBATCH -n 24
#SBATCH --mem=500Gb
#SBATCH --time=0-24:00:00
#SBATCH --output=%j.out
#SBATCH --error=%j.err
#SBATCH --exclude=s1cmp006,s1cmp007

export PATH="/cluster/tufts/patralab/rbator01/software/STAR-2.7.11b/source/:$PATH"

genomeDir=/cluster/tufts/patralab/rbator01/perlis_lab/drugseq_mar24/ref/STAR/
barcode_file=/cluster/tufts/patralab/rbator01/perlis_lab/drugseq_mar24/241126-0144_Joshua_Bowen-8326/barcodes.txt
bamdir=/cluster/tufts/patralab/rbator01/perlis_lab/drugseq_mar24/cspr_jan25_analysis/align_trim/
R1=/cluster/tufts/patralab/rbator01/perlis_lab/drugseq_mar24/cspr_jan25_analysis/trim/241017_DrugSeq-MQ_UDI_1_S1_L001_R1_001_val_1.fq.gz
R2=/cluster/tufts/patralab/rbator01/perlis_lab/drugseq_mar24/cspr_jan25_analysis/trim/241017_DrugSeq-MQ_UDI_1_S1_L001_R2_001_val_2.fq.gz

STAR --runMode alignReads \
  --outSAMmapqUnique 60 \
  --runThreadN 8 \
  --outSAMunmapped Within \
  --soloStrand Forward \
  --quantMode GeneCounts \
  --outBAMsortingThreadN 8 \
  --genomeDir $genomeDir \
  --soloType CB_UMI_Simple \
  --soloCBstart 1 \
  --soloCBlen 14 \
  --soloUMIstart 15 \
  --soloUMIlen 14 \
  --soloUMIdedup NoDedup 1MM_All \
  --soloCellFilter None \
  --soloCBwhitelist $barcode_file \
  --soloBarcodeReadLength 0 \
  --soloFeatures Gene \
  --outSAMattributes NH HI nM AS CR UR CB UB GX GN sS sQ sM \
  --outFilterMultimapNmax 1 \
  --readFilesCommand zcat \
  --outSAMtype BAM SortedByCoordinate \
  --outFileNamePrefix $bamdir \
  --readFilesIn $R2 $R1


#!/bin/bash

bam_dir=/cluster/tufts/patralab/rbator01/perlis_lab/cspr_apr25/analysis/align/
input_bam=Aligned.sortedByCoord.out.bam
barcode_info=/cluster/tufts/patralab/rbator01/perlis_lab/cspr_apr25/241126-0144_Joshua_Bowen-8326/sample_barcodes.txt
 
mkdir ${bam_dir}/demultiplex/
 
while IFS=$'\t' read -r -a line;
do
  sample_id=${line[0]}
  tag_value=${line[1]}
  tag_value=$(echo "$tag_value" | tr -d '\r')
 
  echo "sampleid" $sample_id
  echo "tagval" $tag_value
  sbatch 02_sbatch_bam_demultiplex_qualimap.sh ${bam_dir} ${input_bam} ${sample_id}  ${tag_value}
done < $barcode_info

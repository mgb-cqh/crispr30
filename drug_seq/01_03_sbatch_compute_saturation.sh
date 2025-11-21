#!/bin/bash
#SBATCH --job-name=demultiplex
#SBATCH --partition=patralab
#SBATCH -N 1    
#SBATCH -n 10
#SBATCH --mem=64Gb
#SBATCH --time=0-24:00:00
#SBATCH --output=%j.out
#SBATCH --error=%j.err

module load bc
module load samtools/1.9
cd /cluster/tufts/patralab/rbator01/perlis_lab/cspr_apr25/align/demultiplex/

# Output file
output="saturation_results.tsv"
echo -e "Sample\tTotal_Reads\tUnique_UMIs\tSaturation" > "$output"

# Loop over deduplicated BAM files
for dedup_bam in *.srt.dedup.bam; do
    # Extract sample prefix (e.g., B09)
    sample_name="${dedup_bam%%.srt.dedup.bam}"

    # Define corresponding raw BAM
    raw_bam="${sample_name}.srt.bam"

    # Check that both files exist
    if [[ -f "$raw_bam" && -f "$dedup_bam" ]]; then
        # Count total and unique reads
        total_reads=$(samtools view -c "$raw_bam")
        unique_umis=$(samtools view -c "$dedup_bam")

        # Calculate saturation
        if [[ $total_reads -gt 0 ]]; then
            saturation=$(echo "scale=4; 1 - ($unique_umis / $total_reads)" | bc)
        else
            saturation="NA"
        fi

        # Write to output file
        echo -e "${sample_name}\t${total_reads}\t${unique_umis}\t${saturation}" >> "$output"
    else
        echo "Warning: Missing file for sample ${sample_name}" >&2
    fi
done



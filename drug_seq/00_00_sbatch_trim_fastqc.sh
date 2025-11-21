#!/bin/bash
#SBATCH --job-name=fastqc
#SBATCH --partition=patralab,preempt,largemem,batch
#SBATCH -N 1		
#SBATCH -n 10
#SBATCH --mem=64Gb
#SBATCH --time=0-24:00:00
#SBATCH --output=%j.out
#SBATCH --error=%j.err

module load fastqc/0.11.8
module load trim-galore/0.6.10

fastqc --outdir /cluster/tufts/patralab/rbator01/perlis_lab/cspr_apr25/241126-0144_Joshua_Bowen-8326/ /cluster/tufts/patralab/rbator01/perlis_lab/cspr_apr25/241126-0144_Joshua_Bowen-8326/*001.fastq.gz

trim_galore --paired -o /cluster/tufts/patralab/rbator01/perlis_lab/cspr_apr25/analysis/trim/ /cluster/tufts/patralab/rbator01/perlis_lab/cspr_apr25/241126-0144_Joshua_Bowen-8326/*001.fastq.gz

# *val*.fq.gz are the final names for trim galore
fastqc --outdir /cluster/tufts/patralab/rbator01/perlis_lab/cspr_apr25/analysis/trim/ /cluster/tufts/patralab/rbator01/perlis_lab/cspr_apr25/analysis/trim/*val*fq.gz

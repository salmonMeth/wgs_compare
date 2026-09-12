#!/bin/sh

#SBATCH --job-name=align_8038

#SBATCH --account=project_2019524

#SBATCH --partition=longrun

#SBATCH --qos=normal

#SBATCH --cpus-per-task=16

#SBATCH --mem=100G

#SBATCH --time=10-00:00:00

#SBATCH --output=align_8038_%j.out

#SBATCH --error=align_8038_%j.err



set -euo pipefail



SAMPLE="8038"



R1="/scratch/project_2019524/PinkSalmon_lowinput/E250168591_L01_B5GANIdazlRAADA-8038_1.fq.gz"

R2="/scratch/project_2019524/PinkSalmon_lowinput/E250168591_L01_B5GANIdazlRAADA-8038_2.fq.gz"



export PATH="$HOME/bin/pinksalmon/bwa-mem2-2.3_x64-linux:$HOME/bin/pinksalmon/samtools/bin:$PATH"



REF="/scratch/project_2019524/pinksalmon_reference/GCF_021184085.1_OgorEven_v1.0_genomic.fna"



OUTDIR="/scratch/project_2019524/pinksalmon_alignments/lowinput"

OUT="${OUTDIR}/${SAMPLE}"



mkdir -p "$OUTDIR"



echo "Pink salmon alignment"

echo "Sample: $SAMPLE"

echo "R1: $R1"

echo "R2: $R2"

echo "Output: ${OUT}.sorted.bam"

echo "Start: $(date)"



if [[ ! -f "$REF" ]]; then

    echo "ERROR: Reference not found: $REF"

    exit 1

fi



if [[ ! -f "$R1" ]]; then

    echo "ERROR: R1 not found: $R1"

    exit 1

fi



if [[ ! -f "$R2" ]]; then

    echo "ERROR: R2 not found: $R2"

    exit 1

fi



echo "BWA-MEM2:"

bwa-mem2 version



echo "Samtools:"

samtools --version | head -1



echo "Starting BWA-MEM2 alignment..."



bwa-mem2 mem -t 16 -R "@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tLB:LowInput\tPL:DNBSEQ" "$REF" "$R1" "$R2" | samtools sort -@ 4 -m 4G -o "${OUT}.sorted.bam" -



echo "Indexing BAM..."



samtools index "${OUT}.sorted.bam"



echo "Generating flagstat..."



samtools flagstat "${OUT}.sorted.bam" > "${OUT}.flagstat.txt"



echo "Generating stats..."



samtools stats "${OUT}.sorted.bam" > "${OUT}.stats.txt"



echo "Alignment complete"

echo "Finished: $(date)"



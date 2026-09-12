#!/bin/bash
#SBATCH --job-name=coverage_normal
#SBATCH --account=project_2019524
#SBATCH --partition=longrun
#SBATCH --qos=normal
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=2-00:00:00
#SBATCH --output=coverage_%j.out
#SBATCH --error=coverage_%j.err

set -euo pipefail


#change this to different subdirectories to accomodate the different libraries
BAMDIR="/scratch/project_2019524/pinksalmon_alignments/normalinput"

#change to the path of the samtools binary if needed
export PATH="$HOME/bin/pinksalmon/samtools/bin:$PATH"



#change the sample IDs according to the IDs present in each library
for SAMPLE in 8167 8168
do

    SAMPLE_DIR="${BAMDIR}/${SAMPLE}_align"
    BAM="${SAMPLE_DIR}/${SAMPLE}.sorted.bam"


    echo "=========================================="
    echo "Sample: $SAMPLE"
    echo "BAM: $BAM"
    echo "Start: $(date)"
    echo "=========================================="


    if [[ ! -f "$BAM" ]]; then
        echo "ERROR: BAM not found: $BAM"
        exit 1
    fi

    if [[ ! -f "${BAM}.bai" ]]; then
        echo "ERROR: BAM index not found: ${BAM}.bai"
        exit 1
    fi

    echo "Running samtools coverage..."
    samtools coverage "$BAM" > "${BAMDIR}/${SAMPLE}.coverage.txt"

    echo "Running samtools depth..."
    samtools depth -a "$BAM" > "${BAMDIR}/${SAMPLE}.depth.txt"

    echo "Finished $SAMPLE"
done

echo "=========================================="
echo "All coverage analyses complete"
echo "Finished: $(date)"
echo "=========================================="

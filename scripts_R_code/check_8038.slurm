#!/bin/bash

#SBATCH --job-name=depth_windows_normal_new
#SBATCH --account=project_2019524
#SBATCH --partition=longrun
#SBATCH --qos=normal
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=2-00:00:00
#SBATCH --output=depth_windows_new_%j.out
#SBATCH --error=depth_windows_new_%j.err

set -euo pipefail

BAMDIR="/scratch/project_2019524/pinksalmon_alignments/lowinput"

export PATH="$HOME/bin/pinksalmon/samtools/bin:$PATH"

for SAMPLE in 8038
do

    SAMPLE_DIR="${BAMDIR}/${SAMPLE}_align"
    BAM="${SAMPLE_DIR}/${SAMPLE}.sorted.bam"

    WINDOW=100000
    OUTDIR="${BAMDIR}/depth_windows"

    mkdir -p "$OUTDIR"

    DEPTH="${BAMDIR}/${SAMPLE}.depth_new.txt"
    OUT="${OUTDIR}/${SAMPLE}.depth_${WINDOW}bp_new.txt"

    echo "=========================================="
    echo "Sample: $SAMPLE"
    echo "BAM: $BAM"
    echo "Depth: $DEPTH"
    echo "Windows: $OUT"
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

    echo "Running samtools depth..."

    samtools depth -a "$BAM" > "$DEPTH"

    echo "Generating ${WINDOW} bp windows..."

    awk -v window="$WINDOW" '
    BEGIN {
        OFS="\t"
    }

    {
        chrom = $1
        pos   = $2
        depth = $3

        window_start = int((pos - 1) / window) * window + 1

        key = chrom "\t" window_start

        count[key]++
        sum[key] += depth
    }

    END {
        for (key in count) {
            split(key, a, "\t")

            chrom = a[1]
            start = a[2]
            end = start + window - 1

            mean_depth = sum[key] / count[key]

            print chrom, start, end, mean_depth
        }
    }
    ' "$DEPTH" |
    sort -k1,1 -k2,2n |
    awk '
    BEGIN {
        OFS="\t"
        print "chromosome", "window_start", "window_end", "mean_depth"
    }
    {
        print
    }
    ' > "$OUT"

    echo "Finished $SAMPLE"
    echo

done

echo "=========================================="
echo "All NEW depth and window analyses complete"
echo "Finished: $(date)"
echo "=========================================="

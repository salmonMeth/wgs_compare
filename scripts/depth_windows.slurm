#!/bin/bash

#SBATCH --job-name=depth_windows
#SBATCH --account=project_2019524
#SBATCH --partition=longrun
#SBATCH --qos=normal
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=2-00:00:00
#SBATCH --output=depth_windows_low_inp%j.out
#SBATCH --error=depth_windows_low_inp%j.err

set -euo pipefail


# change these accordingly to the path olf the folder containing the depth.txt files
BAMDIR="/scratch/project_2019524/pinksalmon_alignments/normalinput/coverage_depth"

SAMPLES=(
    8167
    8168
)

# Window size
WINDOW=100000

OUTDIR="${BAMDIR}/depth_windows"

mkdir -p "$OUTDIR"

for SAMPLE in "${SAMPLES[@]}"
do

    DEPTH="${BAMDIR}/${SAMPLE}.depth.txt"
    OUT="${OUTDIR}/${SAMPLE}.depth_${WINDOW}bp.txt"

    echo "=========================================="
    echo "Sample: $SAMPLE"
    echo "Input: $DEPTH"
    echo "Output: $OUT"
    echo "Start: $(date)"
    echo "=========================================="


    if [[ ! -f "$DEPTH" ]]; then
        echo "ERROR: depth file not found:"
        echo "$DEPTH"
        exit 1
    fi


    awk -v window="$WINDOW" '
    BEGIN {
        OFS="\t"
    }

    NR == 1 && $1 == "chromosome" {
        next
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
echo "All samples processed"
echo "Finished: $(date)"
echo "=========================================="
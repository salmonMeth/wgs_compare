#!/bin/bash

#SBATCH --job-name=depth_stats_8038
#SBATCH --account=project_2019524
#SBATCH --partition=longrun
#SBATCH --qos=normal
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=1-00:00:00
#SBATCH --output=depth_stats_8038_lowi%j.out
#SBATCH --error=depth_stats_8038_lowi%j.err

set -euo pipefail

DEPTHDIR="/scratch/project_2019524/pinksalmon_alignments/lowinput/coverage_depth"

SAMPLE=8038

OUTDIR="${DEPTHDIR}/depth_stats"

mkdir -p "$OUTDIR"

DEPTH="${DEPTHDIR}/${SAMPLE}.depth.txt"

GENOME_OUT="${OUTDIR}/${SAMPLE}.depth_stats.txt"
CHR_OUT="${OUTDIR}/${SAMPLE}.chromosome_depth_stats.txt"

echo "============================================"
echo "Processing $SAMPLE"
echo "Input: $DEPTH"
echo "============================================"

if [[ ! -f "$DEPTH" ]]; then
    echo "ERROR: depth file not found:"
    echo "$DEPTH"
    exit 1
fi

# --------------------------------------------------
# Genome-wide statistics
# --------------------------------------------------

awk '
BEGIN {
    total = 0
    d5 = 0
    d8 = 0
    d10 = 0
    d12 = 0
    d15 = 0
    sumdepth = 0
}

NR == 1 && $1 == "chromosome" {
    next
}

{
    depth = $3

    total++

    if (depth >= 5)  d5++
    if (depth >= 8)  d8++
    if (depth >= 10) d10++
    if (depth >= 12) d12++
    if (depth >= 15) d15++

    sumdepth += depth
}

END {
    print "Sample\tTotal_bases\tMean_depth\tPercent_5x\tPercent_8x\tPercent_10x\tPercent_12x\tPercent_15x"

    printf "%s\t%d\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\n",
           "'"$SAMPLE"'",
           total,
           sumdepth / total,
           100 * d5 / total,
           100 * d8 / total,
           100 * d10 / total,
           100 * d12 / total,
           100 * d15 / total
}
' "$DEPTH" > "$GENOME_OUT"


# --------------------------------------------------
# Per-chromosome statistics
# --------------------------------------------------

awk '
BEGIN {
    OFS = "\t"
}

NR == 1 && $1 == "chromosome" {
    next
}

{
    chrom = $1
    depth = $3

    total[chrom]++
    sumdepth[chrom] += depth

    if (depth >= 5)  d5[chrom]++
    if (depth >= 8)  d8[chrom]++
    if (depth >= 10) d10[chrom]++
    if (depth >= 12) d12[chrom]++
    if (depth >= 15) d15[chrom]++
}

END {

    print "Sample",
          "Chromosome",
          "Total_bases",
          "Mean_depth",
          "Percent_5x",
          "Percent_8x",
          "Percent_10x",
          "Percent_12x",
          "Percent_15x"

    for (chrom in total) {

        printf "%s\t%s\t%d\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\n",
               "'"$SAMPLE"'",
               chrom,
               total[chrom],
               sumdepth[chrom] / total[chrom],
               100 * d5[chrom] / total[chrom],
               100 * d8[chrom] / total[chrom],
               100 * d10[chrom] / total[chrom],
               100 * d12[chrom] / total[chrom],
               100 * d15[chrom] / total[chrom]
    }
}
' "$DEPTH" |
awk 'NR == 1 {print; next} {print | "sort -k2,2"}' \
> "$CHR_OUT"


echo "Genome-wide output:"
echo "$GENOME_OUT"

echo "Chromosome-level output:"
echo "$CHR_OUT"

echo
echo "============================================"
echo "8038 finished"
echo "Finished: $(date)"
echo "============================================"

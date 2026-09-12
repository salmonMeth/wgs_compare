#!/bin/bash
#SBATCH --job-name=depth_evenness_short
#SBATCH --output=depth_evenness_shortlow%j.out
#SBATCH --error=depth_evenness_shortlow%j.err
#SBATCH --account=project_2019524
#SBATCH --partition=longrun
#SBATCH --qos=normal
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=1-00:00:00

# ============================================================
# Calculate fast depth evenness statistics
#
# Statistics:
#   Mean
#   SD
#   CV
#   Min
#   Max
#
# Calculated:
#   - genome-wide
#   - per chromosome
#
# No large AWK arrays
# No sorting
# One streaming pass through each depth file
# ============================================================

METHOD="lowinput"

ALIGNMENT_BASE="/scratch/project_2019524/pinksalmon_alignments"
DEPTHDIR="${ALIGNMENT_BASE}/${METHOD}/coverage_depth"

OUTDIR="/scratch/project_2019524/pinksalmon_depth_evenness_short"
mkdir -p "${OUTDIR}"

SAMPLES=(8035 8036 8038)

# Chromosome names
CHR_NAMES=(
  NC_060173.1
  NC_060174.1
  NC_060175.1
  NC_060176.1
  NC_060177.1
  NC_060178.1
  NC_060179.1
  NC_060180.1
  NC_060181.1
  NC_060182.1
  NC_060183.1
  NC_060184.1
  NC_060185.1
  NC_060186.1
  NC_060187.1
  NC_060188.1
  NC_060189.1
  NC_060190.1
  NC_060191.1
  NC_060192.1
  NC_060193.1
  NC_060194.1
  NC_060195.1
  NC_060196.1
  NC_060197.1
  NC_060198.1
  NC_060199.1
)

for SAMPLE in "${SAMPLES[@]}"; do

DEPTH_FILE="${DEPTHDIR}/${SAMPLE}.depth.txt"
OUTFILE="${OUTDIR}/${SAMPLE}_${METHOD}_depth_evenness.txt"

if [[ ! -f "${DEPTH_FILE}" ]]; then
echo "WARNING: Missing depth file:"
echo "  ${DEPTH_FILE}"
continue
fi

echo "Processing ${SAMPLE} ${METHOD}"
echo "Input:  ${DEPTH_FILE}"
echo "Output: ${OUTFILE}"

# --------------------------------------------------------
# Write header
# --------------------------------------------------------

echo -e \
"Sample\tMethod\tLevel\tchr_number\tchromosome\tMean_depth\tSD_depth\tCV_depth\tMin_depth\tMax_depth" \
> "${OUTFILE}"

# --------------------------------------------------------
# One streaming AWK pass
#
# Assumes depth.txt has:
# chromosome position depth
#
# i.e.
# NC_060173.1  1  12
# NC_060173.1  2  14
# ...
# --------------------------------------------------------

awk \
    -v sample="${SAMPLE}" \
    -v method="${METHOD}" \
'
BEGIN {
    OFS="\t"

    chr[1]  = "NC_060173.1"
    chr[2]  = "NC_060174.1"
    chr[3]  = "NC_060175.1"
    chr[4]  = "NC_060176.1"
    chr[5]  = "NC_060177.1"
    chr[6]  = "NC_060178.1"
    chr[7]  = "NC_060179.1"
    chr[8]  = "NC_060180.1"
    chr[9]  = "NC_060181.1"
    chr[10] = "NC_060182.1"
    chr[11] = "NC_060183.1"
    chr[12] = "NC_060184.1"
    chr[13] = "NC_060185.1"
    chr[14] = "NC_060186.1"
    chr[15] = "NC_060187.1"
    chr[16] = "NC_060188.1"
    chr[17] = "NC_060189.1"
    chr[18] = "NC_060190.1"
    chr[19] = "NC_060191.1"
    chr[20] = "NC_060192.1"
    chr[21] = "NC_060193.1"
    chr[22] = "NC_060194.1"
    chr[23] = "NC_060195.1"
    chr[24] = "NC_060196.1"
    chr[25] = "NC_060197.1"
    chr[26] = "NC_060198.1"
    chr[27] = "NC_060199.1"

    for (i = 1; i <= 27; i++) {
        chr_number[chr[i]] = i
    }
}

{
    chromosome = $1
    depth = $3

    if (depth == "" || depth !~ /^[0-9.]+$/)
        next

    # Genome-wide
    n_all++
    sum_all += depth
    sumsq_all += depth * depth

    if (n_all == 1 || depth < min_all)
        min_all = depth

    if (n_all == 1 || depth > max_all)
        max_all = depth

    # Per chromosome
    n[chromosome]++
    sum[chromosome] += depth
    sumsq[chromosome] += depth * depth

    if (n[chromosome] == 1 || depth < min[chromosome])
        min[chromosome] = depth

    if (n[chromosome] == 1 || depth > max[chromosome])
        max[chromosome] = depth
}

END {

    # Genome-wide statistics

    mean_all = sum_all / n_all
    variance_all = (sumsq_all / n_all) - (mean_all * mean_all)

    if (variance_all < 0)
        variance_all = 0

    sd_all = sqrt(variance_all)

    if (mean_all > 0)
        cv_all = sd_all / mean_all
    else
        cv_all = "NA"

    print sample, method, "Genome-wide", "NA", "Genome-wide", mean_all, sd_all, cv_all, min_all, max_all

    # Per chromosome statistics

    for (i = 1; i <= 27; i++) {

        chromosome = chr[i]

        if (!(chromosome in n))
            continue

        mean_chr = sum[chromosome] / n[chromosome]
        variance_chr = (sumsq[chromosome] / n[chromosome]) - (mean_chr * mean_chr)

        if (variance_chr < 0)
            variance_chr = 0

        sd_chr = sqrt(variance_chr)

        if (mean_chr > 0)
            cv_chr = sd_chr / mean_chr
        else
            cv_chr = "NA"

        print sample, method, "Chromosome", i, chromosome, mean_chr, sd_chr, cv_chr, min[chromosome], max[chromosome]
    }
}
' "${DEPTH_FILE}" >> "${OUTFILE}"


echo "Finished ${SAMPLE}"
echo

done

echo "All samples completed."

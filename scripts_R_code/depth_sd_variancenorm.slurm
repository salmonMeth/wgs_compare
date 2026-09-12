#!/bin/bash
#SBATCH --job-name=depth_evenness
#SBATCH --output=depth_evenness_norm%j.out
#SBATCH --error=depth_evenness_norm%j.err
#SBATCH --account=project_2019524
#SBATCH --partition=longrun
#SBATCH --qos=normal
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=1-00:00:00

# ============================================================
# Calculate depth distribution statistics
#
# Statistics:
#   Mean
#   Standard deviation
#   Coefficient of variation (CV)
#   Median
#   25th percentile
#   75th percentile
#   Minimum
#   Maximum
#
# Calculated:
#   - Genome-wide across chromosomes 1-27
#   - Separately for each chromosome
#
# Only the 27 chromosomes are processed.
# Other contigs/scaffolds are ignored.
# ============================================================

METHOD="normalinput"

DEPTHDIR="/scratch/project_2019524/pinksalmon_alignments/${METHOD}/coverage_depth"

OUTDIR="${DEPTHDIR}/depth_stats"

mkdir -p "$OUTDIR"

# ------------------------------------------------------------
# Temporary directory
# ------------------------------------------------------------

if [ -n "$SLURM_TMPDIR" ] && [ -d "$SLURM_TMPDIR" ]; then
    TMPDIR="${SLURM_TMPDIR}/depth_evenness_${SLURM_JOB_ID}"
else
    TMPDIR="${OUTDIR}/tmp_depth_evenness_${SLURM_JOB_ID}"
fi

mkdir -p "$TMPDIR"

trap 'rm -rf "$TMPDIR"' EXIT

# ------------------------------------------------------------
# Function to calculate statistics from a sorted depth file
#
# Input:
#   sorted file containing one depth value per line
#
# Output:
#   mean, SD, CV, median, Q25, Q75, min, max
# ------------------------------------------------------------

calculate_stats() {

    SORTED_FILE="$1"

    awk '
    {
        x[NR] = $1

        sum += $1
        sumsq += ($1 * $1)

        if (NR == 1) {
            min = $1
            max = $1
        }

        if ($1 < min)
            min = $1

        if ($1 > max)
            max = $1
    }

    END {

        n = NR

        if (n == 0) {
            print "NA\tNA\tNA\tNA\tNA\tNA\tNA\tNA"
            exit
        }

        mean = sum / n

        # Population variance
        variance = (sumsq / n) - (mean * mean)

        # Protect against tiny negative values caused by rounding
        if (variance < 0)
            variance = 0

        sd = sqrt(variance)

        if (mean > 0)
            cv = sd / mean
        else
            cv = "NA"

        # ----------------------------------------------------
        # Median
        # ----------------------------------------------------

        if (n % 2 == 1) {
            median = x[(n + 1) / 2]
        }
        else {
            median = (x[n / 2] + x[n / 2 + 1]) / 2
        }

        # ----------------------------------------------------
        # Quantile function
        #
        # Uses linear interpolation between ordered values.
        # ----------------------------------------------------

        p25_position = 1 + (n - 1) * 0.25
        p75_position = 1 + (n - 1) * 0.75

        lower = int(p25_position)
        upper = lower + 1
        fraction = p25_position - lower

        if (upper > n)
            q25 = x[n]
        else
            q25 = x[lower] + fraction * (x[upper] - x[lower])

        lower = int(p75_position)
        upper = lower + 1
        fraction = p75_position - lower

        if (upper > n)
            q75 = x[n]
        else
            q75 = x[lower] + fraction * (x[upper] - x[lower])

        printf "%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%d\t%d\n", \
            mean, sd, cv, median, q25, q75, min, max
    }
    ' "$SORTED_FILE"
}


# ============================================================
# Process samples
# ============================================================

for SAMPLE in 8167 8168
do

    # --------------------------------------------------------
    # Set depth file
    # --------------------------------------------------------

    if [ "$SAMPLE" = "8038" ]; then
        DEPTH_FILE="${DEPTHDIR}/8038.depth_new.txt"
    else
        DEPTH_FILE="${DEPTHDIR}/${SAMPLE}.depth.txt"
    fi

    OUT="${OUTDIR}/${SAMPLE}.depth_evenness.txt"

    echo "================================================"
    echo "Processing sample: ${SAMPLE}"
    echo "Method: ${METHOD}"
    echo "Depth file: ${DEPTH_FILE}"
    echo "Output: ${OUT}"
    echo "================================================"

    if [ ! -f "$DEPTH_FILE" ]; then
        echo "WARNING: depth file not found:"
        echo "$DEPTH_FILE"
        echo "Skipping ${SAMPLE}"
        echo
        continue
    fi

    # create temp file for the chromosomes

    printf "Sample\tMethod\tLevel\tchr_number\tchromosome\tMean_depth\tSD_depth\tCV_depth\tMedian_depth\tQ25_depth\tQ75_depth\tMin_depth\tMax_depth\n" \
        > "$OUT"


    #get just the chromosomes

    GENOME_DEPTH="${TMPDIR}/${SAMPLE}.genome.depth"

    awk '
        $1 ~ /^NC_060173\.1$/ ||
        $1 ~ /^NC_060174\.1$/ ||
        $1 ~ /^NC_060175\.1$/ ||
        $1 ~ /^NC_060176\.1$/ ||
        $1 ~ /^NC_060177\.1$/ ||
        $1 ~ /^NC_060178\.1$/ ||
        $1 ~ /^NC_060179\.1$/ ||
        $1 ~ /^NC_060180\.1$/ ||
        $1 ~ /^NC_060181\.1$/ ||
        $1 ~ /^NC_060182\.1$/ ||
        $1 ~ /^NC_060183\.1$/ ||
        $1 ~ /^NC_060184\.1$/ ||
        $1 ~ /^NC_060185\.1$/ ||
        $1 ~ /^NC_060186\.1$/ ||
        $1 ~ /^NC_060187\.1$/ ||
        $1 ~ /^NC_060188\.1$/ ||
        $1 ~ /^NC_060189\.1$/ ||
        $1 ~ /^NC_060190\.1$/ ||
        $1 ~ /^NC_060191\.1$/ ||
        $1 ~ /^NC_060192\.1$/ ||
        $1 ~ /^NC_060193\.1$/ ||
        $1 ~ /^NC_060194\.1$/ ||
        $1 ~ /^NC_060195\.1$/ ||
        $1 ~ /^NC_060196\.1$/ ||
        $1 ~ /^NC_060197\.1$/ ||
        $1 ~ /^NC_060198\.1$/ ||
        $1 ~ /^NC_060199\.1$/ {
            print $3
        }
    ' "$DEPTH_FILE" > "$GENOME_DEPTH"

    # --------------------------------------------------------
    # Sort genome-wide depth values
    # --------------------------------------------------------

    GENOME_SORTED="${TMPDIR}/${SAMPLE}.genome.depth.sorted"

    sort -n "$GENOME_DEPTH" > "$GENOME_SORTED"

    # --------------------------------------------------------
    # Genome-wide statistics
    # --------------------------------------------------------

    GENOME_STATS=$(calculate_stats "$GENOME_SORTED")

    printf "%s\t%s\tgenome\tNA\tgenome\t%s\n" \
        "$SAMPLE" \
        "$METHOD" \
        "$GENOME_STATS" \
        >> "$OUT"


    # ========================================================
    # Per chromosome
    # ========================================================

    for CHR_NUM in $(seq 1 27)
    do

        CHR=$(printf "NC_%06d.1" $((60172 + CHR_NUM)))

        CHR_DEPTH="${TMPDIR}/${SAMPLE}.chr${CHR_NUM}.depth"

        CHR_SORTED="${TMPDIR}/${SAMPLE}.chr${CHR_NUM}.depth.sorted"

        # ----------------------------------------------------
        # Extract this chromosome only
        # ----------------------------------------------------

        awk -v chr="$CHR" '
            $1 == chr {
                print $3
            }
        ' "$DEPTH_FILE" > "$CHR_DEPTH"

        # ----------------------------------------------------
        # Sort depths
        # ----------------------------------------------------

        sort -n "$CHR_DEPTH" > "$CHR_SORTED"

        # ----------------------------------------------------
        # Calculate statistics
        # ----------------------------------------------------

        CHR_STATS=$(calculate_stats "$CHR_SORTED")

        # ----------------------------------------------------
        # Write result
        # ----------------------------------------------------

        printf "%s\t%s\tchromosome\t%d\t%s\t%s\n" \
            "$SAMPLE" \
            "$METHOD" \
            "$CHR_NUM" \
            "$CHR" \
            "$CHR_STATS" \
            >> "$OUT"

        echo "  Finished chromosome ${CHR_NUM}: ${CHR}"

    done

    echo
    echo "Finished ${SAMPLE}"
    echo "Output: ${OUT}"
    echo

done

echo "================================================"
echo "All samples finished."
echo "================================================"

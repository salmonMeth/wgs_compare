#!/bin/bash
#SBATCH --job-name=depth_evenness
#SBATCH --output=depth_evenness_%j.out
#SBATCH --error=depth_evenness_%j.err
#SBATCH --account=project_2019524
#SBATCH --partition=longrun
#SBATCH --qos=normal
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=1-00:00:00


METHOD="lowinput"

DEPTHDIR="/scratch/project_2019524/pinksalmon_alignments/${METHOD}/coverage_depth"

# NEW output folder
OUTDIR="/scratch/project_2019524/pinksalmon_depth_evenness_stats"

mkdir -p "$OUTDIR"


# temp dir

if [ -n "$SLURM_TMPDIR" ] && [ -d "$SLURM_TMPDIR" ]; then
    TMPDIR="${SLURM_TMPDIR}/depth_evenness_${SLURM_JOB_ID}"
else
    TMPDIR="${OUTDIR}/tmp_${SLURM_JOB_ID}"
fi

mkdir -p "$TMPDIR"


# ============================================================
# Cleanup
# ============================================================

cleanup() {
    echo
    echo "Cleaning temporary files..."
    rm -rf "$TMPDIR"
}

trap cleanup EXIT INT TERM HUP


calculate_basic_stats() {

    SORTED_FILE="$1"

    awk '
    {
        x = $1

        n++

        # Welford online mean/variance calculation
        delta = x - mean
        mean += delta / n
        delta2 = x - mean
        M2 += delta * delta2

        if (n == 1) {
            min = x
            max = x
        }

        if (x < min)
            min = x

        if (x > max)
            max = x
    }

    END {

        if (n == 0) {
            print "NA\tNA\tNA\tNA\tNA\t0"
            exit
        }

        variance = M2 / n
        sd = sqrt(variance)

        if (mean > 0)
            cv = sd / mean
        else
            cv = "NA"

        printf "%.6f\t%.6f\t%.6f\t%d\t%d\t%d\n", \
            mean, sd, cv, min, max, n
    }
    ' "$SORTED_FILE"
}


calculate_quantiles() {

    SORTED_FILE="$1"

    awk '
    NR == FNR {
        n++
        next
    }

    END {
        if (n == 0) {
            print "NA\tNA\tNA"
            exit
        }

        # Calculate positions

        p25 = 1 + (n - 1) * 0.25
        p50 = 1 + (n - 1) * 0.50
        p75 = 1 + (n - 1) * 0.75

        p25_low = int(p25)
        p25_high = p25_low + 1
        p25_frac = p25 - p25_low

        p50_low = int(p50)
        p50_high = p50_low + 1
        p50_frac = p50 - p50_low

        p75_low = int(p75)
        p75_high = p75_low + 1
        p75_frac = p75 - p75_low

        # ----------------------------------------------------
        # Store only the few values required.
        # ----------------------------------------------------

        while ((getline line < "'$SORTED_FILE'") > 0) {

            split(line, a, " ")
            value = a[1]

            i++

            if (i == p25_low)
                q25_low_value = value

            if (i == p25_high)
                q25_high_value = value

            if (i == p50_low)
                q50_low_value = value

            if (i == p50_high)
                q50_high_value = value

            if (i == p75_low)
                q75_low_value = value

            if (i == p75_high)
                q75_high_value = value
        }

        # Q25
        if (p25_high > n)
            q25 = q25_low_value
        else
            q25 = q25_low_value + \
                  p25_frac * (q25_high_value - q25_low_value)

        # Median
        if (p50_high > n)
            median = q50_low_value
        else
            median = q50_low_value + \
                     p50_frac * (q50_high_value - q50_low_value)

        # Q75
        if (p75_high > n)
            q75 = q75_low_value
        else
            q75 = q75_low_value + \
                  p75_frac * (q75_high_value - q75_low_value)

        printf "%.6f\t%.6f\t%.6f\n", median, q25, q75
    }
    ' "$SORTED_FILE" "$SORTED_FILE"
}


# ============================================================
# Function:
# Calculate all statistics
# ============================================================

calculate_stats() {

    SORTED_FILE="$1"

    BASIC=$(calculate_basic_stats "$SORTED_FILE")

    QUANTILES=$(calculate_quantiles "$SORTED_FILE")

    # BASIC:
    #   1 mean
    #   2 SD
    #   3 CV
    #   4 min
    #   5 max
    #   6 n
    #
    # QUANTILES:
    #   1 median
    #   2 Q25
    #   3 Q75

    echo -e "${BASIC}\t${QUANTILES}" |
    awk -F '\t' '
    {
        printf "%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%d\t%d\n", \
            $1,       \
            $2,       \
            $3,       \
            $7,       \
            $8,       \
            $9,       \
            $4,       \
            $5
    }
    '
}


# loop

for SAMPLE in 8034 
do

    echo
    echo "========================================================"
    echo "Processing sample: ${SAMPLE}"
    echo "Method: ${METHOD}"
    echo "========================================================"

    DEPTH_FILE="${DEPTHDIR}/${SAMPLE}.depth.txt"

    # New output location
    OUT="${OUTDIR}/${SAMPLE}.depth_evenness.txt"

    echo "Depth file:"
    echo "${DEPTH_FILE}"

    echo
    echo "Output:"
    echo "${OUT}"


    if [ ! -f "$DEPTH_FILE" ]; then

        echo
        echo "WARNING: depth file not found:"
        echo "$DEPTH_FILE"
        echo "Skipping ${SAMPLE}"
        continue

    fi



    printf "Sample\tMethod\tLevel\tchr_number\tchromosome\tMean_depth\tSD_depth\tCV_depth\tMedian_depth\tQ25_depth\tQ75_depth\tMin_depth\tMax_depth\n" \
        > "$OUT"


    # ========================================================
    # Create genome-wide depth file
    # ========================================================

    GENOME_DEPTH="${TMPDIR}/${SAMPLE}.genome.depth"

    awk '
        $1 == "NC_060173.1" ||
        $1 == "NC_060174.1" ||
        $1 == "NC_060175.1" ||
        $1 == "NC_060176.1" ||
        $1 == "NC_060177.1" ||
        $1 == "NC_060178.1" ||
        $1 == "NC_060179.1" ||
        $1 == "NC_060180.1" ||
        $1 == "NC_060181.1" ||
        $1 == "NC_060182.1" ||
        $1 == "NC_060183.1" ||
        $1 == "NC_060184.1" ||
        $1 == "NC_060185.1" ||
        $1 == "NC_060186.1" ||
        $1 == "NC_060187.1" ||
        $1 == "NC_060188.1" ||
        $1 == "NC_060189.1" ||
        $1 == "NC_060190.1" ||
        $1 == "NC_060191.1" ||
        $1 == "NC_060192.1" ||
        $1 == "NC_060193.1" ||
        $1 == "NC_060194.1" ||
        $1 == "NC_060195.1" ||
        $1 == "NC_060196.1" ||
        $1 == "NC_060197.1" ||
        $1 == "NC_060198.1" ||
        $1 == "NC_060199.1" {
            print $3
        }
    ' "$DEPTH_FILE" > "$GENOME_DEPTH"


    GENOME_LINES=$(wc -l < "$GENOME_DEPTH")

    echo "Depth values from chromosomes 1-27: ${GENOME_LINES}"


    if [ "$GENOME_LINES" -eq 0 ]; then

        echo "ERROR: No chromosome 1-27 depth values found."
        echo "Check the format of the depth file."

        continue

    fi


    # ========================================================
    # Sort genome-wide depth
    # ========================================================

    GENOME_SORTED="${TMPDIR}/${SAMPLE}.genome.sorted"

    echo
    echo "Sorting genome-wide depths..."

    sort -n "$GENOME_DEPTH" > "$GENOME_SORTED"


    # ========================================================
    # Genome-wide statistics
    # ========================================================

    echo "Calculating genome-wide statistics..."

    GENOME_STATS=$(calculate_stats "$GENOME_SORTED")


    if [ -z "$GENOME_STATS" ]; then

        echo "ERROR: Genome statistics calculation failed."

    else

        printf "%s\t%s\tgenome\tNA\tgenome\t%s\n" \
            "$SAMPLE" \
            "$METHOD" \
            "$GENOME_STATS" \
            >> "$OUT"

        echo "Genome-wide statistics completed."

    fi


    # ========================================================
    # Per chromosome
    # ========================================================

    for CHR_NUM in $(seq 1 27)
    do

        # Reference:
        # chromosome 1 = NC_060173.1
        # chromosome 27 = NC_060199.1

        CHR=$(printf "NC_%06d.1" $((60172 + CHR_NUM)))

        CHR_DEPTH="${TMPDIR}/${SAMPLE}.chr${CHR_NUM}.depth"

        CHR_SORTED="${TMPDIR}/${SAMPLE}.chr${CHR_NUM}.sorted"


        # ----------------------------------------------------
        # Extract chromosome
        # ----------------------------------------------------

        awk -v chr="$CHR" '
            $1 == chr {
                print $3
            }
        ' "$DEPTH_FILE" > "$CHR_DEPTH"


        CHR_LINES=$(wc -l < "$CHR_DEPTH")


        # ----------------------------------------------------
        # Check chromosome
        # ----------------------------------------------------

        if [ "$CHR_LINES" -eq 0 ]; then

            echo "WARNING: ${CHR} contains no depth values."

            printf "%s\t%s\tchromosome\t%d\t%s\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\n" \
                "$SAMPLE" \
                "$METHOD" \
                "$CHR_NUM" \
                "$CHR" \
                >> "$OUT"

            continue

        fi


        # ----------------------------------------------------
        # Sort
        # ----------------------------------------------------

        sort -n "$CHR_DEPTH" > "$CHR_SORTED"


        # ----------------------------------------------------
        # Calculate statistics
        # ----------------------------------------------------

        CHR_STATS=$(calculate_stats "$CHR_SORTED")


        # ----------------------------------------------------
        # Write result
        # ----------------------------------------------------

        if [ -z "$CHR_STATS" ]; then

            echo "ERROR: Statistics failed for ${CHR}"

            printf "%s\t%s\tchromosome\t%d\t%s\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\n" \
                "$SAMPLE" \
                "$METHOD" \
                "$CHR_NUM" \
                "$CHR" \
                >> "$OUT"

        else

            printf "%s\t%s\tchromosome\t%d\t%s\t%s\n" \
                "$SAMPLE" \
                "$METHOD" \
                "$CHR_NUM" \
                "$CHR" \
                "$CHR_STATS" \
                >> "$OUT"

            echo "  Chromosome ${CHR_NUM}: ${CHR_LINES} depth values"

        fi

    done


    # ========================================================
    # Finished sample
    # ========================================================

    echo
    echo "--------------------------------------------------------"
    echo "Finished sample: ${SAMPLE}"
    echo "Output: ${OUT}"
    echo "--------------------------------------------------------"

    # Show first few output lines
    echo
    echo "Output preview:"
    head -n 5 "$OUT"

done


echo
echo "========================================================"
echo "All samples finished."
echo "========================================================"

echo
echo "Results are in:"
echo "${OUTDIR}"

echo
echo "Files created:"
ls -lh "${OUTDIR}"

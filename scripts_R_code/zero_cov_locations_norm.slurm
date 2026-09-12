#!/bin/bash
#SBATCH --job-name=zero_depth_all
#SBATCH --output=zero_depth_all_norm%j.out
#SBATCH --error=zero_depth_all_norm%j.err
#SBATCH --account=project_2019524
#SBATCH --partition=longrun
#SBATCH --qos=normal
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=1-00:00:00

METHOD="normalinput"

DEPTHDIR="/scratch/project_2019524/pinksalmon_alignments/${METHOD}/coverage_depth"

FAI_FILE="/scratch/project_2019524/pinksalmon_reference/GCF_021184085.1_OgorEven_v1.0_genomic.fna.fai"

OUTDIR="${DEPTHDIR}/depth_stats"

mkdir -p "$OUTDIR"

#loop

for SAMPLE in 8167 8168
do

    # Use the new depth file for 8038
    if [ "$SAMPLE" = "8038" ]; then
        DEPTH_FILE="${DEPTHDIR}/8038.depth_new.txt"
    else
        DEPTH_FILE="${DEPTHDIR}/${SAMPLE}.depth.txt"
    fi

    OUT="${OUTDIR}/${SAMPLE}.zero_depth_summary.txt"

    echo "Processing ${SAMPLE}..."
    echo "Depth file: ${DEPTH_FILE}"


    awk -v OFS="\t" '

    BEGIN {

        # Map chromosome names to chromosome numbers
        for (i = 1; i <= 27; i++) {

            chr = sprintf("NC_%06d.1", 60172 + i)

            chr_number[chr] = i

            # Explicitly initialise counts
            total_zero[i] = 0
            first15[i] = 0
            middle70[i] = 0
            last15[i] = 0
        }
    }


    # -----------------------------------------------------
    # Read chromosome lengths from .fai
    # -----------------------------------------------------

    FNR == NR {

        if ($1 in chr_number) {
            chr_length[$1] = $2
        }

        next
    }


    {

        chr = $1
        pos = $2
        depth = $3

        if (depth == 0 && chr in chr_number) {

            n = chr_number[chr]

            chr_len = chr_length[chr]

            total_zero[n]++

            relative_position = pos / chr_len


            if (relative_position <= 0.15) {

                first15[n]++

            }
            else if (relative_position < 0.85) {

                middle70[n]++

            }
            else {

                last15[n]++
            }
        }
    }


    END {

        print "chr_number",
              "chromosome",
              "total_not_covered",
              "first_15_percent",
              "middle_70_percent",
              "last_15_percent"


        for (i = 1; i <= 27; i++) {

            chr = sprintf("NC_%06d.1", 60172 + i)

            total = total_zero[i] + 0


            if (total > 0) {

                first_pct = 100 * (first15[i] + 0) / total

                middle_pct = 100 * (middle70[i] + 0) / total

                last_pct = 100 * (last15[i] + 0) / total

            }
            else {

                first_pct = 0

                middle_pct = 0

                last_pct = 0
            }


            print i,
                  chr,
                  total,
                  first_pct,
                  middle_pct,
                  last_pct
        }
    }

    ' "$FAI_FILE" "$DEPTH_FILE" > "$OUT"


    echo "Finished ${SAMPLE}"
    echo "Output: ${OUT}"
    echo

done


echo "All samples finished."


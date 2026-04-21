#!/usr/bin/env bash

set -euo pipefail

# ----------------------------
# Usage:
# map_and_fix_header.sh <ID> <index_prefix> <reads.fq.gz> <threads> <k_multimapper>
# ----------------------------

ID="$1"
idx="$2"
reads="$3"
threads="$4"
n_allow_multimapper="$5"

# ----------------------------
# Try direct mapping first
# ----------------------------
if bowtie2 -p "${threads}" -k "${n_allow_multimapper}" -x "${idx}" -U "${reads}" \
    --un "${ID}_${idx}_unclass.fq" \
    | samtools view -@ "${threads}" -Sb -q 0 -F 4 - > "${ID}.${idx}.bam"; then

   echo "$(date -u) mapping directly" >&2

else
    echo "$(date -u) fallback: rebuilding header" >&2

    # ----------------------------
    # Generate headers
    # ----------------------------
    bowtie2-inspect "${idx}" > "${idx}.fa"

    samtools dict "${idx}.fa" \
        | LC_ALL=C grep "^@SQ" \
        | cut -f1-3 \
        > "${idx}.headers"

    # ----------------------------
    # Generate SAM without SQ header
    # ----------------------------
    bowtie2 -x "${idx}" -U "${reads}" \
        --threads "${threads}" \
        -k "${n_allow_multimapper}" \
        --no-sq --no-unal \
        -S "${ID}_${idx}_noheader.sam" \
        --un "${ID}_${idx}_unclass.fq"


    echo "$(date -u) 1. Headerless SAM created" >&2

    # ----------------------------
    # Extract reference names (robust)
    # ----------------------------
    awk '!/^@/ && !seen[$3]++ {print $3}' "${ID}_${idx}_noheader.sam" > "${ID}_refnames.txt"

    # ----------------------------
    # Match headers (robust awk)
    # ----------------------------
    awk 'NR==FNR { pat[$1]; next }
    {
        for(i=1;i<=NF;i++){
            if($i ~ /^SN:/){
                split($i,a,":");
                if(a[2] in pat){
                    print;
                    break;
                }
            }
        }
    }' "${ID}_refnames.txt" "${idx}.headers" \
    > "${ID}_header_matches.txt"

    # ----------------------------
    # Rebuild SAM safely
    # ----------------------------
    {
        sed -n '1p' "${ID}_${idx}_noheader.sam"
        cat "${ID}_header_matches.txt"
        sed -n '2,$p' "${ID}_${idx}_noheader.sam"
    } > "${ID}_${idx}_aln.sam"

    echo "$(date -u) 2. Headers added to SAM" >&2

    # ----------------------------
    # Convert to BAM
    # ----------------------------
    samtools view -@ "${threads}" -b "${ID}_${idx}_aln.sam" \
        | samtools sort -n -@ "${threads}" \
        -o "${ID}_${idx}_aln.bam"

    echo "$(date -u) 3. BAM created and sorted" >&2
fi

pigz "${ID}_${idx}_unclass.fq"

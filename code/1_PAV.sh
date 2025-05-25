#!/bin/bash

# ================================================================
# Universal scanPAV Pipeline Script
#
# This script performs pairwise Presence-Absence Variation (PAV)
# analysis between genome assemblies using the scanPAV tool.
#
# For each pair A vs B:
#   - Identifies sequences present in A and absent in B
#   - Then sequences absent in A and present in B

# ---- Parameters (edit if needed) ----
cpus=32
aligner="bwa"
score=550
scanPAV_version=$(realpath "$(pwd)/../scanPAV/src")

# ---- Write here assemblies and names ----

assemblies=(
    "/path/to/assembly1.fa"
    "/path/to/assembly2.fa"
    "/path/to/assembly3.fa"
)

names=(
    sample1
    sample2
    sample3
)

# ---- Create unique working directory ----
timestamp=$(date +%Y%m%d_%H%M%S)
workdir="scanpav_run_$timestamp"
mkdir "$workdir"
cd "$workdir" || exit 1

# ----  Create symbolic links to assemblies in the working directory ----

for i in "${!assemblies[@]}"; do
    ln -s "${assemblies[$i]}" "${names[$i]}.fa"
done

# ---- Log file ----

touch log.txt

# ---- Run all pairwise permutations ----

for i in "${!names[@]}"; do
    for j in "${!names[@]}"; do
        if [[ "$i" != "$j" ]]; then
            a="${names[$i]}"
            b="${names[$j]}"

            echo; echo ">>> Presence PAVs: $a vs $b" | tee -a log.txt; echo
            "$scanPAV_version/scanPAV" -nodes "$cpus" -score "$score" -align "$aligner" \
                "$a.fa" "$b.fa" "pavs_present_in_${a}_vs_${b}" | tee -a log.txt

            echo; echo ">>> Absence PAVs: $b vs $a" | tee -a log.txt; echo
            "$scanPAV_version/scanPAV" -nodes "$cpus" -score "$score" -align "$aligner" \
                "$b.fa" "$a.fa" "pavs_absent_in_${a}_vs_${b}" | tee -a log.txt
        fi
    done
done
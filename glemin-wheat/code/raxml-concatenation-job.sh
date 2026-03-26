# !/bin/bash

set -x

subset="$1"

DATADIR="Concatenation10Mb_subset_${subset}"

mkdir concatenation_results_${subset}

for file in "$DATADIR"/*; do
    raxml-ng --msa $file --model GTR+G4 --threads 8
    mv ${file}**raxml** concatenation_results_${subset}
done
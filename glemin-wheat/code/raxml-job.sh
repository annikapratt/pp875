# !/bin/bash

set -x

subset="$1"

DATADIR="OneCopyGenes_${subset}"

mkdir RAxML_results_${subset}

for file in "$DATADIR"/*; do
    raxml-ng --msa $file --model GTR+G4
    mv ${file}**raxml** RAxML_results_${subset}
done

ls -FR
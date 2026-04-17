#!/bin/bash

set -x

ENVNAME=pp875
CONDATAR='/staging/pratt7/conda_envs/'

cp "$CONDATAR/$ENVNAME.tar.gz" .

export ENVDIR=$ENVNAME

export PATH
mkdir $ENVDIR
tar -xzf $ENVNAME.tar.gz -C $ENVDIR
. $ENVDIR/bin/activate

conda-unpack

DATADIR="10concatenation10Mb_OneCopy-phylip-ch3"
mkdir HyDe/
mkdir HyDe/10Mb-concatenation-ch3

for file in "$DATADIR"/*; do
	echo "Processing file: $file"
    
    # 1. Extract the filename from the path
    filename=$(basename "$file")
    
    # 2. Remove the extension
    base="${filename%.*}"

    run_hyde.py -i "$file" -m 07-species_mapping.txt -o H_vulgare -n 47 -t 17 -s 11354214 --prefix "HyDe/10Mb-concatenation-ch3/$base"
done

ls -lh
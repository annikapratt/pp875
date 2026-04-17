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

run_hyde.py -i 10-triticeae_allindividuals_OneCopyGenes.phylip \
-m 07-species_mapping.txt \
-o H_vulgare -n 47 -t 17 -s 11354214 --prefix 10-hyde-job

ls -lh
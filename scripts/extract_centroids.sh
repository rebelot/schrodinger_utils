#!/usr/bin/env bash

CLUSTLOG="$1"
NCLUST="$2"
CMS="$3"
TRJ="$4"
OUT="$5"

indexes=($(grep -E -o 'medoid is frame\s*\d*' "$CLUSTLOG" | awk '{print $4}' | head -n $NCLUST))

for i in $(seq 1 "$NCLUST"); do
    $SCHRODINGER/run frame2cms.py $CMS $TRJ ${OUT}_$i -i ${indexes[$((i -1 ))]}
done

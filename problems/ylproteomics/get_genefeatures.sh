#!/bin/sh

# Mutated genes are not intersting here, but other features (transcription factors) are
grep -v Mutated_Gene ../ylipo/1_WTIL_WTnoIL/input/gene_features.csv >> /tmp/gene_features.csv

cp  /tmp/gene_features.csv P3_00vsP3_400/input/
cp  /tmp/gene_features.csv P3_00vsWT_00/input/
cp  /tmp/gene_features.csv WT_00vsWT_400/input/


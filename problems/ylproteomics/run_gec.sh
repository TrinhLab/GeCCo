#!/bin/sh

gec P3_00vsP3_400/ -p param.csv -n network/network.csv  --write_scores true
gec P3_00vsWT_00/ -p param.csv -n network/network.csv  --write_scores true
gec WT_00vsWT_400/ -p param.csv -n network/network.csv --write_scores true

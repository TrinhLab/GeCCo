#!/bin/sh

gec P3_00vsP3_400/ -p param.csv -n network/network.csv
gec P3_00vsWT_00/ -p param.csv -n network/network.csv
gec WT_00vsWT_400/ -p param.csv -n network/network.csv

# Without network
#gec P3_00vsP3_400/ -p param.csv
#gec P3_00vsWT_00/ -p param.csv
#gec WT_00vsWT_400/ -p param.csv
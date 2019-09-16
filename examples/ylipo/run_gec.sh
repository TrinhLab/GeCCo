#!/bin/sh

gec 1_WTIL_WTnoIL -n network.csv --write_scores true
gec 2_MTIL_MTnoIL -n network.csv --write_scores true
gec 3_MTnoIL_WTnoIL -n network.csv --write_scores true
gec 4_WTIL_MTIL -n network.csv --write_scores true


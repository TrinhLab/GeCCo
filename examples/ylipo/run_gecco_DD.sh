#!/bin/sh

gecco 1_WTIL_WTnoIL_DD -n network.csv --correlation_cutoff 0.95 --write_scores true
gecco 2_MTIL_MTnoIL_DD -n network.csv --correlation_cutoff 0.95 --write_scores true
gecco 3_MTnoIL_WTnoIL_DD -n network.csv --correlation_cutoff 0.95 --write_scores true
gecco 4_WTIL_MTIL_DD -n network.csv --correlation_cutoff 0.95 --write_scores true

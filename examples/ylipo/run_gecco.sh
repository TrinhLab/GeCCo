#!/bin/sh

gecco 1_WTIL_WTnoIL -n network.csv --write_scores true
gecco 2_MTIL_MTnoIL -n network.csv --write_scores true
gecco 3_MTnoIL_WTnoIL -n network.csv --write_scores true
gecco 4_WTIL_MTIL -n network.csv --write_scores true


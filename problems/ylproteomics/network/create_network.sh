#!/bin/sh

# Creates co-expression network, requires coexpy: https://github.com/TrinhLab/coexpy

coexpy ../input_for_proteomics2.csv ./ --is_log_transformed true --param_file_path ./nparam.csv

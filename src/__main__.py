"""
Gene Expression Classifier (GEC)

Command line interface
"""


import argparse
from src.classifier import find_classes
import src.settings as settings
import os
import pandas as pd
import sys

# Default parameters
d_parameters_path = settings.DEFAULT_PARAMETERS
d_write_scores = 'false'
d_verbose_level = 1


# Command line interface
def main(args_=None):
    if args_ is None:
        args_ = sys.argv[1:]

    # Parse input
    parser = argparse.ArgumentParser(description='Gene Expression Classifier (GEC)')
    parser.add_argument('problem_dir', help='Path to problem directory')
    parser.add_argument('-p','--parameters_path', help='Path to .csv parameter file', required=False, default=d_parameters_path)
    parser.add_argument('--write_scores', choices=['false','true'], help='Writes folds changes and Z-scores to a file names scores.csv in the problem directory', default=d_write_scores)
    parser.add_argument('-v','--verbose_level', choices=[0,1,2], default=d_verbose_level,
                    help='Verbose options: 0, remove all output; 1, basic output (default); 2, most descriptive output')
    args = parser.parse_args(args_)
    core(**vars(args))


# core method
def core(problem_dir, parameters_path=d_parameters_path, write_scores=d_write_scores,verbose_level=d_verbose_level):
    # Default paths
    input_dir = os.path.join(problem_dir,'input')

    # Apply gec
    if verbose_level >= 1:
        print('Classifying genes...', end='')
    tpm_file_path = os.path.join(input_dir, 'tpm.csv')
    classified_df, cdf = find_classes(tpm_file_path, param_file_path=parameters_path)
    if verbose_level >= 1:
        print('Done')

    # Include gene features if applicable
    gf_path = os.path.join(input_dir, 'gene_features.csv')
    if os.path.exists(gf_path):
        gf = pd.read_csv(gf_path)
        gf = gf.set_index('Gene')
        classified_df = classified_df.join(gf)
        if verbose_level >= 1:
            print('Gene features included')
    else:
        if verbose_level >= 1:
            print('Gene feature file not found')

    # Create cooexpresion network if applicable
    if verbose_level >=1:
        print('Coexpression network not implemented yet')

    # Write output
    node_output_path = os.path.abspath(os.path.join(problem_dir, 'node_attributes.csv'))
    classified_df.to_csv(node_output_path)
    if verbose_level >=1:
        print('Output files written to:')
        print('\t', node_output_path)

    if write_scores == 'true':
        scores_output_path = os.path.abspath(os.path.join(problem_dir, 'scores.csv'))
        cdf.to_csv(scores_output_path)
        if verbose_level >=1:
            print('\t', scores_output_path)


if __name__ == "__main__":
    main()
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


def main(args_=None):
    if args_ is None:
        args_ = sys.argv[1:]

    # Parse input
    parser = argparse.ArgumentParser(description='Gene Expression Classifier (GEC)')
    parser.add_argument('problem_dir', help='Path to problem directory')
    parser.add_argument('-p','--parameters_path', help='Path to .csv parameter file', required=False)
    parser.add_argument('--write_scores', choices=['false','true'], help='Writes folds changes and Z-scores to a file names scores.csv in the problem directory', default='false')
    parser.add_argument('-v','--verbose_level', choices=[0,1,2], default=1,
                    help='Verbose options: 0, remove all output; 1, basic output (default); 2, most descriptive output')
    args = parser.parse_args(args_)

        # Default paths
    if args.parameters_path is None:
        args.parameters_path = settings.DEFAULT_PARAMETERS
    input_dir = os.path.join(args.problem_dir,'input')

    # Apply gec
    if args.verbose_level >= 1:
        print('Classifying genes...', end='')
    tpm_file_path = os.path.join(input_dir, 'tpm.csv')
    classified_df, cdf = find_classes(tpm_file_path, param_file_path=args.parameters_path)
    if args.verbose_level >= 1:
        print('Done')

    # Include gene features if applicable
    gf_path = os.path.join(input_dir, 'gene_features.csv')
    if os.path.exists(gf_path):
        gf = pd.read_csv(gf_path)
        gf = gf.set_index('Gene')
        classified_df = classified_df.join(gf)
        if args.verbose_level >= 1:
            print('Gene features included')
    else:
        if args.verbose_level >= 1:
            print('Gene feature file not found')

    # Create cooexpresion network if applicable
    if args.verbose_level >=1:
        print('Coexpression network not implemented yet')

    # Write output
    node_output_path = os.path.abspath(os.path.join(args.problem_dir, 'node_attributes.csv'))
    classified_df.to_csv(node_output_path)
    if args.verbose_level >=1:
        print('Output files written to:')
        print('\t',node_output_path)

    if args.write_scores == 'true':
        scores_output_path = os.path.abspath(os.path.join(args.problem_dir, 'scores.csv'))
        cdf.to_csv(scores_output_path)
        if args.verbose_level >=1:
            print('\t',scores_output_path)

if __name__ == "__main__":
    main()
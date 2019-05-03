"""
Gene Expression Classifier (GEC)

Command line interface
"""


import argparse
from gec.classifier import find_classes
import gec.settings as settings
from gec.coexpression import *
import os
import sys
import pandas as pd
import networkx as nx

# Default parameters
d_parameters_path = settings.DEFAULT_PARAMETERS
d_write_scores = 'false'
d_verbose_level = 1
d_network_file_path = None
d_calc_centrality = 'true'

# Command line interface
def cli(args_=None):
    if args_ is None:
        args_ = sys.argv[1:]

    # Parse input
    parser = argparse.ArgumentParser(description='Gene Expression Classifier (GEC)')
    parser.add_argument('problem_dir', help='Path to problem directory')
    parser.add_argument('-p','--parameters_path', help='Path to .csv parameter file', required=False, default=d_parameters_path)
    parser.add_argument('--write_scores', choices=['false','true'], help='Writes folds changes and Z-scores to a file names scores.csv in the problem directory', default=d_write_scores)
    parser.add_argument('-v','--verbose_level', choices=[0,1,2], default=d_verbose_level,
                    help='Verbose options: 0, remove all output; 1, basic output (default); 2, most descriptive output')
    parser.add_argument('-n', '--network_file_path', help='If provided, coexpression analysis is performed.', default=d_network_file_path)
    parser.add_argument('-c','--calc_centrality', choices=['false','true'], help='Determines if centrality metrics are calcualted, can be slow for large graphs', default=d_calc_centrality)
    args = parser.parse_args(args_)
    core(**vars(args))


# core method
def core(problem_dir, parameters_path=d_parameters_path, write_scores=d_write_scores,network_file_path=d_network_file_path, calc_centrality=d_calc_centrality, verbose_level=d_verbose_level):
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
        rename_classes(classified_df)
        if verbose_level >= 1:
            print('Gene features included')
    else:
        if verbose_level >= 1:
            print('Gene feature file not found')

    node_output_path = os.path.abspath(os.path.join(problem_dir, 'gene_attributes.csv'))

    # Coexpresion network analysis if applicable
    if not network_file_path:
        if verbose_level > 0:
            print('Coexpression analysis not performed, only gene classification will be reported')

        # Write output
        classified_df.to_csv(node_output_path)
        if verbose_level > 0:
            print('Output files written to:')
            print('\t', node_output_path)

    else:
        if verbose_level > 0:
            print('Coexpression network analysis in progress...', end='')

        classified_subgraph = get_classified_subgraph(network_file_path, classified_df)
        if calc_centrality.lower() == 'true':
            add_centrality_metrics(classified_subgraph)

        # Add formatting features to node attribute file for simpler cytoscape usage
        add_node_visualization_attributes(classified_subgraph)

        if verbose_level > 0:
            print('Done')

        # Write gene attribute file ( # Also contains useful info for cytoscape visualization)
        write_node_attributes(classified_subgraph, node_output_path)
        if verbose_level > 0:
            print('Output files written to:')
            print('\t', node_output_path)

        # Write two edge lists for "positive" and "negative" networks
        pos_g, neg_g = split_g(classified_subgraph)
        pos_output_path = os.path.abspath(os.path.join(problem_dir, 'positive_network.csv'))
        neg_output_path = os.path.abspath(os.path.join(problem_dir, 'negative_network.csv'))

        nx.write_edgelist(pos_g, path=pos_output_path, delimiter=',', data=False)
        nx.write_edgelist(neg_g, path=neg_output_path, delimiter=',', data=False)
        if verbose_level > 0:
            print('Output files written to:')
            print('\t', pos_output_path)
            print('\t', neg_output_path)

    # Write  scores
    if write_scores == 'true':
        scores_output_path = os.path.abspath(os.path.join(problem_dir, 'scores.csv'))
        cdf.to_csv(scores_output_path)
        if verbose_level >=1:
            print('\t', scores_output_path)


def rename_classes(classified_df):
    """ More precise nomenclature that is easier to understand
    """
    old2new = {
             "highly_expressed"   : "control_overexpressed",
             "lowly_expressed"    : "case_overexpressed",
             "upregulated"        : "control_upregulated",
             "downregulated"      : "case_upregulated",
             "changed_regulation" : "changed_regulation",
             "no_change"          : "no_change"
             }
    classified_df["Class"].replace(old2new, inplace=True)


if __name__ == "__main__":
    cli()

"""
Gene Expression Classifier (GeCCo)

Command line interface
"""


import argparse
from gecco.classifier import find_classes, to_numeric
import gecco.settings as settings
from gecco.coexpression import *
import os
import sys
import pandas as pd
import networkx as nx

# Default parameters
d_parameters_path = settings.DEFAULT_PARAMETERS
d_write_scores = False
d_verbose_level = 1
d_network_file_path = None
d_calc_centrality = True
d_include_no_change = False

def cli(args_=None):
    """
    Command line interface
    """
    if args_ is None:
        args_ = sys.argv[1:]

    # Parse input
    parser = argparse.ArgumentParser(description='Gene Expression Classifier (GeCCo)', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('problem_dir', help='Path to problem directory')
    parser.add_argument('--write_scores', help='Writes folds changes and Z-scores to a file names scores.csv in the problem directory', default=d_write_scores, choices=[True, False], type=lambda x: (str(x).lower() == 'true'))
    parser.add_argument('-v','--verbose_level', choices=[0,1,2], default=d_verbose_level,
                    help='Verbose options: 0, remove all output; 1, basic output (default); 2, most descriptive output')
    parser.add_argument('-n', '--network_file_path', help='If provided, coexpression analysis is performed.', default=d_network_file_path)
    parser.add_argument('-c','--calc_centrality', help='Determines if centrality metrics are calcualted, can be slow for large graphs', default=d_calc_centrality, choices=[True, False], type=lambda x: (str(x).lower() == 'true'))
    parser.add_argument('--include_no_change',  help='Includes genes in the "no_change" category as part of the output' , default=d_include_no_change, choices=[True, False], type=lambda x: (str(x).lower() == 'true'))

    parser.add_argument('-p','--parameters_path', help='Path to .csv parameter file. LEGACY: Paramter input through command line is recommended instead of parameter file. However, if a parameter file is provided this will take precedence over command line versions', required=False, default=d_parameters_path)
    # Command line paramters also present in parameter file:
    parser.add_argument('--z_cutoff',  help='wip' , default=1.5, type=float)
    parser.add_argument('--fc_cutoff',  help='wip' , default=1, type=float)
    parser.add_argument('--correlation_type', choices=["pearson","spearman"], default="pearson", help='wip', type=str)
    parser.add_argument('--correlation_cutoff', help='wip' , default=0.85, type=float)
    parser.add_argument('--seudovariance', help='wip' , default=0.25, type=float)
    parser.add_argument('--floor_and_logtransform',help='wip' , default=True, choices=[True, False], type=lambda x: (str(x).lower() == 'true'))
    parser.add_argument('-mtp', '--min_tpm',help='Used for flooring (any value below this value will be set to this value)' , default=5, type=float)

    args = parser.parse_args(args_)
    core(args.problem_dir, parameters_path=args.parameters_path, write_scores=args.write_scores,network_file_path=args.network_file_path, calc_centrality=args.calc_centrality, verbose_level=args.verbose_level, include_no_change=args.include_no_change, cli_args=args)

def load_params(parameters_path, cli_args):
    """ Load parameters. Non-default parameter file superseed CLI parameters.
    """
    if parameters_path != d_parameters_path:
        prm = pd.read_csv(param_file_path)
        param_dict = dict(zip(prm['param_id'], prm['value'].apply(to_numeric)))
        # Compatibility: (convert strings to logical)
        param_dict['floor_and_logtransform'] = param_dict['floor_and_logtransform'].lower() == 'true'
    else:
        param_dict = vars(cli_args) # Some are not parameters but it does not hurt
    return param_dict

# core method
def core(problem_dir, parameters_path=d_parameters_path, write_scores=d_write_scores,network_file_path=d_network_file_path, calc_centrality=d_calc_centrality, verbose_level=d_verbose_level, include_no_change=d_include_no_change, cli_args=None):
    # Default paths
    input_dir = os.path.join(problem_dir,'input')

    param_dict = load_params(parameters_path, cli_args)
    # Apply gecco
    if verbose_level >= 1:
        print('Classifying genes...', end='')
    tpm_file_path = os.path.join(input_dir, 'tpm.csv')
    classified_df, cdf = find_classes(tpm_file_path, param_dict=param_dict, remove_no_change=(not include_no_change))
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
        if calc_centrality == True:
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
    if write_scores:
        scores_output_path = os.path.abspath(os.path.join(problem_dir, 'scores.csv'))
        out_scores_df = cdf.join(classified_df)
        out_scores_df["Class"] = out_scores_df["Class"].fillna("no_change") # Just for clarity
        out_scores_df.to_csv(scores_output_path)
        if verbose_level >=1:
            print('\t', scores_output_path)


def rename_classes(classified_df):
    """ More precise nomenclature that is easier to understand
    """
    old2new = {
             "highly_expressed"   : "case_overexpressed",
             "lowly_expressed"    : "control_overexpressed",
             "upregulated"        : "case_upregulated",
             "downregulated"      : "control_upregulated",
             "changed_regulation" : "changed_regulation",
             "no_change"          : "no_change"
             }
    classified_df["Class"].replace(old2new, inplace=True)


if __name__ == "__main__":
    cli()

"""
Methods to connects gene classification with co-expression network.
"""

import networkx as nx
import gec.settings
import pandas as pd
import os

def get_classified_subgraph(network_file_path, classified_df):
    g = nx.read_weighted_edgelist(network_file_path, delimiter=',', nodetype=str)
    gs = nx.subgraph(g, list(classified_df.index))  # Only classified genes appear in classified_df
    classes = dict(zip(classified_df.index, classified_df['Class']))
    features = dict(zip(classified_df.index, classified_df['Feature']))
    nx.set_node_attributes(gs, classes, 'class')
    nx.set_node_attributes(gs, features, 'feature')
    return gs


def add_centrality_metrics(g):
    bc = nx.betweenness_centrality(g)
    ec = nx.eigenvector_centrality(g)
    dc = nx.degree_centrality(g)
    nx.set_node_attributes(g, dc, 'centrality_degree')
    nx.set_node_attributes(g, bc, 'centrality_betweenness')
    nx.set_node_attributes(g, ec, 'centrality_eigenvector_weighted')


def add_node_visualization_attributes(g, param_dir=gec.settings.DEFAULT_V_PARAM_DIR):
    """
    Adds node attributes that can be used with cytoscape "passthrough" mapping to quickly generate visually appealing networks.

    Parameters
    ----------
    g : networkx.Graph
         Graph object to be modified
    class_color_path : str
    feature_shape_path : str
    """

    df = pd.read_csv(os.path.join(param_dir, 'class_colors.csv'))
    class_color_d = dict(zip(df['class'], df['color']))
    df = pd.read_csv(os.path.join(param_dir, 'feature_shapes.csv'))
    shape_d = dict(zip(df['feature_id'], df['shape']))
    df = pd.read_csv(os.path.join(param_dir, 'border.csv'))
    border_d = dict(zip(df['type'], df['value']))
    df = pd.read_csv(os.path.join(param_dir, 'size.csv'))
    size_d = dict(zip(df['type'], df['value']))

    available_features = list(shape_d.keys())

    for node in g.nodes(data=True):
        node_id = node[0]
        data = dict(node[1])
        # Color
        g.node[node_id]['color'] = class_color_d[data['class']]
        if not (data['feature'] in available_features):
            g.node[node_id]['shape'] = shape_d['otherwise']
            g.node[node_id]['size'] = size_d['otherwise']
            g.node[node_id]['border'] = border_d['otherwise']
        else:
            g.node[node_id]['shape'] = shape_d[data['feature']]
            g.node[node_id]['size'] = size_d['has_feature']
            g.node[node_id]['border'] = border_d['has_feature']


def split_g(g):
    """
    Returns two networks:
        ‘positive network’ (= upregulated & highly expressed)
        ‘negative network' (= downregulated & lowly expressed)

    Parameters
    ----------
    g : networkx.Graph
    """

    pos_nodes = [x for x,y in g.nodes(data=True) if y['class'] in ['upregulated', 'highly_expressed']]
    neg_nodes = [x for x,y in g.nodes(data=True) if y['class'] in ['downregulated', 'lowly_expressed']]
    p_g = nx.subgraph(g, pos_nodes)
    n_g = nx.subgraph(g, neg_nodes)
    return p_g, n_g


def write_node_attributes(g, output_file_path):
    df = pd.DataFrame(dict(g.nodes.data())).transpose()
    df.index.name = 'Gene'
    df.to_csv(output_file_path)

import pandas as pd
import numpy as np
import re

# Load parameters

def find_classes(tpm_file_path, param_dict, remove_no_change=True):
    """ Assign classes to each gene based on decision tree.
    :param tpm_file_path: Path to csv file formated for gec.
    :param param_dict: Paramter dictionary
    :param remove_no_change: If true removes genes classified as no_change from the output (Default is True).
    :return: classified_df, cdf. Where classified_df contains gene ids and their associated categories and cdf contains the classificaiton variables (z-scores, and fold changes)
    """

    tpm, sample_groups = parse_input(tpm_file_path)
    cdf = calc_classification_variables(tpm, sample_groups, param_dict)
    classified_df = classify(cdf, param_dict)
    if remove_no_change:
        classified_df = classified_df[classified_df['Class'].str.contains("no_change") == False]
    return classified_df, cdf

def to_numeric(x):
    """ Convert numeric values to float type and leave strings alone"""
    try:
        return float(x)
    except:
        return x

def parse_input(tpm_file_path):
    """ loads tpm data table and identifiers replicates
    """
    tpm = pd.read_csv(tpm_file_path)
    check_format(tpm)
    tpm = tpm.set_index('Gene')
    # Identify replicate groups
    no_rep_samples = list(set([re.sub('_rep\d$','',header) for header in  tpm.columns]))
    sample_groups = {}
    for sample_id in no_rep_samples:
        sample_groups[sample_id] = [h for h in tpm.columns if sample_id in h]
    return tpm, sample_groups


def check_format(tpm):
    """ Ensure the input format is correct
    WIP:
    - Check column labels
    - Check that if data is log2
    - ... """
    assert tpm.columns[0] == 'Gene', 'First column of tpm table must correspond to genes and be labeled \"Gene\"'
    return None


def calc_classification_variables(tpm, sample_groups, param_dict):
    """ Calculate z-scores and fold changes, output in a new table"""
    seudovar = param_dict['seudovariance']
    if param_dict['floor_and_logtransform']:
        min_tpm = param_dict['min_tpm']
        # Replace tpm  by minimum relevant value and convert tpm to log2
        tpm[tpm < min_tpm] = min_tpm
        tpm = np.log2(tpm)

    # Calc means and variance/n
    calcdf = pd.DataFrame(index=tpm.index)
    for k, v in sample_groups.items():
        calcdf[k] = tpm[v].mean(axis=1)
        calcdf[k+'_vn'] = tpm[v].var(axis=1)/len(v) # variance divided by the number of samples
    calcdf = calcdf.fillna(0) # e.g. variance for one sample = NaN

    calcdf['varsum_Z'] = calcdf.MT_t1_vn + calcdf.MT_t2_vn + calcdf.WT_t1_vn + calcdf.WT_t2_vn
    calcdf['varsum_fct1'] = calcdf.MT_t1_vn + calcdf.WT_t1_vn
    calcdf['varsum_fct2'] = calcdf.MT_t2_vn + calcdf.WT_t2_vn
    def calc_Z(r):
        diff = (r.MT_t2 - r.MT_t1) - (r.WT_t2 - r.WT_t1)
        return diff/np.sqrt(seudovar + r.varsum_Z)
    def calc_FC_t1(r):
        diff = r.MT_t1 - r.WT_t1
        return diff/np.sqrt(seudovar + r.varsum_fct1)
    def calc_FC_t2(r):
        diff = r.MT_t2 - r.WT_t2
        return diff/np.sqrt(seudovar + r.varsum_fct2)

    cdf = pd.DataFrame(index=tpm.index, columns=['FC_t1', 'FC_t2', 'Z'])
    cdf['Z'] = calcdf.apply(calc_Z, axis=1)
    cdf['FC_t1'] = calcdf.apply(calc_FC_t1, axis=1)
    cdf['FC_t2'] = calcdf.apply(calc_FC_t2, axis=1)
    return cdf


def classify(cdf, param_dict):
    classified_df = pd.DataFrame(index=cdf.index, columns=['Class'])

    z_cutoff = float(param_dict['z_cutoff'])
    fc_cutoff = float(param_dict['fc_cutoff'])

    def class_map(r):
        if (r.FC_t1 >= fc_cutoff) and (r.FC_t2 >= fc_cutoff):
            return "highly_expressed" #aka case_overexpressed
        elif (-fc_cutoff < r.FC_t1) and (r.FC_t1 < fc_cutoff) and (r.FC_t2 >= fc_cutoff) and (r.Z >= z_cutoff):
            return "upregulated" #aka case_upregulated
        elif (r.FC_t1 <= -fc_cutoff) and (r.FC_t2 <= -fc_cutoff):
            return "lowly_expressed" #aka control_overexpressed
        elif (-fc_cutoff < r.FC_t1) and (r.FC_t1 < fc_cutoff) and (r.FC_t2 <= -fc_cutoff) and (r.Z <= -z_cutoff):
            return "downregulated" #aka control_upregulated
        elif abs(r.Z) > z_cutoff:
            return "changed_regulation"
        else:
            return "no_change"

    classified_df['Class'] = cdf.apply(class_map, axis=1)

    return classified_df



#!/usr/bin/env python3

# Create inputs for the following cases:
#   1.   P3_00 (up) vs. WT_00 (down)
#   2.   WT_00 (up) vs. WT_400 (down)
#   2.   P3_00 (up) vs. P3_400 (down)
#12h = time point 1
#16h = time point 2
#Samples:  Strain_[thiamine]_time_replicate

# GEC naming convetion: `Gene|WT_t1_rep1|WT_t1_rep2|WT_t1_rep3|WT_t2_rep1|WT_t2_rep2|MT_t1_rep1|MT_t1_rep2|MT_t2_rep1|MT_t2_rep2` Any number of replicates is acceptable. Note that WT corresponds to _case_ and MT to _control_.

import pandas as pd
import re
import os

df = pd.read_csv('./input_for_proteomics2.csv')

def create_input(wt, mut):
    df1 = df.filter(regex='Gene|{}|{}'.format(wt,mut))
    df1 = df1.rename(columns=lambda x: re.sub('{}_12_0([0-9])'.format(wt),r'WT_t1_rep\1',x))
    df1 = df1.rename(columns=lambda x: re.sub('{}_16_0([0-9])'.format(wt),r'WT_t2_rep\1',x))
    df1 = df1.rename(columns=lambda x: re.sub('{}_12_0([0-9])'.format(mut),r'MT_t1_rep\1',x))
    df1 = df1.rename(columns=lambda x: re.sub('{}_16_0([0-9])'.format(mut),r'MT_t2_rep\1',x))
    dirname =  '{}vs{}/input'.format(wt,mut)
    if not os.path.exists(dirname):
        os.makedirs(dirname)
    df1.to_csv('./{}/tpm.csv'.format(dirname), index=False)

create_input('P3_00', 'WT_00')
create_input('WT_00', 'WT_400')
create_input('P3_00', 'P3_400')

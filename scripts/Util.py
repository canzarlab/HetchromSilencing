#!/usr/bin/env python

import logging
from math import log
import numpy as np
import os
import pandas as pd

# pd.options.mode.chain_assignment = None

logger = logging.getLogger(__name__)


header = ['gene-id', 'gene-name', 'length', 'type', 'category', 'bio_type']


def to_log2_tpm(tpm_df, cols=None, shift=1):

    h = list(set(tpm_df.columns).intersection({'gene-id', 'gene-name', 'type', 'category', 'bio_type'}))
    df = tpm_df[h]

    if cols:
        columns = cols
    else:
        columns = tpm_df.columns

    zero_cols = 0

    for sample in columns:

        if sample not in header:

            new_df, zero_col = to_log2_col(tpm_df[['gene-id', sample]], sample, shift=shift)
            df = pd.merge(df, new_df, on='gene-id', how='left')
            df = df.drop_duplicates()
            zero_cols += zero_col

    if zero_cols:
        print('%d columns contained zero values. Their log-transformed results are NaNs' % zero_cols)

    return df


def to_log2_col(df, col, shift=1):

    zero_shift = 0

    if shift == 0:
        new_df = df[df[col] != 0]
        if len(new_df) < len(df):
            zero_shift = 1
    else:
        new_df = df.copy()

    new_df[col] = new_df[col].apply(lambda x: log((x + shift), 2))

    return new_df, zero_shift


def zero_adjust(x, shift=0.1):
    if x[0] == 0.0 or x[1] == 0.0:
        x[0] += shift
        x[1] += shift
    return x[0], x[1]


def zero_adjust_df(df, pairs):

    new_header = [x for x in header if x in list(df.columns)]
    za_df = df[new_header]

    for pair in pairs:
        t_s = df[list(pair)].apply(zero_adjust, axis=1)
        t_d = t_s.apply(pd.Series)
        t_d.columns = list(pair)
        za_df = pd.concat([za_df, t_d], axis=1)

    return za_df


def ratio_df(df, pairs, new_cols=None):

    new_header = [x for x in header if x in list(df.columns)]

    if new_cols:
        assert len(pairs) == len(new_cols)
    else:
        new_cols = [x + ',' + y for x, y in pairs]

    for index, pair in enumerate(pairs):

        x, y = pair
        name = new_cols[index]
        df[name] = df[x] / df[y]

    df = df.replace([np.inf, -np.inf], np.nan)
    gxp_df = df[new_header + new_cols]

    return gxp_df


def ratio_table(df, pairs, new_cols=None, out_dir='.', file_name=None):

    if not file_name:
        file_name = 'te_table.csv'

    out_path = os.path.join(out_dir, file_name)
    gxp_df = ratio_df(df, pairs, new_cols=new_cols)
    gxp_df.to_csv(out_path, sep='\t', index=None)

    return gxp_df

if __name__ == "__main__":

    logging.basicConfig(level=logging.INFO)
    logger.info('Call to RNAdeg utility module.')
    logger.info('')
    logger.info('-' * 40)






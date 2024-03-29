#!/usr/bin/env python


import os
import logging
import argparse
import itertools

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats.stats import pearsonr

import Util
pd.options.mode.chained_assignment = None


logger = logging.getLogger(__name__)

header = Util.header


def get_args():

    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--in_file', '-f', type=str, help='Replicates gene expression data file.', required=True)
    parser.add_argument(
        '--sample_list', '-l', nargs='+', help='Samples list (None if single sample is passed)')
    parser.add_argument(
        '--out_dir', '-o', type=str, default='', help='Output folder (default: source directory).')
    parser.add_argument('--prefix', '-x', type=str, default='', help='Name prefix for result files. (default:'')')

    args = parser.parse_args()

    return args


def ma_plot(X, Y, sample1, sample2, out_file=None, out_dir='.', prefix=''):

    if not out_file:
        out_file = prefix + 'MAPlot_' + sample1 + '_' + sample2 + '.png'

    XY = [(i, j) for i, j in zip(X, Y)]
    M = [i - j for i, j in XY]
    A = [0.5 * (i + j) for i, j in XY]

    fig, ax = plt.subplots(figsize=(15, 10))

    plt.scatter(A, M, alpha=0.4)

    plt.xlabel('A-Mean', fontsize=15)
    plt.ylabel('M-Mean', fontsize=15)

    plt.rc('grid', linestyle="--", color='lavender')
    ax.grid(True)

    plt.suptitle('MA-plot : ' + sample1 + ' - ' + sample2, fontsize=20, y=0.92)
    plt.savefig(os.path.join(out_dir, out_file))


def correlation_plot(X, Y, sample1, sample2, names=None, alpha=1, out_file=None, out_dir='.', prefix=''):

    if not out_file:
        out_file = prefix + 'CorrPlot_' + sample1 + '_' + sample2 + '.png'

    fig, ax = plt.subplots(figsize=(15, 10))

    if names:

        for a, b, c in list(zip(X, Y, names)):

            ax.annotate(c, (a, b))
            plt.scatter(a, b, c='red', alpha=alpha)
    else:
        plt.scatter(X, Y, c = 'red', alpha=alpha)

    plt.xlabel(sample1, fontsize=15)
    plt.ylabel(sample2, fontsize=15)

    plt.rc('grid', linestyle="--", color='lavender')
    ax.grid(True)

    ax.plot(ax.get_xlim(), ax.get_ylim(), ls="--", c='gray')

    plt.suptitle(prefix + 'Correlation-plot : ' + sample1 + ' - ' + sample2, fontsize=20, y=0.92)
    plt.savefig(os.path.join(out_dir, out_file))
    plt.close()


def run_ma_checks(df, samples=None, out_dir='.', prefix=''):

    if not samples:
        samples = [s for s in list(df.columns) if s not in header]

    sample_comb = [(i, j) for i, j in itertools.combinations(samples, 2)]
    for item in sample_comb:
        x, y = item
        X = list(df[x])
        Y = list(df[y])

        # Prepare and save MA plot.
        ma_plot(X, Y, x, y, out_dir=out_dir, prefix=prefix)

    return


def run_corr_checks(df, samples=None, out_dir='.', prefix=''):

    if not samples:
        samples = [s for s in list(df.columns) if s not in header]

    corr_list = []
    corr_df = pd.DataFrame()

    sample_comb = [(i, j) for i, j in itertools.combinations(samples, 2)]
    for item in sample_comb:
        x, y = item
        X = list(df[x])
        Y = list(df[y])

        # Calculate correlation.
        corr_list.append((x, y) + pearsonr(X, Y))

        # Prepare and save Correlation plot.
        correlation_plot(X, Y, x, y, out_dir=out_dir)

        corr_df = pd.DataFrame(corr_list, columns=['Sample1', 'Sample2', 'Correlation', 'P-value'])

    corr_df.to_csv(os.path.join(out_dir, prefix + 'correlations.csv'), sep='\t', index=None)

    return


def labeled_corr_plots(df, samples=None, out_dir='.', prefix=''):

    if not samples:
        samples = [s for s in list(df.columns) if s not in header]

    sample_comb = [(i, j) for i, j in itertools.combinations(samples, 2)]
    genes = list(df['gene-name'])

    for item in sample_comb:
        x, y = item
        X = list(df[x])
        Y = list(df[y])

        # Prepare and save Correlation plot.
        correlation_plot(X, Y, x, y, names=genes, out_dir=out_dir, prefix=prefix)

    return


def report_corr(df, samples):

    corr_list = []
    corr_df = pd.DataFrame()

    sample_comb = [(i, j) for i, j in itertools.combinations(samples, 2)]
    sample_comb_reverse = [(j, i) for i, j in sample_comb]
    pairs = sample_comb + sample_comb_reverse
    pairs_df = pd.DataFrame(pairs, columns = ['Sample1', 'Sample2'])

    report = df.merge(pairs_df, on=['Sample1', 'Sample2'])

    return report


def repli_merge(df, sub_samples, new_cols=None, out_dir='.', out_file=None):

    new_header = [x for x in header if x in list(df.columns)]
    new_df = df[new_header]

    if new_cols:
        assert len(sub_samples) == len(new_cols)
    else:
        new_cols = [','.join(x) for x in sub_samples]

    for index, sub_sample in enumerate(sub_samples):

        sub_mtx = df[sub_sample]
        new_val = pd.Series(np.mean(sub_mtx, axis=1))

        sample_name = new_cols[index]
        new_df[sample_name] = new_val

    if not out_file:
        out_file = 'merged.csv'
    out_path = os.path.join(out_dir, out_file)

    new_df.to_csv(out_path, sep='\t', index=None)

    return new_df


if __name__ == "__main__":

    logging.basicConfig(level=logging.INFO)
    logger.info('Call to Replicates Check module.')
    logger.info('')
    logger.info('This module runs quality checks on replicates.')
    logger.info('-' * 40)

    params = vars(get_args())

    param_in_file = params['in_file']
    param_samples = params['sample_list']
    param_out_dir = params['out_dir']
    param_prefix = params['prefix']

    if param_out_dir:
        out_dir = param_out_dir
    else:
        out_dir = '.'

    # Load the gene expression table into a data frame.
    g_df = pd.read_csv(param_in_file, sep='\t')
    g_df = Util.to_log2_tpm(g_df)

    samples = [i for i in list(g_df.columns) if i not in header]

    if param_samples:
        sub_samples = []
        for sample in param_samples:
            sub_samples.append([i for i in samples if sample in i])
    else:
        sub_samples = [samples]

    for sub_sample in sub_samples:

        print(sub_sample)

        # Run the checks.
        run_ma_checks(g_df, sub_sample, out_dir=out_dir, prefix=param_prefix)
        run_corr_checks(g_df, sub_sample, out_dir=out_dir, prefix=param_prefix)
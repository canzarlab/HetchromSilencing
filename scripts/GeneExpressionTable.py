#!/usr/bin/env python

import os
import logging
import argparse

import pandas as pd
import numpy as np
import pysam

import Util

pd.options.mode.chained_assignment = None

header = Util.header


logger = logging.getLogger(__name__)


def get_args():

    parser = argparse.ArgumentParser()

    parser.add_argument(
        '--in_dir', '-d', type=str, help='Directory for alignment files(.bam or .sam)', required=True)
    parser.add_argument(
        '--in_gdf', '-g', type=str, help='Gene information table', required=True)
    parser.add_argument(
        '--out_dir', '-o', type=str, default='', help='Save results in (default: current working directory).')
    parser.add_argument('--prefix', '-x', type=str, default='', help='Name prefix for count tables. (default:'')')
    parser.add_argument('-r', action='store_true')

    args = parser.parse_args()

    return args


def out_file_name(in_file, out_dir='.', ext='', prefix=''):

    dir_name = out_dir
    base_name = os.path.splitext(os.path.basename(in_file))[0]

    if not out_dir:
        dir_name = os.path.dirname(in_file)

    out_file = os.path.join(dir_name, prefix + base_name + ext)

    return out_file


def sample_count_column(in_file, h_df, repeat_ids):

    gene_count = []

    sample = os.path.basename(in_file).split('.')[0].rstrip('Aligned')
    st = pysam.AlignmentFile(in_file, 'rb')

    for read in st.fetch(until_eof=True):

        qn = read.query_name
        if read.has_tag('NH'):
            nh = read.get_tag('NH')
        else:
            nh = 1
        if read.has_tag('GE'):
            gene = read.get_tag('GE')
        else:
            continue

        if (gene in repeat_ids and nh <= 15) or (gene not in repeat_ids and nh == 1):
            # Pre-processing step
            gene_count.append((qn, gene))

    count_column = []

    df = pd.DataFrame(gene_count, columns=['QN', 'gene-id'])
    df = df.drop_duplicates()

    for gene, group in df.groupby(['gene-id']):
        count_column.append((gene, len(group)))

    column_df = pd.DataFrame(count_column, columns=['gene-id', sample])
    column_df = h_df.merge(column_df, on='gene-id', how='left')
    column = list(column_df[sample])

    return sample, column


def domain_count_table(in_dir, in_gdf):

    gene_df = pd.read_csv(in_gdf, sep='\t')
    h_df = gene_df[['gene-id']]
    repeat_ids = list(set(gene_df[gene_df['category'] == 'repeat']['gene-id']))

    new_cols = []
    count_mtx = np.empty([len(h_df), 1])

    for in_file in os.listdir(in_dir):

        if in_file.endswith('.bam') or in_file.endswith('.sam'):

            in_path = os.path.join(in_dir, in_file)
            logger.info('Input bam: %s' % in_path)

            sample, column = sample_count_column(in_path, h_df, repeat_ids)
            new_cols.append(sample)
            count_mtx = np.c_[count_mtx, np.array(column)]

    assert np.delete(count_mtx, 0, 1).shape == (len(h_df), len(new_cols))

    gene_count_df = pd.DataFrame(np.delete(count_mtx, 0, 1), columns=new_cols)
    count_df = pd.concat([h_df, gene_count_df], axis=1)
    count_df = pd.merge(gene_df[header], count_df, on='gene-id')
    count_df = count_df.fillna(0)
    count_df = count_df.drop_duplicates()

    return count_df


def domain_tpm_table(count_df):

    tpm_df = count_df.copy()

    # Divide the read counts by the length of each gene in kilobases. This gives you reads per kilobase (RPK)
    # Count up all the RPK values in a sample and divide this number by 1,000,000. ('per million' scaling factor)
    # Divide the RPK values by the 'per million' scaling factor. This gives you TPM.

    tpm_df['length'] /= 1000

    columns = [item for item in tpm_df.columns if item not in header]

    gene_lengths = list(tpm_df['length'])

    for column in columns:

        rpk = [i / float(j) for i, j in list(zip(list(tpm_df[column]), gene_lengths))]
        per_million = sum(rpk) / 1000000

        tpm = [i / per_million for i in rpk]
        tpm_df[column] = pd.Series(tpm).values

    tpm_df = tpm_df.drop(['length'], axis=1)

    return tpm_df


def generate_tpm_csv(in_dir, in_gdf, out_dir='.', prefix=''):

    count_df = domain_count_table(in_dir, in_gdf)
    count_file = out_file_name(prefix + 'pombe_gene_count_matrix', out_dir, ext='.csv')
    count_df.to_csv(count_file, sep='\t', index=None)

    tpm_df = domain_tpm_table(count_df)
    tpm_file = out_file_name(prefix + 'pombe_tpm_matrix', out_dir, ext='.csv')
    tpm_df.to_csv(tpm_file, sep='\t', index=None)

    return count_file, tpm_file


if __name__ == "__main__":

    params = vars(get_args())

    param_in_dir = params['in_dir']
    param_in_gdf = params['in_gdf']
    param_out_dir = params['out_dir']
    param_prefix = params['prefix']
    param_report = params['r']

    if param_out_dir:
        out_dir = param_out_dir
    else:
        out_dir = param_in_dir

    # logging.basicConfig(level=logging.INFO)
    report_file = os.path.join(out_dir, 'GeneExpressionTable.report.txt')
    logging.basicConfig(filename=report_file, level=logging.INFO)
    logging.getLogger().addHandler(logging.StreamHandler())

    logger.info('Call to GeneExpressionTable module....')
    logger.info('')
    logger.info('This module calculates gene expression in given alignment file(s).')
    logger.info('Input: folder containing .bam file(s)')
    logger.info('Output: raw and tpm gene counts data (.csv) files')
    logger.info('-' * 40)

    count_file, tpm_file = generate_tpm_csv(param_in_dir, param_in_gdf, out_dir, param_prefix)

    logger.info('Gene count matrix saved in: %s' % count_file)
    logger.info('TPM expression matrix saved in: %s' % tpm_file)

    logger.info('Finished successfully!')

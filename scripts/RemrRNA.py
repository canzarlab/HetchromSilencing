#!/usr/bin/env python

import os
import logging
import argparse

import pandas as pd
import pysam


logger = logging.getLogger(__name__)


def get_args():

    parser = argparse.ArgumentParser()

    parser.add_argument(
        '--in_dir', '-d', type=str, help='Directory for alignment files(.bam or .sam).', required=True)
    parser.add_argument(
        '--in_gdf', '-g', type=str, help='Gene annotation table (tab-delimited file).', required=True)
    parser.add_argument(
        '--verbose', '-v', type=bool, default=True, help='Verbose (default: True).')
    parser.add_argument(
        '--out_dir', '-o', type=str, default='.', help='Save results in (default: current working directory).')
    parser.add_argument('-s', action='store_true', help='Save report to file.')

    args = parser.parse_args()

    return args


def out_file_name(in_file, out_dir='.', ext='', prefix=''):

    dir_name = out_dir
    base_name = os.path.splitext(os.path.basename(in_file))[0]

    if not out_dir:
        dir_name = os.path.dirname(in_file)

    out_file = os.path.join(dir_name, prefix + base_name + ext)

    return out_file


def bam_short_report(in_file):

    st = pysam.AlignmentFile(in_file, 'rb')
    total = []
    unmapped = []
    for r in st.fetch(until_eof=True):
        if r.is_unmapped:
            unmapped.append(r.query_name)
        total.append(r.query_name)

    st.close()

    return len(set(total)), len(set(unmapped))


def split_rna_iters(gdf):

    mrna_df = gdf[~gdf['type'].str.contains('rRNA')]
    rrna_df = gdf[gdf['type'].str.contains('rRNA')]

    assert len(gdf) == len(rrna_df) + len(mrna_df)

    mrna_iter = list(zip(mrna_df['chr'].tolist(), mrna_df['start'].tolist(), mrna_df['end'].tolist(),
                    mrna_df['gene-id'].tolist()))
    rrna_iter = list(zip(rrna_df['chr'].tolist(), rrna_df['start'].tolist(), rrna_df['end'].tolist()))

    return mrna_iter, rrna_iter


def count_rrna(in_file, rna_iter):

    st = pysam.AlignmentFile(in_file, 'rb')
    rrna = []

    for chro, start, end in rna_iter:

        st.reset()
        for r in st.fetch(chro, start, end):
            rrna.append(r.query_name)

    st.close()

    return len(set(rrna))


def mrna_tagged_file(in_file, rna_iter, rrna_iter, out_dir=''):

    st = pysam.AlignmentFile(in_file, 'rb')

    total, unmapped = bam_short_report(in_file)
    rrna = count_rrna(in_file, rrna_iter)

    out_file = out_file_name(in_file, out_dir, '.tagged.bam')
    tmp_file = out_file_name(in_file, out_dir, '.tmp.bam')

    tagged_reads = pysam.AlignmentFile(tmp_file, 'wb', template=st)

    spliced_list = []
    new_total_list = []

    for chro, start, end, gene in rna_iter:

        st.reset()
        for r in st.fetch(chro, start, end):
            qn = r.query_name
            if 'N' in r.cigarstring:
                spliced_list.append(qn)
            else:
                if r.has_tag('RG'):
                    r.set_tag('RG', None)
                r.tags += [('GE', gene)]
                tagged_reads.write(r)
                new_total_list.append(qn)

    tagged_reads.close()

    pysam.sort("-o", out_file, tmp_file)
    pysam.index(out_file)

    os.remove(tmp_file)

    tagged, _ = bam_short_report(out_file)
    spliced = len(set(spliced_list))
    new_total = len(set(new_total_list))
    assert new_total == tagged

    return out_file, rrna, total, spliced, tagged


def mrna_tagged_files(in_dir, in_gdf, out_dir=''):

    tagged_files = []

    gene_df = pd.read_csv(in_gdf, sep='\t')
    rna_iter, rrna_iter = split_rna_iters(gene_df)

    logger.info('Filtering ribosomal-RNA reads and tagging started.')
    logger.info('Source directory:\t%s' %in_dir)
    logger.info('-' * 50)

    for in_file in os.listdir(in_dir):

        if in_file.endswith('.bam') or in_file.endswith('.sam'):
            logger.info('Sample name:\t%s' % in_file)
            tagged_file, rrna, total, spliced, tagged = mrna_tagged_file(os.path.join(in_dir, in_file), rna_iter, rrna_iter, out_dir)
            tagged_files.append(tagged_file)
            if param_verbose:
                logger.info('-' * 25)
                logger.info('Total reads:\t\t\t%s' % format(total, ','))
                logger.info('Spliced reads ignored:\t\t%s' %format(spliced, ','))
                logger.info('Total reads mapped to genes:\t%s' % format(tagged, ','))
                logger.info('Output file:\t%s' % tagged_file)
                logger.info('')
            logger.info('')

    return tagged_files


if __name__ == "__main__":

    params = vars(get_args())

    param_in_dir = params['in_dir']
    param_in_gdf = params['in_gdf']
    param_verbose = params['verbose']
    param_out_dir = params['out_dir']
    param_save_report = params['s']

    out_directory = '.'
    if param_out_dir:
        out_directory = param_out_dir

    if param_save_report:
        report_file = os.path.join(out_directory, 'RemRNA.report.txt')
        logging.basicConfig(filename=report_file, level=logging.INFO)
    else:
        logging.basicConfig(level=logging.INFO)

    logging.getLogger().addHandler(logging.StreamHandler())

    logger.info('Call to RemrRNA module.')
    logger.info('')
    logger.info('This module removes ribosomal-RNA related reads from given .bam files')
    logger.info('It produces a gene-tagged .bam file as output.')
    logger.info('-' * 40)

    tagged_files = mrna_tagged_files(param_in_dir, param_in_gdf, out_dir=param_out_dir)

    logger.info('Finished successfully!')

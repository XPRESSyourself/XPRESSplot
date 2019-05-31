"""
XPRESSplot
A toolkit for navigating and analyzing gene expression datasets
alias: xpressplot

Copyright (C) 2019  Jordan A. Berg
jordan <dot> berg <at> biochem <dot> utah <dot> edu

This program is free software: you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free Software
Foundation, either version 3 of the License, or (at your option) any later
version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with
this program.  If not, see <https://www.gnu.org/licenses/>.
"""
from __future__ import print_function

"""IMPORT DEPENDENCIES"""
import os
import sys
import pandas as pd
import numpy as np
from sklearn import preprocessing
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

"""IMPORT INTERNAL DEPENDENCIES"""
from .utils_analyze import parallelize, count_threshold_util

"""INITIALIZATION PARAMETERS"""
# Retrieve path for scripts used in this pipeline, appended to argument dictionary for every function
__path__, xpressplot_arguments = os.path.split(__file__)

"""Create gene dictionary"""
def gene_length_dictionary(
    gtf,
    gene_name_prefix='gene_id \"',
    gene_name_location=0,
    sep='\t'):

    # Process gtf data for gene_name and gene_length
    gtf = pd.read_csv(
        str(gtf),
        sep = sep,
        header = None,
        comment = '#',
        low_memory = False)

    # Parse out relevant information
    gtf_genes = gtf.loc[gtf[2] == 'gene']
    gtf_genes['gene_name'] = gtf[8].str.split(';').str[gene_name_location]
    gtf_genes['gene_name'] = gtf_genes['gene_name'].map(lambda x: x.lstrip(gene_name_prefix).rstrip('\"').rstrip(' '))
    gtf_genes['length'] = abs((gtf[4]) - (gtf[3]))

    # Create dictionary
    length_df = gtf_genes[['gene_name','length']].copy()
    length_df = length_df.set_index('gene_name')
    del length_df.index.name
    length_df = length_df[~length_df.index.duplicated()]
    length_df.length = length_df.length / 1e3

    return length_df

"""Perform gene kilobase normalization"""
def rpk(
    data,
    length_df):

    data_c = data.copy()

    # Only accept genes in both length_df and data
    length_df = length_df[length_df.index.isin(data_c.index.values.tolist())]

    # Calculate
    data_rpk = data_c.div(length_df.length, axis=0)
    data_rpk = data_rpk.dropna(axis=0)

    return data_rpk

"""Perform reads per million sample normalization on RNAseq data"""
def rpm(
    data):

    data_c = data.copy()
    data_rpm = data_c / \
        (data_c.sum() / 1e6)

    return data_rpm

"""Perform transcripts per million normalization on RNAseq data"""
def tpm(
    data,
    gtf,
    gene_name_prefix='gene_id \"',
    gene_name_location=0,
    sep='\t'):

    length_df = gene_length_dictionary(
        gtf,
        sep = sep,
        gene_name_prefix = gene_name_prefix,
        gene_name_location = gene_name_location)

    data_rpk = rpk(
        data,
        length_df)
    data_tpm = rpm(data_rpk)

    return data_tpm

"""Perform reads/fragments per kilobase million sample normalization on RNAseq data"""
def r_fpkm(
    data,
    gtf,
    gene_name_prefix='gene_id \"',
    gene_name_location=0,
    sep='\t'):

    length_df = gene_length_dictionary(
        gtf,
        sep = sep,
        gene_name_prefix = gene_name_prefix,
        gene_name_location = gene_name_location)

    data_rpm = rpm(data)
    data_rpkm = rpk(
        data_rpm,
        length_df)

    return data_rpkm

"""Normalize out batch effects from RNAseq data"""
def batch_normalize(
    input_file,
    batch_file):

    # Get output file name
    if input_file.endswith('.txt') or input_file.endswith('.tsv'):
        output_file = str(input_file[:-4]) + '_batched.tsv'
    else:
        raise Exception('Unrecognized input_file delimiter type. Files must be tab-delimited')

    # Run sva combat in R
    os.system('rscript' \
        + ' ' + str(__path__) + '/batch_normalize.r' \
        + ' ' + str(input_file) \
        + ' ' + str(batch_file) \
        + ' ' + str(output_file))

"""Check sample means and medians"""
def check_samples(
    data):

    wid = len(list(data))
    ax = data.boxplot(
        column = list(data),
        figsize = (wid, (wid / 3)))
    ax.set_xlabel('Samples')
    ax.set_ylabel('Expression')

"""Cleans axis of NULL values"""
def clean_df(
    data,
    axis=0):

    data = data.dropna(axis=axis)
    data = data[~data.index.duplicated(keep=False)]

    return data

"""Remove genes from analysis where sequence coverage does not meet minimum"""
def threshold(
    data,
    minimum=None,
    maximum=None):

    data_c = data.copy()
    data_c = parallelize(
        count_threshold_util,
        data_c,
        minimum,
        maximum)

    return data_c

"""Prepare dataframes for analysis plotting functions found within analyze.py"""
def prep_data(
    data,
    info,
    gene_scale=True,
    print_means=False):

    # Convert data to float and drop bad values
    data_c = data.copy()
    data_scaled = data_c.astype(dtype='float')
    data_scaled = clean_df(data_scaled)
    data_c = clean_df(data_c)

    # Gene normalization
    if gene_scale == True:
        data_scaled[data_scaled.columns] = preprocessing.scale(data_scaled[data_scaled.columns], axis=1)

    if print_means == True:
        print(data_scaled.mean(axis=1))

    # Map labels to samples
    labels = pd.Series(info[1].values,index=info[0]).to_dict()
    data_c.loc['label'] = data_c.columns.map(labels.get)

    # Output collapsed dataframe
    newIndex = ['label'] + [ind for ind in data_c.index if ind != 'label']
    data_c = data_c.reindex(index=newIndex)

    return data_scaled, data_c

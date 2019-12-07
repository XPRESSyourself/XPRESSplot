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
if str(matplotlib.get_backend()).lower() != 'agg':
    import matplotlib.pyplot as plt
    plt.switch_backend('agg')
else:
    import matplotlib.pyplot as plt

"""INITIALIZATION PARAMETERS"""
# Retrieve path for scripts used in this pipeline, appended to argument dictionary for every function
__path__, xpressplot_arguments = os.path.split(__file__)
gtf_type_column = 2
gtf_leftCoordinate_column = 3
gtf_rightCoordinate_column = 4
gtf_annotation_column = 8
gtf_transcript_column = 9
gtf_gene_id_column = 10
gtf_gene_name_column = 11
id_search = r'transcript_id \"(.*?)\"; '
gene_id_search = r'gene_id \"(.*?)\"; '
gene_name_search = r'gene_name \"(.*?)\"; '

"""Create gene dictionary"""
# Provide Longest GTF to get longest canonical record
def gene_length_dictionary(
    gtf,
    feature_type='exon', # or CDS
    identifier='gene_name', # or gene_id or transcript_id
    sep='\t'):

    if not str(gtf).endswith('.gtf'):
        print(str(gtf) + ' does not appear to be a GTF file')

    if not str(feature_type).lower() in ['cds', 'exon']:
        print('Must provide CDS or exon as feature_type')

    # Process gtf data for gene_name and gene_length
    gtf = pd.read_csv(
        str(gtf),
        sep = '\t',
        header = None,
        comment = '#',
        low_memory = False)

    # Get relevant metadata
    gtf = gtf[(gtf[gtf_type_column] == feature_type)]
    gtf[gtf_transcript_column] = gtf[gtf_annotation_column].str.extract(id_search)
    gtf[gtf_gene_id_column] = gtf[gtf_annotation_column].str.extract(gene_id_search)
    gtf[gtf_gene_name_column] = gtf[gtf_annotation_column].str.extract(gene_name_search)
    gtf['length'] = abs(gtf[gtf_rightCoordinate_column] - gtf[gtf_leftCoordinate_column]) + 1

    # Return a transcript id mapping dictionary for isoform normalization
    if identifier == 'transcript_id':
        gtf = gtf[[gtf_transcript_column, 'length']]
        gtf.columns = ['transcript', 'length']
        gtf = gtf.dropna()
        length_dict = pd.Series(gtf['length'].values,index=gtf['transcript']).to_dict()
        return length_dict

    else:
        if identifier == 'gene_id':
            col = gtf_gene_id_column
        else:
            col = gtf_gene_name_column

        gtf = gtf[[col, gtf_transcript_column, 'length']]
        gtf.columns = ['gene', 'transcript', 'length']

        # Make transcript gene mapping dictionary
        ref = gtf[['transcript', 'gene']]
        reference = pd.Series(ref['gene'].values,index=ref['transcript']).to_dict()

        # Group by transcript ID and get feature length sum
        records = gtf[['transcript', 'length']]
        records = records.sort_values('transcript')
        records = records.groupby('transcript').sum()
        records = records.reset_index()

        # Map exon space to each bam record based on its transcript ID
        records['gene'] = records['transcript'].map(reference)

        # Get longest transcript based on max exon or CDS space size along (not Ensembl canonical)
        records = records.sort_values('gene')
        records = records.groupby('gene').max()
        records = records.reset_index()

        length_dict = pd.Series(records['length'].values,index=records['gene']).to_dict()
        return length_dict

"""Perform gene kilobase normalization"""
def rpk(
    data,
    length_dict):

    data_c = data.copy()

    data_rpk = data_c.div(data_c.index.map(length_dict) / 1000, axis = 0)
    data_rpk['index'] = data_rpk.index
    data_rpk = data_rpk.drop_duplicates(subset='index', keep='first')
    #data_rpk = data_rpk.reindex(list(length_dict.keys()), axis=0)
    data_rpk = data_rpk.drop('index', axis=1)
    data_rpk = data_rpk.fillna(0)

    return data_rpk

"""Perform reads per million sample normalization on RNAseq data"""
def rpm(
    data):

    data_rpm = data / \
        (data.sum() / 1e6)
    data_rpm = data_rpm.fillna(0)

    return data_rpm

"""Perform reads/fragments per kilobase million sample normalization on RNAseq data"""
def r_fpkm(
    data,
    gtf,
    feature_type='exon', # other option -> CDS
    identifier='gene_name', # other options -> gene_id, transcript_id
    sep='\t'):

    length_dict = gene_length_dictionary(
        gtf,
        feature_type = feature_type,
        identifier = identifier,
        sep = sep)

    data_rpm = rpm(data)
    data_rpkm = rpk(
        data_rpm,
        length_dict)

    return data_rpkm

def rpkm(
    data,
    gtf,
    feature_type='exon', # other option -> CDS
    identifier='gene_name', # other options -> gene_id, transcript_id
    sep='\t'):

    r_fpkm(
        data,
        gtf,
        feature_type = feature_type,
        identifier = identifier,
        sep = sep)

def fpkm(
    data,
    gtf,
    feature_type='exon', # other option -> CDS
    identifier='gene_name', # other options -> gene_id, transcript_id
    sep='\t'):

    r_fpkm(
        data,
        gtf,
        feature_type = feature_type,
        identifier = identifier,
        sep = sep)

"""Perform transcripts per million normalization on RNAseq data"""
def tpm(
    data,
    gtf,
    feature_type='exon', # other option -> CDS
    identifier='gene_name', # other options -> gene_id, transcript_id
    sep='\t'):

    length_dict = gene_length_dictionary(
        gtf,
        feature_type = feature_type,
        identifier = identifier,
        sep = sep)

    data_rpk = rpk(
        data,
        length_dict)
    data_tpm = rpm(data_rpk)

    return data_tpm

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

    if minimum != None:
        data = data[data.min(axis=0) > minimum]

        return data

    if maximum != None:
        data = data[data.max(axis=0) > maximum]

        return data

    return data

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

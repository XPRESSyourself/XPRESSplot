"""
XPRESStools
A toolkit for navigating and analyzing gene expression datasets
alias: xpresstools

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

"""
IMPORT DEPENDENCIES
"""
import os, sys
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from sklearn import preprocessing
from .utils_analyze import parallelize, count_threshold_util

"""
INITIALIZATION PARAMETERS
"""
#Retrieve path for scripts used in this pipeline, appended to argument dictionary for every function
__path__, xpresstools_arguments = os.path.split(__file__)

"""
DESCRIPTION: Perform reads per million sample normalization on RNAseq data
"""
def rpm(data):

    data_c = data.copy()
    data_rpm = data_c / (data_c.sum() / 1e6)

    return data_rpm

"""
DESCRIPTION: Perform reads/fragments per kilobase million sample normalization on RNAseq data
"""
def r_fpkm(data, gtf, gene_name_prefix='gene_id \"', gene_name_location=0, sep='\t'):

    #process gtf data for gene_name and gene_length
    gtf = pd.read_csv(str(gtf),sep=sep,comment='#', low_memory=False, header=None)
    gtf_genes = gtf.loc[gtf[2] == 'gene']
    gtf_genes['gene_name'] = gtf[8].str.split(';').str[gene_name_location]
    gtf_genes['length'] = abs((gtf[4]) - (gtf[3]))
    gtf_genes['gene_name'] = gtf_genes['gene_name'].map(lambda x: x.lstrip(gene_name_prefix).rstrip('\"').rstrip(' '))

    #Create dictionary
    length_df = gtf_genes[['gene_name','length']].copy()
    length_df = length_df.set_index('gene_name')
    del length_df.index.name
    length_df = length_df[length_df.index.isin(data.index.values.tolist())]
    length_df.length = length_df.length / 1e3

    #Perform per-million calculations
    data_c = data.copy()
    data_rpm = rpm(data_c)

    #Perform kilobase calculations
    data_rpkm = data_rpm.div(length_df.length, axis=0)
    data_rpkm = data_rpkm.dropna(axis=0)

    return data_rpkm

"""
DESCRIPTION: Perform log2(TE) normalization on ribosome profiling values
ASSUMPTIONS: Table is such that the paired RPF and RNA files are next to each other, in that order
                Values have been properly normalized already (ie. RPM, RPKM)
"""
def te(data, samples=None, log2=True):

    data_c = data.copy()
    data_c += .1

    #Perform translation efficiency calculations
    y = 0
    z = 1
    if samples == None:
        samples = []
        for x in range(int(len(data_c.columns)/2)):
            name = str(data_c.columns[y]) + '_te'
            samples.append(name)
            data_c[name] = data_c[data_c.columns[y]]/data_c[data_c.columns[z]]
            #Move counters for next set of samples
            y = y + 2
            z = z + 2
    else:
        for x in samples:
            data_c[x] = data_c[data_c.columns[y]]/data_c[data_c.columns[z]]
            #Move counters for next set of samples
            y = y + 2
            z = z + 2

    #Get TE_normalized columns
    df_te = data_c[samples]

    #Perform log2 scaling of data
    if log2 == True:
        df_logte = np.log2(df_te)

    return df_logte

"""
DESCRIPTION: Log-scale a sample-normalized dataframe

VARIABLES:
data= Sample normalized, XPRESStools-formated dataframe
log_base= Log-scale to normalize data with (Options: 2 or 10)

ASSUMPTIONS:
Requires a properly formatted dataframe for XPRESStools usage where samples are normalized
"""
def log_scale(data, log_base=10):

    if log_base == 10:
        data_log = np.log10(data + 1)
    elif log_base == 2:
        data_log = np.log2(data + 1)
    else:
        raise Exception('Invalid log_base option provided')

    return data_log

"""
DESCRIPTION: Normalize out batch effects from RNAseq data

VARIABLES:
input_file= Path and file name for sequence dataframe to batch normalize. Output file will be placed in same location and with same prefix, but will include the '_batch.csv' suffix
batch_info= Path and file name for batch effect dataframe. First column should be sample names as in input_file, and second column are the batches. This file should not have a header of any kind
input_sep= Input dataframe delimiter
batch_sep= Batch file delimiter

ASSUMPTIONS:
Data has already been sample normalized
"""
def batch_normalize(input_file, batch_file):

    #Get output file name
    if input_file.endswith('.txt') or input_file.endswith('.tsv'):
        output_file = str(input_file[:-4]) + '_batched.tsv'
    else:
        raise Exception('Unrecognized input_file delimiter type. Files must be tab-delimited')

    #Run sva combat in R
    os.system('rscript ' + str(__path__) + '/batch_normalize.r ' + str(input_file) + ' ' + str(batch_file) + ' ' + str(output_file))

"""
DESCRIPTION: Check sample means and medians
METHODS: Output density plot for dataframe, barplot for each sample
VARIABLES:
USAGE:
ASSUMPTIONS:
Dataframe has been properly formatted so that probes or genes are rows and samples are columns
"""
def check_samples(data):

    wid = len(list(data))
    ax = data.boxplot(column=list(data), figsize=(wid,wid/3))
    ax.set_xlabel('Samples')
    ax.set_ylabel('Expression')

"""
DESCRIPTION: Cleans axis of NULL values
VARIABLES:
USAGE:
ASSUMPTIONS:
If dataframe has been properly formatted previously and genes are in rows, the default parameters will remove along the gene axis
"""
def clean_df(data, axis=0):

    data = data.dropna(axis=axis)
    data = data[~data.index.duplicated(keep=False)]

    return data

"""
DESCRIPTION: Remove genes from analysis where sequence coverage does not meet minimum

VARIABLES:
data= XPRESStools formatted data
minimum= Float or int of minimum count/read value to accept per gene (all samples need to meet this requirement to keep)
"""
def threshold(data, minimum=None, maximum=None):

    data_c = data.copy()
    data_c = parallelize(count_threshold_util, data_c, minimum, maximum)

    return data_c

"""
DESCRIPTION: Prepare dataframes for analysis plotting functions found within analyze.py

METHODS:
Original dataframe is unformatted besides adding labels from info to the first row of the dataframe
Formatted dataframe is scaled if option provided and dataframe is converted to float

VARIABLES:
data= XPRESStools formatted dataframe of expression values
info= XPRESStools formatted sample info dataframe
gene_scale= Scale genes (rows) of data
print_means= Print appropriate means that were scaled for verification

USAGE:
import XPRESStools as mat
df_scaled, df_collapsed = xp.prep_df(df_collapsed, info)

ASSUMPTIONS:
Requires properly formatted df and info dataframes for XPRESStools usage
"""
def prep_data(data, info, gene_scale=True, print_means=False):

    #Convert data to float and drop bad values
    data_c = data.copy()
    data_scaled = data_c.astype(dtype='float')
    data_scaled = clean_df(data_scaled)
    data_c = clean_df(data_c)

    #gene normalization
    if gene_scale == True:
        data_scaled[data_scaled.columns] = preprocessing.scale(data_scaled[data_scaled.columns], axis=1)

    if print_means == True:
        print(data_scaled.mean(axis=1))

    #Map labels to samples
    labels = pd.Series(info[1].values,index=info[0]).to_dict()
    data_c.loc['label'] = data_c.columns.map(labels.get)

    #Output collapsed dataframe
    newIndex = ['label'] + [ind for ind in data_c.index if ind != 'label']
    data_c = data_c.reindex(index=newIndex)

    return data_scaled, data_c

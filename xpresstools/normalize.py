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
from sklearn import preprocessing
from .utils import parallelize

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
    gtf_genes[gene_name_prefix] = gtf[8].str.split(';').str[gene_name_location]
    gtf_genes['length'] = abs((gtf[4]) - (gtf[3]))
    gtf_genes[gene_name_prefix] = gtf_genes[gene_name_prefix].map(lambda x: x.lstrip(gene_name_prefix).rstrip('\"').rstrip(' '))

    #Create dictionary
    length_df = gtf_genes[[gene_name_prefix,'length']].copy()
    dict_df.columns = ['gene', 'length']
    length_df = length_df.set_index('gene')
    del length_df.index.name

    #Perform per-million calculations
    data_c = data.copy()
    data_rpm = rpm(data_c)

    #Perform kilobase calculations
    data_rpkm = data_rpm.div(length_df.length, axis=0)
    data_rpkm = data_rpkm.dropna(axis=0)

    return data_rpkm

"""
DESCRIPTION: Log-scale a sample-normalized dataframe

VARIABLES:
data= Sample normalized, MICARtools-formated dataframe
log_base= Log-scale to normalize data with (Options: 2 or 10)

ASSUMPTIONS:
Requires a properly formatted dataframe for MICARtools usage where samples are normalized
"""
def log_scale(data, log_base=10):

    if log_base == 10:
        data_log = np.log10(data)
    elif log_base == 2:
        data_log = np.log2(data)
    else:
        print('Invalid log_base option provided')

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
def batch_normalize(input_file, batch_file, input_sep=',', batch_sep=','):

    #Get output file name
    if intput_sep == ',':
        output_file = str(input_file[:-4]) + '_batch.csv'
    elif intput_sep == '\t':
        output_file = str(input_file[:-4]) + '_batch.tsv'
    else:
        raise Exception('Unrecognized input_file delimiter type')

    #Run sva combat in R
    os.system('rscript ' + str(__path__) + '/batch_normalize.r ' + str(input_file) + ' ' + str(batch_file) + str(input_sep) + ' ' + str(batch_sep) + ' ' + str(output_file))

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
    data = data[~data.index.duplicated()]

    return data

"""
DESCRIPTION: Remove genes from analysis where sequence coverage does not meet minimum

VARIABLES:
data= MICARtools formatted data
minimum= Float or int of minimum count/read value to accept per gene (all samples need to meet this requirement to keep)
"""
def threshold_util(data, minimum, maximum):

    data = data.T

    if minimum != None:
        data = data[data.columns[data.min() > minimum]]

    if maximum != None:
        data = data[data.columns[data.max() < maximum]]

    data = data.T

    return data

def threshold(data, minimum=None, maximum=None):

    data_c = data.copy()
    data_c = parallelize(threshold_util, data_c, minimum, maximum)

    return data_c

"""
DESCRIPTION: Normalize samples, prints sample axis means for verification

METHODS: For each sample axis, divide each cell by the sum the axis divided by the factor provided (default: 1e6)

VARIABLES:
data= Dataframe of microarray probe data
axis= Axis where samples are found in the dataframe
factor= Numeric value to scale samples by
print_means= Print appropriate means that were scaled for verification

USAGE:
import micartools as mat
df_norm = mat.sample_norm(df)
"""
def sample_norm(data, axis=1, factor=1e6, print_means=False):

    #Initialize axis variables based on user input
    if axis == 0:
        axis_2 = 1
    elif axis == 1:
        axis_2 = 0
    else:
        pass

    #Perform normalization
    data_norm = data.divide((data.sum(axis=axis_2) / float(factor)),axis=axis)

    if print_means == True:
        print(data_norm.mean(axis=axis_2))

    return data_norm

"""
DESCRIPTION: Prepare dataframes for analysis plotting functions found within analyze.py

METHODS:
Original dataframe is unformatted besides adding labels from info to the first row of the dataframe
Formatted dataframe is scaled if option provided and dataframe is converted to float

VARIABLES:
data= MICARtools formatted dataframe of expression values
info= MICARtools formatted sample info dataframe
gene_scale= Scale genes (rows) of data
print_means= Print appropriate means that were scaled for verification

USAGE:
import micartools as mat
df_scaled, df_collapsed = mat.prep_df(df_collapsed, info)

ASSUMPTIONS:
Requires properly formatted df and info dataframes for MICARtools usage
"""
def prep_data(data, info, gene_scale=True, print_means=False):

    #Convert data to float and drop bad values
    data_scaled = data.astype(dtype='float')
    data_scaled = data_scaled.dropna(axis=0)
    data = data.dropna(axis=0)

    #gene normalization
    if gene_scale == True:
        data_scaled[data_scaled.columns] = preprocessing.scale(data_scaled[data_scaled.columns], axis=1)

    if print_means == True:
        print(data_scaled.mean(axis=1))

    #Map labels to samples
    labels = pd.Series(info[1].values,index=info[0]).to_dict()
    data.loc['label'] = data.columns.map(labels.get)

    #Output collapsed dataframe
    newIndex = ['label'] + [ind for ind in data.index if ind != 'label']
    data = data.reindex(index=newIndex)

    return data_scaled, data

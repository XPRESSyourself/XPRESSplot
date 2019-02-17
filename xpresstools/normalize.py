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
DESCRIPTION: Normalize out batch effects from RNAseq data

VARIABLES:
input_file= Path and file name for sequence dataframe to batch normalize. Output file will be placed in same location and with same prefix, but will include the '_batch.csv' suffix
batch_info= Path and file name for batch effect dataframe. First column should be sample names as in input_file, and second column are the batches. This file should not have a header of any kind
input_sep= Input dataframe delimiter
batch_sep= Batch file delimiter

ASSUMPTIONS:
Data has already been sample normalized
"""
def batch_normalize(input_file, batch_file, input_sep=',', batch_sep='\t'):

    #Get output file name
    if intput_sep == ',':
        output_file = str(input_file[:-4]) + '_batch.csv'
    elif intput_sep == '\t':
        output_file = str(input_file[:-4]) + '_batch.tsv'
    else:
        raise Exception('Unrecognized input_file delimiter type')

    #Run sva combat in R
    os.system('rscript ' + str(__path__) + '/batch_normalize.r ' + str(input_file) + ' ' + str(batch_file) + str(input_sep) + ' ' + str(batch_sep) + ' ' + str(output_file))

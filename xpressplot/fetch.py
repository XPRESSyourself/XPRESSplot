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
import re
import pandas as pd

"""IMPORT INTERNAL DEPENDENCIES"""
from .normalize import clean_df
from .utils import check_directories

"""Get dataframe from user file"""
def get_df(
    file_name,
    delimiter="\t", low_memory=False, gene_axis='row'):

    # Read in file
    data = pd.read_csv(str(file_name), sep=delimiter, index_col=0, header=0, low_memory=low_memory)

    # Check data orientation
    if str(gene_axis).lower() == 'row':
        data = data
    elif str(gene_axis).lower() == 'col':
        data = data.T
    else:
        print("Incorrect gene_axis option specified")

    data = clean_df(data)

    return data

"""Get GEO dataframe and metadata"""
"""DEPRECATED
import GEOparse
def get_geo(
    geo_id,
    output_info=False, output_path="./"):

    # Get data
    gse = GEOparse.get_GEO(geo=str(geo_id).upper(), destdir=output_path) # Import GSE dataset

    data = gse.pivot_samples('VALUE')
    data = clean_df(data)

    # Get metadata
    # Write data to output file
    if output_info != False:
        with open(str(geo_id).upper() + '.txt', 'w+') as f: # Save all information as text file for reference
            for gsm_name, gsm in gse.gsms.items():
                f.write(gsm_name + '\n')
                for key, value in gsm.metadata.items():
                    f.write(" - %s : %s" % (key, ", ".join(value)) + '\n')

    # Populate metadata with sample ids and names
    metadata = pd.DataFrame(columns=['gsm', 'title']) # Create dataframe
    gsm_list, title_list, data_processing_list = [], [], []
    for gsm_name, gsm in gse.gsms.items():
        for key, value in gsm.metadata.items():
            if key == 'title':
                title_list.append(''.join(value))
            if key == 'geo_accession':
                gsm_list.append(''.join(value))
            if key == 'data_processing':
                data_processing_list.append(''.join(value))

    metadata['gsm'], metadata['title'] = gsm_list, title_list
    metadata.columns = range(metadata.shape[1])

    # Output processing style
    print('Data processing summary:\n' + str(set(data_processing_list))) # To determine if all samples have undergone the sample data processing

    # Clean data
    del data.columns.name
    del data.index.name

    # Clean metadata
    metadata[1] = metadata[1].apply(lambda x: x[0:(re.search("\d", x).start()) - 1])

    return data, metadata
"""


"""Get user file with pertinent sample information"""
def get_info(
    file_name,
    delimiter="\t", axis="col",
    sample_ids=0, labels=1):

    # Read in file
    info = pd.read_csv(str(file_name), sep=delimiter, header=None)

    # Reorganize dataframe as necessary so that data is column-wise and
    # sample_ids are the first column, labels are the second column
    if str(axis).lower() == 'col':
        if sample_ids == 0 and labels == 1:
            info = info
        else:
            info = info[[sample_ids, labels]] # Reorder dataframe columns
            info.columns = [0, 1]
    elif str(axis).lower() == 'row':
        info = info.T # Rotate dataframe to make it column-wise
        if sample_ids == 0 and labels == 1:
            info = info
        else:
            info = info[[sample_ids, labels]] # Reorder dataframe columns
            info.columns = [0, 1]
    else:
        print("Incorrect axis option specified")

    return info

"""Drop samples by sample IDs -- pass in a list of names"""
def drop_samples(
    data, ids):

    # Check file formats
    if type(ids) is not list:

        return

    # Drop samples in list
    else:
        data_dropped = data.drop(ids, axis=1)

        return data_dropped

"""Drop samples by label group name"""
def drop_label(
    data, info, label):

    # Check file formats
    if type(label) is not str:

        return

    # Drop samples by name (will grab from info df)
    else:
        # Create list of sample_ids based on name provided
        drop_ids = info[info[1] == str(label)]
        drop_ids_list = list(drop_ids[0])

        # Remove these samples from data
        data_dropped = data.drop(drop_ids_list, axis=1)

        return data_dropped

"""Keep samples by list of label names"""
def keep_labels(
    data, info,
    label_list=None):

    if label_list == None:
        label_list = list(set(info[1]))

        keep_ids = info[info[1].isin(label_list)]
        keep_ids_list = list(keep_ids[0])

        # Drop samples not given in list to keep
        data_dropped = data[keep_ids_list]

    else:
        # Check file formats
        if type(label_list) is not list:

            return

        # Keep samples by name (will grab from info df)
        else:

            # Create list of sample_ids based on what is not provided in keep list
            drop_ids = info[~info[1].isin(label_list)]
            drop_ids_list = list(drop_ids[0])

            # Drop samples not given in list to keep
            data_dropped = data.drop(drop_ids_list, axis=1)

    return data_dropped

"""Rename column names using dictionary"""
def rename_cols(
    data, converters):

    data_c = data.copy()
    dictionary = pd.Series(converters[1].values,index=converters[0]).to_dict()
    data_set = data_c.rename(columns=dictionary, inplace=False)
    return data_set

"""Rename values in a column (selected by providing column name) with a dictionary of keys and sort_values"""
def rename_rows(
    data, converters,
    label='index'):

    data_c = data.copy()

    if label == 'index':
        data_c['index'] = data_c.index
        dictionary = pd.Series(converters[1].values,index=converters[0]).to_dict()
        data_c[label] = data_c['index'].replace(dictionary)
        data_c = data_c.set_index('index')
        data_c.index.name = None

    else:
        dictionary = pd.Series(converters[1].values,index=converters[0]).to_dict()
        data_c[label] = data_c[label].replace(dictionary)

    return data_c


"""Compiles expression counts from multiple files into one table"""
def catenate_files(
    directory,
    file_suffix='txt', save_file=None, delimiter='\t',
    drop_rows=0):

    # Walk through raw data files within given directory
    if directory[:-1] != '/':
        directory = directory + '/'

    file_list = []
    for subdir, dirs, files in os.walk(directory):
        for f in files:
            if f.endswith(file_suffix): # Ignore hidden files and other uninterested files (particular for some submodules)
                file_list.append(f)
            else:
                pass

    # Sort files in alphabetical order (helps in formatting the count tables correctly)
    file_list = sorted(file_list)

    # Get gene list from the first, as well as length(#rows)
    with open(str(directory) + str(file_list[0])) as f:
        gene_names = pd.read_csv(f, header=None, usecols=[0], dtype=str, sep=delimiter)
        row_num = len(gene_names) - drop_rows
        if drop_rows > 0:
            gene_names = gene_names[:-drop_rows]

    # Populate dataframe with expression values
    data = pd.DataFrame(index=range(row_num))
    data['gene_names'] = gene_names

    for x in file_list:

        with open(str(directory) + x) as f:
            reader = pd.read_csv(f, header=None, usecols=[1], sep=delimiter)
            data[x] = pd.Series(data.index)
            length = len(reader)
            data[x] = reader

    # Remove gene_names label
    data = data.set_index('gene_names')
    data.index.name = None

    if save_file != None:
        data.to_csv(str(save_file),sep=delimiter)

    return data

"""Collate counts files"""
def count_table(
    file_list,
    gene_column=0, sample_column=1,
    sep='\t', drop_rows=5):

    # Read in first count file to get gene names
    df = pd.read_csv(str(file_list[0]), sep=sep, comment='#', header=None)
    pos_starter = [gene_column,sample_column]
    colname = df.columns[pos_starter]
    df = df[colname]

    # For the rest of the files in the file list, add the counts for that sample only
    for f in file_list[1:]:
        df_pull = pd.read_csv(str(f), sep=sep, comment='#', header=None)
        df = pd.concat([df, df_pull[df_pull.columns[sample_column]]], axis=1)
        df_pull = None

    # Final formatting clean up of table
    df_counts = df.copy()
    df = None
    df_counts = df_counts.set_index(0)
    df_counts.index.name = None

    # Remove path and file suffix from each file's name before adding as column names to table
    c = 0
    for x in file_list:
        file_list[c] = x[(x.rfind('/')+1):(x.find('.'))]
        c += 1
    df_counts.columns = file_list
    df_counts = df_counts[:-drop_rows]

    return df_counts

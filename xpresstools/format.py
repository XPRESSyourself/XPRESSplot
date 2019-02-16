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

"""
DESCRIPTION: Collate HTseq counts files

VARIABLES:
file_list= List of files with the path names appended to each file to be collated into a single count table
gene_column= Column location in all count files of gene names
gene_column= Column location in all count files of samples
sep= Separator of counts files

ASSUMPTIONS:
No headers are included in the count files
"""
def count_table(file_list, gene_column=0, sample_column=1, sep=','):

    #Read in first count file to get gene names
    df = pd.read_csv(str(file_list[gene_column]), sep=sep, comment='#', header=None)
    pos_starter = [gene_column,sample_column]
    colname = df.columns[pos_starter]
    df = df[colname]

    #For the rest of the files in the file list, add the counts for that sample only
    for f in file_list[1:]:
        df_pull = pd.read_csv(str(f), sep=sep, comment='#', header=None)
        df = pd.concat([df, df_pull[df_pull.columns[sample_column]]], axis=1)
        del df_pull

    #Final formatting clean up of table
    df_counts = df.copy()
    del df
    df_counts = df_counts.set_index(0)
    del df_counts.index.name

    #Remove path and file suffix from each file's name before adding as column names to table
    c = 0
    for x in file_list:
        file_list[c] = x[(x.rfind('/')+1):(x.find('.'))]
        c += 1
    df_counts.columns = file_list

    return df_counts

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

"""IMPORT DEPENDENCIES"""
import csv
import pandas as pd
pd.options.mode.chained_assignment = None

"""IMPORT INTERNAL DEPENDENCIES"""
from .utils import check_directories

"""Convert row names (genes) of dataframe using GTF as reference for new name"""
def convert_names(
    data,
    gtf,
    orig_name_label='gene_id',
    orig_name_location=0,
    new_name_label='gene_name',
    new_name_location=2,
    refill=None,
    add_space=True,
    sep='\t'):

    if add_space == True:
        orig_name_label = orig_name_label + ' \"'
        new_name_label = new_name_label + ' \"'

    # Import reference GTF
    gtf = pd.read_csv(str(gtf),sep=sep,comment='#', low_memory=False, header=None)

    # Parse out old and new names from GTF
    gtf_genes = gtf.loc[gtf[2] == 'gene']
    gtf_genes['original'] = gtf[8].str.split(';').str[orig_name_location]
    gtf_genes['new'] = gtf[8].str.split(';').str[new_name_location]
    gtf_genes['original'] = gtf_genes['original'].map(lambda x: x.lstrip(str(orig_name_label)).rstrip('\"').rstrip(' '))
    gtf_genes['new'] = gtf_genes['new'].map(lambda x: x.lstrip(str(new_name_label)).rstrip('\"').rstrip(' '))
    gtf_genes = gtf_genes[['original','new']].copy()

    # Create dictionary
    if refill != None:
        gene_dict = {}
        for index, row in gtf_genes.iterrows():
            if row[1] == str(refill):
                gene_dict[row[0]] = row[0]
            else:
                gene_dict[row[0]] = row[1]
    else:
        gene_dict = pd.Series(gtf_genes['new'].values,index=gtf_genes['original']).to_dict()

    # Replace old gene names/ids with new
    data_names = data.copy()
    data_names['new'] = data_names.index.to_series().map(gene_dict).fillna(data_names.index.to_series())
    data_names = data_names.set_index('new')
    data_names.index.name = None

    return data_names

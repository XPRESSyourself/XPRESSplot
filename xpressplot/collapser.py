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
import pandas as pd
import numpy as np

"""Create dictionary of probes and their gene names for the microarray probe collapser"""
def prep_collapser(
    reference, gene_list=None, no_multimappers=True):

    # Get probe/gene_name reference file and import into a dataframe
    df = pd.read_csv(str(reference), sep="\t", low_memory=False, comment='#') # No index
    df = df[['ID','Gene Symbol']]
    df = df.dropna()

    # If a custom gene list is given so that only certain probes are looked at and collapsed downstream
    # Take the gene list file, join genes for regex search
    # Create dictionary for just the probes corresponding with the genes of interest
    if gene_list != None:
        if type(gene_list) == list:

            gene_list = custom_list(gene_list)

            search = '|'.join(gene_list)
            search = search.upper()

            df = df[df['Gene Symbol'].str.contains(search)] # Only df where probes of interest are
            df_dict = df.set_index('ID')['Gene Symbol'].to_dict()
    else: # Take full probe set and create dictionary
        df_dict = df.set_index('ID')['Gene Symbol'].to_dict()

    if no_multimappers == False: # If allowing for multimappers, return dictionary as is
        return df_dict
    else: # Remove any probes that map to several genes
        df_dict_mod = {}
        for key, value in df_dict.items():
            if '///' in value:
                continue
            else:
                df_dict_mod[key] = value

        return df_dict_mod

"""Collapse probes of microarray dataset using previously prepared probe collapser dictionary (see prep_collapser function))"""
def run_collapser(
    data, dict):

    # This appears to fix the SettingwithCopyWarning error
    data_c = data.copy()

    # Get list of probes to find in df
    dict_list = list(dict.keys())

    # Only keep df rows where probes of interest are (in cases where only looking at certain genes, default: all genes)
    # Column 'name' header for probes in these files
    try:
        data_c = data_c[data_c['name'].isin(dict_list)]
    except:
        data_c['name'] = data_c.index
        data_c = data_c[data_c['name'].isin(dict_list)]

    # Map gene names in place of probes
    # Will set off SettingwithCopyWarning -- should be fine as we are replacing values in same column and it appears to change as expected
    data_c['name'] = data_c['name'].map(dict)

    # Set gene names as indices (allows for multiple indices with same name)
    # Needed to remove strings from df to allow for next step
    data_c = data_c.set_index('name', drop=True)

    # Force data to float
    data_numeric = data_c.apply(pd.to_numeric)

    # Reset indices to its own column to allow for sorting in next step
    data_numeric['name_sort'] =  data_numeric.index

    # groupby index (gene name) and collapse, taking the mean of rows with same name in 'name_sort'
    # Sets name_sort column names (post-collapse) as indices
    data_collapsed = data_numeric.groupby('name_sort').mean()

    # Remove double header for indices
    data_collapsed.index.name = None

    return data_collapsed

"""Ties prep_collapser and probe_collapse functions together"""
def probe_collapse(
    data, reference,
    gene_list=None, no_multimappers=True):

    dict = prep_collapser(reference, gene_list=gene_list, no_multimappers=no_multimappers)
    data_collapsed = run_collapser(data, dict)

    return data_collapsed

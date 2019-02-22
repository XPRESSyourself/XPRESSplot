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
import pandas as pd
import numpy as np
import scipy.stats as stats

"""
DESCRIPTION: Default data prep for all analysis functions
"""
def analysis_prep(data):
    data_c = data.copy()
    data_c = data_c.dropna(axis=0)
    return data_c

"""
DESCRIPTION: Axis-agnostic list-reader from dataframe
VARIABLES:
USAGE:
ASSUMPTIONS:
"""
def custom_list(file, delimiter=','):

    gene_df = pd.read_csv(str(file), sep=delimiter, index_col=False, header=None)
    genes = gene_df.values.tolist()
    gene_list = [val for sublist in genes for val in sublist]

    return gene_list

"""
DESCRIPTION
"""
def calculate_fc(data, label_comp, label_base):

    # Average every by cell line
    data['log$_2$(Fold Change)'] = np.log2((data.filter(regex=str(label_comp)).mean(axis=1)) / \
                                      (data.filter(regex=str(label_base)).mean(axis=1)))
    data['-log$_1$$_0$(P-Value)'] = ''

    return data

"""
DESCRIPTION
"""
def calculate_p(data, label_comp, label_base, drop_index):

    # Calculate p-value using 1-way ANOVA with replicates and append to df_oxsm_volc
    for row in data.iterrows():
        index, row_data = row
        comp_row = data.loc[index].filter(regex=str(label_comp)).values.tolist()
        base_row = data.loc[index].filter(regex=str(label_base)).values.tolist()

        # Append p_value to df_oxsm_volc
        try:
            statistic, p_value = stats.ttest_ind(comp_row, base_row)
            data.loc[index,'-log$_1$$_0$(P-Value)'] = float(-1 * (np.log10(p_value)))
        except:
            drop_index.append(index)

    data = data.drop(labels=drop_index, axis=0)

    return data

"""
DESCRIPTION: Reset plotting object to avoid bleed through
"""
def reset_plot(whitegrid, ax=False):

    if ax == True:
        del ax

    plt.close()
    plt.clf()

    if whitegrid == True:
        sns.set_style("whitegrid")
    else:
        sns.set_style("darkgrid")

"""
DESCRIPTION: Return the subsetted dataframe using a gene list
"""
def data_subset(data, gene_list):

    data_c = data.copy()

    #Check file formats
    if type(gene_list) is list:
        data_sub = data_c.reindex(labels=gene_list, axis=0)
    elif type(gene_list) is str:
        genes = custom_list(str(gene_list))
        data_sub = data_c.reindex(labels=genes, axis=0)
    else:
        print('Incorrect gene_list type provided')
        return

    return data_sub

"""
DESCRIPTION: Prepare a color map from palette dictionary
"""
def prep_palette(info, palette):

    if type(palette) is dict:
        info = info.T
        info.columns = info.iloc[0]
        info = info.reindex(info.index.drop(0))
        info = info.rename({1: 'samples'})
        labels = info.iloc[0]
        color_map = labels.map(palette)
        return color_map
    else:
        print('Error: a dictionary was not provided as palette')
        return

"""
DESCRIPTION: Add labels to dataframe for each sample
"""
def label(data, info):

    data_c = data.copy()

    #Prep data_scaled by adding labels from info
    labels = pd.Series(info[1].values,index=info[0]).to_dict()
    data_c.loc['label'] = data_c.columns.map(labels.get)

    #Output collapsed dataframe
    newIndex = ['label'] + [ind for ind in data_c.index if ind != 'label']
    data_c = data_c.reindex(index=newIndex)
    data_c = data_c.T.set_index('label', drop=True)
    data_c = data_c.T.dropna(axis=0)

    return data_c

"""
DESCRIPTION: Unstack dataframe for use with some seaborns functions
"""
def unstack_data(data):

    data_c = data.copy()
    data_unstacked = data_c.unstack().reset_index()
    data_unstacked = data_unstacked.rename(columns={data_unstacked.columns[0]: 'type', data_unstacked.columns[1]: 'gene', data_unstacked.columns[2]: 'expr'})

    return data_unstacked

"""
DESCRIPTION: Convert gene in dataframe to array
"""
def get_array(data, gene):

    data_c = data.copy()

    gene = data_c[str(gene)].values.tolist()
    gene = np.array(gene).astype(np.float)
    gene = np.ndarray.tolist(gene)

    return gene

"""
"""
def make_linreg(data, gene1, gene2):

    data_c = data.copy()
    gene_a = get_array(data_c, gene1)
    gene_b = get_array(data_c, gene2)

    slope, intercept, r_value, p_value, std_err = linregress(gene_a, gene_b)
    x = np.linspace(data_c[str(gene1)].min(), data_c[str(gene1)].max(), 100)
    y = (slope * x) + intercept

    title = 'r = ' + "%.2f" % round(r_value,4)

    return x, y, r_value, title

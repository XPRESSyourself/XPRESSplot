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
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from mpl_toolkits.mplot3d import Axes3D

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

"""
Initial variable checks
"""
def init_pca(principle_components, _3d_pca, plotly_login):

    if len(principle_components) != 2 and _3d_pca == False:
        print('Incompatible options provided for principle_components and _3d_pca')
        return

    if plotly_login != None and _3d_pca == False:
        print('Only provide plotly login for _3d_pca')
        return

    if _3d_pca == True:
        if principle_components == [1,2]:
            principle_components = [1,2,3]

        elif len(principle_components) != 3:
            print('Incompatible options provided for principle_components and _3d_pca')
            return

        else:
            print('Not sure if this will ever be triggered, but if it is need to check case')
            return

        if plotly_login != None and type(plotly_login) is not list:
            print('Provide proper plotly login information in form of list')
            print("['userid','api key']")
            return

    return principle_components

"""
DESCRIPTION: Generate scree plot for principle components
"""
def make_scree(pca, n_components, save_fig, dpi, bbox_inches, scree_only, grid, whitegrid):

    vari = 'Explained variation per principal component: {}'.format(np.round(pca.explained_variance_ratio_, decimals=4)*100)
    scree = np.round(pca.explained_variance_ratio_, decimals=4)*100

    #Plot scree
    sing_vals = np.arange(n_components) + 1

    ax = sns.lineplot(x=sing_vals, y=scree, color="red")
    ax.set(xlabel='Principal Component', ylabel='Proportion of Variance Explained', title='Scree Plot')
    plt.savefig(str(save_fig[:-4]) + '_scree.pdf', dpi=dpi, bbox_inches=bbox_inches)

    if scree_only == True:
        plt.show()

    if grid == False:
        ax.grid(False)

    #Remove scree from memory to prevent plot bleeding
    reset_plot(whitegrid, ax=True)

    return scree

"""
DESCRIPTION: Add confidence intervals to scatterplot
"""
def set_confidence(df_pca, pca_plot, unique_labels, palette):

    for x in unique_labels:
        #slice df into label specific datasets
        df_slice = df_pca[df_pca['label'] == x]

        #make numpy array from label-specific dfs
        x_slice = df_slice.PCa.values
        y_slice = df_slice.PCb.values

        #pca maths
        cov = np.cov(x_slice, y_slice)
        lambda_, v = np.linalg.eig(cov)
        theta = np.degrees(np.arctan2(*v[:,0][::-1]))
        lambda_ = np.sqrt(lambda_)

        #plot
        pca_plot.add_patch(patches.Ellipse(xy=(np.mean(x_slice), np.mean(y_slice)),
                          width=lambda_[0]*ci*2, height=lambda_[1]*ci*2,
                          angle=theta,
                          alpha=0.3, facecolor=palette[x], edgecolor='black', linewidth=1, linestyle='solid')
                          )

"""
DESCRIPTION: 2D Non-interactive PCA scatterplot
"""
def pca2(df_pca, unique_labels, palette, principle_components, scree, order_legend, save_fig, dpi, bbox_inches):

    pca_plot = sns.scatterplot(df_pca.PCa, df_pca.PCb, hue=df_pca['label'], palette=palette)
    set_confidence(df_pca, pca_plot, unique_labels, palette)

    # Put the legend out of the figure
    handles,labels = pca_plot.get_legend_handles_labels()

    if order_legend != None:
        if type(order_legend) is list:
            plt.legend([handles[idx] for idx in order_legend],[labels[idx] for idx in order_legend], bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
        else:
            plt.legend(handles, labels, bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
            print('order_legend datatype is invalid -- plotting samples in default order...')
    else:
        plt.legend(handles, labels, bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)

    plt.xlabel('PC' + str(principle_components[0]) + ' (' + str(round(scree[(principle_components[0] - 1)],2)) + '%)')
    plt.ylabel('PC' + str(principle_components[1]) + ' (' + str(round(scree[(principle_components[1] - 1)],2)) + '%)')

    if grid == False:
        plt.grid(False)

    if save_fig != None:
        #Save plot
        plt.title(str(title))
        plt.savefig(str(save_fig), dpi=dpi, bbox_inches=bbox_inches)

    plt.show()

"""
DESCRIPTION: 3D Non-interactive PCA scatterplot
"""
def pca3(df_pca, palette, save_fig, dpi, bbox_inches):

    fig = plt.figure()
    ax = Axes3D(fig)

    df_pca.columns = ['PCa', 'PCb', 'PCc', 'label']
    unique_labels = df_pca['label'].unique() #Gather unique labels

    pca0 = df_pca.loc[df_pca['label'] == unique_labels[0]]
    pca1 = df_pca.loc[df_pca['label'] == unique_labels[1]]
    pca2 = df_pca.loc[df_pca['label'] == unique_labels[2]]

    x0 = pca0.PCa.values
    y0 = pca0.PCb.values
    z0 = pca0.PCc.values
    ax.scatter(x0, y0, z0, c=palette[unique_labels[0]], label=str(unique_labels[0]))

    x1 = pca1.PCa.values
    y1 = pca1.PCb.values
    z1 = pca1.PCc.values
    ax.scatter(x1, y1, z1, c=palette[unique_labels[1]], label=str(unique_labels[1]))

    x2 = pca2.PCa.values
    y2 = pca2.PCb.values
    z2 = pca2.PCc.values
    ax.scatter(x2, y2, z2, c=palette[unique_labels[2]], label=str(unique_labels[2]))

    ax.set_xlabel(str(pc_list[0]))
    ax.set_ylabel(str(pc_list[1]))
    ax.set_zlabel(str(pc_list[2]))
    ax.legend()

    plt.show()

    if save_fig != None:
        plt.savefig(str(save_fig), dpi=dpi, bbox_inches=bbox_inches)

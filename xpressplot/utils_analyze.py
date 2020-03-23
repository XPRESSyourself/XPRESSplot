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
import pandas as pd
from functools import partial
from multiprocessing import cpu_count, Pool
import numpy as np
import scipy.stats as stats
import matplotlib
if str(matplotlib.get_backend()).lower() != 'agg':
    import matplotlib.pyplot as plt
    plt.switch_backend('agg')
else:
    import matplotlib.pyplot as plt
import matplotlib.patches as patches
from mpl_toolkits.mplot3d import Axes3D
import seaborn as sns

"""Default data prep for all analysis functions"""
def analysis_prep(
    data):

    data_c = data.copy()
    data_c = data_c.dropna(axis=0)

    return data_c

"""Axis-agnostic list-reader from dataframe"""
def custom_list(
    file, delimiter=','):

    gene_df = pd.read_csv(str(file), sep=delimiter, index_col=False, header=None)
    genes = gene_df.values.tolist()
    gene_list = [val for sublist in genes for val in sublist]

    return gene_list

"""Reset plotting object to avoid bleed through"""
def reset_plot(
    whitegrid, ax=False):

    if ax == True:
        ax = None

    plt.close()
    plt.clf()

    if whitegrid == True:
        sns.set_style("whitegrid")
    else:
        sns.set_style("darkgrid")

"""Return the subsetted dataframe using a gene list"""
def data_subset(
    data, gene_list):

    data_c = data.copy()

    # Check file formats
    if type(gene_list) is list:
        data_sub = data_c.reindex(labels=gene_list, axis=0)
    elif type(gene_list) is str:
        genes = custom_list(str(gene_list))
        data_sub = data_c.reindex(labels=genes, axis=0)
    else:
        print('Incorrect gene_list type provided')
        return

    return data_sub

"""Prepare a color map from palette dictionary"""
def prep_palette(
    info, palette):

    if type(palette) is dict:
        info = info.T
        info.columns = info.iloc[0]
        info = info.reindex(info.index.drop(0))
        info = info.rename({1: 'samples'})
        labels = info.iloc[0]
        color_map = labels.map(palette)

        return color_map

    else:
        return print('Error: a dictionary was not provided as palette')

"""Add labels to dataframe for each sample"""
def label(
    data, info):

    data_c = data.copy()

    # Prep data_scaled by adding labels from info
    labels = pd.Series(info[1].values,index=info[0]).to_dict()
    data_c.loc['label'] = data_c.columns.map(labels.get)

    # Output collapsed dataframe
    newIndex = ['label'] + [ind for ind in data_c.index if ind != 'label']
    data_c = data_c.reindex(index=newIndex)
    data_c = data_c.T.set_index('label', drop=True)
    data_c = data_c.T.dropna(axis=0)

    return data_c

"""Unstack dataframe for use with some seaborns functions"""
def unstack_data(data):

    data_c = data.copy()
    data_unstacked = data_c.unstack().reset_index()
    data_unstacked = data_unstacked.rename(columns={data_unstacked.columns[0]: 'type', data_unstacked.columns[1]: 'gene', data_unstacked.columns[2]: 'expr'})

    return data_unstacked

"""Convert gene in dataframe to array"""
def get_array(data, gene):

    data_c = data.copy()

    gene = data_c.loc[str(gene)].values.tolist()
    gene = np.array(gene).astype(np.float)
    gene = np.ndarray.tolist(gene)

    return gene

"""Calculate linreg values for two genes (rows) of dataframe"""
def make_linreg(
    data, gene1, gene2):

    data_c = data.copy()
    gene_a = get_array(data_c, gene1)
    gene_b = get_array(data_c, gene2)

    slope, intercept, r_value, p_value, std_err = stats.linregress(gene_a, gene_b)
    x = np.linspace(data_c.loc[str(gene1)].min(), data_c.loc[str(gene1)].max(), 100)
    y = (slope * x) + intercept

    title = 'r = ' + "%.2f" % round(r_value,4)

    return x, y, r_value, title

"""Add X-axis threshold lines to plot"""
def add_x_threshold(
    threshold, threshold_color, ax):

    for x in threshold:
        if type(x) is int or type(x) is float:
            ax.axvline(x, ls='--', color=str(threshold_color))
        else:
            print('Invalid X threshold provided')

"""Add Y-axis threshold lines to plot"""
def add_y_threshold(
    threshold, threshold_color, ax):

    for x in threshold:
        if type(x) is int or type(x) is float:
            ax.axhline(x, ls='--', color=str(threshold_color))
        else:
            print('Invalid Y threshold provided')

"""Initial variable checks"""
def init_pca(
    principle_components, _3d_pca, plotly_login):

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

"""Generate scree plot for principle components"""
def make_scree(
    pca, n_components,
    save_fig, dpi, bbox_inches,
    scree_only, save_scree,
    grid, whitegrid):

    vari = 'Explained variation per principal component: {}'.format(np.round(pca.explained_variance_ratio_, decimals=4)*100)
    scree = np.round(pca.explained_variance_ratio_, decimals=4)*100

    # Plot scree
    sing_vals = np.arange(n_components) + 1

    ax = sns.lineplot(
            x = sing_vals,
            y = scree,
            color = "red")

    ax.set(xlabel='Principal Component', ylabel='Proportion of Variance Explained', title='Scree Plot')

    if save_scree == True:
        plt.savefig(str(save_fig[:-4]) + '_scree.png', dpi=dpi, bbox_inches=bbox_inches)

    if scree_only == True:
        plt.show()

    if grid == False:
        ax.grid(False)

    # Remove scree from memory to prevent plot bleeding
    reset_plot(whitegrid, ax=True)

    return scree

"""Add confidence intervals to scatterplot
NOTE: This code is adapted from Jaime (https://stackoverflow.com/a/20127387/9571488) and Ben (https://stackoverflow.com/a/25022642/9571488) on Stack Overflow
"""
def set_confidence(
    df_pca, pca_plot,
    unique_labels, palette, ci):

    for x in unique_labels:
        # Slice df into label specific datasets
        df_slice = df_pca[df_pca['label'] == x]

        # Make numpy array from label-specific dfs
        x_slice = df_slice.PCa.values
        y_slice = df_slice.PCb.values

        # Confidence interval math
        cov = np.cov(x_slice, y_slice)
        lambda_, v = np.linalg.eig(cov)
        theta = np.degrees(np.arctan2(*v[:, 0][::-1]))
        lambda_ = np.sqrt(lambda_)

        # Plot
        pca_plot.add_patch(
            patches.Ellipse(
                xy = (np.mean(x_slice), np.mean(y_slice)),
                width = lambda_[0] * ci * 2,
                height = lambda_[1] * ci * 2,
                angle = theta,
                alpha = 0.3,
                facecolor = palette[x],
                edgecolor = 'black',
                linewidth = 1,
                linestyle = 'solid'
                )
            )

"""2D Non-interactive PCA scatterplot"""
def pca2(
    df_pca,
    unique_labels,
    palette,
    principle_components,
    scree,
    order_legend,
    save_fig,
    dpi,
    bbox_inches,
    ci,
    grid,
    title,
    size,
    highlight_color,
    highlight_names):

    ax = sns.scatterplot(
            df_pca.PCa,
            df_pca.PCb,
            hue = df_pca['label'],
            palette = palette,
            s = size)

    set_confidence(df_pca, ax, unique_labels, palette, ci)

    # Put the legend out of the figure
    make_legend(ax, order_legend, highlight_color, highlight_names)

    plt.xlabel('PC' + str(principle_components[0]) + ' (' + str(round(scree[(principle_components[0] - 1)],2)) + '%)')
    plt.ylabel('PC' + str(principle_components[1]) + ' (' + str(round(scree[(principle_components[1] - 1)],2)) + '%)')

    if grid == False:
        plt.grid(False)

    if save_fig != None:
        # Save plot
        plt.title(str(title))
        plt.savefig(str(save_fig), dpi=dpi, bbox_inches=bbox_inches)

    plt.show()

"""3D Non-interactive PCA scatterplot"""
def pca3(
    df_pca, palette,
    save_fig, dpi, bbox_inches,
    size, order_legend):

    fig = plt.figure()
    ax = Axes3D(fig)

    pc_list = df_pca.columns.tolist()
    while 'label' in pc_list:
        pc_list.remove('label')

    df_pca.columns = ['PCa', 'PCb', 'PCc', 'label']
    unique_labels = df_pca['label'].unique() # Gather unique labels

    # Convert color tuples to arrays
    for x in unique_labels:
        pca0 = df_pca.loc[df_pca['label'] == str(x)]
        x0 = pca0.PCa.values
        y0 = pca0.PCb.values
        z0 = pca0.PCc.values
        ax.scatter(x0, y0, z0, c=np.array([palette[str(x)]]), label=str(x), s=size)

    ax.set_xlabel(str(pc_list[0]))
    ax.set_ylabel(str(pc_list[1]))
    ax.set_zlabel(str(pc_list[2]))

    plt.show()

"""Parallelize function on a chunk of a dataframe"""
def parallelize(
    func, *args):

    cores = cpu_count() # Number of CPU cores on your system
    partitions = cpu_count() # Define as many partitions as you want

    data_split = np.array_split(args[0], partitions)
    pool = Pool(cores)

    if func == calculate_fc:
        func = partial(calculate_fc, label_comp=args[1], label_base=args[2])
    elif func == calculate_p:
        func = partial(calculate_p, label_comp=args[1], label_base=args[2])
    else:
        return

    data = pd.concat(pool.map(func, data_split))

    pool.close()
    pool.join()

    return data

"""Calculate log2 fold changes between two group names"""
def calculate_fc(
    data, label_comp, label_base):

    # Average every by cell line
    data['log$_2$(Fold Change)'] = np.log2(
                                    (data.filter(regex=str(label_comp)).mean(axis=1)) \
                                    / (data.filter(regex=str(label_base)).mean(axis=1))
                                    )

    return data

"""Calculate p values between two group names"""
def calculate_p(
    data, label_comp, label_base):

    drop_index = []

    data['-log$_1$$_0$(P-Value)'] = '' # Initialize column to store p-values

    # Calculate p-value using T-test for two independent groups
    for row in data.iterrows():
        index, row_data = row
        comp_row = data.loc[index].filter(regex=str(label_comp)).values.tolist()
        base_row = data.loc[index].filter(regex=str(label_base)).values.tolist()

        if len(comp_row) < 2 or len(base_row) < 2:
            raise Exception('Calculating a P-value requires replicates for each sample type')

        # Append p_value to dataframe
        try:
            statistic, p_value = stats.ttest_ind(comp_row, base_row)
            data.loc[index,'-log$_1$$_0$(P-Value)'] = float(-1 * (np.log10(p_value)))
        except:
            drop_index.append(index)

    data = data.drop(labels=drop_index, axis=0)

    return data

"""Format threshold as list"""
def make_threshold_list(threshold):

    if threshold != None and type(threshold) != list:
        threshold = [threshold]

    return threshold

"""Output matrix of genes meeting threshold criteria"""
def output_threshold(
    data_c, x_threshold, y_threshold,
    save_threshold_hits, save_threshold_hits_delimiter):

    x_threshold = make_threshold_list(x_threshold)
    y_threshold = make_threshold_list(y_threshold)

    if len(x_threshold) == 2 and len(y_threshold) == 1:

        data_c = data_c.rename(columns = {'log$_2$(Fold Change)':'log2 Fold Change', '-log$_1$$_0$(P-Value)':'-log10 P-Value'})
        df_c = data_c[['log2 Fold Change', '-log10 P-Value']].copy()

        # Get up- and down-regulated hits
        df_up = df_c.loc[(df_c['log2 Fold Change'] > max(x_threshold)) & (df_c['-log10 P-Value'] > y_threshold[0])]
        df_down = df_c.loc[(df_c['log2 Fold Change'] < min(x_threshold)) & (df_c['-log10 P-Value'] > y_threshold[0])]

        # Append hits tables and output
        thresh_hits = df_up.append(df_down)
        thresh_hits.to_csv(str(save_threshold_hits), sep=save_threshold_hits_delimiter)

"""Generate legend based on plotted values"""
def make_legend(ax, order_legend, highlight_color, highlight_names):

    handles, labels = ax.get_legend_handles_labels()

    if len(labels) > 1:
        if order_legend != None:
            if type(order_legend) is list:
                plt.legend(
                        [handles[idx] for idx in order_legend],
                        [labels[idx] for idx in order_legend],
                        bbox_to_anchor = (1.05, 1),
                        loc = 2,
                        borderaxespad = 0.0)
            else:
                plt.legend(
                        handles,
                        labels,
                        bbox_to_anchor = (1.05, 1),
                        loc = 2,
                        borderaxespad = 0.0)
                print('order_legend datatype is invalid -- plotting samples in default order...')
        else:
            plt.legend(
                    handles,
                    labels,
                    bbox_to_anchor = (1.05, 1),
                    loc = 2,
                    borderaxespad = 0.0)

    elif len(labels) == 0:
        if highlight_names != None and type(highlight_names) is list:
            if type(highlight_color) is str:
                highlight_color = [highlight_color]

            if len(highlight_names) == len(highlight_color):
                # Make a custom legend using the provided info
                # Adapted from https://stackoverflow.com/a/47749903
                f = lambda m,c: plt.plot(
                                        [], [],
                                        marker = 'o',
                                        color = c,
                                        ls = 'none')[0]
                handles = [f("s", highlight_color[i]) for i in range(len(highlight_color))]
                plt.legend(
                        handles,
                        highlight_names,
                        bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
            else:
                raise Exception('Highlight colors and names lists are not of equal size')
        else:
            pass

    else:
        pass

    return ax

"""Highlight specific points on scatter plot"""
def highlight_markers(
    data, x, y,
    highlight_points, highlight_color,
    highlight_names, alpha_highlights,
    highlight_size, ax):

    if all(isinstance(z, list) for z in highlight_points):
        p = 0
        for h in highlight_points:
            data_sub = data_subset(data, highlight_points[p]).T
            data_genes = data_sub.dropna(axis=1)
            if type(alpha_highlights) is list:
                ax = sns.scatterplot(
                        x = data_genes.loc[str(x)],
                        y = data_genes.loc[str(y)],
                        color = str(highlight_color[p]),
                        alpha = alpha_highlights[p],
                        s = highlight_size)
            else:
                ax = sns.scatterplot(
                        x = data_genes.loc[str(x)],
                        y = data_genes.loc[str(y)],
                        color = str(highlight_color[p]),
                        alpha = alpha_highlights,
                        s = highlight_size)

            p += 1

    else:
        data_sub = data_subset(data, highlight_points).T
        data_genes = data_sub.dropna(axis=1)
        ax = sns.scatterplot(
                x = data_genes.loc[str(x)],
                y = data_genes.loc[str(y)],
                color = str(highlight_color),
                alpha = alpha_highlights,
                s = highlight_size)

    return ax

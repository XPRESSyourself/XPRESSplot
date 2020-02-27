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
import pandas as pd
import numpy as np
from scipy.stats import linregress
from sklearn.decomposition import PCA
import matplotlib
if str(matplotlib.get_backend()).lower() != 'agg':
    import matplotlib.pyplot as plt
    plt.switch_backend('agg')
else:
    import matplotlib.pyplot as plt
import seaborn as sns
sns.set(font='arial')
jakes_cmap = sns.diverging_palette(212, 61, s=99, l=77, sep=1, n=16, center='dark') #Custom aesthetics
import plotly
import plotly.offline as py
import plotly_express as px

"""IMPORT INTERNAL DEPENDENCIES"""
from .utils_analyze import *

"""INITIALIZATION PARAMETERS"""
#Retrieve path for scripts used in this pipeline, appended to argument dictionary for every function
__path__, xpressplot_arguments = os.path.split(__file__)

"""Plots heatmap for prep_data formatted data"""
def heatmap(
    data,
    info,
    sample_palette=None,
    gene_info=None,
    gene_palette=None,
    gene_list=None,
    col_cluster=True,
    row_cluster=False,
    metric='euclidean',
    method='centroid',
    font_scale=.8,
    cmap=jakes_cmap,
    center=0,
    xticklabels=True,
    yticklabels=True,
    linewidths=0,
    linecolor='#DCDCDC',
    cbar_kws=None,
    figsize=(16,6.5),
    save_fig=None,
    dpi=600,
    bbox_inches='tight'):

    reset_plot(True)
    data_c = analysis_prep(data)

    # Set colors for plotting
    if sample_palette != None:
        sample_palette = prep_palette(info, sample_palette)

    if gene_palette != None and isinstance(gene_info, pd.DataFrame):
        gene_palette = prep_palette(gene_info, gene_palette)

    # Custom panel heatmap
    if gene_list != None:
        data_sub = data_subset(data_c, gene_list)
        plot_data = data_sub.dropna(axis=0)
    else:
        plot_data = data_c.dropna(axis=0)

    # Generate clustermap
    sns.set(font_scale=float(font_scale))
    ax = sns.clustermap(
            plot_data,
            cmap = cmap,
            center = float(center),
            metric = str(metric),
            method = str(method),
            xticklabels = xticklabels,
            yticklabels = yticklabels,
            linewidths = float(linewidths),
            linecolor = linecolor,
            col_cluster = col_cluster,
            row_cluster = row_cluster,
            col_colors = sample_palette,
            row_colors = gene_palette,
            figsize = figsize,
            cbar_kws = cbar_kws)

    # Save figure
    if save_fig is not None:
        plt.savefig(str(save_fig), dpi=dpi, bbox_inches=str(bbox_inches))

"""Create violin plots of subset of gene expressions or all gene expression by sample"""
def multigene_overview(
    data, info,
    palette=None, gene_list=None, order=None,
    scale='area', grid=False, whitegrid=False,
    title=None,
    save_fig=None, dpi=600, bbox_inches='tight'):

    reset_plot(whitegrid)
    data_c = data.copy()

    # Custom panel heatmap
    if gene_list != None:
        data_sub = data_subset(data_c, gene_list)
        plot_data = data_sub.dropna(axis=0)
    else:
        plot_data = data_c.dropna(axis=0)

    plot_data_c = analysis_prep(plot_data)
    plot_data_c = label(plot_data_c, info)
    plot_data_c = unstack_data(plot_data_c)

    # Final formatting check
    try:
        plot_data_c[["expr"]] = plot_data_c[["expr"]].apply(pd.to_numeric)
    except:
        raise Exception('Data is not properly formatted for downstream plotting')

    # Plot
    if not order != None and type(order) is list or palette != None and type(palette) is list:
        return

    ax = sns.violinplot(
            x = plot_data_c['type'],
            y = plot_data_c['expr'],
            data = plot_data_c,
            order = order,
            palette = palette,
            scale = scale)

    ax.set_xlabel('')
    ax.set_ylabel('Expression')

    if grid == False:
        ax.grid(False)

    # Save fig
    if save_fig is not None:
        if title is None:
            plt.savefig(str(save_fig), dpi=dpi, bbox_inches=bbox_inches)
        else:
            ax.set_title(str(title))
            plt.savefig(str(save_fig), dpi=dpi, bbox_inches=bbox_inches)

    plt.show()

"""Plot boxplot overlaid with swarmplot of each sample type's gene expression for the given gene"""
def gene_overview(
    data, info, gene_name, palette,
    order=None,
    grid=False, whitegrid=False,
    save_fig=None, dpi=600, bbox_inches='tight'):

    reset_plot(whitegrid)
    data_c = analysis_prep(data).T

    # Prep data_scaled by adding labels from info to column (samples are rows)
    labels = pd.Series(info[1].values,index=info[0]).to_dict()
    data_c['label'] = data_c.index.to_series().map(labels)
    gene_df = data_c[[str(gene_name), 'label']]
    gene_df[[str(gene_name)]] = gene_df[[str(gene_name)]].apply(pd.to_numeric, errors='coerce')

    # Swarm plot
    ax = sns.catplot(
            x = 'label',
            y = str(gene_name),
            data = gene_df,
            color = 'black',
            order = order,
            kind = 'swarm')

    #Boxplot, fliersize=0 removes outlier diamonds from sns
    ax = sns.boxplot(
            x = 'label',
            y = str(gene_name),
            data = gene_df,
            width = 0.3,
            fliersize = 0,
            order = order,
            palette = palette)

    plt.setp(ax.collections, sizes=[12]) #Resize markers for catplot

    if grid == False:
        ax.grid(False)

    fig = ax.get_figure()

    if save_fig != None:
        fig.savefig(str(save_fig), dpi=dpi, bbox_inches=bbox_inches)

    plt.show()

"""Create scatterplot comparing gene expression of two genes across samples"""
def scatter(
    data,
    info,
    x,
    y,
    palette=None,
    add_linreg=False,
    order_legend=None,
    title=None,
    alpha=1,
    highlight_points=None,
    highlight_color='DarkRed',
    highlight_names=None,
    alpha_highlights=1,
    size=30,
    highlight_size=30,
    y_threshold=None,
    x_threshold=None,
    threshold_color='b',
    label_points=None,
    grid=False,
    whitegrid=False,
    figsize=(10,10),
    save_fig=None,
    dpi=600,
    bbox_inches='tight'):

    reset_plot(whitegrid)
    plt.figure(figsize=figsize)

    data_c = analysis_prep(data)

    #Prep data_scaled by adding labels from info
    labels = pd.Series(info[1].values,index=info[0]).to_dict()
    data_c.loc['label'] = data_c.columns.map(labels.get)

    if palette == None:
        ax = sns.scatterplot(data_c[str(x)], data_c[str(y)], color='black', alpha=alpha, s=size)
    else:
        # Transpose dataframe so whatever you want to plot is in columns
        if x in data.columns or y in data.columns:
            ax = sns.scatterplot(
                    data_c[str(x)],
                    data_c[str(y)],
                    hue = data_c['label'],
                    palette = palette,
                    alpha = alpha,
                    s = size)
        else:
            ax = sns.scatterplot(
                    data_c.loc[str(x)],
                    data_c.loc[str(y)],
                    hue = data_c.loc['label'],
                    palette = palette,
                    alpha = alpha,
                    s = size)

    if add_linreg == True:
        _x, _y, r_value, title = make_linreg(data_c, x, y)
        ax.plot(_x, _y, '-k')

    # Plot thresholds
    y_threshold = make_threshold_list(y_threshold)
    if type(y_threshold) != list and y_threshold != None or type(y_threshold) == list and any(x is None for x in y_threshold) == False:
        add_y_threshold(y_threshold, threshold_color, ax)

    x_threshold = make_threshold_list(x_threshold)
    if type(x_threshold) != list and x_threshold != None or type(x_threshold) == list and any(x is None for x in x_threshold) == False:
        add_x_threshold(x_threshold, threshold_color, ax)

    # Plot selected genes if user-specified
    if highlight_points != None:
        ax = highlight_markers(data_c, x, y, highlight_points, highlight_color, highlight_names, alpha_highlights, highlight_size, ax)

    # Label points
    if label_points != None and type(label_points) is dict:
        for key, value in label_points.items():
            ax.text(value[0], value[1], str(key), horizontalalignment='left', size='medium', color='black', weight='semibold')

    # Put the legend out of the figure
    ax = make_legend(ax, order_legend, highlight_color, highlight_names)

    plt.xlabel(str(x))
    plt.ylabel(str(y))

    if title != None:
        plt.title(str(title))

    if grid == False:
        plt.grid(False)

    if save_fig != None:
        plt.savefig(str(save_fig), dpi=dpi, bbox_inches=bbox_inches)

    plt.show()

"""Create scatterplot with r value and jointplot density distributions for axes"""
def jointplot(
    data, info, x, y,
    kind='reg',
    palette=None, order=None,
    title_pad=0, title_pos='right',
    grid=False, whitegrid=False, figsize=(10,10),
    save_fig=None, dpi=600, bbox_inches='tight'):

    reset_plot(whitegrid)
    plt.figure(figsize=figsize)

    data_c = data.copy()
    data_c = analysis_prep(data_c)

    # Prep data_scaled by adding labels from info
    if 'label' not in data_c.index:
        labels = pd.Series(info[1].values,index=info[0]).to_dict()
        data_c.loc['label'] = data_c.columns.map(labels.get).T

    # Get data as array
    gene_a = get_array(data_c, x)
    gene_b = get_array(data_c, y)

    # Plot
    ax = sns.jointplot(
            x = gene_a,
            y = gene_b,
            kind = kind)

    ax.ax_joint.collections[0].set_visible(False)

    if kind.lower() == 'reg':
        r_value = stats.pearsonr(gene_a, gene_b)[0]

        ax = sns.scatterplot(
                x = str(x),
                y = str(y),
                data = data_c.T,
                hue = 'label',
                palette = palette,
                hue_order = order)

        plt.xlabel(str(x))
        plt.ylabel(str(y))
        ax.set_title('r: ' + str(round(r_value,2)), pad=float(title_pad), loc=title_pos.lower())

    if grid == False:
        plt.grid(False)

    if kind.lower() == 'reg':
        fig = ax.get_figure()
    plt.show()

    if save_fig != None:
        plt.savefig(str(save_fig), dpi=dpi, bbox_inches=bbox_inches)

"""Calculates r, r^2 values, and p-values for every gene against target gene for given dataset"""
def linreg(
    data, gene_name,
    save_file, delimiter=','):

    data_c = data.copy()
    data_c = analysis_prep(data_c)

    if 'label' in data_c.index:
        data_c = data_c.drop(labels='label',axis=0)

    # Gene expression array for gene of interest
    interest = get_array(data_c, gene_name)

    # Get expression arrays for all genes and run linear regression, save values to matrix
    lm_interest=[]
    for row in data_c.iterrows():
        index, data = row
        gene = get_array(data_c, index)

        # For genes of interest, calculate pearson r and output stats to table
        if len(gene) is not len(interest):
            continue
        else:
            slope, intercept, r_value, p_value, std_err = linregress(gene, interest)
            lm_interest.append([index, slope, intercept, r_value, (r_value**2), p_value, std_err])

    df_lm_interest = pd.DataFrame(lm_interest, columns=['gene', 'slope', 'intercept', 'r_value', 'r_squared', 'p_value', 'std_err'])

    # Save linear modeling metrics to .csv
    df_lm_interest.to_csv(str(save_file), sep=delimiter)

"""
Provide DESeq2 output
@param data: DESeq2 table
@param highlight_names: list of gene names to highlight
@param label_points: list of gene names to label

"""
fold_change = 'log2FoldChange'
fdr = 'padj'
highlight_position = 6

def rna_volcano(
    file,
    title=None,
    order_legend=None,
    alpha=1,
    highlight_points=None,
    highlight_color='DarkRed',
    highlight_names=None,
    alpha_highlights=1,
    size=30,
    highlight_size=30,
    y_threshold=None,
    x_threshold=None,
    threshold_color='b',
    label_points=False,
    grid=False,
    whitegrid=False,
    figsize=(10,10),
    interactive=False,
    save_fig=None,
    dpi=600,
    bbox_inches='tight'):

    reset_plot(whitegrid)

    data = pd.read_csv(
        str(file),
        sep = '\t',
        index_col=0)

    data = data.fillna(1)
    data[fdr] = np.log10(data[fdr]) * -1

    # Plot all genes
    if interactive == False:
        data = data.rename(columns = {fold_change:'log$_2$(Fold Change)', fdr:'-log$_1$$_0$(FDR)'})
        plt.figure(figsize=figsize)
        ax = sns.scatterplot(
            data['log$_2$(Fold Change)'],
            data['-log$_1$$_0$(FDR)'],
            linewidth=0,
            palette = highlight_color,
            alpha = alpha,
            s = size,
            color = 'black')

        # remove duplicate labels for labeling
        # in most cases, seems these are rows with NAs
        data = data.loc[~data.index.duplicated()]

        # Plot thresholds
        y_threshold = make_threshold_list(y_threshold)
        if type(y_threshold) != list and y_threshold != None or type(y_threshold) == list and any(x is None for x in y_threshold) == False:
            add_y_threshold(y_threshold, threshold_color, ax)

        x_threshold = make_threshold_list(x_threshold)
        if type(x_threshold) != list and x_threshold != None or type(x_threshold) == list and any(x is None for x in x_threshold) == False:
            add_x_threshold(x_threshold, threshold_color, ax)

        # Plot selected genes if user-specified
        if highlight_points != None:
            ax = highlight_markers(data, 'log$_2$(Fold Change)', '-log$_1$$_0$(FDR)', highlight_points, highlight_color, highlight_names, alpha_highlights, highlight_size, ax)

        # Label points
        if label_points != None:
            if type(label_points) != list:
                label_points = [label_points]
            data_label = data.loc[label_points]

            for index, row in data_label.iterrows():
                ax.text(row[1] + 0.2, row[5] + 0, str(index), horizontalalignment='left', size='large', color=highlight_color, weight='semibold')

        # Put the legend out of the figure
        ax = make_legend(ax, order_legend, highlight_color, highlight_names)

        plt.xlabel('log$_2$(Fold Change)')
        plt.ylabel('-log$_1$$_0$(FDR)')

        if title != None:
            plt.title(str(title))

        if grid == False:
            plt.grid(False)

        if save_fig != None:
            plt.savefig(str(save_fig), dpi=dpi, bbox_inches=bbox_inches)

        plt.show()

    else:
        # Get genes column
        data['genes'] = data.index.tolist()

        # Make labels
        data['label'] = 'All'
        if all(isinstance(z, list) for z in highlight_points):
            p = 0
            for h in highlight_points:
                data.loc[h, 'label'] = highlight_names[p]
                p += 1

        else:
            data.loc[highlight_points, 'label'] = highlight_names[0]

        sc = px.scatter(
            data,
            x=fold_change,
            y=fdr,
            color='label',
            hover_name='genes',
            log_x=False,
            log_y=False,
            opacity=alpha,
            width=1400,
            height=1000,
            title=str(title))

        if save_fig != None:
            py.offline.plot(sc, filename=str(save_fig))





"""Plot volcano plot for dataframe, can highlight subset of genes"""
def volcano(
    data,
    info,
    label_comp,
    label_base,
    order_legend=None,
    title=None,
    alpha=1,
    highlight_points=None,
    highlight_color='DarkRed',
    highlight_names=None,
    alpha_highlights=1,
    size=30,
    y_threshold=None,
    x_threshold=None,
    threshold_color='b',
    save_threshold_hits=None,
    save_threshold_hits_delimiter=',',
    label_points=None,
    grid=False,
    whitegrid=False,
    figsize=(10,10),
    return_data=False,
    interactive=False,
    save_fig=None,
    dpi=600,
    bbox_inches='tight'):

    reset_plot(whitegrid)

    data_c = analysis_prep(data)
    info_c = info.copy()

    # Add labels to column name
    if 'label' in data_c.index:
        data_c = data_c.drop(labels='label',axis=0)

    # Convert sample names to include their category (ie: Sample1_Normal, Sample2_Treatment)
    info_c["id"] = info_c[0] + '_' + info_c[1]
    labels = pd.Series(info_c['id'].values,index=info_c[0]).to_dict()
    data_c = data_c.rename(labels, axis='columns')

    data_c += 1e-7 # Avoid divide by 0 error

    # Calculate parameters
    data_c = parallelize(calculate_fc, data_c, label_comp, label_base)
    data_c = parallelize(calculate_p, data_c, label_comp, label_base)

    # Plot all genes
    if interactive == False:
        scatter(
            data_c,
            info_c,
            'log$_2$(Fold Change)',
            '-log$_1$$_0$(P-Value)',
            palette = None,
            add_linreg = False,
            order_legend = order_legend,
            title = title,
            save_fig = save_fig,
            dpi = dpi,
            bbox_inches = bbox_inches,
            grid = grid,
            whitegrid = whitegrid,
            alpha = alpha,
            highlight_points = highlight_points,
            highlight_color = highlight_color,
            highlight_names = highlight_names,
            alpha_highlights = alpha_highlights,
            y_threshold = y_threshold,
            x_threshold = x_threshold,
            threshold_color = threshold_color,
            label_points = label_points,
            size = size,
            figsize = figsize)

    else:
        # Get genes column
        data_c['genes'] = data_c.index.tolist()

        # Make labels
        data_c['label'] = 'All'
        if all(isinstance(z, list) for z in highlight_points):
            p = 0
            for h in highlight_points:
                data_c.loc[h, 'label'] = highlight_names[p]
                p += 1

        else:
            data_c.loc[highlight_points, 'label'] = highlight_names

        sc = px.scatter(
            data_c,
            x = 'log$_2$(Fold Change)',
            y = '-log$_1$$_0$(P-Value)',
            color='label',
            hover_name='genes',
            log_x=False,
            log_y=False,
            opacity=alpha,
            width=1400,
            height=1000,
            title=str(title))

        if save_fig != None:
            py.offline.plot(sc, filename=str(save_fig))

    # Output threshold hits outer bounds, specific to volcano plot
    if save_threshold_hits != None:
        output_threshold(data_c, x_threshold, y_threshold, save_threshold_hits, save_threshold_hits_delimiter)

    # Save hits if user-specified
    if return_data == True:
        df_c = data_c[['log$_2$(Fold Change)', '-log$_1$$_0$(P-Value)']].copy()
        return df_c

"""Plot a 2-D PCA with confidence intervals or a 3-D PCA with no confidence intervals"""
def pca(
    data,
    info,
    palette,
    grouping='samples',
    gene_list=None,
    gene_labels=False,
    _3d_pca=False,
    ci=2,
    size=30,
    principle_components=[1,2],
    n_components=10,
    scree_only=False,
    save_scree=False,
    order_legend=None,
    title=None,
    grid=False,
    whitegrid=False,
    figsize=(10,10),
    save_fig=None,
    dpi=600,
    bbox_inches='tight',
    return_pca=False,
    plotly_login=None):

    highlight_color = None
    highlight_names = None

    reset_plot(whitegrid)
    principle_components = init_pca(principle_components, _3d_pca, plotly_login)

    data_c = data.copy()
    data_c = analysis_prep(data)

    # Add labels to column name
    if 'label' in data_c.index:
        data_c = data_c.drop(labels='label',axis=0)

    # Custom panel heatmap
    if gene_list != None:
        data_sub = data_subset(data_c, gene_list)
        plot_data = data_sub.dropna(axis=0)
    else:
        plot_data = data_c.dropna(axis=0)

    # Format data for sample-wise or gene-wise PCA analysis
    if str(grouping) == 'samples':
        plot_data = plot_data.T
    elif str(grouping) == 'genes':
        plot_data = plot_data
        print('PCA by genes is not yet implemented...')
        return
    else:
        print('Incorrect grouping variable provided')
        return

    # Prep PCA
    pca = PCA(n_components=n_components)
    pca_result = pca.fit_transform(plot_data.values)

    # Record PCs
    x = 1
    while x <= n_components:
        pc = 'PC' + str(x)
        component = x - 1
        plot_data[pc] = pca_result[:,component]
        x += 1

    # Scree
    if save_scree == True or scree_only == True:
        scree = make_scree(
                    pca,
                    n_components,
                    save_fig,
                    dpi,
                    bbox_inches,
                    scree_only,
                    save_scree,
                    grid,
                    whitegrid)

        if scree_only == True:
            return
    else:
        scree = np.round(pca.explained_variance_ratio_, decimals=4)*100

    if gene_labels == False:

        # Prep data_scaled by adding labels from info to column (samples are rows)
        labels = pd.Series(info[1].values,index=info[0]).to_dict()
        plot_data['label'] = plot_data.index.to_series().map(labels)

        # Plot PCa & PCb
        pc_list = []
        for p in principle_components:
            pc_list.append('PC' + str(p))

        pc_list.append('label')

        df_pca = plot_data[pc_list] #Prepare pca data

        # 2D PCA
        if _3d_pca == False:
            df_pca.columns = ['PCa', 'PCb', 'label']
            unique_labels = df_pca['label'].unique() #Gather unique labels

            # Non-interactive
            if plotly_login == None:
                pca2(
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
                    highlight_names)

            # Plotly
            else:
                interactive_scatter(
                    df_pca,
                    info,
                    'PCa',
                    'PCb',
                    plotly_login,
                    palette,
                    highlight = 'sample',
                    save_fig = save_fig)

        # 3D PCA
        elif _3d_pca == True:

            #Plotly
            if plotly_login != None:
                interactive_3D(
                    df_pca,
                    size,
                    palette,
                    plotly_login,
                    save_fig)

            # Non-interactive
            else:
                pca3(
                    df_pca,
                    palette,
                    save_fig,
                    dpi,
                    bbox_inches,
                    size,
                    order_legend)

        else:
            return

    # Gene scatter
    else:
        print('This feature has not been implemented yet')

    if return_pca == True:
        return df_pca

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
from scipy.stats import linregress
from sklearn.decomposition import PCA
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
sns.set(font='arial')
jakes_cmap = sns.diverging_palette(212, 61, s=99, l=77, sep=1, n=16, center='dark') #Custom aesthetics
from .utils_analyze import *
from .interactive import *

"""
INITIALIZATION PARAMETERS
"""
#Retrieve path for scripts used in this pipeline, appended to argument dictionary for every function
__path__, xpresstools_arguments = os.path.split(__file__)

"""
DESCRIPTION: Run Differential Expression analysis using DESeq2
"""
def diff_xpress(input_file, sample_file, equation):

    #Get output file name
    if input_file.endswith('.txt') or input_file.endswith('.tsv'):
        output_file = str(input_file[:-4]) + '_diffx.tsv'
    else:
        raise Exception('Unrecognized input_file delimiter type. Files must be tab-delimited')

    if sample_file.endswith('.txt') or sample_file.endswith('.tsv'):
        pass
    else:
        raise Exception('Unrecognized sample_file delimiter type. Files must be tab-delimited')

    if equation.startswith('~'):
        raise Exception('Tilde should not be included in design formula, script will automatically add this syntax.')

    #Run deseq2 in R
    os.system('rscript ' + str(__path__) + '/diffxpress.r ' + str(input_file) + ' ' + str(sample_file) + ' ' + str(output_file) + ' ' + str(equation))

"""
DESCRIPTION: Plots heatmap for prep_data formatted data

METHODS: We use centroid method as default as it looks at the gene set for each sample holistically to determine a 'centroid' then used for clustering. Other valid options are available (see VARIABLES)

VARIABLES:
data_scaled= Scaled (generally) dataframe as created with the XPRESStools prep_data function
info= XPRESStools formatted info matrix
sample_palette= Dictionary of labels and colors for plotting, or valid seaborns.clustermap col_colors option
gene_info= XPRESStools formatted info matrix for genes (column0) and groupings (column1)
gene_palette= Dictionary of labels and colors for plotting, or valid seaborns.clustermap col_colors option
gene_list= List of genes (either as list or as .csv file path and name with list of genes) to plot (IMPORTANT: Gene names are case-sensitive)
save_fig= If not None, provide full file path, name, and extension to save the file as
dpi= Set dpi of saved figure
bbox_inches= Format saved figure (often useful for making sure no text is cut off)
font_scale= Scaling factor for font
cmap= A valid seaborns.clustermap cmap option (default: Rutter Lab colorblind-friendly scale)
Please see seaborns.heatmap documentation for descriptions of the following options:
center, metric, method, xticklabels, linewidths, linecolor, col_cluster, row_cluster, figsize
figsize= (Width, Height)

USAGE:
import xpresstools as xp
col_colors = {'adenocarcinoma': 'grey',
        'adenoma': (0.8705882352941177, 0.5607843137254902, 0.0196078431372549),
        'normal': (0.00784313725490196, 0.6196078431372549, 0.45098039215686275)}
cbar_kws= {'label': 'z-score',"shrink": 0.2, 'aspect': 10}

xp.heatmap(data_scaled, data_labeled, color_dict=color_dict, gene_list='/path/to/gene_list.csv', figsize=(40,8))

ASSUMPTIONS:
Data has been scaled and labeled with the XPRESStools prep_data function
"""
def heatmap(data, info, sample_palette=None, gene_info=None, gene_palette=None, gene_list=None, col_cluster=True, row_cluster=False, metric='euclidean', method='centroid', font_scale=.8, cmap=jakes_cmap, center=0, xticklabels=True, yticklabels=True, linewidths=0, linecolor='#DCDCDC', cbar_kws=None, figsize=(16,6.5), save_fig=None, dpi=600, bbox_inches='tight'):

    reset_plot(True)
    data_c = analysis_prep(data)

    #Set colors for plotting
    if sample_palette != None:
        sample_palette = prep_palette(info, sample_palette)

    if gene_palette != None and isinstance(gene_info, pd.DataFrame):
        gene_palette = prep_palette(gene_info, gene_palette)

    #Custom panel heatmap
    if gene_list != None:
        data_sub = data_subset(data_c, gene_list)
        plot_data = data_sub.dropna(axis=0)
    else:
        plot_data = data_c.dropna(axis=0)

    #Generate clustermap
    sns.set(font_scale=float(font_scale))
    ax = sns.clustermap(plot_data,
                    cmap=cmap,
                    center=float(center),
                    metric=str(metric),
                    method=str(method),
                    xticklabels=xticklabels,
                    yticklabels=yticklabels,
                    linewidths=float(linewidths),
                    linecolor=linecolor,
                    col_cluster=col_cluster,
                    row_cluster=row_cluster,
                    col_colors=sample_palette,
                    row_colors=gene_palette,
                    figsize=figsize,
                    cbar_kws=cbar_kws
                   )

    #Save figure
    if save_fig is not None:
        plt.savefig(str(save_fig), dpi=dpi, bbox_inches=str(bbox_inches))

"""
DESCRIPTION: Create violin plots of subset of gene expressions or all gene expression by sample

VARIABLES:
data_labeled= Unscaled dataframe with sample labels as created with the XPRESStools prep_data function
gene_list= List of genes (either as list or as .csv file path and name with list of genes) to plot (IMPORTANT: Gene names are case-sensitive)
order= List of samples in order to plot
palette= Dictionary of matplotlib compatible colors for samples
save_fig= If not None, provide full file path, name, and extension to save the file as
dpi= Set dpi of saved figure
bbox_inches= Format saved figure (often useful for making sure no text is cut off)
title= Provide title for figure and saved file if save_fig option used

ASSUMPTIONS:
Data has been scaled and labeled with the XPRESStools prep_data function
"""
def multigene_overview(data, info, palette=None, gene_list=None, order=None, scale='area', save_fig=None, dpi=600, bbox_inches='tight', title=None, grid=False, whitegrid=False):

    reset_plot(whitegrid)
    data_c = data.copy()

    #Custom panel heatmap
    if gene_list != None:
        data_sub = data_subset(data_c, gene_list)
        plot_data = data_sub.dropna(axis=0)
    else:
        plot_data = data_c.dropna(axis=0)

    plot_data_c = analysis_prep(plot_data)
    plot_data_c = label(plot_data_c, info)
    plot_data_c = unstack_data(plot_data_c)

    #Final formatting check
    try:
        plot_data_c[["expr"]] = plot_data_c[["expr"]].apply(pd.to_numeric)
    except:
        raise Exception('Data is not properly formatted for downstream plotting')

    #Plot
    if not order != None and type(order) is list or palette != None and type(palette) is list:
        return

    ax = sns.violinplot(x=plot_data_c['type'], y=plot_data_c['expr'], data=plot_data_c, order=order, palette=palette, scale=scale)
    ax.set_xlabel('')
    ax.set_ylabel('Expression')

    if grid == False:
        ax.grid(False)

    #Save fig
    if save_fig is not None:
        if title is None:
            plt.savefig(str(save_fig), dpi=dpi, bbox_inches=bbox_inches)
        else:
            ax.set_title(str(title))
            plt.savefig(str(save_fig), dpi=dpi, bbox_inches=bbox_inches)

    plt.show()

"""
DESCRIPTION: Plot boxplot overlaid with swarmplot of each sample type's gene expression for the given gene

VARIABLES:
data= Dataframe (can be sample-normalized, prep_data() scaled, or log_scale() scaled)
info= XPRESStools formatted sample info dataframe
gene_name= Name of gene to plot
palette= Dictionary of matplotlib compatible colors for samples
order= List of samples in order to plot
save_fig= If not None, provide full file path, name, and extension to save the file as
dpi= Set dpi of saved figure
bbox_inches= Format saved figure (often useful for making sure no text is cut off)
grid= Control plot gridlines (default: False)

ASSUMPTIONS:
data and info dataframes are properly formatted for XPRESStools and any appropriate sample/gene normalizations have been performed
"""
def gene_overview(data, info, gene_name, palette, order=None, grid=False, whitegrid=False, save_fig=None, dpi=600, bbox_inches='tight'):

    reset_plot(whitegrid)
    data_c = analysis_prep(data).T

    #Prep data_scaled by adding labels from info to column (samples are rows)
    labels = pd.Series(info[1].values,index=info[0]).to_dict()
    data_c['label'] = data_c.index.to_series().map(labels)
    gene_df = data_c[[str(gene_name), 'label']]
    gene_df[[str(gene_name)]] = gene_df[[str(gene_name)]].apply(pd.to_numeric, errors='coerce')

    ax = sns.catplot(x='label', y=str(gene_name), data=gene_df, color='black', order=order, kind='swarm') #Swarm plot
    ax = sns.boxplot(x='label', y=str(gene_name), data=gene_df, width =0.3, fliersize=0, order=order, palette=palette) #Boxplot, fliersize=0 removes outlier diamonds from sns
    plt.setp(ax.collections, sizes=[12]) #Resize markers for catplot

    if grid == False:
        ax.grid(False)

    fig = ax.get_figure()

    if save_fig != None:
        fig.savefig(str(save_fig), dpi=dpi, bbox_inches=bbox_inches)

    plt.show()

"""
DESCRIPTION: Create scatterplot comparing gene expression of two genes across samples

VARIABLES:
data= Dataframe of expression data (samples are columns, genes are rows)
info= XPRESStools formatted sample info dataframe
gene1= Name of gene to plot on x_axis
gene1= Name of gene to plot on y_axis
palette= Dictionary of matplotlib compatible colors for samples
add_linreg= Add linear regression line and print r value in plot title
order_legend= Order of samples to display in legend
title= Plot title to display; if add_linreg used, title variable is void
save_fig= If not None, provide full file path, name, and extension to save the file as
dpi= Set dpi of saved figure
bbox_inches= Format saved figure (often useful for making sure no text is cut off)
grid= Control plot gridlines (default: False)
whitegrid= Use whitegrid background in plot
alpha= Control opacity of points on plot (useful for getting an idea of density in large datasets)

highlight_points= If plotting multiple highlight colors, make sure this and highlight_colors are lists of lists

thresholds can be list of ints or int or None

ASSUMPTIONS:
XPRESStools formatted data and info dataframes, palette is a dictionary of labels and colors to plot points with
"""
def scatter(data, info, x, y, palette=None, add_linreg=False, order_legend=None, title=None, alpha=1, highlight_points=None, highlight_color='DarkRed', highlight_names=None, alpha_highlights=1, size=30, y_threshold=None, x_threshold=None, threshold_color='b', label_points=None, grid=False, whitegrid=False, save_fig=None, dpi=600, bbox_inches='tight', figsize=(10,10)):

    reset_plot(whitegrid)
    plt.figure(figsize=figsize)

    data_c = analysis_prep(data)

    #Prep data_scaled by adding labels from info
    labels = pd.Series(info[1].values,index=info[0]).to_dict()
    data_c.loc['label'] = data_c.columns.map(labels.get)

    if palette == None:
        ax = sns.scatterplot(data_c[str(x)], data_c[str(y)], color='black', alpha=alpha, s=size)
    else:
        #Transpose dataframe so whatever you want to plot is in columns
        if x in data.columns or y in data.columns:
            ax = sns.scatterplot(data_c[str(x)], data_c[str(y)], hue=data_c['label'], palette=palette, alpha=alpha, s=size)
        else:
            ax = sns.scatterplot(data_c.loc[str(x)], data_c.loc[str(y)], hue=data_c.loc['label'], palette=palette, alpha=alpha, s=size)

    if add_linreg == True:
        _x, _y, r_value, title = make_linreg(data_c, x, y)
        ax.plot(_x, _y, '-k')

    #Plot thresholds
    y_threshold = make_threshold_list(y_threshold)
    if type(y_threshold) != list and y_threshold != None or type(y_threshold) == list and any(x is None for x in y_threshold) == False:
        add_y_threshold(y_threshold, threshold_color, ax)

    x_threshold = make_threshold_list(x_threshold)
    if type(x_threshold) != list and x_threshold != None or type(x_threshold) == list and any(x is None for x in x_threshold) == False:
        add_x_threshold(x_threshold, threshold_color, ax)

    #Plot selected genes if user-specified
    if highlight_points != None:
        ax = highlight_markers(data_c, x, y, highlight_points, highlight_color, highlight_names, alpha_highlights, ax)

    #Label points
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

"""
DESCRIPTION: Create scatterplot with r value and jointplot density distributions for axes

VARIABLES:
data= Dataframe of expression data (samples are columns, genes are rows)
info= XPRESStools formatted sample info dataframe
gene1= Name of gene to plot on x_axis
gene1= Name of gene to plot on y_axis
kind= See seaborn jointplot documentation for options (https://seaborn.pydata.org/generated/seaborn.jointplot.html)
palette= Dictionary of matplotlib compatible colors for samples
order= Order of samples to display in plot
save_fig= If not None, provide full file path, name, and extension to save the file as
dpi= Set dpi of saved figure
bbox_inches= Format saved figure (often useful for making sure no text is cut off)
grid= Control plot gridlines (default: False)
whitegrid= Use whitegrid background in plot

Full functionality only supported for kind='reg' currently

ASSUMPTIONS:
XPRESStools formatted data and info dataframes, palette (if used) is a dictionary of labels and colors to plot points with
"""
def jointplot(data, info, x, y, kind='reg', palette=None, order=None, title_pad=0, title_pos='right', grid=False, whitegrid=False, save_fig=None, dpi=600, bbox_inches='tight', figsize=(10,10)):

    reset_plot(whitegrid)
    plt.figure(figsize=figsize)

    data_c = data.copy()
    data_c = analysis_prep(data_c)

    #Prep data_scaled by adding labels from info
    if 'label' not in data_c.index:
        labels = pd.Series(info[1].values,index=info[0]).to_dict()
        data_c.loc['label'] = data_c.columns.map(labels.get).T

    #Get r
    gene_a = get_array(data_c, x)
    gene_b = get_array(data_c, y)

    #Plot
    ax = sns.jointplot(x=gene_a, y=gene_b, kind=kind)
    ax.ax_joint.collections[0].set_visible(False)

    if kind.lower() == 'reg':
        r_value = stats.pearsonr(gene_a, gene_b)[0]

        ax = sns.scatterplot(x=str(x), y=str(y), data=data_c.T, hue='label', palette=palette, hue_order=order)

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

"""
DESCRIPTION: Calculates r, r^2 values, and p-values for every gene against target gene for given dataset

VARIABLES:
data= Dataframe of expression data (samples are columns, genes are rows)
gene_name= Gene to run correlations against
save_file= File path and name to save correlations matrix with
delimiter= Delimiter to use for save_file

ASSUMPTIONS:
data dataframes is properly formatted for XPRESStools, sample scaling required but gene scaling is not necessary
data should not contain labels
"""
def linreg(data, gene_name, save_file, delimiter=','):

    data_c = data.copy()
    data_c = analysis_prep(data_c)

    if 'label' in data_c.index:
        data_c = data_c.drop(labels='label',axis=0)

    #Gene expression array for gene of interest
    interest = get_array(data_c, gene_name)

    #Get expression arrays for all genes and run linear regression, save values to matrix
    lm_interest=[]
    for row in data_c.iterrows():
        index, data = row
        gene = get_array(data_c, index)

        if len(gene) is not len(interest):
            continue
        else:
            slope, intercept, r_value, p_value, std_err = linregress(gene, interest)
            lm_interest.append([index, slope, intercept, r_value, (r_value**2), p_value, std_err])

    df_lm_interest = pd.DataFrame(lm_interest, columns=['gene', 'slope', 'intercept', 'r_value', 'r_squared', 'p_value', 'std_err'])

    # Save linear modeling metrics to .csv
    df_lm_interest.to_csv(str(save_file), sep=delimiter)

"""
DESCRIPTION: Plot volcano plot for dataframe, can highlight subset of genes

VARIABLES:
data= Sample-normalized, XPRESStools-formatted data -- should NOT be gene-normalized
info= XPRESStools formatted sample info dataframe
label_comp= Sample label to compare against base
label_base= Sample label to use as base for comparison
highlight_points= If provided with a list, or a file path and name to a .csv-type series of gene names, will highlight genes of interest in a different color (Gene names are case-sensitive). Must be a list of lists to allow for multiple highlighting groups
highlight_color= Color to use for highlighted genes
y_threshold= -log10(P-Value) threshold to use to identify significant hits. Will create a dotted line on the plot and use to pull significant hits if return_threshold_hits is not None
x_threshold= log2(Fold Change) threshold to use to identify significant hits. Will take the positive and negative of the number. Will create a dotted line on the plot and use to pull significant hits if return_threshold_hits is not None
save_threshold_hits= Will save a .csv-type matrix of significant genes/hits to the file path and name specified
save_threshold_hits_delimiter= Delimiter to use in exporting return_threshold_hits
save_fig= If not None, provide full file path, name, and extension to save the file as
dpi= Set dpi of saved figure
bbox_inches= Format saved figure (often useful for making sure no text is cut off)
return_data= Return dataframe with PCs (will not print plot to stout)
label_genes= Dictionary for add labels to non-plotly volcano plot

Need to output threshold to figure out coordinates to label

LIMITATIONS:
Will only perform comparison between 2 sample types

ASSUMPTIONS:
data should ONLY be sample normalized. If using a previous function that returned a modified original
y_threshold must be a postive integer or float
"""
def volcano(data, info, label_comp, label_base, order_legend=None, title=None, alpha=1, highlight_points=None, highlight_color='DarkRed', highlight_names=None, alpha_highlights=1, size=30, y_threshold=None, x_threshold=None, threshold_color='b', save_threshold_hits=None, save_threshold_hits_delimiter=',', label_points=None, grid=False, whitegrid=False, return_data=False, plotly_login=False, save_fig=None, dpi=600, bbox_inches='tight', figsize=(10,10)):

    reset_plot(whitegrid)

    data_c = analysis_prep(data)
    info_c = info.copy()

    #Add labels to column name
    if 'label' in data_c.index:
        data_c = data_c.drop(labels='label',axis=0)

    #Convert sample names to include their category (ie: Sample1_Normal, Sample2_Treatment)
    info_c["id"] = info_c[0] + '_' + info_c[1]
    labels = pd.Series(info_c['id'].values,index=info_c[0]).to_dict()
    data_c = data_c.rename(labels, axis='columns')

    data_c += 1e-7 #Avoid divide by 0 error

    #Calculate parameters
    data_c = parallelize(calculate_fc, data_c, label_comp, label_base)
    data_c = parallelize(calculate_p, data_c, label_comp, label_base)

    #Plot all genes
    if plotly_login == False:
        scatter(data_c, info_c, 'log$_2$(Fold Change)', '-log$_1$$_0$(P-Value)', palette=None, add_linreg=False, order_legend=order_legend, title=title, save_fig=save_fig, dpi=dpi, bbox_inches=bbox_inches, grid=grid, whitegrid=whitegrid, alpha=alpha, highlight_points=highlight_points, highlight_color=highlight_color, highlight_names=highlight_names,  alpha_highlights=alpha_highlights, y_threshold=y_threshold, x_threshold=x_threshold, threshold_color=threshold_color, label_points=label_points, size=size, figsize=figsize)

    else:
        data_c = data_c.rename(columns = {'log$_2$(Fold Change)':'log2 Fold Change', '-log$_1$$_0$(P-Value)':'-log10 P-Value'})
        interactive_scatter(data_c, info_c, x='log2 Fold Change', y='-log10 P-Value', plotly_login=plotly_login, file_name=save_fig, highlight='sample', palette=None)

    #Output threshold hits outer bounds, specific to volcano plot
    if save_threshold_hits != None:
        output_threshold(data_c, x_threshold, y_threshold, save_threshold_hits, save_threshold_hits_delimiter)

    #Save hits if user-specified
    if return_data == True:
        df_c = data_c[['log2 Fold Change', '-log10 P-Value']].copy()
        return df_c

"""
DESCRIPTION: Plot a 2-D PCA with confidence intervals or a 3-D PCA with no confidence intervals

RETURNS: Dataframe with PCs calculated

METHODS: 3-D PCA -- option to output as interactive plot by providing plotly credentials, or as a static plot without these credentials (default)

VARIABLES:
data_scaled= Scaled dataframe as created with the XPRESStools prep_data and scaling function -- can be a dataframe prepared using the prep_data() function or log_scale() function
info= XPRESStools formatted sample info dataframe
palette= Dictionary of matplotlib compatible colors for samples; for plotly 3-D PCA,
grouping= Perform PCA sample-wise (default) or gene-wise (grouping='genes')
gene_list= List of genes (either as list or as .csv file path and name with list of genes) to plot (IMPORTANT: Gene names are case-sensitive) (only functional when grouping='samples')
gene_labels= Option for grouping='genes', not currently implemented
ci= For 2-D PCA, confidence interval to display for each sample type (i.e. 1 == CI1 == 68%, 2 == CI2 == 95%, 3 == CI3 == 99%)
principle_components= Principle components to plot (default: [1,2] for 2-D PCA, or [1,2,3] for 3-D PCA)
n_components= Number of components to calculate in dataframe (more applicable if you want to perform a deeper survey of the principle components)
_3d_pca= Plot 3 principle components (default: False)
plotly_login= ['userid','api key'], usage of this option creates a plotly interactive plot that can be viewed online using your login credentials
title= Provide title for figure and saved file if save_fig option used
save_fig= If not None, provide full file path, name, and extension to save the file as. if saving scree using this, don't need extension
dpi= Set dpi of saved figure
bbox_inches= Format saved figure (often useful for making sure no text is cut off)
order_legend= List of integers to reorder samples in legend (i.e. if samples are displayed 1:Sample_A, 2:Sample_C, 3:Sample_B, provide the list [1,3,2]) (Not currently compatible with 3-D PCA options)
grid= For non-plotly options, remove gridlines from plot
fig_size= Option not used in function currently
size= Marker size

palette examples:
colors = {'adenocarcinoma': (0.5725490196078431, 0.5843137254901961, 0.5686274509803921),
        'adenoma': (0.8705882352941177, 0.5607843137254902, 0.0196078431372549),
        'normal': (0.00784313725490196, 0.6196078431372549, 0.45098039215686275)}
colors = {'adenocarcinoma': 'gray',
        'adenoma': 'orange',
        'normal': 'green'}
colors = {'adenocarcinoma': '#59656d',
        'adenoma': '#f0833a',
        'normal': '#0a5f38'}

USAGE:


ASSUMPTIONS:
Plotly and api key generation, addition
Data has been scaled and labeled with the XPRESStools prep_data function

FEATURES TO ADD:
Allow for compatibility with adding labels for gene classes for plotting
Option to order legend and rearrange in 3d static plot
Add options to vary marker size and opacity
"""
def pca(data, info, palette, grouping='samples', gene_list=None, gene_labels=False, _3d_pca=False, principle_components=[1,2], n_components=10, ci=2, scree_only=False, save_scree=False, size=30, order_legend=None, title=None, figsize=(10,10), grid=False, whitegrid=False, save_fig=None, dpi=600, bbox_inches='tight', return_pca=False, plotly_login=None):

    highlight_color = None
    highlight_names = None

    reset_plot(whitegrid)
    principle_components = init_pca(principle_components, _3d_pca, plotly_login)

    data_c = data.copy()
    data_c = analysis_prep(data)

    #Add labels to column name
    if 'label' in data_c.index:
        data_c = data_c.drop(labels='label',axis=0)

    #Custom panel heatmap
    if gene_list != None:
        data_sub = data_subset(data_c, gene_list)
        plot_data = data_sub.dropna(axis=0)
    else:
        plot_data = data_c.dropna(axis=0)

    #Format data for sample-wise or gene-wise PCA analysis
    if str(grouping) == 'samples':
        plot_data = plot_data.T
    elif str(grouping) == 'genes':
        plot_data = plot_data
        print('PCA by genes is not yet implemented...')
        return
    else:
        print('Incorrect grouping variable provided')
        return

    #Prep PCA
    pca = PCA(n_components=n_components)
    pca_result = pca.fit_transform(plot_data.values)

    #Record PCs
    x = 1
    while x <= n_components:
        pc = 'PC' + str(x)
        component = x - 1
        plot_data[pc] = pca_result[:,component]
        x += 1

    #Scree
    if save_scree == True or scree_only == True:
        scree = make_scree(pca, n_components, save_fig, dpi, bbox_inches, scree_only, save_scree, grid, whitegrid)
        if scree_only == True:
            return
    else:
        scree = np.round(pca.explained_variance_ratio_, decimals=4)*100

    if gene_labels == False:

        #Prep data_scaled by adding labels from info to column (samples are rows)
        labels = pd.Series(info[1].values,index=info[0]).to_dict()
        plot_data['label'] = plot_data.index.to_series().map(labels)

        #Plot PCa & PCb
        pc_list = []
        for p in principle_components:
            pc_list.append('PC' + str(p))

        pc_list.append('label')

        df_pca = plot_data[pc_list] #Prepare pca data

        #2D PCA
        if _3d_pca == False:
            df_pca.columns = ['PCa', 'PCb', 'label']
            unique_labels = df_pca['label'].unique() #Gather unique labels

            #Non-interactive
            if plotly_login == None:
                pca2(df_pca, unique_labels, palette, principle_components, scree, order_legend, save_fig, dpi, bbox_inches, ci, grid, title, size, highlight_color, highlight_names)

            #Plotly
            else:
                interactive_scatter(df_pca, info, 'PCa', 'PCb', plotly_login, palette, highlight='sample', save_fig=save_fig)

        #3D PCA
        elif _3d_pca == True:

            #Plotly
            if plotly_login != None:
                interactive_3D(df_pca, size, palette, plotly_login, save_fig)

            #Non-interactive
            else:
                pca3(df_pca, palette, save_fig, dpi, bbox_inches, size, order_legend)

        else:
            return

    #Gene scatter
    else:
        print('This feature has not been implemented yet')

    if return_pca == True:
        return df_pca

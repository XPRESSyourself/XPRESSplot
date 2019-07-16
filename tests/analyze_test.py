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

"""
IMPORT DEPENDENCIES AND DATASETS
"""
import pandas as pd
import xpressplot as xp
import numpy as np

import os
__path__  =  os.path.dirname(os.path.realpath(__file__)) + '/'
#__path__ = '/Users/jordan/scripts/XPRESSyourself/XPRESSplot/tests/'

data_loc = str(__path__ ) + 'large_test.csv'
meta_loc = str(__path__ ) + 'sample_info_test.csv'
save_threshold = str(__path__ ) + 'threshold.csv'
linreg_file = str(__path__ ) + 'linreg_results.csv'
pca_file = str(__path__ ) + 'pca_test.pdf'
probe_loc = str(__path__ ) + 'GPL570.txt'

"""
%matplotlib inline
data_loc = '/Users/jordan/scripts/XPRESSyourself/XPRESSplot/tests/large_test.csv'
meta_loc = '/Users/jordan/scripts/XPRESSyourself/XPRESSplot/tests/sample_info_test.csv'
save_threshold = '/Users/jordan/scripts/XPRESSyourself/XPRESSplot/tests/threshold.csv'
linreg_file = '/Users/jordan/scripts/XPRESSyourself/XPRESSplot/tests/linreg_results.csv'
pca_file = '/Users/jordan/scripts/XPRESSyourself/XPRESSplot/tests/pca_test.pdf'
probe_loc = '/Users/jordan/scripts/XPRESSyourself/XPRESSplot/tests/GPL570.txt'
"""

"""
Get test data
"""
geo = xp.get_df(data_loc, delimiter=',')
meta = xp.get_info(meta_loc, delimiter=',')

meta[1] = meta[1].str.capitalize() #Make sample types look nice
meta = meta.replace('Normal_colon', 'Normal')
geo = xp.keep_labels(geo, meta, label_list=['Normal','Adenoma','Adenocarcinoma'])
geo_clean = xp.clean_df(geo)

#Collapse multi-mapping probes
geo_collapsed = xp.probe_collapse(geo_clean, probe_loc)

#Scale sorted dataset
geo_scaled, geo_labeled = xp.prep_data(geo_collapsed, meta)

geo_colors = {'Adenocarcinoma': (0.5725490196078431, 0.5843137254901961, 0.5686274509803921),
        'Adenoma': (0.8705882352941177, 0.5607843137254902, 0.0196078431372549),
        'Normal': (0.00784313725490196, 0.6196078431372549, 0.45098039215686275)}

gene_colors = {'Group1': (0.5725490196078431, 0.5843137254901961, 0.5686274509803921),
        'Group2': (0.8705882352941177, 0.5607843137254902, 0.0196078431372549)}

"""
NON-INTERACTIVE TESTS
"""
#Single-gene summary
xp.gene_overview(geo_labeled, meta, gene_name='SEC62', palette=geo_colors, order=['Normal','Adenoma','Adenocarcinoma'])
xp.gene_overview(geo_labeled, meta, 'CCL5', geo_colors, order=None, dpi=600, bbox_inches='tight', grid=True, whitegrid=True)

#Multi-gene summary
"""/anaconda3/lib/python3.6/site-packages/scipy/stats/stats.py:1713: FutureWarning:

Using a non-tuple sequence for multidimensional indexing is deprecated; use `arr[tuple(seq)]` instead of `arr[seq]`. In the future this will be interpreted as an array index, `arr[np.array(seq)]`, which will result either in an error or a different result.

This appears to be fine and not of concern (https://stackoverflow.com/a/52595447/9571488)
"""

xp.multigene_overview(geo_labeled, meta, palette=geo_colors, gene_list=['SEC62','CCL5','STX6'], order=None, dpi=600, bbox_inches='tight', title=None, grid=False, whitegrid=False)

xp.multigene_overview(geo_labeled, meta, palette=geo_colors, gene_list=['STX6'], order=['Normal','Adenoma','Adenocarcinoma'], dpi=600, bbox_inches='tight', title=None, grid=False, whitegrid=False)

xp.multigene_overview(geo_labeled, meta, palette=geo_colors, order=['Normal','Adenoma','Adenocarcinoma'], dpi=600, bbox_inches='tight', title=None, grid=False, whitegrid=False)

try:
    xp.multigene_overview(geo_labeled, palette=geo_colors, order=['Normal','Adenoma','Adenocarcinoma'], save_fig=None, dpi=600, bbox_inches='tight', title=None, grid=False, whitegrid=False)
except:
    pass
else:
    fail('multigene_overview() failed to catch an error in not providing metadata')

#Heatmap
xp.heatmap(geo_scaled, meta, sample_palette=geo_colors, gene_list=['SEC62','STX6','CCL5'], cbar_kws={'label':'z-score'}, figsize=(20,2))

genes = np.array([['SEC62','Group1'],['STX6','Group2'],['CCL5','Group1']])
gene_info2 = pd.DataFrame({0:genes[:,0],1:genes[:,1]})
xp.heatmap(geo_scaled, meta, sample_palette=geo_colors, gene_palette=gene_colors, gene_info=gene_info2, gene_list=['SEC62','STX6','CCL5'], figsize=(20,2), row_cluster=True)

try:
    xp.heatmap(geo_scaled, meta, sample_palette=geo_colors, gene_palette=gene_colors, xticklabels=True, linewidths=.5, linecolor='red', gene_list=['SEC62','STX6','CCL5'], figsize=(10,2))
except:
    pass
else:
    fail('heatmap() failed to catch error when only providing part of gene labeling variables')

xp.heatmap(geo_scaled, meta, sample_palette=geo_colors, xticklabels=True, linewidths=.5, linecolor='black', gene_list=['SEC62','STX6','CCL5'], figsize=(20,2))

#Scatterplot basic
xp.scatter(geo_labeled, meta, 'SEC62', 'STX6', palette=geo_colors, add_linreg=True, order_legend=[1,3,2], alpha=.7)

try:
    xp.scatter(geo_labeled, meta, 'SEC62', palette=geo_colors, add_linreg=True, order_legend=[1,3,2], alpha=.7)
except:
    pass
else:
    fail('scatter() failed to catch error when only providing one gene/axis to plot')

xp.scatter(geo_labeled, meta, 'SEC62', 'STX6', palette=geo_colors, add_linreg=False, alpha=.7)

try:
    xp.scatter(geo_labeled, meta, 'SEC62', 'STX6', highlight_points=['Normal'], highlight_color='DarkRed', palette=geo_colors, add_linreg=False, alpha=.7)
except:
    pass
else:
    fail('scatter() failed to catch error when only providing one gene/axis to plot')

xp.scatter(geo_labeled, meta, 'SEC62', 'STX6', palette=geo_colors, add_linreg=True, alpha=.2, title='this is a title', y_threshold=5, x_threshold=[7])

#Volcano Plot
xp.volcano(geo_labeled, meta, 'Adenoma', 'Normal', highlight_points=['STX6','SCARB1','CCL5'])

xp.volcano(geo_labeled, meta, 'Adenoma', 'Normal', highlight_points=['STX6','SCARB1','CCL5'], y_threshold=2, x_threshold=[-1,1], save_threshold_hits=save_threshold)

xp.volcano(geo_labeled, meta, 'Adenoma', 'Normal', highlight_points=['STX6','SCARB1','CCL5'], highlight_names=['gene_set'], y_threshold=2, x_threshold=[-1,1], save_threshold_hits=save_threshold)

xp.volcano(geo_labeled, meta, 'Adenoma', 'Normal', highlight_points=[['STX6','SCARB1','CCL5'],['BEST4']], highlight_color=['blue','red'], alpha=.3, y_threshold=2, x_threshold=[-1,1], label_points={'BEST4':[-1.24288077425345,21.782377963035827]})

try:
    xp.volcano(geo_labeled, meta, 'Adenoma', highlight_points=[['STX6','SCARB1','CCL5'],['BEST4']], highlight_color=['blue','red'], alpha=.3, y_threshold=2, x_threshold=[-1,1], label_points={'BEST4':[-1.24288077425345,21.782377963035827]})
except:
    pass
else:
    fail('scatter() failed to catch error when only providing one condition to plot')

#Jointplot
xp.jointplot(geo_labeled, meta, 'STX6', 'STX6', kind='reg')

xp.jointplot(geo_labeled, meta, 'STX6', 'CCL5', kind='reg', palette=geo_colors, order=['Normal','Adenoma','Adenocarcinoma'], title_pad=-305, title_pos='center')

xp.jointplot(geo_labeled, meta, 'STX6', 'CCL5', kind='kde', palette=geo_colors, order=['Normal','Adenoma','Adenocarcinoma'])

#Linear Regression
xp.linreg(geo_labeled, 'STX6', linreg_file, delimiter=',')

#2D-PCA
xp.pca(geo_labeled, meta, geo_colors, grouping='samples', gene_list=None, gene_labels=False, ci=2, principle_components=[1,2], n_components=10, _3d_pca=False, scree_only=False, save_scree=None, size=10)

xp.pca(geo_labeled, meta, geo_colors, grouping='samples', gene_list=['STX6','SCARB1','CCL5'], gene_labels=False, ci=2, principle_components=[1,2], n_components=2, _3d_pca=False, scree_only=False, save_scree=None, size=20, order_legend=[1,3,2])

xp.pca(geo_labeled, meta, geo_colors, grouping='samples', gene_labels=False, ci=1, principle_components=[1,2], n_components=10, _3d_pca=False, scree_only=False, save_scree=None, size=20)

df_pca = xp.pca(geo_labeled, meta, geo_colors, grouping='samples', gene_labels=False, ci=1, principle_components=[2,3], n_components=10, _3d_pca=False, scree_only=False, save_scree=None, size=30, return_pca=True)
df_pca.head()

xp.pca(geo_labeled, meta, geo_colors, _3d_pca=False, scree_only=True, save_scree=True, save_fig=pca_file)

#3D-PCA
xp.pca(geo_labeled, meta, geo_colors, _3d_pca=True, order_legend=[1,3,2])

"""
INTERACTIVE TESTS
"""
#2D-plotly




#3D-plotly




#Scatter



#Volcano

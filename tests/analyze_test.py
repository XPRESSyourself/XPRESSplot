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
IMPORT DEPENDENCIES AND DATASETS
"""
import pandas as pd
import xpresstools as xp
%matplotlib inline
import numpy as np

data_loc = './large_test.csv'
meta_loc = './sample_info_test.csv'

data_loc = '/Users/jordan/scripts/XPRESSyourself/XPRESStools/tests/large_test.csv'
meta_loc = '/Users/jordan/scripts/XPRESSyourself/XPRESStools/tests/sample_info_test.csv'

"""
Get test data
"""
geo = xp.get_df(data_loc)
meta = xp.get_info(meta_loc)

meta[1] = meta[1].str.capitalize() #Make sample types look nice
meta = meta.replace('Normal_colon', 'Normal')
geo = xp.keep_labels(geo, meta, label_list=['Normal','Adenoma','Adenocarcinoma'])
geo_clean = xp.clean_df(geo)

#Collapse multi-mapping probes
geo_collapsed = xp.probe_collapse(geo_clean,"/Users/jordan/scripts/XPRESSyourself/XPRESStools/tests/GPL570.txt")

#Scale sorted dataset
geo_scaled, geo_labeled = xp.prep_data(geo_collapsed, meta)

geo_colors = {'Adenocarcinoma': (0.5725490196078431, 0.5843137254901961, 0.5686274509803921),
        'Adenoma': (0.8705882352941177, 0.5607843137254902, 0.0196078431372549),
        'Normal': (0.00784313725490196, 0.6196078431372549, 0.45098039215686275)}

gene_colors = {'Group1': (0.5725490196078431, 0.5843137254901961, 0.5686274509803921),
        'Group2': (0.8705882352941177, 0.5607843137254902, 0.0196078431372549)}

"""
Unit tests of intermediate functions for analysis
"""







"""
Run a couple of test instances of each analysis function to make sure there are no errors or exits
"""
#Single-gene summary
xp.gene_overview(geo_labeled, meta, gene_name='SEC62', palette=geo_colors, order=['Normal','Adenoma','Adenocarcinoma'])
xp.gene_overview(geo_labeled, meta, 'CCL5', geo_colors, order=None, save_fig=None, dpi=600, bbox_inches='tight', grid=True, whitegrid=True)

#Multi-gene summary
"""FIX: /anaconda3/lib/python3.6/site-packages/scipy/stats/stats.py:1713: FutureWarning:

Using a non-tuple sequence for multidimensional indexing is deprecated; use `arr[tuple(seq)]` instead of `arr[seq]`. In the future this will be interpreted as an array index, `arr[np.array(seq)]`, which will result either in an error or a different result.
"""

xp. multigene_overview(geo_labeled, meta, palette=geo_colors, gene_list=['SEC62','CCL5','STX6'], order=None, save_fig=None, dpi=600, bbox_inches='tight', title=None, grid=False, whitegrid=False)




#Scatterplot basic
xp.scatter(geo_labeled, meta, 'SEC62', 'STX6', palette=geo_colors, add_linreg=True, save_fig='/Users/jordan/Desktop/figure.pdf', order_legend=[1,3,2], alpha=.7)

xp.heatmap(geo_scaled, meta, sample_palette=geo_colors, gene_list=['SEC62','STX6','CCL5'], cbar_kws={'label':'z-score'}, figsize=(10,2), save_fig='/Users/jordan/Desktop/heatmap.pdf')

genes = np.array([['MPC1','Group1'],['MYC','Group2'],['MPC2','Group1']])
gene_info2 = pd.DataFrame({0:genes[:,0],1:genes[:,1]})
xp.heatmap(geo_scaled, meta, sample_palette=geo_colors, gene_palette=gene_colors, gene_info=gene_info2, gene_list=['SEC62','STX6','CCL5'], save_fig='/Users/jordan/Desktop/heatmap.pdf')

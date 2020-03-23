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
import os
import sys
import pandas as pd
import numpy as np
import xpressplot as xp

__path__  =  os.path.dirname(os.path.realpath(__file__)) + '/'

geo = 'GSE20916'
data_small_file = str(__path__ ) + 'small_test.csv'
metadata_file = str(__path__ ) + 'sample_info_test.csv'
count_dir = str(__path__ )

data_small = pd.read_csv(data_small_file, index_col=0)
metadata = pd.read_csv(metadata_file, header=None)

"""
Make test data_set
"""
def create_test():
    df_test = pd.DataFrame(columns=['fGSM523242','fGSM523243','fGSM523244','fGSM523245','fGSM523246'], index=['1007_s_at','1053_at','121_at','1294_at','1405_i_at'])
    df_test.loc['1007_s_at'] = pd.Series({'fGSM523242':66,'fGSM523243':59,'fGSM523244':1,'fGSM523245':82,'fGSM523246':45})
    df_test.loc['1053_at'] = pd.Series({'fGSM523242':35,'fGSM523243':0,'fGSM523244':7,'fGSM523245':72,'fGSM523246':2})
    df_test.loc['121_at'] = pd.Series({'fGSM523242':20,'fGSM523243':70,'fGSM523244':85,'fGSM523245':78,'fGSM523246':36})
    df_test.loc['1294_at'] = pd.Series({'fGSM523242':96,'fGSM523243':7,'fGSM523244':93,'fGSM523245':38,'fGSM523246':85})
    df_test.loc['1405_i_at'] = pd.Series({'fGSM523242':73,'fGSM523243':41,'fGSM523244':92,'fGSM523245':77,'fGSM523246':26})

    meta_data = np.array([['fGSM523242','normal'],['fGSM523243','adenoma'],['fGSM523244','adenoma'],['fGSM523245','normal'],['fGSM523246','adenocarcinoma']])
    meta_test = pd.DataFrame({0:meta_data[:,0],1:meta_data[:,1]})

    return df_test, meta_test

"""
Get dataframe
"""
df = xp.get_df(data_small_file, delimiter=',')
assert df.equals(data_small), 'get_df() failed'

"""
Get metadata
"""
#meta = xp.get_info(metadata_file)
#assert meta.equals(metadata), 'get_info() failed'

"""
Get geo dataset
Not running currently due to connection time out in testing
"""
#df, meta = xp.get_geo(geo)
#assert df.shape == data_large.shape, 'get_geo() failed at accessing values dataframe'
#assert meta.equals(metadata), 'get_geo() failed'

"""
drop_sample
"""
df_test, meta_test = create_test()
df_test = xp.drop_samples(df_test, ['fGSM523246'])
df_truth = pd.DataFrame(columns=['fGSM523242','fGSM523243','fGSM523244','fGSM523245'], index=['1007_s_at','1053_at','121_at','1294_at','1405_i_at'])
df_truth.loc['1007_s_at'] = pd.Series({'fGSM523242':66,'fGSM523243':59,'fGSM523244':1,'fGSM523245':82})
df_truth.loc['1053_at'] = pd.Series({'fGSM523242':35,'fGSM523243':0,'fGSM523244':7,'fGSM523245':72})
df_truth.loc['121_at'] = pd.Series({'fGSM523242':20,'fGSM523243':70,'fGSM523244':85,'fGSM523245':78})
df_truth.loc['1294_at'] = pd.Series({'fGSM523242':96,'fGSM523243':7,'fGSM523244':93,'fGSM523245':38})
df_truth.loc['1405_i_at'] = pd.Series({'fGSM523242':73,'fGSM523243':41,'fGSM523244':92,'fGSM523245':77})
assert df_test.equals(df_truth), 'drop_sample() failed'

"""
drop_label
"""
df_test, meta_test = create_test()
df_test = xp.drop_label(df_test, meta_test, 'adenocarcinoma')
df_truth = pd.DataFrame(columns=['fGSM523242','fGSM523243','fGSM523244','fGSM523245'], index=['1007_s_at','1053_at','121_at','1294_at','1405_i_at'])
df_truth.loc['1007_s_at'] = pd.Series({'fGSM523242':66,'fGSM523243':59,'fGSM523244':1,'fGSM523245':82})
df_truth.loc['1053_at'] = pd.Series({'fGSM523242':35,'fGSM523243':0,'fGSM523244':7,'fGSM523245':72})
df_truth.loc['121_at'] = pd.Series({'fGSM523242':20,'fGSM523243':70,'fGSM523244':85,'fGSM523245':78})
df_truth.loc['1294_at'] = pd.Series({'fGSM523242':96,'fGSM523243':7,'fGSM523244':93,'fGSM523245':38})
df_truth.loc['1405_i_at'] = pd.Series({'fGSM523242':73,'fGSM523243':41,'fGSM523244':92,'fGSM523245':77})
assert df_test.equals(df_truth), 'drop_label() failed'

"""
keep_labels
"""
df_test, meta_data = create_test()
df_test = xp.keep_labels(df_test, meta_test, ['normal','adenoma'])
df_truth = pd.DataFrame(columns=['fGSM523242','fGSM523243','fGSM523244','fGSM523245'], index=['1007_s_at','1053_at','121_at','1294_at','1405_i_at'])
df_truth.loc['1007_s_at'] = pd.Series({'fGSM523242':66,'fGSM523243':59,'fGSM523244':1,'fGSM523245':82})
df_truth.loc['1053_at'] = pd.Series({'fGSM523242':35,'fGSM523243':0,'fGSM523244':7,'fGSM523245':72})
df_truth.loc['121_at'] = pd.Series({'fGSM523242':20,'fGSM523243':70,'fGSM523244':85,'fGSM523245':78})
df_truth.loc['1294_at'] = pd.Series({'fGSM523242':96,'fGSM523243':7,'fGSM523244':93,'fGSM523245':38})
df_truth.loc['1405_i_at'] = pd.Series({'fGSM523242':73,'fGSM523243':41,'fGSM523244':92,'fGSM523245':77})
assert df_test.equals(df_truth), 'keep_labels() failed'

"""
catenate_files
"""
counts_truth = pd.DataFrame(columns=['fGSM523242_counts.txt','fGSM523243_counts.txt','fGSM523244_counts.txt','fGSM523245_counts.txt'], index=['1007_s_at','1053_at','121_at','1294_at','1405_i_at'])
counts_truth.loc['1007_s_at'] = pd.Series({'fGSM523242_counts.txt':66,'fGSM523243_counts.txt':59,'fGSM523244_counts.txt':1,'fGSM523245_counts.txt':82})
counts_truth.loc['1053_at'] = pd.Series({'fGSM523242_counts.txt':35,'fGSM523243_counts.txt':0,'fGSM523244_counts.txt':7,'fGSM523245_counts.txt':72})
counts_truth.loc['121_at'] = pd.Series({'fGSM523242_counts.txt':20,'fGSM523243_counts.txt':70,'fGSM523244_counts.txt':85,'fGSM523245_counts.txt':78})
counts_truth.loc['1294_at'] = pd.Series({'fGSM523242_counts.txt':96,'fGSM523243_counts.txt':7,'fGSM523244_counts.txt':93,'fGSM523245_counts.txt':38})
counts_truth.loc['1405_i_at'] = pd.Series({'fGSM523242_counts.txt':73,'fGSM523243_counts.txt':41,'fGSM523244_counts.txt':92,'fGSM523245_counts.txt':77})
#Default
counts = xp.catenate_files(count_dir, file_suffix='counts.txt', save_file=None, delimiter='\t', drop_rows=0)

assert (counts == counts_truth).all().all(), 'catenate_files() failed'
#drop_rows
counts_truth_drop = counts_truth[:-2]
counts = xp.catenate_files(count_dir, file_suffix='counts.txt', save_file=None, delimiter='\t', drop_rows=2)
assert (counts == counts_truth_drop).all().all(), 'catenate_files() failed when dropping tailing rows'

"""
rename_cols
"""
df_test, meta_data = create_test()
convert_data = np.array([['fGSM523242','normal1'],['fGSM523243','adenoma2'],['fGSM523244','adenoma3'],['fGSM523245','normal4'],['fGSM523246','adenocarcinoma5']])
convert_test = pd.DataFrame({0:convert_data[:,0],1:convert_data[:,1]})
df_test = xp.rename_cols(df_test, convert_test)
df_truth = pd.DataFrame(columns=['normal1','adenoma2','adenoma3','normal4','adenocarcinoma5'], index=['1007_s_at','1053_at','121_at','1294_at','1405_i_at'])
df_truth.loc['1007_s_at'] = pd.Series({'normal1':66,'adenoma2':59,'adenoma3':1,'normal4':82,'adenocarcinoma5':45})
df_truth.loc['1053_at'] = pd.Series({'normal1':35,'adenoma2':0,'adenoma3':7,'normal4':72,'adenocarcinoma5':2})
df_truth.loc['121_at'] = pd.Series({'normal1':20,'adenoma2':70,'adenoma3':85,'normal4':78,'adenocarcinoma5':36})
df_truth.loc['1294_at'] = pd.Series({'normal1':96,'adenoma2':7,'adenoma3':93,'normal4':38,'adenocarcinoma5':85})
df_truth.loc['1405_i_at'] = pd.Series({'normal1':73,'adenoma2':41,'adenoma3':92,'normal4':77,'adenocarcinoma5':26})
assert df_test.equals(df_truth), 'rename_cols() failed'

"""
rename_rows
"""
#Test renaming index
df_test, meta_data = create_test()
convert_data = np.array([['1007_s_at','Gene1'],['1294_at','Gene2']])
convert_test = pd.DataFrame({0:convert_data[:,0],1:convert_data[:,1]})
df_test = xp.rename_rows(df_test, convert_test, label='index')
df_truth = pd.DataFrame(columns=['fGSM523242','fGSM523243','fGSM523244','fGSM523245','fGSM523246'], index=['Gene1','1053_at','121_at','Gene2','1405_i_at'])
df_truth.loc['Gene1'] = pd.Series({'fGSM523242':66,'fGSM523243':59,'fGSM523244':1,'fGSM523245':82,'fGSM523246':45})
df_truth.loc['1053_at'] = pd.Series({'fGSM523242':35,'fGSM523243':0,'fGSM523244':7,'fGSM523245':72,'fGSM523246':2})
df_truth.loc['121_at'] = pd.Series({'fGSM523242':20,'fGSM523243':70,'fGSM523244':85,'fGSM523245':78,'fGSM523246':36})
df_truth.loc['Gene2'] = pd.Series({'fGSM523242':96,'fGSM523243':7,'fGSM523244':93,'fGSM523245':38,'fGSM523246':85})
df_truth.loc['1405_i_at'] = pd.Series({'fGSM523242':73,'fGSM523243':41,'fGSM523244':92,'fGSM523245':77,'fGSM523246':26})
assert df_test.equals(df_truth), 'rename_rows() failed'

#Test renaming a column
df_test['some_col'] = df_test.index
df_test = xp.rename_rows(df_test, convert_test, label='some_col')
df_truth['some_col'] = ['Gene1','1053_at','121_at','Gene2','1405_i_at']
assert df_test.equals(df_truth), 'rename_rows() failed'

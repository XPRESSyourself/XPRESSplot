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
gtf = str(__path__) + 'transcripts.gtf'
gpl_ref = str(__path__) + './GPL570.txt'

#gtf = '/Users/jordan/scripts/XPRESSyourself/XPRESSplot/tests/transcripts.gtf'
#gpl_ref = '/Users/jordan/scripts/XPRESSyourself/XPRESSplot/tests/GPL570.txt'

"""
Make test data_set
"""
def create_data():
    df_test = pd.DataFrame(columns=['fGSM523242','fGSM523243','fGSM523244','fGSM523245','fGSM523246'], index=['ENSG00000227232','ENSG00000240361','ENSG00000238009','ENSG00000241860','ENSG00000187634'])
    df_test.loc['ENSG00000227232'] = pd.Series({'fGSM523242':66,'fGSM523243':59,'fGSM523244':1,'fGSM523245':82,'fGSM523246':45})
    df_test.loc['ENSG00000240361'] = pd.Series({'fGSM523242':35,'fGSM523243':0,'fGSM523244':7,'fGSM523245':72,'fGSM523246':2})
    df_test.loc['ENSG00000238009'] = pd.Series({'fGSM523242':20,'fGSM523243':70,'fGSM523244':85,'fGSM523245':78,'fGSM523246':36})
    df_test.loc['ENSG00000241860'] = pd.Series({'fGSM523242':96,'fGSM523243':7,'fGSM523244':93,'fGSM523245':38,'fGSM523246':85})
    df_test.loc['ENSG00000187634'] = pd.Series({'fGSM523242':73,'fGSM523243':41,'fGSM523244':92,'fGSM523245':77,'fGSM523246':26})

    return df_test

def create_full():
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
rpm()
"""
data_rpm = create_data()
data_rpm = xp.rpm(data_rpm)
rpm_truth = pd.DataFrame(columns=['fGSM523242','fGSM523243','fGSM523244','fGSM523245','fGSM523246'], index=['ENSG00000227232','ENSG00000240361','ENSG00000238009','ENSG00000241860','ENSG00000187634'], dtype='float')
rpm_truth.loc['ENSG00000227232'] = pd.Series({'fGSM523242':float(227586.206897), 'fGSM523243':float(333333.333333), 'fGSM523244':float(3597.12230216), 'fGSM523245':float(236311.239193), 'fGSM523246':float(231958.762887)})
rpm_truth.loc['ENSG00000240361'] = pd.Series({'fGSM523242':float(120689.655172), 'fGSM523243':float(0.000000), 'fGSM523244':float(25179.8561151), 'fGSM523245':float(207492.795389), 'fGSM523246':float(10309.2783505)})
rpm_truth.loc['ENSG00000238009'] = pd.Series({'fGSM523242':float(68965.5172414), 'fGSM523243':float(395480.225989), 'fGSM523244':float(305755.395683), 'fGSM523245':float(224783.861671), 'fGSM523246':float(185567.010309)})
rpm_truth.loc['ENSG00000241860'] = pd.Series({'fGSM523242':float(331034.482759), 'fGSM523243':float(39548.0225989), 'fGSM523244':float(334532.374101), 'fGSM523245':float(109510.086455), 'fGSM523246':float(438144.329897)})
rpm_truth.loc['ENSG00000187634'] = pd.Series({'fGSM523242':float(251724.137931), 'fGSM523243':float(231638.418079), 'fGSM523244':float(330935.251799), 'fGSM523245':float(221902.017291), 'fGSM523246':float(134020.618557)})
data_rpm = data_rpm.round(decimals=4)
rpm_truth = rpm_truth.round(decimals=4)
assert data_rpm.equals(rpm_truth), 'rpm() failed'

"""
Other Normalization Tests
"""
from xpressplot.normalize import rpm, rpk

# Covers all normalization methods
# Edge cases covered:
# - Row with 0s
# - Column with 0s
# - Missing gene in data
# - Gene in dictionary not used by normalization method

length_dict = {
    'MPC2':2223, # -> 2.223
    'SLC1A1':3698, # -> 3.698
    'ATF4':2019 # -> 2.019
}

df = pd.DataFrame()
df['sample1'] = [100000,222300,605700,0] # -> .928000
df['sample2'] = [150000,444600,201900,0] # -> .796500
df['sample3'] = [0,0,0,0] # -> 0
df.index = ['FAKEGENE','MPC2','ATF4','SLC1A1']

# Test RPK
rpk_df = rpk(
    df,
    length_dict)

rpk_truth = pd.DataFrame()
rpk_truth['sample1'] = [100000.0,300000.0,0.0]
rpk_truth['sample2'] = [200000.0,100000.0,0.0]
rpk_truth['sample3'] = [0.0,0.0,0.0]
rpk_truth.index = ['MPC2','ATF4','SLC1A1']
rpk_df = rpk_df.round(decimals=4)
rpk_truth = rpk_truth.round(decimals=4)

assert rpk_df.loc[['MPC2','ATF4','SLC1A1']].equals(rpk_truth.loc[['MPC2','ATF4','SLC1A1']]), 'RPK failed'

# Test RPM
rpm_df = rpm(df)

rpm_truth = pd.DataFrame()
rpm_truth['sample1'] = [107758.620690,239547.413793,652693.965517,0.000000]
rpm_truth['sample2'] = [188323.917137,558192.090395,253483.992467,0.000000]
rpm_truth['sample3'] = [0.0,0.0,0.0,0.0]
rpm_truth.index = ['FAKEGENE','MPC2','ATF4','SLC1A1']
rpm_df = rpm_df.round(decimals=4)
rpm_truth = rpm_truth.round(decimals=4)

assert rpm_df.equals(rpm_truth), 'RPM failed'

# Test RPKM/FPKM
rpm_rpm = rpm(df)
r_fpkm = rpk(
    rpm_rpm,
    length_dict)

r_fpkm_truth = pd.DataFrame()
r_fpkm_truth['sample1'] = [107758.6207,323275.8621,0.0000]
r_fpkm_truth['sample2'] = [251098.5562,125549.2781,0.0000]
r_fpkm_truth['sample3'] = [0.0,0.0,0.0]
r_fpkm_truth.index = ['MPC2','ATF4','SLC1A1']
r_fpkm = r_fpkm.round(decimals=4)
r_fpkm_truth = r_fpkm_truth.round(decimals=4)

assert r_fpkm.loc[['MPC2','ATF4','SLC1A1']].equals(r_fpkm_truth.loc[['MPC2','ATF4','SLC1A1']]), 'RPKM/FPKM failed'

# Test TPM
tpm_rpk = rpk(
    df,
    length_dict)
tpm_df = rpm(tpm_rpk)

tpm_truth = pd.DataFrame()
tpm_truth['sample1'] = [250000.0,750000.0,0.0]
tpm_truth['sample2'] = [666666.666667,333333.333333,0.0]
tpm_truth['sample3'] = [0.0,0.0,0.0]
tpm_truth.index = ['MPC2','ATF4','SLC1A1']
tpm_df = tpm_df.round(decimals=4)
tpm_truth = tpm_truth.round(decimals=4)

assert tpm_df.loc[['MPC2','ATF4','SLC1A1']].equals(tpm_truth.loc[['MPC2','ATF4','SLC1A1']]), 'TPM failed'



"""
batch_normalize()
"""


"""
check_samples()
"""
check_test = pd.DataFrame(columns=['fGSM523242_rpf','fGSM523243_rna','fGSM523244_rpf','fGSM523245_rna'], index=['ENSG00000227232','ENSG00000240361','ENSG00000238009','ENSG00000241860','ENSG00000187634'], dtype='float')
check_test.loc['ENSG00000227232'] = pd.Series({'fGSM523242_rpf':66.34,'fGSM523243_rna':59.13,'fGSM523244_rpf':1.90,'fGSM523245_rna':82.49})
check_test.loc['ENSG00000240361'] = pd.Series({'fGSM523242_rpf':35.73,'fGSM523243_rna':0.00,'fGSM523244_rpf':7.38,'fGSM523245_rna':72.94})
check_test.loc['ENSG00000238009'] = pd.Series({'fGSM523242_rpf':20.02,'fGSM523243_rna':70.21,'fGSM523244_rpf':85.10,'fGSM523245_rna':78.87})
check_test.loc['ENSG00000241860'] = pd.Series({'fGSM523242_rpf':96.23,'fGSM523243_rna':7.49,'fGSM523244_rpf':93.49,'fGSM523245_rna':38.39})
check_test.loc['ENSG00000187634'] = pd.Series({'fGSM523242_rpf':73.91,'fGSM523243_rna':41.28,'fGSM523244_rpf':92.27,'fGSM523245_rna':77.93})
xp.check_samples(check_test)

"""
clean_df()
"""
clean_test = pd.DataFrame(columns=['fGSM523242_rpf','fGSM523243_rna','fGSM523244_rpf','fGSM523245_rna'], index=['ENSG00000227232','ENSG00000240361','ENSG00000238009','ENSG00000241860','ENSG00000187634','ENSG00000227200','ENSG00000227201','ENSG00000241899'], dtype='float')
clean_test.loc['ENSG00000227232'] = pd.Series({'fGSM523242_rpf':66.34,'fGSM523243_rna':59.13,'fGSM523244_rpf':1.90,'fGSM523245_rna':82.49})
clean_test.loc['ENSG00000240361'] = pd.Series({'fGSM523242_rpf':35.73,'fGSM523243_rna':0.00,'fGSM523244_rpf':7.38,'fGSM523245_rna':72.94})
clean_test.loc['ENSG00000238009'] = pd.Series({'fGSM523242_rpf':20.02,'fGSM523243_rna':70.21,'fGSM523244_rpf':85.10,'fGSM523245_rna':78.87})
clean_test.loc['ENSG00000241860'] = pd.Series({'fGSM523242_rpf':96.23,'fGSM523243_rna':7.49,'fGSM523244_rpf':93.49,'fGSM523245_rna':38.39})
clean_test.loc['ENSG00000187634'] = pd.Series({'fGSM523242_rpf':73.91,'fGSM523243_rna':41.28,'fGSM523244_rpf':92.27,'fGSM523245_rna':77.93})
clean_test.loc['ENSG00000227200'] = pd.Series({'fGSM523242_rpf':66.34,'fGSM523243_rna':59.13,'fGSM523244_rpf':1.90,'fGSM523245_rna':82.49})
clean_test.loc['ENSG00000227201'] = pd.Series({'fGSM523242_rpf':35.73,'fGSM523243_rna':0.00,'fGSM523244_rpf':7.38,'fGSM523245_rna':72.94})
clean_test.loc['ENSG00000241899'] = pd.Series({'fGSM523242_rpf':np.NaN,'fGSM523243_rna':7.49,'fGSM523244_rpf':93.49,'fGSM523245_rna':38.39})
clean_test = clean_test.rename(index={'ENSG00000227201':'ENSG00000227200'})

clean_test = xp.clean_df(clean_test)
assert clean_test.equals(check_test), 'clean_df() failed'

"""
prep_data()
"""
data, info = create_full()

scale_truth = pd.DataFrame(columns=['fGSM523242','fGSM523243','fGSM523244','fGSM523245','fGSM523246'], index=['1007_s_at','1053_at','121_at','1294_at','1405_i_at'])
scale_truth.loc['1007_s_at'] = pd.Series({'fGSM523242':0.559708,'fGSM523243':0.305295,'fGSM523244':-1.802695,'fGSM523245':1.141222,'fGSM523246':-0.203530})
scale_truth.loc['1053_at'] = pd.Series({'fGSM523242':0.429685,'fGSM523243':-0.844805,'fGSM523244':-0.589907,'fGSM523245':1.777003,'fGSM523246':-0.771977})
scale_truth.loc['121_at'] = pd.Series({'fGSM523242':-1.493989,'fGSM523243':0.482187,'fGSM523244':1.075040,'fGSM523245':0.798375,'fGSM523246':-0.861613})
scale_truth.loc['1294_at'] = pd.Series({'fGSM523242':0.912156,'fGSM523243':-1.609020,'fGSM523244':0.827172,'fGSM523245':-0.730858,'fGSM523246':0.600550})
scale_truth.loc['1405_i_at'] = pd.Series({'fGSM523242':0.458554,'fGSM523243':-0.851601,'fGSM523244':1.236459,'fGSM523245':0.622324,'fGSM523246':-1.465737})

label_truth = pd.DataFrame(columns=['fGSM523242','fGSM523243','fGSM523244','fGSM523245','fGSM523246'], index=['label','1007_s_at','1053_at','121_at','1294_at','1405_i_at'])
label_truth.loc['label'] = pd.Series({'fGSM523242':'normal','fGSM523243':'adenoma','fGSM523244':'adenoma','fGSM523245':'normal','fGSM523246':'adenocarcinoma'})
label_truth.loc['1007_s_at'] = pd.Series({'fGSM523242':66,'fGSM523243':59,'fGSM523244':1,'fGSM523245':82,'fGSM523246':45})
label_truth.loc['1053_at'] = pd.Series({'fGSM523242':35,'fGSM523243':0,'fGSM523244':7,'fGSM523245':72,'fGSM523246':2})
label_truth.loc['121_at'] = pd.Series({'fGSM523242':20,'fGSM523243':70,'fGSM523244':85,'fGSM523245':78,'fGSM523246':36})
label_truth.loc['1294_at'] = pd.Series({'fGSM523242':96,'fGSM523243':7,'fGSM523244':93,'fGSM523245':38,'fGSM523246':85})
label_truth.loc['1405_i_at'] = pd.Series({'fGSM523242':73,'fGSM523243':41,'fGSM523244':92,'fGSM523245':77,'fGSM523246':26})

data_noscale, data_label = xp.prep_data(data, info, gene_scale=False, print_means=True)
assert data_noscale.astype(int).equals(data), 'prep_data() failed'
assert data_label.equals(label_truth), 'prep_data() failed'

data_scale, data_label_s = xp.prep_data(data, info, gene_scale=True, print_means=False)
data_scale = np.around(data_scale.astype(float), decimals=3)
scale_truth = np.around(scale_truth.astype(float), decimals=3)
assert data_scale.equals(scale_truth), 'prep_data() failed'
assert data_label_s.equals(label_truth), 'prep_data() failed'

"""
probe_collapse()
"""
probe_test = pd.DataFrame(columns=['fGSM523242','fGSM523243','fGSM523244','fGSM523245','fGSM523246'], index=['1007_s_at','1053_at','121_at','218024_at','240362_at'])
probe_test.loc['1007_s_at'] = pd.Series({'fGSM523242':66,'fGSM523243':59,'fGSM523244':1,'fGSM523245':82,'fGSM523246':45})
probe_test.loc['1053_at'] = pd.Series({'fGSM523242':35,'fGSM523243':0,'fGSM523244':7,'fGSM523245':72,'fGSM523246':2})
probe_test.loc['121_at'] = pd.Series({'fGSM523242':20,'fGSM523243':70,'fGSM523244':85,'fGSM523245':78,'fGSM523246':36})
probe_test.loc['218024_at'] = pd.Series({'fGSM523242':96,'fGSM523243':7,'fGSM523244':93,'fGSM523245':38,'fGSM523246':85})
probe_test.loc['240362_at'] = pd.Series({'fGSM523242':73,'fGSM523243':41,'fGSM523244':92,'fGSM523245':77,'fGSM523246':26})

probe_truth = pd.DataFrame(columns=['fGSM523242','fGSM523243','fGSM523244','fGSM523245','fGSM523246'], index=['RFC2','PAX8','MPC1'])
probe_truth.loc['RFC2'] = pd.Series({'fGSM523242':35,'fGSM523243':0,'fGSM523244':7,'fGSM523245':72,'fGSM523246':2})
probe_truth.loc['PAX8'] = pd.Series({'fGSM523242':20,'fGSM523243':70,'fGSM523244':85,'fGSM523245':78,'fGSM523246':36})
probe_truth.loc['MPC1'] = pd.Series({'fGSM523242':84.5,'fGSM523243':24,'fGSM523244':92.5,'fGSM523245':57.5,'fGSM523246':55.5})

multimap_truth = pd.DataFrame(columns=['fGSM523242','fGSM523243','fGSM523244','fGSM523245','fGSM523246'], index=['DDR1 /// MIR4640','RFC2','PAX8','MPC1'])
multimap_truth.loc['DDR1 /// MIR4640'] = pd.Series({'fGSM523242':66,'fGSM523243':59,'fGSM523244':1,'fGSM523245':82,'fGSM523246':45})
multimap_truth.loc['RFC2'] = pd.Series({'fGSM523242':35,'fGSM523243':0,'fGSM523244':7,'fGSM523245':72,'fGSM523246':2})
multimap_truth.loc['PAX8'] = pd.Series({'fGSM523242':20,'fGSM523243':70,'fGSM523244':85,'fGSM523245':78,'fGSM523246':36})
multimap_truth.loc['MPC1'] = pd.Series({'fGSM523242':84.5,'fGSM523243':24,'fGSM523244':92.5,'fGSM523245':57.5,'fGSM523246':55.5})

probe_collapse = xp.probe_collapse(probe_test, gpl_ref)
probe_collapse = probe_collapse.sort_index().astype(float)
probe_truth = probe_truth.sort_index().astype(float)
assert probe_collapse.equals(probe_truth), 'probe_collapse() failed during no_multimappers test'

multimap_collapse = xp.probe_collapse(probe_test, gpl_ref, no_multimappers=False)
multimap_collapse = multimap_collapse.sort_index().astype(float)
multimap_truth = multimap_truth.sort_index().astype(float)
assert multimap_collapse.equals(multimap_truth), 'probe_collapse() failed during multimappers test'

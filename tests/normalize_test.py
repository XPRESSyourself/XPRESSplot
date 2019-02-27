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
import os, sys
import pandas as pd
import numpy as np
import xpresstools as xp

gtf = './transcripts.gtf'

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
r_fpkm()
"""
data_rpkm = create_data()
data_rpkm = xp.r_fpkm(data_rpkm, gtf, gene_name_location=0)
rpkm_truth = pd.DataFrame(columns=['fGSM523242','fGSM523243','fGSM523244','fGSM523245','fGSM523246'], index=['ENSG00000227232','ENSG00000240361','ENSG00000238009','ENSG00000241860','ENSG00000187634'], dtype='float')
rpkm_truth.loc['ENSG00000227232'] = pd.Series({'fGSM523242':float(15006.343591), 'fGSM523243':float(21978.988087), 'fGSM523244':float(237.183325), 'fGSM523245':float(15581.645733), 'fGSM523246':float(15294.656659)})
rpkm_truth.loc['ENSG00000240361'] = pd.Series({'fGSM523242':float(18516.363175), 'fGSM523243':float(0.000000), 'fGSM523244':float(3863.126130), 'fGSM523245':float(31833.813346), 'fGSM523246':float(1581.662834)})
rpkm_truth.loc['ENSG00000238009'] = pd.Series({'fGSM523242':float(1552.298488), 'fGSM523243':float(8901.598676), 'fGSM523244':float(6882.042759), 'fGSM523245':float(5059.508906), 'fGSM523246':float(4176.803149)})
rpkm_truth.loc['ENSG00000241860'] = pd.Series({'fGSM523242':float(10220.899184), 'fGSM523243':float(1221.070230), 'fGSM523244':float(10328.898793), 'fGSM523245':float(3381.193234), 'fGSM523246':float(13527.983509)})
rpkm_truth.loc['ENSG00000187634'] = pd.Series({'fGSM523242':float(12188.260201), 'fGSM523243':float(11215.727404), 'fGSM523244':float(16023.592301), 'fGSM523245':float(10744.299486), 'fGSM523246':float(6489.159858)})
data_rpkm = data_rpkm.round(decimals=4)
rpkm_truth = rpkm_truth.round(decimals=4)
assert data_rpkm.equals(rpkm_truth), 'r_fpkm() failed'

"""
te()
"""
te_test = pd.DataFrame(columns=['fGSM523242_rpf','fGSM523243_rna','fGSM523244_rpf','fGSM523245_rna'], index=['ENSG00000227232','ENSG00000240361','ENSG00000238009','ENSG00000241860','ENSG00000187634'], dtype='float')
te_test.loc['ENSG00000227232'] = pd.Series({'fGSM523242_rpf':66.34,'fGSM523243_rna':59.13,'fGSM523244_rpf':1.90,'fGSM523245_rna':82.49})
te_test.loc['ENSG00000240361'] = pd.Series({'fGSM523242_rpf':35.73,'fGSM523243_rna':0.00,'fGSM523244_rpf':7.38,'fGSM523245_rna':72.94})
te_test.loc['ENSG00000238009'] = pd.Series({'fGSM523242_rpf':20.02,'fGSM523243_rna':70.21,'fGSM523244_rpf':85.10,'fGSM523245_rna':78.87})
te_test.loc['ENSG00000241860'] = pd.Series({'fGSM523242_rpf':96.23,'fGSM523243_rna':7.49,'fGSM523244_rpf':93.49,'fGSM523245_rna':38.39})
te_test.loc['ENSG00000187634'] = pd.Series({'fGSM523242_rpf':73.91,'fGSM523243_rna':41.28,'fGSM523244_rpf':92.27,'fGSM523245_rna':77.93})

te_truth = pd.DataFrame(columns=['fGSM523242_rpf_te','fGSM523244_rpf_te'], index=['ENSG00000227232','ENSG00000240361','ENSG00000238009','ENSG00000241860','ENSG00000187634'], dtype='float')
te_truth.loc['ENSG00000227232'] = pd.Series({'fGSM523242_rpf_te':0.165724,'fGSM523244_rpf_te':-5.367895})
te_truth.loc['ENSG00000240361'] = pd.Series({'fGSM523242_rpf_te':8.485024,'fGSM523244_rpf_te':-3.287577})
te_truth.loc['ENSG00000238009'] = pd.Series({'fGSM523242_rpf_te':-1.805100,'fGSM523244_rpf_te':0.109549})
te_truth.loc['ENSG00000241860'] = pd.Series({'fGSM523242_rpf_te':3.665813,'fGSM523244_rpf_te':1.281871})
te_truth.loc['ENSG00000187634'] = pd.Series({'fGSM523242_rpf_te':0.838787,'fGSM523244_rpf_te':0.243395})

te_test1 = xp.te(te_test, samples=None, log2=True)
te_test1 = te_test1.round(decimals=4)
te_truth = te_truth.round(decimals=4)
assert te_test1.equals(te_truth), 'te() failed'

te_test2 = xp.te(te_test, samples=['sample1','sample2'], log2=True)
te_test2 = te_test2.round(decimals=4)
te_truth.columns = ['sample1','sample2']
te_truth = te_truth.round(decimals=4)
assert te_test2.equals(te_truth), 'te() failed'

"""
log_scale()
"""
log_test = pd.DataFrame(columns=['fGSM523242_rpf','fGSM523243_rna','fGSM523244_rpf','fGSM523245_rna'], index=['ENSG00000227232','ENSG00000240361','ENSG00000238009','ENSG00000241860','ENSG00000187634'], dtype='float')
log_test.loc['ENSG00000227232'] = pd.Series({'fGSM523242_rpf':66.34,'fGSM523243_rna':59.13,'fGSM523244_rpf':1.90,'fGSM523245_rna':82.49})
log_test.loc['ENSG00000240361'] = pd.Series({'fGSM523242_rpf':35.73,'fGSM523243_rna':0.00,'fGSM523244_rpf':7.38,'fGSM523245_rna':72.94})
log_test.loc['ENSG00000238009'] = pd.Series({'fGSM523242_rpf':20.02,'fGSM523243_rna':70.21,'fGSM523244_rpf':85.10,'fGSM523245_rna':78.87})
log_test.loc['ENSG00000241860'] = pd.Series({'fGSM523242_rpf':96.23,'fGSM523243_rna':7.49,'fGSM523244_rpf':93.49,'fGSM523245_rna':38.39})
log_test.loc['ENSG00000187634'] = pd.Series({'fGSM523242_rpf':73.91,'fGSM523243_rna':41.28,'fGSM523244_rpf':92.27,'fGSM523245_rna':77.93})

log_truth = pd.DataFrame(columns=['fGSM523242_rpf','fGSM523243_rna','fGSM523244_rpf','fGSM523245_rna'], index=['ENSG00000227232','ENSG00000240361','ENSG00000238009','ENSG00000241860','ENSG00000187634'], dtype='float')
log_truth.loc['ENSG00000227232'] = pd.Series({'fGSM523242_rpf':1.822430,'fGSM523243_rna':1.772542,'fGSM523244_rpf':0.301030,'fGSM523245_rna':1.916927})
log_truth.loc['ENSG00000240361'] = pd.Series({'fGSM523242_rpf':1.554247,'fGSM523243_rna':-1.000000,'fGSM523244_rpf':0.873902,'fGSM523245_rna':1.863561})
log_truth.loc['ENSG00000238009'] = pd.Series({'fGSM523242_rpf':1.303628,'fGSM523243_rna':1.847017,'fGSM523244_rpf':1.930440,'fGSM523245_rna':1.897462})
log_truth.loc['ENSG00000241860'] = pd.Series({'fGSM523242_rpf':1.983762,'fGSM523243_rna':0.880242,'fGSM523244_rpf':1.971229,'fGSM523245_rna':1.585348})
log_truth.loc['ENSG00000187634'] = pd.Series({'fGSM523242_rpf':1.869290,'fGSM523243_rna':1.616790,'fGSM523244_rpf':1.965531,'fGSM523245_rna':1.892262})
#Check output
log_test1 = xp.log_scale(log_test)
log_test1 = log_test1.round(decimals=4)
log_truth = log_truth.round(decimals=4)
assert log_test1.equals(log_truth), 'log_scale() failed'
#Check that it runs
log_test2 = xp.log_scale(log_test, log_base=2)
#Check that it catches error inputs
try:
    xp.log_scale(log_test, log_base=3)
except:
    pass
else:
    fail('log_scale() failed')

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
threshold()
"""
threshold_test1 = pd.DataFrame(columns=['fGSM523242_rpf','fGSM523243_rna','fGSM523244_rpf','fGSM523245_rna'], index=['ENSG00000227232','ENSG00000240361','ENSG00000238009','ENSG00000241860','ENSG00000187634'], dtype='float')
threshold_test1.loc['ENSG00000227232'] = pd.Series({'fGSM523242_rpf':66.34,'fGSM523243_rna':59.13,'fGSM523244_rpf':1.90,'fGSM523245_rna':82.49})
threshold_test1.loc['ENSG00000240361'] = pd.Series({'fGSM523242_rpf':35.73,'fGSM523243_rna':0.00,'fGSM523244_rpf':7.38,'fGSM523245_rna':72.94})
threshold_test1.loc['ENSG00000238009'] = pd.Series({'fGSM523242_rpf':20.02,'fGSM523243_rna':70.21,'fGSM523244_rpf':85.10,'fGSM523245_rna':78.87})
threshold_test1.loc['ENSG00000241860'] = pd.Series({'fGSM523242_rpf':96.23,'fGSM523243_rna':7.49,'fGSM523244_rpf':93.49,'fGSM523245_rna':38.39})
threshold_test1.loc['ENSG00000187634'] = pd.Series({'fGSM523242_rpf':73.91,'fGSM523243_rna':41.28,'fGSM523244_rpf':92.27,'fGSM523245_rna':77.93})

threshold_truth1 = pd.DataFrame(columns=['fGSM523242_rpf','fGSM523243_rna','fGSM523244_rpf','fGSM523245_rna'], index=['ENSG00000238009','ENSG00000241860','ENSG00000187634'], dtype='float')
threshold_truth1.loc['ENSG00000238009'] = pd.Series({'fGSM523242_rpf':20.02,'fGSM523243_rna':70.21,'fGSM523244_rpf':85.10,'fGSM523245_rna':78.87})
threshold_truth1.loc['ENSG00000241860'] = pd.Series({'fGSM523242_rpf':96.23,'fGSM523243_rna':7.49,'fGSM523244_rpf':93.49,'fGSM523245_rna':38.39})
threshold_truth1.loc['ENSG00000187634'] = pd.Series({'fGSM523242_rpf':73.91,'fGSM523243_rna':41.28,'fGSM523244_rpf':92.27,'fGSM523245_rna':77.93})

threshold_truth2 = pd.DataFrame(columns=['fGSM523242_rpf','fGSM523243_rna','fGSM523244_rpf','fGSM523245_rna'], index=['ENSG00000238009','ENSG00000187634'], dtype='float')
threshold_truth2.loc['ENSG00000238009'] = pd.Series({'fGSM523242_rpf':20.02,'fGSM523243_rna':70.21,'fGSM523244_rpf':85.10,'fGSM523245_rna':78.87})
threshold_truth2.loc['ENSG00000187634'] = pd.Series({'fGSM523242_rpf':73.91,'fGSM523243_rna':41.28,'fGSM523244_rpf':92.27,'fGSM523245_rna':77.93})

threshold_test_util = xp.threshold_util(threshold_test1, 5, None)
threshold_test1 = xp.threshold(threshold_test1, 5, None)
assert threshold_test1.equals(threshold_test_util), 'threshold() failed at setting thresholding minimum'
assert threshold_test1.equals(threshold_truth1), 'threshold() failed at setting thresholding minimum'

threshold_test2 = xp.threshold(threshold_test1,minimum=None, maximum=93)
assert threshold_test2.equals(threshold_truth2), 'threshold() failed at setting maximum threshold'

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

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

#gtf = './transcripts.gtf'
gtf = '/Users/jordan/scripts/XPRESSyourself/XPRESStools/tests/transcripts.gtf'

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

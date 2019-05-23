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
import os, sys
import pandas as pd
import numpy as np
import xpressplot as xp

__path__  =  os.path.dirname(os.path.realpath(__file__)) + '/'
gtf = str(__path__) + 'transcripts.gtf'
gtf_output = str(__path__)

#gtf = '/Users/jordan/scripts/XPRESSyourself/XPRESSplot/tests/transcripts.gtf'
#gtf_output = '/Users/jordan/scripts/XPRESSyourself/XPRESSplot/tests/'

"""
Test dataset
"""
def create_data():
    df_test = pd.DataFrame(columns=['fGSM523242','fGSM523243','fGSM523244','fGSM523245','fGSM523246'], index=['ENSG00000227232','ENSG00000240361','ENSG00000238009','ENSG00000241860','ENSG00000187634'])
    df_test.loc['ENSG00000227232'] = pd.Series({'fGSM523242':66,'fGSM523243':59,'fGSM523244':1,'fGSM523245':82,'fGSM523246':45})
    df_test.loc['ENSG00000240361'] = pd.Series({'fGSM523242':35,'fGSM523243':0,'fGSM523244':7,'fGSM523245':72,'fGSM523246':2})
    df_test.loc['ENSG00000238009'] = pd.Series({'fGSM523242':20,'fGSM523243':70,'fGSM523244':85,'fGSM523245':78,'fGSM523246':36})
    df_test.loc['ENSG00000241860'] = pd.Series({'fGSM523242':96,'fGSM523243':7,'fGSM523244':93,'fGSM523245':38,'fGSM523246':85})
    df_test.loc['ENSG00000187634'] = pd.Series({'fGSM523242':73,'fGSM523243':41,'fGSM523244':92,'fGSM523245':77,'fGSM523246':26})

    df_truth = pd.DataFrame(columns=['fGSM523242','fGSM523243','fGSM523244','fGSM523245','fGSM523246'], index=['WASH7P','OR4G11P','AL627309.1','AL627309.5','SAMD11'])
    df_truth.loc['WASH7P'] = pd.Series({'fGSM523242':66,'fGSM523243':59,'fGSM523244':1,'fGSM523245':82,'fGSM523246':45})
    df_truth.loc['OR4G11P'] = pd.Series({'fGSM523242':35,'fGSM523243':0,'fGSM523244':7,'fGSM523245':72,'fGSM523246':2})
    df_truth.loc['AL627309.1'] = pd.Series({'fGSM523242':20,'fGSM523243':70,'fGSM523244':85,'fGSM523245':78,'fGSM523246':36})
    df_truth.loc['AL627309.5'] = pd.Series({'fGSM523242':96,'fGSM523243':7,'fGSM523244':93,'fGSM523245':38,'fGSM523246':85})
    df_truth.loc['SAMD11'] = pd.Series({'fGSM523242':73,'fGSM523243':41,'fGSM523244':92,'fGSM523245':77,'fGSM523246':26})

    return df_test, df_truth

"""
convert_names_gtf
"""
data, data_truth = create_data()
data = xp.convert_names(data, gtf, orig_name_label='gene_id \"', orig_name_location=0, new_name_label='gene_name \"', new_name_location=2, refill=None, sep='\t')
assert data.equals(data_truth), 'convert_names() failed'

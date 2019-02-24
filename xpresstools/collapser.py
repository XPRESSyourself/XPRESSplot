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
import pandas as pd
from .utils_collapser import prep_collapser, run_collapser

"""
DESCRIPTION: Ties prep_collapser and probe_collapse functions together

VARIABLES:
data= Dataframe of microarray probe data to be collapsed
reference= Full path and file name for GPL reference file, accessed from NCBI (should be a .txt file)
gene_list= Full path and file name to .csv file listing gene names to get probes for. Resulting function dictionary will only contain these genes and their probes (default: None)
no_multimappers= Do not allow ambiguous probes in the probe-gene dictionary (default: True)

USAGE:
import micartools as mat
df_collapsed = mat.auto_collapse(df, "~/Desktop/GPL570.csv")

ASSUMPTIONS:
See assumptions for prep_collapser and probe_collapse functions
"""
def probe_collapse(data, reference, gene_list=None, no_multimappers=True):

    dict = prep_collapser(reference, gene_list=gene_list, no_multimappers=no_multimappers)
    data_collapsed = run_collapser(data, dict)

    return data_collapsed

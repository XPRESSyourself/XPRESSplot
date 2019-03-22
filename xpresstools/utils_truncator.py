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
import csv
import pandas as pd
pd.options.mode.chained_assignment = None
import numpy as np
from multiprocessing import cpu_count, Pool
from functools import partial
from .utils import check_directories

"""
DESCRIPTION: Create a GTF reference file with only protein coding genes and the first n nucleotides of each first exon

RETURNS: Return the truncated, coding only table when running the function

METHODS: Considers strandedness
Multiprocesses chunks of a dataframe on cores of computer

VARIABLES:
input_gtf= GTF reference file path and name
truncate_amount= Number of nucleotides to truncate from the 5' end of each exon 1 (considers strandedness)
save_coding_path= Location to save coding-only GTF, set to None to avoid outputting file
save_truncated_path= Location to save coding-only truncated GTF
sep= Separator type of the GTF file (generally always tab-delimited)

ASSUMPTIONS:
Input file is a properly formatted GTF file
Protein coding transcripts are denoted by 'protein_coding in the final column of the GTF'
"""
def truncate(input_gtf, truncate_amount=45, save_coding_path='./', save_truncated_path='./', sep='\t', return_files=False):

    if save_coding_path != None:
        save_coding_path = check_directories(save_coding_path)
    if save_truncated_path != None:
        save_truncated_path = check_directories(save_truncated_path)

    #Import gtf reference file to
    if str(input_gtf).endswith('.gtf'):
        gtf = pd.read_csv(str(input_gtf), sep=sep, header=None, comment='#', low_memory=False)
    else:
        raise Exception('Error: A GTF-formatted file was not provided')

    #Get only protein_coding coordinates
    gtf_coding = gtf[gtf.iloc[:, 8].str.contains('protein_coding') == True]

    #Save to .gtf file (tsv)
    if save_coding_path != None:
        gtf_coding.to_csv(str(save_coding_path) + 'transcripts_coding.gtf', sep='\t', header=None, index=False, quoting=csv.QUOTE_NONE)

    if truncate_amount != None:
        gtf_coding_c = gtf_coding.copy()

        print("Multiprocessing reference chunks -- this may take a while...")
        gtf_truncated = parallelize_truncator(execute_truncator, gtf_coding_c, truncate_amount)

        if save_truncated_path != None:
            gtf_truncated.to_csv(str(save_truncated_path) + 'transcripts_coding_truncated.gtf', sep='\t', header=None, index=False, quoting=csv.QUOTE_NONE)

    if return_files != False:
        if truncate_amount == None:
            return gtf_coding
        else:
            return gtf_coding, gtf_truncated

"""
DESCRIPTION: GTF truncation function
"""
def execute_truncator(gtf, truncate_amount):

    #gtf[8] = gtf[8].replace({'\"':''}, regex=True)

    gtf['plus'] = gtf[[2,3,4,6,8]].apply(lambda x:
        (x[3] + truncate_amount) if x[2] == "exon" and x[3] + truncate_amount <= x[4] and x[6] == "+" and "exon_number \"1;" in x[8] else (
        "delete_this" if x[2] == "exon" and x[3] + truncate_amount > x[4] and x[6] == "+" and "exon_number \"1;" in x[8] else x[3]), axis=1)

    gtf['minus'] = gtf[[2,3,4,6,8]].apply(lambda x:
        (x[4] - truncate_amount) if x[2] == "exon" and x[3] <= x[4] - truncate_amount and x[6] == "-" and "exon_number \"1;" in x[8] else (
        "delete_this" if x[2] == "exon" and x[3] > x[4] - truncate_amount and x[6] == "-" and "exon_number \"1;" in x[8] else x[4]), axis=1)

    #remove exon1s that are too short
    gtf = gtf[~gtf['plus'].isin(['delete_this'])]
    gtf = gtf[~gtf['minus'].isin(['delete_this'])]

    #copy new coordinates back to original columns
    gtf[3] = gtf['plus']
    gtf[4] = gtf['minus']

    #remove placeholder columns
    gtf = gtf.drop(columns=['plus','minus'])
    return gtf

"""
DESCRIPTION: Parallelization scheme for truncation via chunking dataframe and parallel processing truncation on chunks
"""
def parallelize_truncator(func, *args):

    cores = cpu_count() #Number of CPU cores on your system

    #Split dataframe into chunks for processing by number of cores
    data_split = np.array_split(args[0], cores)
    pool = Pool(cores)

    func = partial(execute_truncator, truncate_amount=args[1])
    data = pd.concat(pool.map(func, data_split))

    pool.close()
    pool.join()

    return data

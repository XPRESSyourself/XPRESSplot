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
pd.options.mode.chained_assignment = None
import numpy as np
from multiprocessing import cpu_count, Pool
from functools import partial

"""
DESCRIPTION: GTF truncation function
"""
def execute_truncator(gtf, truncate_amount):

    gtf[8] = gtf[8].replace({'\"':''}, regex=True)

    gtf['plus'] = gtf[[2,3,4,6,8]].apply(lambda x:
        (x[3] + truncate_amount) if x[2] == "exon" and x[3] + truncate_amount <= x[4] and x[6] == "+" and "exon_number 1;" in x[8] else (
        "delete_this" if x[2] == "exon" and x[3] + truncate_amount > x[4] and x[6] == "+" and "exon_number 1;" in x[8] else x[3]),axis=1)

    gtf['minus'] = gtf[[2,3,4,6,8]].apply(lambda x:
        (x[4] - truncate_amount) if x[2] == "exon" and x[3] <= x[4] - truncate_amount and x[6] == "-" and "exon_number 1;" in x[8] else (
        "delete_this" if x[2] == "exon" and x[3] > x[4] - truncate_amount and x[6] == "-" and "exon_number 1;" in x[8] else x[4]),axis=1)

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

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
from multiprocessing import cpu_count, Pool

"""
DESCRIPTION: Check directory formatting
"""
#Check directory formatting
def check_directories(directory):

    #Check input directory name is formatted correctly and fix if necessary
    if directory.endswith('/'):
        pass
    else:
        directory += '/'

    return directory

"""
DESCRIPTION: Parallelize function on a chunk of a dataframe
"""
def parallelize(func, *args):

    cores = cpu_count() #Number of CPU cores on your system
    partitions = cpu_count() #Define as many partitions as you want

    data_split = np.array_split(args[0], partitions)
    pool = Pool(cores)

    if func == calculate_fc:
        func = partial(calculate_fc, label_comp=args[1], label_base=args[2])
    elif func == calculate_p:
        func = partial(calculate_p, label_comp=args[1], label_base=args[2], drop_index=args[3])
    elif func == threshold_util:
        func = partial(threshold_util, minimum=args[1], maximum=args[2])
    elif func == truncate:
        func = partial(truncate, truncate_amount=args[1])
    else:
        return

    data = pd.concat(pool.map(func, data_split))

    pool.close()
    pool.join()

    return data

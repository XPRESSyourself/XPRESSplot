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
import xpresstools as xp
%matplotlib inline

"""
utils.py tests
"""
#check_directories()
dir1 = '/Users/scripts/XPRESSyourself/XPRESStools/tests/'
dir2 = '/Users/scripts/XPRESSyourself/XPRESStools/tests'
truth_dir = '/Users/scripts/XPRESSyourself/XPRESStools/tests/'

dir1 = xp.check_directories(dir1)
dir2 = xp.check_directories(dir2)

assert dir1 == truth_dir, 'check_directories() failed'
assert dir2 == truth_dir, 'check_directories() failed'

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

"""IMPORT DEPENDENCIES"""
from setuptools import setup
import re
import os

__path__  =  os.path.dirname(os.path.realpath(__file__)) + '/'

"""Get version"""
with open(str(__path__) + 'xpressplot/__init__.py', 'r') as fd:
    version = re.search(r'^__version__\s*=\s*[\'"]([^\'"]*)[\'"]',
                        fd.read(), re.MULTILINE).group(1)

"""Setup arguments"""
setup(
    name = 'XPRESSplot',
    version = version,
    description = 'A toolkit for navigating and analyzing gene expression datasets',
    long_description = open('README.md').read(),
    long_description_content_type='text/markdown',
    author = 'Jordan Berg',
    author_email = 'jordan.berg@biochem.utah.edu',
    url = 'https://github.com/XPRESSyourself/XPRESSplot',
    packages = ['xpressplot'],
    exclude= ['tests','docs','recipes'],
    package_dir = {'xpressplot': 'xpressplot'},
    license = 'GPL-3.0',
    zip_safe = False,

    install_requires = [
        'pandas',
        'numpy',
        'scipy',
        'scikit-learn',
        'matplotlib<3.0.0,>=2.1.1',
        'seaborn',
        'plotly',
        'plotly_express'

    ],

    classifiers=[
        'Development Status :: 5 - Production/Stable',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics'
        ]
)

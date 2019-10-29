# ![XPRESSplot](https://raw.githubusercontent.com/XPRESSyourself/XPRESSplot/master/docs/content/xpressplot.png)


### A toolkit for navigating and analyzing gene expression datasets

[![Build Status](https://travis-ci.org/XPRESSyourself/XPRESSplot.svg?branch=master)](https://travis-ci.org/XPRESSyourself/XPRESSplot)
[![Documentation Status](https://readthedocs.org/projects/xpressplot/badge/?version=latest)](https://xpressplot.readthedocs.io/en/latest/?badge=latest)
[![codecov.io](https://codecov.io/gh/XPRESSyourself/XPRESSplot/XPRESSplot.svg?branch=master)](https://codecov.io/gh/XPRESSyourself/XPRESSplot)
[![PyPi Status](https://img.shields.io/pypi/v/XPRESSplot.svg)](https://pypi.org/project/XPRESSplot/)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/xpressplot/badges/version.svg)](https://anaconda.org/bioconda/xpressplot)
[![DOI](https://zenodo.org/badge/170940002.svg)](https://zenodo.org/badge/latestdoi/170940002)

##### Find documentation [here](https://xpressplot.readthedocs.io/en/latest/)

-----

### Development Notes:     
XPRESSplot supports Python 2.7 and >=3.5   

### Citation:       
```
Berg JA, et. al. (2019). XPRESSyourself: Enhancing and Automating the Ribosome
Profiling and RNA-Seq Analysis Toolkit. bioRxiv 704320; doi: https://doi.org/10.1101/704320
```

### Installation:   
```
pip install xpressplot
```

### Other Requirements:
- Tested on 64-bit Linux, compatible with Mac OS X
- Python and R are required
- If [PyPi](https://pip.pypa.io/en/stable/installing/) and [Conda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html?highlight=conda) are not already installed, these should be installed
- If using the interactive notebook provided, Jupyter needs to be [installed](https://jupyter.org/install) if not already

### QuickStart:   
- Download the repository and modify the [interactive Jupyter notebook](https://github.com/XPRESSyourself/XPRESSplot/blob/master/example_notebook.ipynb) to get started fast!   
- Read the instructions as you navigate through the code blocks for a guide on how to use the example code   
- Code blocks are run by selecting the block and pressing Shift + Enter   
- See [documentation](https://xpressplot.readthedocs.io/en/latest/) for more detailed instructions   

### Important Notes:    
* If working with XPRESSplot within an interactive notebook (i.e. Jupyter Notebook, Atom Hydrogen, etc), you must include the following line of code after importing XPRESSplot

```
import XPRESSplot as xp
%matplotlib inline
```

* Assumes all dataframes are i * j (or genes * samples, except in certain cases, see [documentation](https://xpressplot.readthedocs.io/en/latest/content/general-usage.html) for help)    

```
>>> geo.head()
name       GSM523242  GSM523243  GSM523244  GSM523245  GSM523246  GSM523247    ...     
GeneA    8.98104    8.59941    8.25395    8.72981    8.70794    8.10693    ...       
GeneB      5.84313    6.59168    8.27881    6.64005    4.65107    7.19090    ...       
GeneC       6.17189    5.73603    5.55673    5.69374    6.77618    5.84524    ...       
GeneD      6.97009    6.80003    5.56620    7.43816    7.36375    5.85687    ...       
GeneE   10.24611   10.13807    8.84743    9.72365   10.42940    9.17510    ...   

[5 rows x 145 columns]   
```

### Updates
Information on updates to the software can be found [here](https://github.com/XPRESSyourself/XPRESSplot/releases).

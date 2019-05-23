# ![XPRESSplot](https://raw.githubusercontent.com/XPRESSyourself/XPRESSplot/master/docs/content/xpressplot.png)


### A toolkit for navigating and analyzing gene expression datasets

[![Build Status](https://travis-ci.org/XPRESSyourself/XPRESSplot.svg?branch=master)](https://travis-ci.org/XPRESSyourself/XPRESSplot)
[![codecov.io](https://codecov.io/gh/XPRESSyourself/XPRESSplot/XPRESSplot.svg?branch=master)](https://codecov.io/gh/XPRESSyourself/XPRESSplot)
[![Documentation Status](https://readthedocs.org/projects/xpressplot/badge/?version=latest)](https://xpressplot.readthedocs.io/en/latest/?badge=latest)
[![PyPi Status](https://img.shields.io/pypi/v/XPRESSplot.svg)](https://pypi.org/project/XPRESSplot/)
[![DOI](https://zenodo.org/badge/170940002.svg)](https://zenodo.org/badge/latestdoi/170940002)

##### Find documentation [here](https://xpressplot.readthedocs.io/en/latest/)

-----

### Development Notes:
<b><i>XPRESSplot is still in beta production</i></b>  
Interactive scatter plotting functions are not currently tested    
XPRESSplot supports Python 2.7 and >=3.5   

### Citation:    
```
Berg, JA (2019). XPRESSyourself suite: Gene expression processing and analysis made easy. https://github.com/XPRESSyourself. DOI: 10.5281/zenodo.2581692.
```

### Installation:   
```
pip install xpressplot
```

### Other Requirements:
- Tested on 64-bit Linux, compatible with Mac OS X
- Python3 is recommended   
  - Current test cases build to:
    - Python2.7
    - Python3.5
    - Python3.6
    - Python3.7
- If [PyPi](https://pip.pypa.io/en/stable/installing/) and [Conda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html?highlight=conda) are not already installed, these should be installed
- If using this package to perform batch effect normalization or differential expression analysis, you must install [R](https://www.r-project.org/)   
- If using the interactive notebook provided, Jupyter needs to be [installed](https://jupyter.org/install) if not already

### QuickStart:   
Download the repository and modify the [interactive Jupyter notebook](https://github.com/XPRESSyourself/XPRESSplot/blob/master/example_notebook.ipynb) to get started quick!   
Read the instructions as you navigate through the code blocks for a guide on how to use the example code   
Code blocks are run by selecting the block and pressing Shift + Enter   
See [documentation](https://xpressplot.readthedocs.io/en/latest/) for more detailed instructions   

### Important Notes:    
* If working with XPRESSplot within an interactive notebook (i.e. Jupyter Notebook, Atom Hydrogen, etc), you must include the following line of code after importing XPRESSplot

```
import XPRESSplot as xp
%matplotlib inline
```

* Assumes all dataframes are columns=samples and rows=genes (except in certain cases, see [documentation](https://xpressplot.readthedocs.io/en/latest/content/general-usage.html) for help)    

```
>>> geo.head()
name       GSM523242  GSM523243  GSM523244  GSM523245  GSM523246  GSM523247    ...     
1007_s_at    8.98104    8.59941    8.25395    8.72981    8.70794    8.10693    ...       
1053_at      5.84313    6.59168    8.27881    6.64005    4.65107    7.19090    ...       
121_at       6.17189    5.73603    5.55673    5.69374    6.77618    5.84524    ...       
1294_at      6.97009    6.80003    5.56620    7.43816    7.36375    5.85687    ...       
1405_i_at   10.24611   10.13807    8.84743    9.72365   10.42940    9.17510    ...   

[5 rows x 145 columns]   
```

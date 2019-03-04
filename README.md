# ![XPRESStools](https://raw.githubusercontent.com/XPRESSyourself/XPRESStools/master/docs/content/xpresstools.png)


### A toolkit for navigating and analyzing gene expression datasets

[![Build Status](https://travis-ci.org/XPRESSyourself/XPRESStools.svg?branch=master)](https://travis-ci.org/XPRESSyourself/XPRESStools)
[![codecov.io](https://codecov.io/gh/XPRESSyourself/XPRESStools/XPRESStools.svg?branch=master)](https://codecov.io/gh/XPRESSyourself/XPRESStools)
[![Documentation Status](https://readthedocs.org/projects/xpresstools/badge/?version=latest)](https://xpresstools.readthedocs.io/en/latest/?badge=latest)
[![PyPi Status](https://img.shields.io/pypi/v/XPRESStools.svg)](https://pypi.org/project/XPRESStools/)
[![DOI](https://zenodo.org/badge/170940002.svg)](https://zenodo.org/badge/latestdoi/170940002)

##### Find documentation [here](https://xpresstools.readthedocs.io/en/latest/)

-----

### Development Notes:
<b><i>XPRESStools is still in beta production</i></b>  
Interactive scatter plotting functions are not currently tested    
XPRESStools supports Python 2.7 and >=3.5   

### Citation:    
```
Berg, JA (2019). XPRESSyourself suite: Gene expression processing and analysis made easy. https://github.com/XPRESSyourself.
```

### Installation:   
```
pip install xpresstools
```

### QuickStart:   
```
>>> import pandas as pd
>>> import xpresstools as xp   

>>> geo, metadata = xp.get_geo('GSE20916')
>>> geo = xp.rpm(geo)
```

### Important Notes:    
* If working with XPRESStools within an interactive notebook (i.e. Jupyter Notebook, Atom Hydrogen, etc), you must include the following line of code after importing XPRESStools

```
import XPRESStools as xp
%matplotlib inline
```

* Assumes all dataframes are columns=samples and rows=genes (except in certain cases, see documentation for help)    

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

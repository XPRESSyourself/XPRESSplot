[![Build Status](https://travis-ci.org/dwyl/learn-travis.svg?branch=master)](https://travis-ci.org/dwyl/learn-travis)

<p><img src="https://github.com/RiboPipe/ribopipe/blob/master/docs/content/xpressyourself.png" class="center" width="17%" height="17%" align="right">

<span style="color:green">XPRESS</span><span style="color:black">tools</span>   
A toolkit for navigating and analyzing gene expression datasets

<br />

Citation:    
```
Berg, JA (2019). XPRESSyourself suite: Gene expression processing and analysis made easy. https://github.com/XPRESSyourself.
```

<br />

Installation:   
```
conda install -c bioconda xpresstools
```

<br />

QuickStart:   
```
>>> import xpresstools as xp   

>>> geo = xp.get_geo('GSE20917')
>>> geo = xp.rpm(geo)
```

<br />

Important Notes:    
Assumes all dataframes are columns=samples and rows=genes (except in certain cases, see documentation for help)    

<br />

  +--------------+--------------+--------------+--------------+------------+
  |              |    Sample1   |   Sample2    |    Sample3   |    ...     |
  |--------------|---------------------------------------------------------|
  |    Gene1     |    1.07238   |   5.30284    |    3.29814   |    ...     |
  |--------------|---------------------------------------------------------|
  |    Gene2     |    9.32782   |   2.43879    |    6.32782   |    ...     |
  |--------------|---------------------------------------------------------|
  |    Gene3     |    0.23497   |   3.23475    |    4.98723   |    ...     |
  |--------------|---------------------------------------------------------|
  |     ...      |      ...     |     ...      |      ...     |    ...     |
  |--------------+--------------+--------------+--------------+------------+


<br />

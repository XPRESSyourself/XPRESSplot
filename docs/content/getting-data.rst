###############
Retrieving Data
###############
#explain drop rows in catenate and all variables for all


=====================
Import data from file
=====================
| **xpresstools.get_df(file_name, delimiter=',', low_memory=False, gene_axis='row')**
|
| Purpose:
| Get sequence dataframe from user-provided file.
|
| Assumptions:
| 1) Dataset does not contain axis labels (i.e. a column header for 'gene names')
| 2) Dataset only has gene names and sample_ids as column headers and row indices. Orientation is flexible, but needs to be specified in options if genes are not rows
| 3) If orientation is not default, it is then specified or else function will not be able to properly format the dataframe for downstream application
|
| Parameters:
| **file_name**: Full path of file to import into pandas dataframe
| **delimiter**: Delimiter type for importing file, default: ','
| **low_memory**: Specify memory limits for importing large files, default: False (allows for large imports)
| **gene_axis**: Orientiation of the data, where categorical data is either column-wise, (default: 'col') or row-wise ('row'). Case insensitive
|
| Examples:

.. ident with TABs
.. code-block:: python

  > import pandas as pd
  > import xpresstools as xp
  > data = xp.get_df('/path/to/data.csv')
  > data
                GSM523242 GSM523243 GSM523244 GSM523245 ...
  1007_s_at     8.98104   8.59941   8.25395   8.72981   ...
  1053_at       5.84313   6.59168   8.27881   6.64005   ...
  121_at        6.17189   5.73603   5.55673   5.69374   ...
  1294_at       6.97009   6.80003   5.56620   7.43816   ...
  ...           ...       ...       ...       ...       ...


=========================
Import metadata from file
=========================
**xpresstools.get_info(file_name, delimiter=",", axis="col", sample_ids=0, labels=1)**

Purpose:
Get sample metadata from user-provided file

Assumptions:
1) Data categories are not labeled
2) If orientation is not default, it is then specified or else function will not be able to properly format the dataframe for downstream application

Parameters:
**file_name**: full path of file to import into pandas dataframe
**delimiter**: delimiter type for importing file, default: ','
**axis**: Orientiation of the data, where categorical data is either column-wise, (default: 'col') or row-wise ('row'). Case insensitive
**sample_ids**: Column or row number where sample IDs are found (default: 0)
**labels**: Column or row number where categorical label data are found (default: 1)

Examples:

.. code-block:: python

  > import pandas as pd
  > import xpresstools as xp
  > metadata = xp.get_info('/path/to/metadata.csv')
  > metadata
    0	        1
  0	GSM523242	mucosa_normal_colon
  1	GSM523243	mucosa_normal_colon
  2	GSM523244	mucosa_adenoma
  3	GSM523245	colonic_crypt_epithelial_cells_normal_colon
  4	GSM523246	mucosa_normal_colon


===================
Import data from GEO dataset
===================
**xpresstools.get_geo(geo_id, output_info=False)**


.. code-block:: python

  > import pandas as pd
  > import xpresstools as xp
  > data, metadata = xp.get_geo('GSE20916')
  > data
                GSM523242	GSM523243	GSM523244	GSM523245	GSM523246	GSM523247	GSM523248	GSM523249	GSM523250
  1007_s_at	    8.98104	  8.59941	  8.25395	  8.72981	  8.70794	  8.10693	  8.96594	  9.35994	  8.40191
  1053_at	      5.84313	  6.59168	  8.27881	  6.64005	  4.65107	  7.19090	  6.24983	  6.98251	  7.41631
  121_at	      6.17189	  5.73603	  5.55673	  5.69374	  6.77618	  5.84524	  6.02640	  6.49465	  6.18855
  1294_at	      6.97009	  6.80003	  5.56620	  7.43816	  7.36375	  5.85687	  6.26649	  6.44538	  6.52518
  1405_i_at	    10.24611	10.13807	8.84743	  9.72365	  10.42940	9.17510	  7.89429	  7.81446	  7.57219
  1438_at	      2.18618	  2.88067	  7.08605	  2.19053	  2.09528	  6.60998	  2.37048	  2.47629	  2.18834
  1487_at	      9.19325	  9.59890	  9.27702	  8.79801	  8.94495	  9.05543	  8.89486	  9.29561	  9.25947
  1552256_a_at	10.27871	9.94561	  10.58918	9.51612	  10.27936	10.69953	11.47816	11.33376	11.47500
  1552257_a_at	9.82578	  9.67882	  9.60929	  9.95520	  9.15464	  9.40873	  9.45647	  9.72854	  9.46014

  > metadata
    0	        1
  0	GSM523242	mucosa_normal_colon
  1	GSM523243	mucosa_normal_colon
  2	GSM523244	mucosa_adenoma
  3	GSM523245	colonic_crypt_epithelial_cells_normal_colon
  4	GSM523246	mucosa_normal_colon

###############
Retrieving Data
###############
#explain drop rows in catenate and all variables for all


=====================
Import data from file
=====================
| **xpresstools.get_df ( file_name, delimiter=',', low_memory=False, gene_axis='row' )**
|
| Purpose:
| Get sequence dataframe from user-provided file.
|
| Assumptions:
|   - Dataset does not contain axis labels (i.e. a column header for 'gene names')
|   - Dataset only has gene names and sample_ids as column headers and row indices. Orientation is flexible, but needs to be specified in options if genes are not rows
|   - If orientation is not default, it is then specified or else function will not be able to properly format the dataframe for downstream application
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
| **xpresstools.get_info ( file_name, delimiter=",", axis="col", sample_ids=0, labels=1 )**
|
| Purpose:
| Get sample metadata from user-provided file
|
| Assumptions:
| * Data categories are not labeled
| * If orientation is not default, it is then specified or else function will not be able to properly format the dataframe for downstream application
|
| Parameters:
| **file_name**: full path of file to import into pandas dataframe
| **delimiter**: delimiter type for importing file, default: ','
| **axis**: Orientiation of the data, where categorical data is either column-wise, (default: 'col') or row-wise ('row'). Case insensitive
| **sample_ids**: Column or row number where sample IDs are found (default: 0)
| **labels**: Column or row number where categorical label data are found (default: 1)
|
| Examples:

.. ident with TABs
.. code-block:: python

  > import pandas as pd
  > import xpresstools as xp
  > metadata = xp.get_info('/path/to/metadata.csv')
  > metadata
      0         1
  0   GSM523242 mucosa_normal_colon
  1   GSM523243 mucosa_normal_colon
  2   GSM523244 mucosa_adenoma
  3   GSM523245 colonic_crypt_epithelial_cells_normal_colon
  ... ...       ...

===================
Import data from GEO dataset
===================
| **xpresstools.get_geo ( geo_id, output_info=False )**
|
| Purpose:
| Get sample data and metadata from a GEO database
|
| Parameters:
| **geo_id**: GEO ID for dataset of interest, input is case insensitive (ex: GSE20716)
| **output_info**: Output long-form metadata to txt file if True (default: False)
|
| Examples:

.. ident with TABs
.. code-block:: python

  > import pandas as pd
  > import xpresstools as xp
  > data, metadata = xp.get_geo('GSE20916')
  > data
                GSM523242 GSM523243 GSM523244 GSM523245 ...
  1007_s_at     8.98104   8.59941   8.25395   8.72981   ...
  1053_at       5.84313   6.59168   8.27881   6.64005   ...
  121_at        6.17189   5.73603   5.55673   5.69374   ...
  1294_at       6.97009   6.80003   5.56620   7.43816   ...
  ...           ...       ...       ...       ...       ...
  > metadata
      0         1
  0   GSM523242 mucosa_normal_colon
  1   GSM523243 mucosa_normal_colon
  2   GSM523244 mucosa_adenoma
  3   GSM523245 colonic_crypt_epithelial_cells_normal_colon
  ... ...       ...

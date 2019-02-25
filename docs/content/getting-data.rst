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
| Returns:
| **data**: Pandas dataframe with data matrix
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
| Returns:
| **metadata**: Pandas dataframe with metadata
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

============================
Import data from GEO dataset
============================
| **xpresstools.get_geo ( geo_id, output_info=False )**
|
| Purpose:
| Get sample data and metadata from a GEO database
|
| Parameters:
| **geo_id**: GEO ID for dataset of interest, input is case insensitive (ex: GSE20716)
| **output_info**: Output long-form metadata to txt file if True (default: False)
|
| Returns:
| **data**: Pandas dataframe with data matrix
| **metadata**: Pandas dataframe with metadata
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

===========================
Catenate raw counts files
===========================
| **xpresstools.catenate_files ( directory, file_suffix='txt', save_file=None, delimiter='\t', drop_rows=0 )**
|
| Purpose:
| Compiles expression counts from multiple files into one table. For example, HTSeq-count outputs each alignment file's counts as a separate count file. This module will collect all single count files and compile them into a single count table.
|
| Assumptions:
|   - File length of each is the same and ordered the same (same genes in the same order)
|   - Files to parse are expected to be header-less and column[0] should be gene identifiers and column[1] should be expression values
|
| Parameters:
| **directory**: Path to directory containing raw counts files (only tested currently with HTSeq-count output files)
| **file_suffix**: Common suffix of all count files (default: 'txt'). This feature is useful for modification if there other files in the directory that are not count files, as if they do not contain the same suffix, they will not be used in the function.
| **save_file**: Include if you want the resulting counts table saved for later use (default: None)
| **delimiter**: Delimiter style for expression files, will also output files if saved in this same format
| **drop_rows**: Number of rows to drop from the end of each count file. HTSeq-count provides 5 lines of summary statistics at the end of each file, so for HTSeq-count files, use drop_rows=5
|
| Returns:
| **count_table**: Pandas dataframe with the catenated counts. Samples are along columns, genes are along rows
|
| Examples:

.. ident with TABs
.. code-block:: python

  > counts = xp.catenate_files(count_dir, file_suffix='counts.txt', drop_rows=5)
  > counts
          S1_counts.txt S2_counts.txt S3_counts.txt S4_counts.txt
  Gene1   66            59            1             82
  Gene2   35            0             7             72
  Gene3   20            70            87            78
  Gene4   96            7             93            38
  ...     ...           ...           ...           ...

============================
Drop samples
============================
| **xpresstools.drop_samples ( data, ids )**
|
| Purpose:
| Drop samples by sample IDs -- pass in a list of names
|
| Assumptions:
|   - Dataframe axes have been properly formatted (samples are columns, genes are rows)
|
| Parameters:
| **data**: Dataframe containing expression data
| **ids**: List of sample IDs to remove from the dataframe
|
| Returns:
| **data**: Pandas dataframe with modified data matrix
|
| Examples:

.. ident with TABs
.. code-block:: python

  > data
                GSM523242 GSM523243 GSM523244 GSM523245 ...
  1007_s_at     8.98104   8.59941   8.25395   8.72981   ...
  1053_at       5.84313   6.59168   8.27881   6.64005   ...
  121_at        6.17189   5.73603   5.55673   5.69374   ...
  1294_at       6.97009   6.80003   5.56620   7.43816   ...
  ...           ...       ...       ...       ...       ...
  > data = xp.drop_samples(data, metadata, ['GSM523244'])
  > data
                GSM523242 GSM523243 GSM523245 ...
  1007_s_at     8.98104   8.59941   8.72981   ...
  1053_at       5.84313   6.59168   6.64005   ...
  121_at        6.17189   5.73603   5.69374   ...
  1294_at       6.97009   6.80003   7.43816   ...
  ...           ...       ...       ...       ...

============================
Drop label
============================
| **xpresstools.drop_label ( data, info, label )**
|
| Purpose:
| Drop samples by label group name
|
| Assumptions:
|   - Dataframe axes have been properly formatted (samples are columns, genes are rows)
|   - Only one string is given to drop per call instance of function
|
| Parameters:
| **data**: Dataframe containing expression data
| **info**: Dataframe containing sample information data
| **label**: Name of sample type to drop (string)
|
| Returns:
| **data**: Pandas dataframe with modified data matrix
|
| Examples:

.. ident with TABs
.. code-block:: python

  > data
                GSM523242 GSM523243 GSM523244 GSM523245 ...
  1007_s_at     8.98104   8.59941   8.25395   8.72981   ...
  1053_at       5.84313   6.59168   8.27881   6.64005   ...
  121_at        6.17189   5.73603   5.55673   5.69374   ...
  1294_at       6.97009   6.80003   5.56620   7.43816   ...
  ...           ...       ...       ...       ...       ...
  > data = xp.drop_label(data, metadata, 'mucosa_adenoma')
  > data
                GSM523242 GSM523243 GSM523245 ...
  1007_s_at     8.98104   8.59941   8.72981   ...
  1053_at       5.84313   6.59168   6.64005   ...
  121_at        6.17189   5.73603   5.69374   ...
  1294_at       6.97009   6.80003   7.43816   ...
  ...           ...       ...       ...       ...

============================
Keep labels
============================
| **xpresstools.keep_labels ( data, info, label_list=None )**
|
| Purpose:
| Keep samples by list of label names
|
| Assumptions:
|   - Dataframe axes have been properly formatted (samples are columns, genes are rows)
|   - Labels provided are in list format
|
| Parameters:
| **data**: Dataframe containing expression data
| **info**: Dataframe containing sample information data
| **labels**: List of sample types to keep
|
| Returns:
| **data**: Pandas dataframe with modified data matrix
|
| Examples:

.. ident with TABs
.. code-block:: python

  > data
                GSM523242 GSM523243 GSM523244 GSM523245 ...
  1007_s_at     8.98104   8.59941   8.25395   8.72981   ...
  1053_at       5.84313   6.59168   8.27881   6.64005   ...
  121_at        6.17189   5.73603   5.55673   5.69374   ...
  1294_at       6.97009   6.80003   5.56620   7.43816   ...
  ...           ...       ...       ...       ...       ...
  > data = xp.keep_labels(data, metadata, ['mucosa_normal_colon', 'mucosa_adenoma'])
  > data
                GSM523242 GSM523243 GSM523244 ...
  1007_s_at     8.98104   8.59941   8.25395   ...
  1053_at       5.84313   6.59168   8.27881   ...
  121_at        6.17189   5.73603   5.55673   ...
  1294_at       6.97009   6.80003   5.56620   ...
  ...           ...       ...       ...       ...

======================================
Rename dataframe column names
======================================
| **xpresstools.rename_cols ( data, converters )**
|
| Purpose:
| Rename column names using dataframe
|
| Parameters:
| **data**: Dataframe to rename column names
| **converters**: Dataframe where column 0 contains old names and column 1 contains new names
|
| Returns:
| **data**: Pandas dataframe with modified data matrix
|
| Examples:

.. ident with TABs
.. code-block:: python

  > data
                GSM523242 GSM523243 GSM523244 GSM523245 ...
  1007_s_at     8.98104   8.59941   8.25395   8.72981   ...
  1053_at       5.84313   6.59168   8.27881   6.64005   ...
  121_at        6.17189   5.73603   5.55673   5.69374   ...
  1294_at       6.97009   6.80003   5.56620   7.43816   ...
  ...           ...       ...       ...       ...       ...
  > conversion_table
      0         1
  0   GSM523242 normal
  1   GSM523244 adenoma
  2   GSM523245 normal
  > data = xp.rename_cols(data, conversion_table)
  > data
                normal    GSM523243 adenoma   normal ...
  1007_s_at     8.98104   8.59941   8.25395   8.72981   ...
  1053_at       5.84313   6.59168   8.27881   6.64005   ...
  121_at        6.17189   5.73603   5.55673   5.69374   ...
  1294_at       6.97009   6.80003   5.56620   7.43816   ...
  ...           ...       ...       ...       ...       ...

======================================
Rename dataframe row names
======================================
| **xpresstools.rename_rows ( data, converters, label='index' )**
|
| Purpose:
| Rename values in an index (row names) or a column
|
| Parameters:
| **data**: Dataframe to rename rows of a column
| **converters**: Dataframe where column 0 contains old names and column 1 contains new names
| **label**: Name of column to convert names; if 'index' is provided, will rename the index of the dataframe
|
| Returns:
| **data**: Pandas dataframe with modified data matrix
|
| Examples:

.. ident with TABs
.. code-block:: python

  > data
                GSM523242 GSM523243 GSM523244 GSM523245 ...
  1007_s_at     8.98104   8.59941   8.25395   8.72981   ...
  1053_at       5.84313   6.59168   8.27881   6.64005   ...
  121_at        6.17189   5.73603   5.55673   5.69374   ...
  1294_at       6.97009   6.80003   5.56620   7.43816   ...
  ...           ...       ...       ...       ...       ...
  > conversion_table
      0         1
  0   1007_s_at Gene1
  1   121_at Gene2
  > data = xp.rename_rows(data, conversion_table)
  > data
                GSM523242 GSM523243 GSM523244 GSM523245 ...
  Gene1         8.98104   8.59941   8.25395   8.72981   ...
  1053_at       5.84313   6.59168   8.27881   6.64005   ...
  Gene2         6.17189   5.73603   5.55673   5.69374   ...
  1294_at       6.97009   6.80003   5.56620   7.43816   ...
  ...           ...       ...       ...       ...       ...

.. ident with TABs
.. code-block:: python

  > data
       gene_names    GSM523242 GSM523243 GSM523244 GSM523245 ...
  0    1007_s_at     8.98104   8.59941   8.25395   8.72981   ...
  1    1053_at       5.84313   6.59168   8.27881   6.64005   ...
  2    121_at        6.17189   5.73603   5.55673   5.69374   ...
  3    1294_at       6.97009   6.80003   5.56620   7.43816   ...
  ...  ...           ...       ...       ...       ...       ...
  > conversion_table
      0         1
  0   1007_s_at Gene1
  1   121_at Gene2
  > data = xp.rename_rows(data, conversion_table, label='gene_names')
  > data
       gene_names    GSM523242 GSM523243 GSM523244 GSM523245 ...
  0    Gene1         8.98104   8.59941   8.25395   8.72981   ...
  1    1053_at       5.84313   6.59168   8.27881   6.64005   ...
  2    Gene2         6.17189   5.73603   5.55673   5.69374   ...
  3    1294_at       6.97009   6.80003   5.56620   7.43816   ...
  ...  ...           ...       ...       ...       ...       ...

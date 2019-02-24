#############
General Usage
#############

XPRESStools is intended as a all-in-one toolkit and interface for analysis of sequencing data

======================================
Sequence data
======================================
Required format for all functions (unless otherwise noted in documentation).

.. code-block:: python

  >>> geo.head()
  name       GSM523242  GSM523243  GSM523244  GSM523245  GSM523246  GSM523247    ...
  1007_s_at    8.98104    8.59941    8.25395    8.72981    8.70794    8.10693    ...
  1053_at      5.84313    6.59168    8.27881    6.64005    4.65107    7.19090    ...
  121_at       6.17189    5.73603    5.55673    5.69374    6.77618    5.84524    ...
  1294_at      6.97009    6.80003    5.56620    7.43816    7.36375    5.85687    ...
  1405_i_at   10.24611   10.13807    8.84743    9.72365   10.42940    9.17510    ...


===========
Metadata
===========
Required format for all functions (unless otherwise noted in documentation).

.. code-block:: python

  >>> geo.head()
      0 	        1
  0 	GSM523242 	mucosa_normal_colon_1 (micro)
  1 	GSM523243 	mucosa_normal_colon_2 (micro)
  2 	GSM523244 	mucosa_adenoma_3 (micro)
  3 	GSM523245 	colonic_crypt_epithelial_cells_normal_colon_4 ...
  4 	GSM523246 	mucosa_normal_colon_5 (micro)

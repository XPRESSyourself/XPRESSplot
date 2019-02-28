###########################
Reference Building
###########################

======================
Truncate Reference GTF
======================
| **xpresstools.truncate ( input_gtf, truncate_amount=45, save_coding_path='./', save_truncated_path='./', sep='\t', return_files=False )**
|
| Purpose:
| Create a GTF reference file with only protein coding genes and the first n nucleotides of each first exon
| By default, will output modified GTFs to the current working directory
|
| Assumptions:
|   - Input file is a properly formatted GTF file
|   - Protein coding transcripts are denoted by 'protein_coding in the final column of the GTF'
|
| Parameters:
| **input_gtf**: GTF reference file path and name
| **truncate_amount**: Number of nucleotides to truncate from the 5' end of each exon 1. Considers strandedness. (default: 45). If value equals None, truncation step will not be performed and only a coding-only reference will be created
| **save_coding_path**: Location to save coding-only GTF, set to None to avoid outputting file
| **save_truncated_path**: Location to save coding-only truncated GTF
| **sep**: Separator type of the GTF file (generally always tab-delimited)
| **return_files**: If not false, will return coding only and truncated/coding only dataframes
|
| Returns:
| **coding**: Pandas dataframe with coding_only reference entries if return_files=True
| **truncated**: Pandas dataframe with coding_only/truncated reference entries if return_files=True
|
| Examples:

.. ident with TABs
.. code-block:: python

  > import pandas as pd
  > import xpresstools as xp
  > gtf_coding, gtf_truncated = xp.truncate ('./transcripts.gtf', truncate_amount=45, save_coding_path=None, save_truncated_path=None, sep='\t', return_files=True)
  > gtf_orig
    0 1       2           3     4     5 6 7 8
  0 1 havana  gene        11869 14409 . + . gene_id "ENSG00000223972"; gene_version "5"; g...
  1 1 havana  transcript  11869 14409 . + . gene_id "ENSG00000223972"; gene_version "5"; t...
  2 1 havana  exon        11869 12227 . + . gene_id "ENSG00000223972"; exon "5"; t...
  3 1 havana  exon        12613 12721 . + . gene_id "ENSG00000223972"; exon "1"; t...
  4 1 havana  CDS         13221 14409 . + . gene_id "ENSG00000223972"; gene_version "5"; t...
  > gtf_truncated
    0 1       2           3     4     5 6 7 8
  0 1 havana  gene        11869 14409 . + . gene_id "ENSG00000223972"; gene_version "5"; g...
  1 1 havana  transcript  11869 14409 . + . gene_id "ENSG00000223972"; gene_version "5"; t...
  2 1 havana  exon        11869 12227 . + . gene_id "ENSG00000223972"; exon "5"; t...
  3 1 havana  exon        12658 12721 . + . gene_id "ENSG00000223972"; exon "1"; t...
  4 1 havana  exon        13221 14409 . + . gene_id "ENSG00000223972"; gene_version "5"; t...

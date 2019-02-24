"""
XPRESStools
A toolkit for navigating and analyzing gene expression datasets
alias: xpresstools

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

"""
IMPORT DEPENDENCIES
"""
import re
import csv
import pandas as pd
pd.options.mode.chained_assignment = None
from .utils import check_directories
from .utils_truncator import execute_truncator, parallelize_truncator

"""
DESCRIPTION: Collate HTseq counts files

VARIABLES:
file_list= List of files with the path names appended to each file to be collated into a single count table
gene_column= Column location in all count files of gene names
gene_column= Column location in all count files of samples
sep= Separator of counts files

ASSUMPTIONS:
No headers are included in the count files
"""
def count_table(file_list, gene_column=0, sample_column=1, sep='\t'):

    #Read in first count file to get gene names
    df = pd.read_csv(str(file_list[gene_column]), sep=sep, comment='#', header=None)
    pos_starter = [gene_column,sample_column]
    colname = df.columns[pos_starter]
    df = df[colname]

    #For the rest of the files in the file list, add the counts for that sample only
    for f in file_list[1:]:
        df_pull = pd.read_csv(str(f), sep=sep, comment='#', header=None)
        df = pd.concat([df, df_pull[df_pull.columns[sample_column]]], axis=1)
        del df_pull

    #Final formatting clean up of table
    df_counts = df.copy()
    del df
    df_counts = df_counts.set_index(0)
    del df_counts.index.name

    #Remove path and file suffix from each file's name before adding as column names to table
    c = 0
    for x in file_list:
        file_list[c] = x[(x.rfind('/')+1):(x.find('.'))]
        c += 1
    df_counts.columns = file_list

    return df_counts

"""
DESCRIPTION: Create a GTF reference file with only protein coding genes and the first n nucleotides of each first exon

RETURNS: Return the truncated, coding only table when running the function

METHODS: Considers strandedness
Multiprocesses chunks of a dataframe on cores of computer

VARIABLES:
input_gtf= GTF reference file path and name
truncate_amount= Number of nucleotides to truncate from the 5' end of each exon 1 (considers strandedness)
save_coding_path= Location to save coding-only GTF, set to None to avoid outputting file
save_truncated_path= Location to save coding-only truncated GTF
sep= Separator type of the GTF file (generally always tab-delimited)

ASSUMPTIONS:
Input file is a properly formatted GTF file
Protein coding transcripts are denoted by 'protein_coding in the final column of the GTF'
"""
def truncate(input_gtf, truncate_amount=45, save_coding_path='./', save_truncated_path='./', sep='\t', return_files=False):

    save_coding_path = check_directories(save_coding_path)
    save_truncated_path = check_directories(save_truncated_path)

    #Import gtf reference file to
    if str(input_gtf).endswith('.gtf'):
        gtf = pd.read_csv(str(input_gtf), sep=sep, header=None, comment='#', low_memory=False)
    else:
        raise Exception('Error: A GTF-formatted file was not provided')

    #Get only protein_coding coordinates
    gtf_coding = gtf[gtf.iloc[:, 8].str.contains('protein_coding') == True]

    #Save to .gtf file (tsv)
    if save_coding_path != None:
        gtf_coding.to_csv(str(save_coding_path) + 'transcripts_coding.gtf', sep='\t', header=None, index=False, quoting=csv.QUOTE_NONE)

    gtf_coding_c = gtf_coding.copy()

    print("Multiprocessing reference chunks -- this may take a while...")
    gtf_truncated = parallelize_truncator(execute_truncator, gtf_coding_c, truncate_amount)

    if save_truncated_path != None:
        gtf_truncated.to_csv(str(save_truncated_path) + 'transcripts_coding_truncated.gtf', sep='\t', header=None, index=False, quoting=csv.QUOTE_NONE)

    if return_files != False:
        return gtf_coding, gtf_truncated

"""
DESCRIPTION: Convert row names (genes) of dataframe using GTF as reference for new name

RETURNS: Return the renamed dataframe

VARIABLES:
data= Dataframe to convert rows names
gtf= Path and name of gtf reference file
orig_name_label= Label of original name (usually a \"gene_id \"')
orig_name_location= Position in last column of GTF where relevant data is found (i.e. 0 would be the first sub-string before the first comma, 3 would be the third sub-string after the second comma before the third comma)
new_name_label= Label of original name (usually \"gene_name \")
new_name_location= Position in last column of GTF where relevant data is found (i.e. 0 would be the first sub-string before the first comma, 3 would be the third sub-string after the second comma before the third comma)
refill= In some cases, where common gene names are unavailable, the dataframe will fill the gene name with the improper field of the GTF. In this case, specify this improper string and these values will be replaced with the original name
sep= GTF delimiter (usually tab-delimited)
"""
def convert_names_gtf(data, gtf, orig_name_label='gene_id \"', orig_name_location=0, new_name_label='gene_name \"', new_name_location=1, refill=None, sep='\t'):

    #Import reference GTF
    gtf = pd.read_csv(str(gtf),sep=sep,comment='#', low_memory=False, header=None)

    #Parse out old and new names from GTF
    gtf_genes = gtf.loc[gtf[2] == 'gene']
    gtf_genes[orig_name_label] = gtf[8].str.split(';').str[orig_name_location]
    gtf_genes[new_name_label] = gtf[8].str.split(';').str[new_name_location]
    gtf_genes[orig_name_label] = gtf_genes[orig_name_label].map(lambda x: x.lstrip(str(orig_name_label)).rstrip('\"').rstrip(' '))
    gtf_genes[new_name_label] = gtf_genes[new_name_label].map(lambda x: x.lstrip(str(new_name_label)).rstrip('\"').rstrip(' '))
    gtf_genes = gtf_genes[[orig_name_label,new_name_label]].copy()

    #Create dictionary
    gene_dict = {}

    if refill != None:
        for index, row in gtf_genes.iterrows():
            if row[1] == str(refill):
                gene_dict[row[0]] = row[0]
            else:
                gene_dict[row[0]] = row[1]
    else:
        gene_dict[gtf_genes[orig_name_label]] = gtf_genes[new_name_label]

    #Replace old gene names/ids with new
    data_names = data.copy()
    data_names[new_name_label] = data_names.index.to_series().map(gene_dict).fillna(data_names.index.to_series())
    data_names = data_names.set_index(new_name_label)
    del data_names.index.name

    return data_names

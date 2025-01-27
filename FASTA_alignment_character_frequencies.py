#!/usr/bin/env python
# coding: utf-8

### Author: Cory Dunn
### Provided courtesy of GRO Biosciences
### Version 1.4

version = '1.4'

import pandas as pd
import argparse

from Bio import SeqIO

if __name__ == "__main__" :

    # Collect input from user

    ap = argparse.ArgumentParser()
    ap.add_argument('-i','--input_file',required=True,type=str,help="Input alignment in FASTA format.\n")
    ap.add_argument('-o','--output_file_prefix',required=True,type=str,help="Output files prefix.\n")
    ap.add_argument('-r','--reference_sequence',required=True,type=str,help="Reference sequence used for comparison and site labelling.\n")
    ap.add_argument('-d','--deletion_types',default='',type=str,help="Deletion of gaps '(-)', Xs ('X'), asterisks ('*'), or reference sequence from final dataframe count and frequency output: Enter codes 'g', 'x', 'a', or 'r', respectively following this flag.\n")
    ap.add_argument('-f','--focus_on_changed',required=False,default='n',type=str,help="Enter code 'y' to also generate analyses focused only on those columns/sites with any frequencies that differ from reference\n")

    args = vars(ap.parse_args())

    reference_acc = args['reference_sequence']
    alignfile = args['input_file']
    output_file_prefix = args['output_file_prefix']
    deletion_types = args['deletion_types']
    focus_flag = args['focus_on_changed']

    print('Running FASTA_alignment_character_frequencies.py ' + version)
    print('Selected reference is ' + reference_acc)
    print('Alignment file is ' + alignfile)
    print('Output file prefix is ' + output_file_prefix)
    if 'g' in deletion_types:
        print('Gaps not included in calculations.')
    if 'x' in deletion_types:
        print('Stop codons not included in calculations.')

    # Obtain sequence data from FASTA file

    alignment_record_name_list = []
    sequence_records = []

    for record in SeqIO.parse(alignfile,"fasta"):

        alignment_record_name_list.append(record.name)
        record_x_toward_seq_dataframe = list(record.seq)
        record_x_toward_seq_dataframe_UPPER = [x.upper() for x in record_x_toward_seq_dataframe]
        record_x_toward_seq_dataframe_ASCII = [ord(ele) for sub in record_x_toward_seq_dataframe_UPPER for ele in sub] # Ensure int to avoid memory problems
    
        if str(record.name) == reference_acc:
            dataframe_columns_from_reference = record_x_toward_seq_dataframe_UPPER
        
        if 'r' in deletion_types and str(record.name) == reference_acc:
            print('Reference sequence not included in calculations.')

        else: sequence_records.append(record_x_toward_seq_dataframe_ASCII)

    numbering_of_columns = list(range(1,len(dataframe_columns_from_reference)+1,1))

    # Build the dataframe of ASCII numbers representing sequence characters

    amino_acid_ASCII_DF = pd.DataFrame(sequence_records, columns = numbering_of_columns)
    amino_acid_ASCII_DF.dropna(inplace = True)
    amino_acid_ASCII_DF = amino_acid_ASCII_DF.astype(int)

    #

    del sequence_records

    # Convert back to characters

    amino_acid_sequences_DF = amino_acid_ASCII_DF.map(chr)
    del amino_acid_ASCII_DF

    # Generate and populate frequency and count dataframes

    index_string_list = ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','X','Y','-','*']
    amino_acid_freq_DF = pd.DataFrame(amino_acid_sequences_DF,index = index_string_list)
    amino_acid_num_DF = amino_acid_freq_DF.copy()

    for (columnName, columnData) in amino_acid_sequences_DF.items():
        amino_acid_freq_DF[columnName] = amino_acid_sequences_DF[columnName].value_counts(normalize=True)

    for (columnName, columnData) in amino_acid_sequences_DF.items():
        amino_acid_num_DF[columnName] = amino_acid_sequences_DF[columnName].value_counts()

    amino_acid_freq_DF = amino_acid_freq_DF.fillna(0)
    amino_acid_freq_DF = amino_acid_freq_DF.round(2)
    amino_acid_num_DF = amino_acid_num_DF.fillna(0)
    amino_acid_num_DF = amino_acid_num_DF.round(2)

    new_column_names_list = []

    for i in range(len(numbering_of_columns)):
        new_column_names_list.append(dataframe_columns_from_reference[i] + str(numbering_of_columns[i]))

    amino_acid_num_DF.columns = new_column_names_list

    # Handle odd characters

    if 'g' in deletion_types:
        amino_acid_num_DF.drop('-', inplace = True)
        
    if 'x' in deletion_types:
        amino_acid_num_DF.drop('X', inplace = True)

    if 'a' in deletion_types:
        amino_acid_num_DF.drop('*', inplace = True)

    amino_acid_perc_DF = 100 * amino_acid_num_DF.div(amino_acid_num_DF.sum())
    amino_acid_perc_DF = amino_acid_perc_DF.round(1)

    # Generate output files

    amino_acid_num_DF.to_csv(output_file_prefix + '_charnum.csv')
    amino_acid_perc_DF.to_csv(output_file_prefix + '_charfreqs.csv')

    if focus_flag == 'y':

        print('Down-selected based upon any difference from reference.')

        list_of_column_indices_that_changed = []

        colcount = 0
        for (colname,colval) in amino_acid_perc_DF.items():
            column_values_to_check = pd.Series(colval.values)
            big_dataframe_index = amino_acid_perc_DF.index
            column_values_to_check.index = big_dataframe_index

            value_to_seek = column_values_to_check[dataframe_columns_from_reference[colcount]]
            
            if value_to_seek != 100: 

                list_of_column_indices_that_changed.append(colcount) 
            
            colcount += 1

        truncate_perc_DF = amino_acid_perc_DF.iloc[:,list_of_column_indices_that_changed]
        truncate_num_DF = amino_acid_num_DF.iloc[:,list_of_column_indices_that_changed]

        truncate_perc_DF.to_csv(output_file_prefix + '_interesting_charfreqs.csv')
        truncate_num_DF.to_csv(output_file_prefix + '_interesting_charnum.csv')
        
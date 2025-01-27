FASTA_alignment_character_frequencies
Version: 1.4

----

Author:

Cory Dunn
Contact: cory.david.dunn@gmail.com
Provided courtesy of GRO Biosciences

----

License:

GPLv3

----


This software computes and reports character counts and frequencies from a FASTA multiple sequence alignment. Counts and frequencies are indexed on a selected reference found within that alignment.


----

Requirements:

This script is implemented using Python 3 (current version tested under Python version 3.12.8) 

Dependencies: 

Biopython (tested under version 1.85),
Pandas (tested under version 2.2.3)

----

Usage: FASTA_alignment_character_frequencies.py [-h] -i INPUT_FILE -o
                                                OUTPUT_FILE_PREFIX -r
                                                REFERENCE_SEQUENCE
                                                [-d DELETION_TYPES]
                                                [-f FOCUS_ON_CHANGED]

Required arguments:

  -i INPUT_FILE, --input_file INPUT_FILE
                        Input alignment in FASTA format.

  -o OUTPUT_FILE_PREFIX, --output_file_prefix OUTPUT_FILE_PREFIX
                        Output files prefix.

  -r REFERENCE_SEQUENCE, --reference_sequence REFERENCE_SEQUENCE
                        Reference sequence used for comparison and site
                        labelling. The reference sequence _must_ be in the alignment being 				analyzed, although counts of reference characters can be masked, as 
			seen below, using '-d r'.

Optional arguments:

  -d DELETION_TYPES, --deletion_types DELETION_TYPES
                        Deletion of gaps '(-)', Xs ('X'), asterisks ('*'), or
                        reference sequence from final dataframe count and
                        frequency output: Enter codes 'g', 'x', 'a', or 'r',
                        respectively following this flag.
  -f FOCUS_ON_CHANGED, --focus_on_changed FOCUS_ON_CHANGED
                        Enter code 'y' to also generate analyses focused only
                        on those columns/sites with any frequencies that


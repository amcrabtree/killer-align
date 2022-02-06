#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May  9 10:04:13 2020

@author: angela
"""
#!/usr/bin/env python

import csv
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import sys

depth_path = sys.argv[1]
consensus_path = sys.argv[2] 
output_filename = sys.argv[3]

# import file containing killer assay data
with open(depth_path, 'r') as f:
    reader = csv.reader(f, delimiter='\t')
    # get all the rows as a list
    depth_data_row_data = list(reader)
    # transform data into numpy array and save start and end nucleotide positions
    depth_data = np.array(depth_data_row_data).astype(str)
    start_pos = int(depth_data[0,1])-1
    end_pos = int(depth_data[-1,1])-1

# parse fasta file
consensus_record = SeqIO.read(consensus_path, "fasta")
consensus_header = str(consensus_record.id)+"_trunc"
consensus_desc = str(consensus_record.description)+" truncated assembly"
consensus_seq = str(consensus_record.seq)

# outupt truncated consensus
trunc_seq = consensus_seq[start_pos:end_pos]
output_record = SeqRecord(Seq(trunc_seq,),
    id=consensus_header,
    description=consensus_desc,)
SeqIO.write(output_record, output_filename, "fasta")

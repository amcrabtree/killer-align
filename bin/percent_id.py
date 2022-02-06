#!/usr/bin/env python3

from Bio import SeqIO
import sys
from Bio import Align

trunc_consensus_path = sys.argv[1] 
trunc_ref_path = sys.argv[2]
og_ref_path = sys.argv[3]
read_depth_file_len = int(sys.argv[4])-1

# parse fasta files
consensus_record = SeqIO.read(trunc_consensus_path, "fasta")
consensus_seq = str(consensus_record.seq)
ref_record = SeqIO.read(trunc_ref_path, "fasta")
ref_seq = str(ref_record.seq)
og_ref_record = SeqIO.read(og_ref_path, "fasta")
og_ref_seq = str(og_ref_record.seq)

# print output
aligner = Align.PairwiseAligner()
aln_len = len(consensus_seq)
ref_len = len(og_ref_seq)
cov = round((100*read_depth_file_len/ref_len), 1)
percent_id = 100*aligner.align(ref_seq, consensus_seq).score/aln_len
percent_id = round(percent_id, 1)
print("Percent Identity:")
print(percent_id)
print("Percent of reference covered with >0 reads:")
print(cov)
print("Consensus sequence length (bp):")
print(aln_len)

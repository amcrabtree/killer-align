# killer-align
Aligns raw Illumina reads to various reference sequences, calls variants, and creates consensus sequence. Optimized for dsRNA elements in yeast. 

More specifically, this command line pipeline aligns raw reads to a reference sequence with BWA-MEM, calls variants and creates a consensus sequence with bcftools, and aligns consensus to reference with MUSCLE. The pipeline generates a consensus sequence fasta file and an html report containing graphs of read depth across the reference sequence and a numeric summary of alignment metrics. 

### **Dependencies**

* Java, BWA, Samtools, Python3, BCFtools, Muscle, (optionally) [fastp](https://github.com/OpenGene/fastp) 

<b>CLI Usage</b>
```
killer-align.sh sample_name -1 for_read.fastq.gz -2 rev_read.fastq.gz -r reference_genome.fasta -o output_folder
```

**Options**

flag | description
------------ | -------------
-s	[ARG]	| (required) sample name, used for file naming purposes
-1	[ARG]	| (required) forward read file
-2	[ARG]	| (required) reverse read file
-r	[ARG]	| (required) reference genome file
-o	[ARG]	| (required) destination of output folder
-t		| tests program with associated test file
-h		| help (print options)

<p>&nbsp;</p>

**Files Required for Program to Run:**

1. **Forward read file** – File contains raw NGS reads in forward orientation; “R1” should be in the filename somewhere (ex: “NCYC190_S3_L001_R1_001.fastq.gz”); can be .fastq or .fastq.gz

2. **Reverse read file** – File contains raw NGS reads in reverse orientation; “R2” should be in the filename somewhere (ex: “NCYC190_S3_L001_R2_001.fastq.gz”); can be .fastq or .fastq.gz

3. **Reference genome file** – File contains the nucleotide sequence of your reference; must be in fasta format; you don’t have to worry about creating index files since the script automatically creates them

**Program Output**

|    Filename or Folder         |     Description                                                                                                |
|-------------------------------|----------------------------------------------------------------------------------------------------------------|
| SAMPLE_REF_aligned.fasta          | A fasta-formatted pairwise alignment using MUSCLE of the sample and the reference.                         |
| SAMPLE_REF_consensus_trunc.fasta | A consensus sequence that has been truncated on the 5' and 3' ends where there is 0 coverage. CAUTION: my python script does not cut out areas of zero coverage in the middle of the sequence.                   |
| SAMPLE_REF_muscle_input.fasta    | This is a fasta-formatted file containing the truncated sample consensus sequence and the reference sequence (before alignment)                                                                                   |
| SAMPLE_REF_read_depth.jpeg       | A graph of sample read depth across the reference genome                                                    |
| SAMPLE_REF_read_depth.txt        | A tab-delimited list containing the read depth at each nucleotide position (if >0 read depth)               |
| SAMPLE_REF_sorted.bam            | A file that contains information on where each read is aligned along the reference genome, sorted from start of reference genome to end of reference genome, in a binary format (not human readable)              |
| SAMPLE_REF_sorted.bam.bai        | An index file for the corresponding bam file; this makes it easier for some programs to utilize the information in the bam file                                                                                   |
| SAMPLE_REF.bam                   | A file that contains information on where each read is aligned along the reference genome, in a binary format (not human readable)                                                                                |
| SAMPLE_REF.sam                   | A file that contains information on where each read is aligned along the reference genome, in a human-readable format                                                                                             |
| troubleshooting                  | A folder for storing troubleshooting information. This will have several text files containing output from different command line apps. If there are errors or the output is messed up, try browsing these files. |

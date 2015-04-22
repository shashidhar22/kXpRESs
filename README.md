#KXpress
=====

KXpress is a k-mer based rna quatification tool, which uses an alignment free methodology to quantify gene expression values from RNASeq data. KXpress is written in python and uses KAnalyze (P.Audano and F.Vannberg, 2014) to count k-mers from RNASeq data and map them to their respective transcripts.

#####Requirements:
  KAnalyze 0.9.7
  Python3.4
  PyFasta
  BioPython


#####Execution:

######1. Indexing

An index of the transcriptome has to be built before the analysis. The index step requires a reference transcript file and the k-mer length. k-mer size can take a maximum value of 31. For any given k-mer length, the indexing step has to be performed only once.
To run the index step:
```python
kexpress.py -r <reference_transcript_file> -m index -d <probe/fasta> -k <kmer length>  -o <index_file_name>
  -r The reference trancripts file, can either be a reference cDNA file or a microarray probe list, for quick analysis. 
  -m The mode of execution (index/express)
  -d The type of referehce file, (Fasta or tab delimited file with probe name and sequence.
  -k k-mer length
  -o Output file path
```
  
  
######2. Expression quantification

Once an index is built, express mode is used to analyze the RNASeq data and quantify the gene expression. 
To run the express step:
```python
  kexpress.py -f  <fastq_file> -o <output_file> -m express -i <index_file_path>
    -f The path to the fastq file
    -o Output file path
    -m The mode of execution (index/express)
    -i The index file path
```

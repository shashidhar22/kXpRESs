#! /usr/local/bin/python3.4
""" KXpress is a k-mer based rna expression quantification tool developed by
Shashidhar Ravishankar, at the Vannberg Lab, Georgia Institute of Technology.
    There are two modules to KXpress, namely the index and the express module.
The index module, takes as input the reference cDNA fasta file or probes list,
and creates a hash based index, for quick k-mer look up. Indexing needs to be
done once for a given k-mer size.
    The express module, takes as input the index file, and a sequence file in
fastq or fastqgz format and using KAnalyze, performs a quick k-mer count and 
calculates gene expression in term of k-mers per kilobase per million k-mers 
mapped."""
import os
import sys
import csv
import time
import math
import glob
import pickle
import profile
import argparse
import subprocess
import statistics as st
import logging
from operator import itemgetter
from pyfasta import Fasta
from Bio import SeqIO
from multiprocessing import Pool
from collections import Counter
from itertools import repeat
import numpy as np
import collections
global kanalyze
kanalyze = os.path.abspath(os.path.dirname(__file__))+'/kanalyze-0.9.7/'


transcript_index = dict()
transcript_rescue = dict()
transcript_count = dict()
transcript_kpkm = dict()

def non_rescue(kc_chunk,index):
    for lines in open(index):
        line = lines.strip().split('\t')
        if int(line[0]) in kc_chunk:
            for genes in line[2].split(','):
                try:
                    transcript_index[int(genes)].append(kc_chunk[int(line[0])])
                except KeyError:
                    transcript_index[int(genes)] = [kc_chunk[int(line[0])]]
    return
def express(kc_chunk,index,transcript_order,transcript_length,kmer_count,klen,seq_len):
    transcript_count = dict()
    for lines in open(index):
        line = lines.strip().split('\t')
        if int(line[1]) == 1:
            try:
                if int(line[0]) in kc_chunk:
                    transcript_index[int(line[2])].append(kc_chunk[int(line[0])])
                    transcript_count[int(line[2])] += int(klen)
                    #transcript_kpkm[int(line[2])] = (sum(transcript_index[int(line[2])])*10**9)/float(( transcript_length[transcript_order[int(line[2])]] - klen + 1) * kmer_count * seq_len)
                    transcript_kpkm[int(line[2])] = (sum(transcript_index[int(line[2])])*10**9)/float((transcript_count[int(line[2])] - klen + 1) * kmer_count * seq_len)
                else:
                    continue
            except KeyError:
                if int(line[0]) in kc_chunk:
                    transcript_index[int(line[2])] = [kc_chunk[int(line[0])]]
                    transcript_count[int(line[2])] = int(klen)
                    #transcript_kpkm[int(line[2])] = (sum(transcript_index[int(line[2])])*10**9)/float((transcript_length[transcript_order[int(line[2])]] - klen + 1) * kmer_count * seq_len)
                    transcript_kpkm[int(line[2])] = (sum(transcript_index[int(line[2])])*10**9)/float((transcript_count[int(line[2])] - klen + 1) * kmer_count * seq_len)
                else:
                    continue
        else:
            if int(line[0]) in kc_chunk:
                res_trans = [int(gene) for gene in line[2].split(',')]
                transcript_rescue[int(line[0])] = [res_trans,kc_chunk[int(line[0])]]
            else:
                continue
    return

def get_next(iterator):
    try:
        return(next(iterator).strip().split('\t'))
    except StopIteration:
        return(None)

def rescue(transcript_index,transcript_rescue,transcript_order,transcript_length,transcript_kpkm,kmer_count,seq_len,klen):
     logging.info('Rescue algorithm started')
     transcript_ec = dict()
     sigma_kpkm = sum(transcript_kpkm.values())
     out = open('log','w')
     print(len(transcript_rescue))
     for count,kmers in enumerate(transcript_rescue):	#loop through all kmers that have more than one transcript index
         if len(transcript_rescue[kmers][0]) <= 10:
             gen = [transcript_order[val] for val in transcript_rescue[kmers][0]]
             out.write(str(kmers)+'\t'+','.join(gen)+'\t'+str(transcript_rescue[kmers][1])+'\n')
             for indices in transcript_rescue[kmers][0]:	#loop through transcript index of kmer
                 try:
                     opt = transcript_rescue[kmers][1] * transcript_kpkm[indices]/sigma_kpkm
                     transcript_index[indices].append(opt)
                     transcript_count[indices] += int(klen)
                     #transcript_kpkm[indices] = (sum(transcript_index[indices])*10**9)/float((transcript_length[transcript_order[indices]] - klen + 1) * kmer_count * seq_len)
                     transcript_kpkm[indices] = (sum(transcript_index[indices])*10**9)/float((transcript_count[indices] - klen + 1) * kmer_count * seq_len)
                     sigma_kpkm = sum(transcript_kpkm.values())
                 except KeyError:
                     #continue
                     opt = transcript_rescue[kmers][1]
                     transcript_index[indices] = [opt]
                     transcript_count[indices] = int(klen)
                     #transcript_kpkm[indices] = (sum(transcript_index[indices])*10**9)/float((transcript_length[transcript_order[indices]] - klen + 1) * kmer_count * seq_len)
                     transcript_kpkm[indices] = (sum(transcript_index[indices])*10**9)/float((transcript_count[indices] - klen + 1) * kmer_count * seq_len)
                     sigma_kpkm = sum(transcript_kpkm.values())
     logging.info('Number of transcripts rescued :' +str(len(transcript_index)))
     logging.info('Rescue completed')
     out.close()
     return(transcript_index)

def merge(file1,file2,database,index,output):
     result = csv.writer(open(output+str(index)+'.mkx','w'),delimiter='\t')
     with open(file1) as f1, open(file2) as f2:
         line1 = get_next(f1) 
         line2 = get_next(f2) 
         while line1 != None and line2 != None:
             if int(line1[0]) < int (line2[0]):
                 result.writerow([line1[0],line1[1],line1[2]])
                 line1 = get_next(f1) 
             elif int(line1[0]) == int(line2[0]):
                 merge_list = line1[2].split(',') + line2[2].split(',')
                 result.writerow([line1[0],str(len(merge_list)),','.join(merge_list)])
                 line1 = get_next(f1) 
                 line2 = get_next(f2) 
             else:
                 result.writerow([line2[0], line2[1], line2[2]])
                 line2 = get_next(f2) 
         while line1 != None:
             result.writerow([line1[0],line1[1],line1[2]])
             line1 = get_next(f1) 
         while line2 != None:
             result.writerow([line2[0],line2[1],line2[2]])
             line2 = get_next(f2) 
     return(output+str(index)+'.mkx')

def lazy_function(fasta_file, chunk_size=1000):
    """Function performs chunking of fasta file"""
    seq = True
    while seq:
        chunk = list()
        while len(chunk) < chunk_size:
            try:
                seq = next(fasta_file)
            except StopIteration:    #account for end of file
                seq = None
            if seq is None:
               break
            chunk.append([seq.id,str(seq.seq).upper()])
        if chunk:
            yield chunk
    return

def kmerize(arguments):
    """Implements a quick k-mer algorithm, to k-merize a given sequence and
    store them as base of four integers. It also maps the k-mers generated
    to the respective transcripts."""
    fasta_file = arguments[0]
    index = arguments[1]
    klen = arguments[2]
    database = arguments[3]
    power = arguments[4]
    output = arguments[5]
    order = index * (10**power)
    transcript_length = dict()
    transcript_kmers = dict()
    transcript_order = dict()
    mask = (1 << (klen*2)) -1         #Create bitmask of 1's
    kmers = 0
    for count, lines in enumerate(fasta_file):
        tindex = count + order
        sequence = lines[1]
        header = lines[0]
        kmer = 0
        bit_counter = 1
        for nuc in sequence:
            kmer =  kmer << 2         #left shift k-kmer 
            if nuc == 'A':            #add the value of the character using bitwise OR
                kmer = kmer | 0x00   
            elif nuc == 'C':
                kmer = kmer | 0x01
            elif nuc == 'G':
                kmer = kmer | 0x02
            elif nuc == 'T':
                kmer = kmer | 0x03
            else:
                bit_counter = 0
            if bit_counter == klen:   #if length equals k-mer length, store k-mer
                kmers += 1
                try:
                    transcript_kmers[kmer & mask].append(tindex)  #k-mer to transcript index mapping
                except KeyError:
                    transcript_kmers[kmer & mask] = [tindex]
            else:
                bit_counter += 1
        logging.debug('Indexing '+header+ '; Line number '+str(tindex))
        transcript_order[tindex] = header         #transcript index to transcript name mapping
        transcript_length[header] = len(sequence) #transcript to transcript length mapping
        count += 1
    temp_file = csv.writer(open(output+str(index)+'.tkx','w'),delimiter='\t')
    for values in sorted(transcript_kmers.keys()):
        #recurrence = Counter(transcript_kmers[values])   #to calculate the number of occurences of k-mers in trasncripts
        #prior = [str(recurrence[gene]) for gene in transcript_kmers[values]]
        temp_file.writerow([str(values),str(len(transcript_kmers[values])),str(','.join(str(x) for x in transcript_kmers[values]))])
    return (transcript_order, transcript_length)

def probe_list(database, line_count, klen):
    """Method to create index from a given probe list"""
    #must be updated 
    logging.basicConfig(level=logging.INFO)
    dbfile = csv.reader(open(database),delimiter='\t')
    transcript_index = dict()
    transcript_length = dict()
    kmers = 0
    perc = [float(i) for i in range(10,110,10)]
    logging.info('Indexing reference probes file')
    for count, lines in enumerate(dbfile):
        if len(lines) > 20:
            kindex = 0
            for i in range(len(lines[17])):
                if len(lines[17][i:i+klen]) == klen:
                    if dec(lines[17][i:i+klen],klen) not in transcript_index:
                        transcript_index[dec(lines[17][i:i+klen],klen)] = [[lines[2],kindex,0]]
                    else:
                        transcript_index[dec(lines[17][i:i+klen])].append([lines[2],kindex,0])
                    kindex += 1
                    kmers += 1
            transcript_length[lines[2]] = [len(lines[17]),kindex]
    kmer_str = str(kmers)
    logging.info('Reference indexing completed; %s kmers indexed to transcripts' %kmers)
    return (transcript_index, transcript_length)

#Create transcript index, maintains refseq id, kmer and kmer index
def transcript_list(database, line_count,klen,output):
    """Method to create index from a fasta file"""
    logging.basicConfig(level=logging.INFO)
    transcript_length = dict()
    transcript_kmers = dict()
    transcript_order = dict()
    logging.info('Indexing reference fasta file')
    fasta_file = SeqIO.parse(open(database),'fasta')
    sequences = 0
    for lines in fasta_file:
        sequences += 1
    fasta_file = SeqIO.parse(open(database),'fasta')
    fasta_handle = True
    fasta_list = list()
    count = 0
    power = int(math.log10(sequences/4)+1)
    for fragments in lazy_function(fasta_file,int(sequences/4)):
        fasta_list.append(fragments)
    pool = Pool(5) 
    results = pool.map(kmerize,zip(fasta_list,range(5),repeat(klen),
                       repeat(database),repeat(power),repeat(output)))
    for i in range(len(results)):            #redundancy to be rectified
        for k, v in results[i][1].items():
            transcript_length[k] = v
        for k, v in results[i][0].items():
            transcript_order[k] = v
    logging.info('Number of transcripts read : %s ' %count)
    logging.info('Reference indexing completed; \
                  %s kmers indexed to transcripts' %len(transcript_order))
    temp_list = glob.glob(output+'*.tkx')
    i = 0
    logging.info('Merging index files')
    while len(temp_list) > 1:
        merged_file = merge(temp_list[0],temp_list[1],database, i,output)
        file = temp_list.pop(0)
        os.remove(file)
        file = temp_list.pop(0)
        os.remove(file)
        temp_list.insert(0,merged_file)
        i += 1
    logging.info('Merging complete')
    return (transcript_length, transcript_order)

def reporter(transcript_index,transcript_length,klen,output_file,seq_number,seq_len,transcript_order):
    logging.info('Generating report')
    output_csv = csv.writer(open(output_file,'w'),delimiter='\t')
    output_log = csv.writer(open(output_file+'.log','w'),delimiter='\t')
    output_csv.writerow(['RefseqID','Length','Number of kmers aligned','KPKM','RPKM'])
    for refseqs in sorted(transcript_index):
        output_log.writerow([transcript_order[refseqs],'\t'.join([str(vals) for vals in transcript_index[refseqs]])])
        if refseqs != '' and transcript_length[transcript_order[refseqs]] > int(klen):
            kpkm = 0
            kpkm = (sum(transcript_index[refseqs])*10**9)/ float((transcript_length[transcript_order[refseqs]]-int(klen)+1) * seq_number * seq_len)
            rpkm = kpkm / float(seq_len - int(klen) +1)
            output_csv.writerow([transcript_order[refseqs], str(transcript_length[transcript_order[refseqs]]),str(len(transcript_index[refseqs])),str(kpkm),str(rpkm)])
        else:
            continue
    logging.info('Report generation complete')
    return

def kxpress(refer, seqfile, loglevel,klen,output,mode,index,db_type,type,seq_len):
    numeric_level = getattr(logging, loglevel.upper(), None)
    if not isinstance(numeric_level, int):
         raise ValueError('Invalid log level: %s' % loglevel)
    logging.basicConfig(level=numeric_level, format='%(levelname)s:%(asctime)s:%(message)s', datefmt='%m/%d/%Y;%I:%M:%S')
    db_count = 0
    seq_file_count = 0
    if mode == 'index':
        if db_type == 'fasta':
            with open(refer) as dbfile:
                for line in SeqIO.parse(dbfile,"fasta"):
                    db_count += 1
            flanking_kmers = transcript_list(refer, db_count, int(klen),output)
            pickle.dump(flanking_kmers, open(output+'.kxp','wb'))
            logging.info('Indexing complete.')
        elif db_type == 'probes':
            flanking_kmers = probe_list(refer, db_count, int(klen))
            pickle.dump(flanking_kmers, open(output+'.kxp','wb'))
            logging.info('Indexing complete.')
    elif mode == 'express':
        logging.info('Parsing input files')
        logging.info('Loading index')
        transcript_length, transcript_order = pickle.load(open(index+'.kxp','rb'))
        if type == 'fastq':
            fastq_parser = open(seqfile) 
            seq_id = fastq_parser.readline()
            seq_len = len(fastq_parser.readline()) -1
            fastq_parser.close()
            logging.info('Kmerizing seqeunce files')
            kmerize = subprocess.Popen([kanalyze+'count','-r','-m','dec','-d','6','-k',klen,
                                        '-f','fastq','-o',os.path.abspath(os.path.dirname(__file__))+'/tmp.kc',seqfile],
                                        stdout=subprocess.PIPE, shell=False)
            kmerize.wait();
            logging.info('Kmers counted')
            kmer_file = os.path.abspath(os.path.dirname(__file__))+'/tmp.kc'
        elif type == 'kc':
            kmer_file = seqfile
        logging.info('Retriving index')
        kmer_index = glob.glob(index+'*.mkx')[0]
        logging.info('Merging transcritps and k-mers')
        kmer_count = 0
        for lines in open(kmer_file):
            if int(lines.split('\t')[1]) > 3:
                kmer_count += 1
        logging.info('Total kmers found : '+str(kmer_count))
        kmer_csv = open(kmer_file)
        kmer = next(kmer_csv).split('\t')
        kmer_count = 0
        kc_chunk = dict()
        while kmer_csv:
            if int(kmer[1]) > 3:
                kmer_count += 1
                if len(kc_chunk) <= 10000000:
                    kc_chunk[int(kmer[0])] = int(kmer[1])
                    try:
                        kmer = next(kmer_csv).split('\t')
                    except StopIteration:
#                        non_rescue(kc_chunk,kmer_index)
                        express(kc_chunk,kmer_index,transcript_order,transcript_length,kmer_count,int(klen),seq_len)
                        break
                else:
                    kc_chunk[int(kmer[0])] = int(kmer[1])
#                    non_rescue(kc_chunk,kmer_index)
                    express(kc_chunk,kmer_index,transcript_order,transcript_length,kmer_count,int(klen),seq_len)
                    kc_chunk = dict()
            else:
                kmer = next(kmer_csv).split('\t')
        logging.info('Merging complete')
        #os.remove(kmer_file)
        rescue_expression = rescue(transcript_index, transcript_rescue,transcript_order,transcript_length,transcript_kpkm,seq_len,kmer_count,int(klen))
        full_transcripts = 0
        incomplete_transcripts = 0
        unknown_kmers = 0
        for transcripts in transcript_index:
            if transcripts != '':
                full_transcripts += 1
            elif transcripts == '':
                unknown_kmers = len(transcript_calls[transcripts])
        reporter(transcript_index,transcript_length,klen,output,kmer_count,seq_len,transcript_order)
        logging.info('Number of complete transcripts found = %s' %full_transcripts)
        logging.info('Number of kmers with no refseq id = %s' %unknown_kmers)
    return

if __name__ == '__main__':
    parser = argparse.ArgumentParser(prog='kXpRESs')
    parser.add_argument('-r','--reference',dest='refer', type=str, help='Transcript reference file')
    parser.add_argument('-f','--seqfile',dest='seqfile', type=str, help='Fastq file')
    parser.add_argument('-k','--kmer_length',dest='klen', type=str, help='Kmer length',default='20')
    parser.add_argument('-o','--output_file',dest='output',type=str, help='Output file path',default='transcripts.kx')
    parser.add_argument('-m','--mode',dest='mode',type=str,choices=['index','express'],help='Enter mode of execution, run index first if the database hasnt been indexed yet')
    parser.add_argument('-d','--db_type',dest='db_type',type=str,choices=['probes','fasta'],help='Database type being used for indexing')
    parser.add_argument('-t','--type', dest='file_type',type=str,choices=['fastq','kc'],help='Input file type') 
    parser.add_argument('-i','--index',dest='index',type=str,help='Path to indexed file')
    parser.add_argument('-l','--log',type=str,default='info',choices=['info','debug','error'],help='Verbosity parameter')
    parser.add_argument('-s','--seqlen',type=int,default=36,help='Read length, must be provided when runnnig using kc file')
    parser.add_argument('-v','--version',action='version',version='%(prog)s 0.9.0')
    args = parser.parse_args()
    kxpress(args.refer, args.seqfile,args.log,args.klen,args.output,args.mode,args.index,args.db_type,args.file_type,args.seqlen)

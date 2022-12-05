#!/usr/bin/env python
import os
import sys
import inspect
import re
import argparse
import random


 # realpath() will make your script run, even if you symlink it :)
cmd_folder = os.path.realpath(os.path.abspath(os.path.split(inspect.getfile( inspect.currentframe() ))[0]))
if cmd_folder not in sys.path:
     sys.path.insert(0, cmd_folder)

 # use this if you want to include modules from a subfolder
cmd_subfolder = os.path.realpath(os.path.abspath(os.path.join(os.path.split(inspect.getfile( inspect.currentframe() ))[0],"bin")))
if cmd_subfolder not in sys.path:
     sys.path.insert(0, cmd_subfolder)


from fastaIO import FastaReader,FastaWriter,SequenceUtility
from PopGenomeDefinitionIO import *

parser = argparse.ArgumentParser(description="""           
Description
-----------
This script generates an empty TE landscape for manual editing
""",formatter_class=argparse.RawDescriptionHelpFormatter,
epilog="""
Authors
-------
    Robert Kofler
""")

    

parser.add_argument("--chassis", type=str, required=True, dest="ref_fasta", default=None, help="the chassis, i.e. the sequence into which TEs will be inserted; a fasta file")
parser.add_argument("--te-seqs", type=str, required=True, dest="te_fasta", default=None, help="TE sequences in a fasta file")
parser.add_argument("--N", type=int, required=True, dest="pop_size", default=None, help="the number of haploid genomes")
parser.add_argument("--insert-count", type=int, required=False, dest="insert_count", default=1, help="number of empty default insertions")
parser.add_argument("--output", type=str, required=True, dest="output", default=None, help="the output file")
parser.add_argument("--min-distance", type=int, required=False, dest="min_distance", default=1000, help="minimum distance between TE insertions")
args = parser.parse_args()

# load the genomes of the chasis and TEs
header,seq=SequenceUtility.load_chasis(args.ref_fasta)
tetuples=FastaReader.readAllTuples(args.te_fasta)


genomesize=len(seq)
tefamilycount=len(tetuples)
insertions=args.insert_count

mindist=args.min_distance
maxdist=int(genomesize/(insertions+1))
if(maxdist<=mindist):
    raise ValueError("Genome of "+str(genomesize)+ " is too small for "+str(insertions)+ " insertions with a min-distance between insertions of "+str(mindist))


pw=PopGenDefinitionWriter(args.output,args.pop_size)
pw.write_chasis_info(header,seq)
pw.write_header(tetuples)

counter=1
for i in range(0,insertions):
     dist=random.randint(mindist,maxdist)
     counter+=dist
     pw.write_empty(counter)
pw.close()
    
    






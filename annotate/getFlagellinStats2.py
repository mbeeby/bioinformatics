#!/usr/bin/python
'''

Gets flagellin stats from a list of proteins regardless of genome content.
Meant for use with a list of proteins from many different organisms.

@author: Morgan Beeby 

'''

import sys
import os
import glob
import numpy
import operator
import re
import argparse
import fasta
from annotateGenes import *
from Bio import SeqIO
             
def extractInserts(outputSeqs, curOrg, genomeFile):
    '''
    Extracts all flagellin insert sequences and appends them to a FASTA file.
    
    curOrg is of type organismClass
    genomeFile is the input FASTA genome file
    rename is bool: if True, rename as a integer for output as a lookup table
    '''
    seqs       = fasta.getProteinSeqs(genomeFile)
    
    
    sorted_proteins = sorted(curOrg.families["FliC"].items(), key=operator.itemgetter(1))[::-1]
    for protein,data in sorted_proteins:
        outputSeqs[protein] = \
            seqs[protein][data.notes['insert_start']-1:data.notes['insert_end']-1]
   
def extractConservedEnds(outputSeqs, curOrg, genomeFile):
    '''
    Extracts all flagellin conserved N- and C-terminal domains, 
    concatenates them both into a single sequence with excised insert,
    and appends them to a FASTA file.
    
    curOrg is of type organismClass
    genomeFile is the input FASTA genome file
    rename is bool: if True, rename as a integer for output as a lookup table
    '''
    
    seqs       = fasta.getProteinSeqs(genomeFile)
    
    sorted_proteins = sorted(curOrg.families["FliC"].items(), key=operator.itemgetter(1))[::-1]
    for protein,data in sorted_proteins:
        outputSeqs[protein] = \
            seqs[protein][0:(data.notes['insert_start']-1)] + \
            seqs[protein][(data.notes['insert_end']-1):-1]
            
def appendLookupTable(curOrg, filename):
    '''
    Outputs a lookup table of protein accession, renames, taxID, and organism
    '''
    f = open(filename,'a')

    sorted_proteins = sorted(curOrg.families["FliC"].items(), key=operator.itemgetter(1))[::-1]
    for protein,data in sorted_proteins:
        f.write(protein)
        f.write("\t")
        f.write(protein)
        f.write("\t")
        f.write(repr(curOrg.organism))
        f.write("\t")
        f.write(repr(curOrg.taxID))
        f.write("\n")
    f.close()


def main():
    global rename
    global renameCount
    parser = argparse.ArgumentParser(description = "Get statistics on flagellins")
    parser.add_argument('proteins', metavar = '.gbff genome dir', type=str, \
        help = "The genome directory (.gbff, etc)")
    parser.add_argument('SVGout', metavar = 'SVGTools output filename', type=str, \
        help = "The SVGTools output filename")
    parser.add_argument('--insertFASTA', type=str, \
        help = "The output FASTA file to store flagellin inserts in")
    parser.add_argument('--conservedFASTA', type=str, \
        help = "FASTA of concatenated conserved termini output filename")
    parser.add_argument('--renameStart', type=int, \
        help = "starting number for renames", default=0)
    parser.add_argument('--hmmdir', \
        help = "Location of HMMs (overrides HMM_DIR environment variable)")
    
    args   = parser.parse_args()
     
    #-----------
    if "HMM_DIR" in os.environ:
        HMMDir = os.environ["HMM_DIR"]
    if args.hmmdir:
        HMMDir = args.hmmdir
    else:
        if not "HMM_DIR" in os.environ:
            sys.exit("ERROR: HMM_DIR not specified on command line or environment var") 
        
    tempFAA = ".temp"+str(random.random()*1000000)+".faa"

    insertSeqs    = {}
    endSeqs       = {}
    organismCount = 0
    
    # Setup the output file to be empty:
    outputFile = open(args.insertFASTA, 'w')
    outputFile.close()
    outputFile = open(args.conservedFASTA, 'w')
    outputFile.close()
  
    curOrg   = organismClass(args.proteins, HMMDir, annotateAll=False)
    curOrg.annotateFlagellins()
    
    # Now output the proteins:
    if len(curOrg.families["FliC"]) > 0:
        sorted_proteins = sorted(curOrg.families["FliC"].items(), key=operator.itemgetter(1))[::-1]
        for protein,data in sorted_proteins:
            sys.stdout.write("FliC\t"+str(protein)+" ("+repr(data.notes["insert_length"])+"; "+str(data.evidence.keys())+")\n")
            
    if len(curOrg.families["FlgL"]) > 0:
        sorted_proteins = sorted(curOrg.families["FlgL"].items(), key=operator.itemgetter(1))[::-1]
        for protein,data in sorted_proteins:
            sys.stdout.write("FlgL\t"+str(protein)+" ("+str(data.evidence.keys())+")\n")
            
    if args.insertFASTA:
        extractInserts(insertSeqs, curOrg, args.proteins)
    if args.conservedFASTA:
        extractConservedEnds(endSeqs, curOrg, args.proteins)
        
    # ----------------
    if args.insertFASTA:
        fasta.appendFASTA(args.insertFASTA, insertSeqs)
    if args.conservedFASTA:
        fasta.appendFASTA(args.conservedFASTA, endSeqs)   
    

if __name__ == '__main__':
    main()

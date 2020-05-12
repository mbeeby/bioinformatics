#!/usr/bin/python
#
# Takes a two-column file (ignoring empty lines)
# with desired name and PFxxxxx number and (1) works
# out current PFxxxxx version, then outputs as
# <desired name>.hmm

import sys
import numpy
import os
import glob
import operator
import re
import argparse
import annotateGenes

def getProteinSeqs(FAAFile):
    # Takes a FASTA file and returns a dictionary with keys as protein IDs
    # and values as protein lengths

    header_re = re.compile('>([^\s]+)')
    seq_re    = re.compile('([A-Z]+)')

    seqs = {}

    f       = open(FAAFile, 'r') 
    lines   = f.readlines()
    curProt = ""
    for line in lines:
        m = header_re.match(line)
        if m:
            curProt = m.group(1)
            seqs[curProt] = ""
        else:
            n = seq_re.match(line)
            if not n:
                sys.exit("ERROR: Unexpected line! "+line+"\n")
                
            seqs[curProt] = seqs[curProt] + n.group(1)

    return seqs

def main():
    parser = argparse.ArgumentParser(description = "Compare hook and filament lengths")
    parser.add_argument('genomesDir', metavar = 'Genome dir', type=str, \
        help = "The FASTA genome directory")
    parser.add_argument('--hmmdir', \
        help = "Location of HMMs (overrides HMM_DIR environment variable)")
    args   = parser.parse_args()
    
    dataFile    = "lengths.csv"
    
    #-----------
    if "HMM_DIR" in os.environ:
        HMMDir = os.environ["HMM_DIR"]
    if args.hmmdir:
        HMMDir = args.hmmdir
    else:
        if not "HMM_DIR" in os.environ:
            sys.exit("ERROR: HMM_DIR not specified on command line or environment var")
            
    # ----------------- 

    genomeFiles = glob.glob(args.genomesDir+"/*.faa")
    organisms   = {}


    familyList = ["FlgE", "FliC"]
    for family in familyList:
        sys.stdout.write("\t"+family)
    sys.stdout.write("\n")
    organismCount = 0
    for genomeFile in genomeFiles:

        # Scrape organism from genome file
        seqs = getProteinSeqs(genomeFile)
        with open(genomeFile, 'r') as f:
            first_line = f.readline()
        m = re.search('\[(.*)\]', first_line)
        if m:
            organism = m.group(1)
        else:
            organism = "COULD NOT IDENTIFY ORGANISM"

        organisms[organism] = annotateGenes.organismClass(genomeFile, HMMDir)
        curOrg = organisms[organism]
        

        # Now output the proteins:
        if (len(curOrg.families["FlgE"]) > 0 and len(curOrg.families["FliC"]) > 0):
            sys.stdout.write(organism)

            for family in familyList:
                if len(curOrg.families[family]) > 0:
                    sys.stdout.write(",")
                    sys.stdout.write(str(len(seqs[curOrg.families[family].keys()[0]])))
                else:
                    sys.stdout.write(",")
            sys.stdout.write("\n")
            organismCount += 1
            #if organismCount > 5:
            #    sys.exit(1)

if __name__ == '__main__':
    main()

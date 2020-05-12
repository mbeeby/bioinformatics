#!/usr/bin/python
#
# Takes a two-column file (ignoring empty lines)
# with desired name and PFxxxxx number and (1) works
# out current PFxxxxx version, then outputs as
# <desired name>.hmm

import sys
import os
import argparse
import annotateGenes

def main():
    parser = argparse.ArgumentParser(\
        description = "Spit out flagellar annotations for a genome")
    parser.add_argument('genomeFile', metavar = 'Genome file', type=str, \
        help = "The FASTA genome directory")
    parser.add_argument('--hmmdir', \
        help = "Location of HMMs (overrides HMM_DIR environment variable)")
    args = parser.parse_args()
    
    # -----------
    if "HMM_DIR" in os.environ:
        HMMDir = os.environ["HMM_DIR"]
    if args.hmmdir:
        HMMDir = args.hmmdir
    else:
        if not "HMM_DIR" in os.environ:
            sys.exit("ERROR: HMM_DIR not specified on command line or environment var") 
    genomeFile  = args.genomeFile
    # ------------

    organism = annotateGenes.organismClass(genomeFile, HMMDir)
   
    for family in sorted(organism.families.keys()):
        if len(organism.families[family]) > 0:           
            for protein in organism.families[family]:
                sys.stdout.write(family+"\t"+protein+"\n")
        else:
            sys.stdout.write(family+"\t-\n")

if __name__ == '__main__':
    main()

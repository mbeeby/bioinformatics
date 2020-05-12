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
sys.path.append(os.environ["HOME"]+"/projects/FlhFG/bioinformatics/bin/")
from optparse import OptionParser
from annotate import *
from workspace.bioinformatics.src import fasta

def main():
    parser = OptionParser("usage: %prog <genomes dir> <HMM dir>")
    (options,args) = parser.parse_args()
    if len(args) != 2:
        parser.error("Incorrect number of arguments!")

    genomesDir  = args[0]
    HMMDir      = args[1]
    dataFile    = "lengths.csv"

    genomeFiles = glob.glob(genomesDir+"/*.faa")

    familyList = ["TonB", "TolA", "RcsF", "IgaA", "LpoB", "ExbD_TolR", "Secretin"]

    
    sys.stdout.write("Organism")
    for family in familyList:
        sys.stdout.write(","+family)
    sys.stdout.write("\n")

    organismCount = 0
    for genomeFile in genomeFiles:

        # Scrape organism from genome file
        seqs = fasta.getProteinSeqs(genomeFile)
        with open(genomeFile, 'r') as f:
            first_line = f.readline()
        m = re.search('\[(.*)\]', first_line)
        if m:
            organism = m.group(1)
        else:
            organism = "COULD NOT IDENTIFY ORGANISM"

        families = {}

        # Get a list of families
        families["TonB"] = getTonB(genomeFile,HMMDir)
        families["TolA"] = getTolA(genomeFile,HMMDir)
        families["RcsF"] = getRcsF(genomeFile,HMMDir)
        families["IgaA"] = getIgaA(genomeFile,HMMDir)
        families["Secretin"]  = getStraightforwardsProtein(genomeFile, HMMDir, ["Secretin_N.hmm", "Secretin.hmm"])
        families["LpoB"]      = getStraightforwardsProtein(genomeFile, HMMDir, ["LpoB.hmm"])
        families["ExbD_TolR"] = getStraightforwardsProtein(genomeFile, HMMDir, ["ExbD_TolR.hmm"])

        # Now output the proteins:
        if (len(families["TonB"]) > 0 or len(families["TolA"]) > 0 or len(families["RcsF"]) > 0):
            sys.stdout.write(organism)

            for family in familyList:
                if len(families[family]) > 0:
                    sys.stdout.write(",")
                    sys.stdout.write(str(len(seqs[families[family].keys()[0]])))
                else:
                    sys.stdout.write(",")
            sys.stdout.write("\n")
            organismCount += 1
            if organismCount > 200:
                sys.exit(1)

main()

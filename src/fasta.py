#!/usr/bin/python
#
# Takes a two-column file (ignoring empty lines)
# with desired name and PFxxxxx number and (1) works
# out current PFxxxxx version, then outputs as
# <desired name>.hmm

import sys
import re

def getProteinSeqs(FAAFile):
    # Takes a FASTA file and returns a dictionary with keys as protein IDs
    # and values as protein lengths

    header_re = re.compile('>([^\s]+)')
    seq_re    = re.compile('([A-Z\-]+)')

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
            if line.rstrip() != "":
                if not n:
                    sys.exit("ERROR in getProteinSeqs: Unexpected line in "+FAAFile+"! ("+line+")\n")
                seqs[curProt] = seqs[curProt] + n.group(1)

    return seqs

def getSequenceOrder(FAAFile):
    # Takes a FASTA file and returns an ordered list of sequence names
    
    header_re = re.compile('>([^\s]+)')
    seq_re    = re.compile('([A-Z]+)')

    order = []

    f       = open(FAAFile, 'r') 
    lines   = f.readlines()
    curProt = ""
    for line in lines:
        m = header_re.match(line)
        if m:
            curProt = m.group(1)
            order.append(curProt)

    return order

def appendFASTA(FAAFilename, proteins):
    lineLen = 60
    
    f = open(FAAFilename,'a')
    for protein in proteins:
        f.write(">"+protein+"\n")
        for i in range(0, len(proteins[protein]), lineLen):
            f.write(proteins[protein][i:i+lineLen]+"\n")
    f.close() 


def writeFASTA(FAAFilename, proteins):
    # Takes a dictionary of keys as protein IDs and values as protein lengthS
    # and outputs as a FASTA file
    lineLen = 60
    
    f = open(FAAFilename,'w')
    for protein in proteins:
        f.write(">"+protein+"\n")
        for i in range(0, len(proteins[protein]), lineLen):
            f.write(proteins[protein][i:i+lineLen]+"\n")
    f.close()

def printFASTA(proteins):
    # Takes a dictionary of keys as protein IDs and values as protein lengthS
    # and prints a FASTA file to stdout
    lineLen = 60
    
    for protein in proteins:
        sys.stdout.write(">"+protein+"\n")
        for i in range(0, len(proteins[protein]), lineLen):
            sys.stdout.write(proteins[protein][i:i+lineLen]+"\n")
        
def importFSA(FSAFile):
    # Takes an FSA alignment file and returns a dictionary with keys as protein IDs
    # and values as protein sequences, including gaps

    header_re  = re.compile('\s+([0-9]+)\s+([0-9]+)')
    seq_re     = re.compile('([^\s]*)\s*(.{10}) (.{10}) (.{10}) (.{10}) (.{10}) (.{10})')

    seqs       = {}
    accessions = []

    f       = open(FSAFile, 'r')
    lines   = f.readlines()
    lineCount = 1
    m = header_re.match(lines[0])
    if not m:
        sys.exit("ERROR: Couldn not parse first line, '"+lines[0]+"'\n")
    curSeqNum  = 0
    firstRound = True
    while lineCount < len(lines):
        m = seq_re.match(lines[lineCount])
        lineCount += 1
        if m:
            if firstRound:
                curProt = m.group(1)
                accessions.append(curProt)
                seqs[curProt] = ""
            else:
                curProt = accessions[curSeqNum]
            seqs[curProt] += m.group(2) + m.group(3) + m.group(4) + m.group(5) + m.group(6) + m.group(7)
            curSeqNum += 1
        else:
            curSeqNum = 0
            firstRound = False

    return seqs
    
def splitFastaFile(FAAFilename, outputDir = '.'):
    # Takes as input a FASTA file and outputs a set of individual protein
    # FASTA files with filenames of accession numbers
    proteins = getProteinSeqs(FAAFilename)
    for protein in proteins:
        writeFASTA(outputDir+"/"+protein+".faa", {protein:proteins[protein]})
    
if __name__ == '__main__':
    sys.stdout.write("Hello, world!")
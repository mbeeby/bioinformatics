#!/usr/bin/python
'''
Gets statistics on different flagellin insert classes

'''

import argparse
import re
import sys
import random
import fasta
import numpy
     
def getClassLengths(classes,seqs,rename):
    classLens = {}
    for protein in seqs:
        accession = rename[protein]
        seq       = seqs[protein]
        
        if accession in classes:
            the_class = classes[accession]
            if not the_class in classLens:
                classLens[the_class] = []
            classLens[the_class].append(len(seq))
    return classLens
        
def getClasses(filename):
    classes = {}
    f = open(filename, 'r')
    lines = f.readlines()
    for line in lines:
        line.rstrip()
        (protein,the_class) = line.rstrip().split("\t")
        classes[protein] = the_class
    return classes
        
def loadRenameFile(filename):
    renames = {}
    f = open(filename, 'r')
    lines = f.readlines()
    for line in lines:
        line.rstrip()
        (leaf,rename) = line.rstrip().split("\t")
        renames[leaf] = rename
    return renames

def main():
    parser = argparse.ArgumentParser(description = "Get class stats")
    parser.add_argument('classFile', metavar = 'Class file', type=str, \
        help = "The file listing different classes")
    parser.add_argument('insertFile', metavar = 'Insert file', type=str, \
        help = "The FASTA file of extracted insert sequences")
    parser.add_argument('renameFile', type=str, \
        help = "File for renaming leaves")
    
    args      = parser.parse_args()
    seqs      = fasta.getProteinSeqs(args.insertFile)
    classes   = getClasses(args.classFile)
    rename    = loadRenameFile(args.renameFile) 
    classLens = getClassLengths(classes,seqs,rename)
    
    for the_class in classLens:
        print the_class+"\t"+str(len(classLens[the_class]))+"\t"+str(numpy.mean(classLens[the_class]))
    

if __name__ == '__main__':
    main()

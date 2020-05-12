#!/usr/bin/python
'''
Created on 20 Apr 2018

@author: mbeeby
'''
from optparse import OptionParser
import fasta

def main():
    parser = OptionParser("usage: %prog <FASTA file> <list of seq IDs>")
    #parser.add_option("-s", "--svg", dest="svg", help="output an SVGTools file of this name")

    (options,args) = parser.parse_args()
    if (len(args) < 2):
        parser.error("Incorrect number of arguments!")

    seqs = fasta.getProteinSeqs(args[0])
    
    count = 1
    while count < len(args):
        newSeqs = {}
        newSeqs[args[count]] = seqs[args[count]]
        fasta.printFASTA(newSeqs)
        count+= 1

main()     
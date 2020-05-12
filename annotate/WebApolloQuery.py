#!/usr/bin/python

import apollo
import sys
import re
from optparse import OptionParser
from subprocess import Popen, PIPE, STDOUT

class ApolloFASTAClass:
    def __init__(self,id, type, length, location, form, name, org):
        self.id       = id
        self.type     = type
        self.length   = length
        self.location = location
        self.form     = form
        self.name     = name
        self.org      = org
        self.seq = ""

    def appendToSeq(self,newSeq):
        self.seq = self.seq + newSeq

    def outputFASTA(self):
        sys.stdout.write("> "+self.name+" ("+self.org+")\n")
        for i in range(0,len(self.seq),60):
            sys.stdout.write(self.seq[i:i+60] +"\n")


def parseApolloFASTA(FASTAText,org):
     # Takes a FASTA file and parses
    header_re = re.compile('>([^\s]+)\s+\((.+)\)\s+([0-9]+)\s+residues\s+\[(.+)\]\s+\[(.+)\]\s+name=(.+)$')
    seq_re    = re.compile('([A-Z]+)')

    seqs = {}

    lines = FASTAText.split('\n')

    curProt = ""
    for line in lines:
        m = header_re.match(line)
        if m:
            curProt = m.group(1)
            seqs[curProt] = ApolloFASTAClass(m.group(1),m.group(2),m.group(3),m.group(4),m.group(5),m.group(6),org)
        else:
            n = seq_re.match(line)
            if (not n):
                if (line != ''):
                    sys.stderr.write("ERROR: Unexpected line! '"+line+"'\n")
                    sys.exit(1)
            else:
                seqs[curProt].appendToSeq(n.group(1))

    return seqs

def main():
    parser = OptionParser("usage: %prog <organism name> <family>")
    (options,args) = parser.parse_args()
    if len(args) != 2:
        sys.exit("ERROR: Most specific organism name and family -- hint: put organism name in quotes!")

    # ------------------
    organism  = args[0] 
    family    = args[1]

    wa   = apollo.ApolloInstance("http://vroomfondel.bc.ic.ac.uk:8080/apollo", "beebylab@gmail.com", "AutoQueryPassword")

    seqs = parseApolloFASTA(wa.io.write_text(organism, export_format='text',seq_type='peptide',export_type="FASTA"),organism)
    for seq in seqs:
        if seqs[seq].name == family:
            seqs[seq].outputFASTA()

main()

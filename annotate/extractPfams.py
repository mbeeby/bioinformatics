#!/usr/bin/python
#
# Takes a two-column file (ignoring empty lines)
# with desired name and PFxxxxx number and (1) works
# out current PFxxxxx version, then outputs as
# <desired name>.hmm

import sys
import os
import re
from optparse import OptionParser
import subprocess

def main():
    parser = OptionParser("usage: %prog <Pfam-A.hmm location> <config file> <outputDir>")
    #parser.add_option("-y", "--yAxis", action="store", dest = "yAxisClusterMethod", help="y-Axis clustering method", choices = clusteringOptions, default = 'corrCoef')
    (options,args) = parser.parse_args()
    if len(args) != 3:
        parser.error("Incorrect number of arguments!")

    Pfam_re = re.compile("(PF[0-9]+)\.([0-9]+)")

    PfamFile      = args[0]
    configFile    = args[1]
    outputDir     = args[2]

    # ---------------------
    # Load desired Pfams from config file
    PfamName      = {}

    file          = open(configFile,'r')
    lines         = file.readlines()
    for line in lines:
        line=line.strip()
        cols = line.split("\t")
        if len(cols) == 2:
            PfamName[cols[1]] = cols[0]
    file.close()

    # ---------------------
    for name in PfamName:
        commands = ["grep", name, PfamFile]
        try:
            process  = subprocess.Popen(commands, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        except:
            sys.stderr.write("Error executing "+" ".join(commands)+"\n\n")
            sys.exit(1)
        stdout, stderr = process.communicate()
        line = stdout.strip()
        m = Pfam_re.search(line)
        if not m:
            sys.stderr.write("ERROR: Could not extract Pfam from "+line+"\n")
            sys.exit(1)

        Pfam_ver = m.group(1)+"."+m.group(2)

        commands = ["hmmfetch", "-o", outputDir+"/"+PfamName[name]+".hmm", PfamFile, Pfam_ver]
        sys.stdout.write(" ".join(commands)+"\n")
        process  = subprocess.Popen(commands, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

main()

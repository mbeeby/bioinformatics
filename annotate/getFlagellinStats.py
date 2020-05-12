#!/usr/bin/python
'''

Takes a two-column file (ignoring empty lines)
with desired name and PFxxxxx number and (1) works
out current PFxxxxx version, then outputs as
<desired name>.hmm

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

rename      = {}
renameCount = -1

def plotSVG(organisms, outputFilename):
    heightPerOrg   = 10
    heightPerProt  = 20
    proteinXOffset = 150
    pixScaling     = 0.5 # Scaling factor from number of residues to pixels
    fontsize       = 12
    rectHeight     = 12
    
    outputLines = []
    
    orgY = 0
    for orgName in organisms:
        curOrg = organisms[orgName]
        if len(curOrg.families["FliC"]) > 0:
                
            #output.write("  <circle cx='150' cy='100' r='80' fill='green' />\n")
            outputLines.append("<text x='0' y='"+str(orgY+fontsize)+"' font-size='"+str(fontsize)+"' text-anchor='start' fill='black'>"+orgName+"</text>\n")
            outputLines.append("<text x='0' y='"+str(orgY+fontsize*2)+"' font-size='"+str(fontsize)+"' text-anchor='start' fill='black'>"+"%.2f"%(curOrg.number_flagellar_systems)+"</text>\n")
            # Now output the proteins:
            protX = proteinXOffset
            protY = fontsize / 2
            sorted_proteins = sorted(curOrg.families["FliC"].items(), key=operator.itemgetter(1))[::-1]
            for protein,data in sorted_proteins:
                # Print a line representing the protein
                outputLines.append('<line \
                    x1="'+str(protX)+'" y1="'+str(orgY+protY)+'" \
                    x2="'+str(protX+data.notes["length"]*pixScaling)+'" y2="'+str(orgY+protY)+'" \
                    style="stroke:rgb(128,128,128);stroke-width:2" />')
                
                outputLines.append("<text x='"+str(protX)+"' y='"+str(orgY+protY)+"' \
                            font-size='"+str(fontsize*0.7)+"' text-anchor='end' \
                            fill='black'>"+str(data.notes['length'])+" aa</text>\n")
                
                insert_label_x = (data.notes['insert_start'] + data.notes['insert_end']) * 0.5 * pixScaling
                outputLines.append("<text x='"+str(protX + insert_label_x)+"' y='"+str(orgY+protY+fontsize/2)+"' \
                            font-size='"+str(fontsize*0.7)+"' text-anchor='middle' \
                            fill='black'>"+str(data.notes['insert_length'])+" aa insert</text>\n")
                
                m = re.search('prot_([A-Z]{2}_[0-9]+\.[0-9])', protein)
                if not m:
                    m = re.search('>([A-Z]{2}_[0-9]+\.[0-9])', protein)
                    
                if m:
                    outputLines.append("<a xlink:href='https://www.ncbi.nlm.nih.gov/protein/"+m.group(1)+"' target='_blank'>")
                    outputLines.append("<text x='"+str(protX+data.notes["length"]*pixScaling)+"' y='"+str(orgY+protY)+"' \
                            font-size='"+str(fontsize*0.7)+"' text-anchor='start' \
                            fill='black'>"+m.group(1)+"</text>\n")
                    outputLines.append("</a>")
                    
                # Print box representing N-terminal domain from top domain
                N = data.notes['FliCn'][0]
                domX = protX + (N["from"]) * pixScaling
                domY = orgY + protY - (rectHeight/2)
                domWidth = (N["to"] - N["from"]) * pixScaling
                
                outputLines.append("  <rect \
                    x='"+str(domX)+"' y='"+str(domY)+"' \
                    width='"+str(domWidth)+"' height='"+str(rectHeight)+"' \
                    style=\" fill-opacity:0.4; \
                    fill:#ff0000\" />\n")
                
                # Print box representing C-terminal domain from top domain
                N = data.notes['FliCc'][0]
                domX = protX + (N["from"]) * pixScaling
                domY = orgY + protY - (rectHeight/2)
                domWidth = (N["to"] - N["from"]) * pixScaling
                
                outputLines.append("  <rect \
                    x='"+str(domX)+"' y='"+str(domY)+"' \
                    width='"+str(domWidth)+"' height='"+str(rectHeight)+"' \
                    style=\" fill-opacity:0.4; \
                    fill:#0000ff\" />\n")
                    
                protY += heightPerProt
                    
            #if len(curOrg.families["FlgL"]) > 0:
            #    sorted_proteins = sorted(curOrg.families["FlgL"].items(), key=operator.itemgetter(1))[::-1]
            #    for protein,data in sorted_proteins:
            #        sys.stdout.write("FlgL\t"+str(protein)+" ("+str(data.evidence.keys())+")\n")
            
            orgY += heightPerOrg + protY
    
    # ------------------------------
    # Output image
    
    imageHeight = orgY + heightPerOrg
    imageWidth  = 600
    output = open(outputFilename, 'w')
    output.write("<svg version='1.1'\n")
    output.write("baseProfile='full'\n")
    output.write("     width='"+str(imageWidth)+"' height='"+str(imageHeight)+"'\n")
    output.write("     xmlns='http://www.w3.org/2000/svg'")
    output.write("     xmlns:xlink='http://www.w3.org/1999/xlink'>\n")
    output.write("  <rect width='100%' height='100%' fill='#f8f8f8' />\n")
    for outputLine in outputLines:
        output.write(outputLine)
    output.write("</svg>\n")
    output.close()
    
def getRenamedProtein(protein):
    global rename
    global renameCount
    
    # If we're not renaming: just return the protein accession
    if renameCount == -1:
        return protein
    if not protein in rename:
        rename[protein] = str(renameCount)
        renameCount += 1
    return rename[protein]      

def outputRenameLookupTable(filename):
    global rename
    global renameCount
    output = open(filename,'w')
    for name in rename:
        output.write(rename[name]+"\t"+name+"\n")
    output.close()
        
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
        outputSeqs[getRenamedProtein(protein)] = \
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
        outputSeqs[getRenamedProtein(protein)] = \
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
        f.write(getRenamedProtein(protein))
        f.write("\t")
        f.write(repr(curOrg.organism))
        f.write("\t")
        f.write(repr(curOrg.taxID))
        f.write("\n")
    f.close()
            
def getFAAfromGenomeFile(genomeFile, outputFasta):
    # Outputs all protein translations from a gbff file as a FASTA file.
    # Also extract TAXID and ORGANISM metadata.
    taxID_re = re.compile('taxon:\s*([0-9]+)')
    taxID    = None
    organism = None
    output = open(outputFasta,'w')
    try:
        contig = SeqIO.parse(open(genomeFile,"r"), "genbank")
    except:
        sys.stderr.write("ERROR: Parsing gbk file "+genomeFile+"!\n")
        sys.exit(1)
    for record in contig:
        for feature in record.features:
            if feature.type == 'source':
                if "db_xref" in feature.qualifiers:
                    for xref in feature.qualifiers["db_xref"]:
                        m = taxID_re.search(xref)
                        if m:
                            if (taxID) and (taxID != int(m.group(1))):
                                sys.exit("ERROR: Multiple taxIDs mentioned in "+genomeFile+"\n")
                            taxID = int(m.group(1))
                if "organism" in feature.qualifiers:
                    if (organism) and (organism != feature.qualifiers["organism"][0]):
                        sys.exit("ERROR: Multiple organisms mentioned in "+genomeFile+"\n")
                    organism = feature.qualifiers["organism"][0]
                else:
                    sys.stderr.write("ERROR: Couldn't find tax_id for record!\n")
                    sys.exit(1)
                # Then break -- just take the first organism, as subsequent
                # organisms are likely phage
                break
                    
  
        for feature in record.features:
            if feature.type == 'CDS':
                if "protein_id" in feature.qualifiers:
                    protein_id  = feature.qualifiers["protein_id"][0]
                    translation = feature.qualifiers["translation"][0]
                    output.write(">"+protein_id+"\n")
                    for i in range(0, len(translation), 60):
                        output.write(translation[i:i+60]+"\n") 
    output.close()
                    
    return (taxID,organism)

def main():
    global rename
    global renameCount
    parser = argparse.ArgumentParser(description = "Get statistics on flagellins")
    parser.add_argument('genomesDir', metavar = '.gbff genome dir', type=str, \
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
    parser.add_argument('--rename', \
        help = "Rename sequences with sequential integers, and output a lookup table specified as filename here, together with organisms")
    parser.add_argument('--orgs', help='Output a lookup table of proteins, renames, taxID, and organism')
    args   = parser.parse_args()
     
    #-----------
    if "HMM_DIR" in os.environ:
        HMMDir = os.environ["HMM_DIR"]
    if args.hmmdir:
        HMMDir = args.hmmdir
    else:
        if not "HMM_DIR" in os.environ:
            sys.exit("ERROR: HMM_DIR not specified on command line or environment var") 
    # -----------------
    if args.rename:
        rename = {}
        renameFile = args.rename
        renameCount = args.renameStart
    # ----------------
    if args.orgs: 
        outputFile = open(args.orgs, 'w')
        outputFile.close()
        
    tempFAA = ".temp"+str(random.random()*1000000)+".faa"

    organisms     = {}
    insertSeqs    = {}
    endSeqs       = {}
    genomeFiles   = glob.glob(args.genomesDir+"/*.gbff")
    organismCount = 0
    
    # Setup the output file to be empty:
    outputFile = open(args.insertFASTA, 'w')
    outputFile.close()
    outputFile = open(args.conservedFASTA, 'w')
    outputFile.close()

    for genomeFile in genomeFiles:
        (taxID, organism) = getFAAfromGenomeFile(genomeFile, tempFAA)
  
        curOrg   = organismClass(tempFAA, HMMDir, annotateAll=False, organism=organism,taxID=taxID)
        organisms[curOrg.organism] = curOrg 
        curOrg.annotateFlagellins()
        curOrg.annotateCoreFlagellarGenes()
        
        if curOrg.number_flagellar_systems > 0:
            sys.stdout.write(curOrg.organism + "\t" + str(curOrg.number_flagellar_systems))
            sys.stdout.write("\n")
            # Now output the proteins:
            if len(curOrg.families["FliC"]) > 0:
                sorted_proteins = sorted(curOrg.families["FliC"].items(), key=operator.itemgetter(1))[::-1]
                for protein,data in sorted_proteins:
                    sys.stdout.write("FliC\t"+str(protein)+" ("+repr(data.notes["insert_length"])+"; "+str(data.evidence.keys())+")\n")
                    
            if len(curOrg.families["FlgL"]) > 0:
                sorted_proteins = sorted(curOrg.families["FlgL"].items(), key=operator.itemgetter(1))[::-1]
                for protein,data in sorted_proteins:
                    sys.stdout.write("FlgL\t"+str(protein)+" ("+str(data.evidence.keys())+")\n")
                    
            else:
                sys.stdout.write("\t0\t-")
            sys.stdout.write("\n")
            
            organismCount += 1

        if args.insertFASTA:
            extractInserts(insertSeqs, curOrg, tempFAA)
        if args.conservedFASTA:
            extractConservedEnds(endSeqs, curOrg, tempFAA)
        if args.orgs: 
            appendLookupTable(curOrg, args.orgs)
    if args.insertFASTA:
        fasta.appendFASTA(args.insertFASTA, insertSeqs)
    if args.conservedFASTA:
        fasta.appendFASTA(args.conservedFASTA, endSeqs)   
    

    if renameCount != -1:
        outputRenameLookupTable(renameFile)

    plotSVG(organisms,args.SVGout)
    os.system("rm "+tempFAA)
    print "DONE"

if __name__ == '__main__':
    main()

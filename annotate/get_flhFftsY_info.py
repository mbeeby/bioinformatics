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
import glob
import operator
import re
import argparse
import fasta
from annotateGenes import *

def plotDomain(domain, colour, protX, protY, orgX, orgY, notes, label, SVGparams):
    outputLines = []
    N    = notes[domain][0]
    domX = orgX+protX + (N["from"]) * SVGparams.pixScaling
    domY = orgY + protY - (SVGparams.rectHeight/2)
    domWidth = (N["to"] - N["from"]) * SVGparams.pixScaling
    
    outputLines.append("  <rect \
        x='"+str(domX)+"' y='"+str(domY)+"' \
        width='"+str(domWidth)+"' height='"+str(SVGparams.rectHeight)+"' \
        style=\" fill-opacity:1; \
        fill:"+colour+"\" />\n")
    outputLines.append("<text x='"+str(domX+(domWidth/2))+"' y='"+str(orgY+protY+SVGparams.rectHeight*0.25)+"' \
                font-size='"+str(SVGparams.smallFontSize())+"' text-anchor='middle' \
                fill='black'>"+label+"</text>\n")
    return outputLines

def plotProtein(protein, protX, protY, orgX, orgY, notes, SVGparams, label=False, plotLength = False):
    outputLines = []
    outputLines.append('<line \
        x1="'+str(orgX+protX)+'" y1="'+str(orgY+protY)+'" \
        x2="'+str(orgX+protX+notes["length"]*SVGparams.pixScaling)+'" y2="'+str(orgY+protY)+'" \
        style="stroke:rgb(128,128,128);stroke-width:2" />')
        
    if plotLength:
        outputLines.append("<text x='"+str(orgX+protX)+"' y='"+str(orgY+protY)+"' \
            font-size='"+str(SVGparams.smallFontSize())+"' text-anchor='end' \
            fill='black'>"+str(notes['length'])+" aa</text>\n")
                
    #insert_label_x = (data.notes['insert_start'] + data.notes['insert_end']) * 0.5 * pixScaling
    #outputLines.append("<text x='"+str(protX + insert_label_x)+"' y='"+str(orgY+protY+fontsize/2)+"' \
    #            font-size='"+str(SVGparams.smallFontSize())+"' text-anchor='middle' \
    #            fill='black'>"+str(data.notes['insert_length'])+" aa insert</text>\n")
    
    if label:
        outputLines.append("<text x='"+str(orgX+protX+notes["length"]*SVGparams.pixScaling)+"' y='"+str(orgY+protY)+"' \
                            font-size='"+str(SVGparams.smallFontSize())+"' text-anchor='start' \
                            fill='black'>"+label+"</text>\n")
    return outputLines

def getLink(protein):
    link = ""        
    m = re.search('prot_([A-Z]{2}_[0-9]+\.[0-9])', protein)
    if not m:
        m = re.search('([A-Z]{2}_[0-9]+\.[0-9])', protein)
    if m:
        link = "https://www.ncbi.nlm.nih.gov/protein/"+m.group(1)
    return link
    
class SVGparamsClass:
    
    def __init__(self):
        self.heightPerOrg     = 30
        self.heightPerProt    = 20
        self.proteinXOffset   = 250
        self.pixScaling       = 0.33 # Scaling factor from number of residues to pixels
        self.fontsize         = 12
        self.rectHeight       = 12
        self.smallFontScaling = 0.7
        self.imageWidth       = 1000
        self.imageHeight      = 900
        self.border           = 30
        
        
    def smallFontSize(self):
        return self.fontsize * self.smallFontScaling
        

def plotSVG(organisms, outputFilename):
    SVGparams = SVGparamsClass()
    
    outputLines = []
    
    orgY = SVGparams.border
    outputLines.append("<text x='"+str(SVGparams.proteinXOffset)+"'     y='"+str(SVGparams.fontsize*1.5)+"' font-size='"+str(SVGparams.fontsize)+"' fill='black'>FlhF</text>\n")
    outputLines.append("<text x='"+str(SVGparams.proteinXOffset+250)+"' y='"+str(SVGparams.fontsize*1.5)+"' font-size='"+str(SVGparams.fontsize)+"' fill='black'>FtsY</text>\n")
    outputLines.append("<text x='"+str(SVGparams.proteinXOffset+500)+"' y='"+str(SVGparams.fontsize*1.5)+"' font-size='"+str(SVGparams.fontsize)+"' fill='black'>SRP</text>\n")
    
    for orgName in organisms:
        curOrg = organisms[orgName]
        if curOrg.number_flagellar_systems > 0 and \
           ((len(curOrg.families["FlhF"]) > 0) or \
            (len(curOrg.families["FtsY"]) > 0) or \
            (len(curOrg.families["Ffh"]) > 0)):
                
            #output.write("  <circle cx='150' cy='100' r='80' fill='green' />\n")
            outputLines.append("<text x='0' y='"+str(orgY+SVGparams.fontsize)+"' font-size='"+str(SVGparams.fontsize)+"' text-anchor='start' fill='black'>"+orgName+"</text>\n")
            #outputLines.append("<text x='0' y='"+str(orgY+SVGparams.fontsize*2)+"' font-size='"+str(SVGparams.fontsize)+"' text-anchor='start' fill='black'>"+"%.2f"%(curOrg.number_flagellar_systems)+"</text>\n")
            
            
            # Now output the proteins:
            protX = SVGparams.proteinXOffset
                        
            orgX     = 0
            protY    = SVGparams.fontsize / 2
            maxProtY = 0

            sorted_proteins = sorted(curOrg.families["FlhF"].items(), key=operator.itemgetter(1))[::-1]
            for protein,data in sorted_proteins:
 
                link = getLink(protein)
                
                outputLines.append("<a xlink:href='"+link+"' target='_blank'>")                           
                outputLines = outputLines + plotProtein(protein, protX, protY, orgX, orgY, data.notes, SVGparams)
                if "Hits SRP54n" in data.evidence:
                    outputLines = outputLines + plotDomain("SRP54n", "#00aa00", protX, protY, orgX, orgY, data.notes, 'S54n', SVGparams)
                if "Hits SRP54" in data.evidence:
                    outputLines = outputLines + plotDomain("SRP54", "#88aa88", protX, protY, orgX, orgY, data.notes, 'S54', SVGparams)
                if "Hits SPB" in data.evidence:
                    outputLines = outputLines + plotDomain("SPB", "#880000", protX, protY, orgX, orgY, data.notes, 'SPB', SVGparams)
                outputLines.append("</a>")
                
                protY += SVGparams.heightPerProt
            if protY > maxProtY:
                maxProtY = protY
            
            
            orgX = 250
            protY = SVGparams.fontsize / 2

            sorted_proteins = sorted(curOrg.families["FtsY"].items(), key=operator.itemgetter(1))[::-1]
            for protein,data in sorted_proteins: 
                #print "FtsY evidence:"
                #print data.evidence    
                link = getLink(protein)
                outputLines.append("<a xlink:href='"+link+"' target='_blank'>")
                outputLines.append("<text x='"+str(orgX+protX)+"' y='"+str(orgY+protY+SVGparams.smallFontSize()*0.5)+"' font-size='"+str(SVGparams.fontsize*SVGparams.smallFontScaling)+"' text-anchor='end' fill='black'>"+str(data.notes["length"])+" aa</text>\n")                           
         
                outputLines = outputLines + plotProtein(protein, protX, protY, orgX, orgY, data.notes, SVGparams)
                if "Hits SRP54n" in data.evidence:
                    outputLines = outputLines + plotDomain("SRP54n", "#00aa00", protX, protY, orgX, orgY, data.notes, 'S54n', SVGparams)
                if "Hits SRP54" in data.evidence:
                    outputLines = outputLines + plotDomain("SRP54", "#88aa88", protX, protY, orgX, orgY, data.notes, 'S54', SVGparams)
                if "Hits SPB" in data.evidence:
                    outputLines = outputLines + plotDomain("SPB", "#880000", protX, protY, orgX, orgY, data.notes, 'SPB', SVGparams)
                if 'A-domain length' in data.notes:
                    outputLines.append("<text x='"+str(orgX+protX+data.notes["A-domain length"]*SVGparams.pixScaling*0.5)+"' y='"+str(orgY+protY-SVGparams.smallFontSize()*0.5)+"' font-size='"+str(SVGparams.fontsize*SVGparams.smallFontScaling)+"' text-anchor='middle' fill='black'>"+str(data.notes["A-domain length"])+" aa</text>\n")

                outputLines.append("</a>")
                protY += SVGparams.heightPerProt
            if protY > maxProtY:
                maxProtY = protY
                
            orgX = 500
            protY = SVGparams.fontsize / 2            
            sorted_proteins = sorted(curOrg.families["Ffh"].items(), key=operator.itemgetter(1))[::-1]
            for protein,data in sorted_proteins:   
                link = getLink(protein)                
                outputLines.append("<a xlink:href='"+link+"' target='_blank'>")                           
           
                outputLines = outputLines + plotProtein(protein, protX, protY, orgX, orgY, data.notes, SVGparams)
                if "Hits SRP54n" in data.evidence:
                    outputLines = outputLines + plotDomain("SRP54n", "#00aa00", protX, protY, orgX, orgY, data.notes, 'S54n', SVGparams)
                if "Hits SRP54" in data.evidence:
                    outputLines = outputLines + plotDomain("SRP54", "#88aa88", protX, protY, orgX, orgY, data.notes, 'S54', SVGparams)
                if "Hits SPB" in data.evidence:
                    outputLines = outputLines + plotDomain("SPB", "#880000", protX, protY, orgX, orgY, data.notes, 'SPB', SVGparams)
                outputLines.append("</a>")
                protY += SVGparams.heightPerProt
            if protY > maxProtY:
                maxProtY = protY                                
            
            orgY += SVGparams.heightPerOrg + maxProtY 
    
    # ------------------------------
    # Output image
    
    SVGparams.imageHeight = orgY + SVGparams.heightPerOrg
    output = open(outputFilename, 'w')
    output.write("<svg version='1.1'\n")
    output.write("baseProfile='full'\n")
    output.write("     width='"+str(SVGparams.imageWidth)+"' height='"+str(SVGparams.imageHeight)+"'\n")
    output.write("     xmlns='http://www.w3.org/2000/svg'")
    output.write("     xmlns:xlink='http://www.w3.org/1999/xlink'>\n")
    output.write("  <rect width='100%' height='100%' fill='#f8f8f8' />\n")
    for outputLine in outputLines:
        output.write(outputLine)
    output.write("</svg>\n")
    output.close()




def main():
    parser = argparse.ArgumentParser(description = "Show FlhF and FtsY information graphically")
    parser.add_argument('genomesDir', metavar = 'Genome dir', type=str, \
        help = "The FASTA genome directory")
    parser.add_argument('SVGout', metavar = 'SVGTools output filename', type=str, \
        help = "The SVGTools output filename")
    parser.add_argument('--maxOrgs', metavar = 'Max orgs plotted', type=int, \
        help = "The maximum number of organisms plotted")
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
    # ----------------- 

    organisms     = {}
    genomeFiles   = glob.glob(args.genomesDir+"/*.faa")
    organismCount = 0
    
    for genomeFile in genomeFiles:
        # Scrape organism from genome file
        with open(genomeFile, 'r') as f:
            first_line = f.readline()
        m = re.search('\[([^\]]+)\]', first_line)
        if m:
            organism = m.group(1)
        else:
            organism = "COULD NOT IDENTIFY ORGANISM"
        
        organisms[organism] = organismClass(genomeFile, HMMDir, annotateAll=False)
        curOrg = organisms[organism]
        curOrg.annotateCoreFlagellarGenes()
        
        if curOrg.number_flagellar_systems > 0:
            curOrg.annotateSRP()
            sys.stdout.write(organism + "\t" + str(curOrg.number_flagellar_systems))
            sys.stdout.write("\n")
            # Now output the proteins:
            organismCount += 1
            if args.maxOrgs:
                if organismCount > args.maxOrgs:
                    sys.stderr.write("Prematurely killed")
                    break

    plotSVG(organisms,args.SVGout)
    print "DONE"

if __name__ == '__main__':
    main()

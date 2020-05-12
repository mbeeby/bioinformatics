#!/usr/bin/python
'''

Created 20th April 2018

Plots HMMs against a provided set of FASTA files for giant flagellins project

@author: Morgan Beeby 

'''

import os
import operator
import re
import sys
import argparse
import fasta
import SVGTools
import annotateGenes

def getLink(protein):
    link = ""        
    m = re.search('prot_([A-Z]{2}_[0-9]+\.[0-9])', protein)
    if not m:
        m = re.search('([A-Z]{2}_[0-9]+\.[0-9])', protein)
    if m:
        link = "https://www.ncbi.nlm.nih.gov/protein/"+m.group(1)
    return link

def recalcSlivers(slivers,From,To):
    # FUNCTION: To split a contiguous gene into slivers according to each
    # new domain annotated to it.
    output = []
    for sliver in slivers:
        if (From <= sliver[1]) and (To >= sliver[0]):
            # Confirmed there is some overlap
            if (From > sliver[0]):
                output.append([sliver[0],From-1])
            if (To < sliver[1]):
                output.append([To+1,sliver[1]])
        else:
            output.append(sliver)
            
    return output
        
    

def main():
    parser = argparse.ArgumentParser(description = "Plots an SVG of domain structures HMMs specified in a directory")
    parser.add_argument('fastaFile', metavar = 'FASTA file', type=str, \
        help = "The FASTA file to plot")
    parser.add_argument('SVGout', metavar = 'SVG output filename', type=str, \
        help = "The SVGTools output filename")
    parser.add_argument('xdomain', \
        help = "Filename of x-domain HMM (not directory -- this is HMM_DIR)")
    parser.add_argument('--hmmdir', \
        help = "Location of HMMs (overrides HMM_DIR environment variable)")
    parser.add_argument('--extractSlivers', help="label, and extract to FAA, sliver sections > X aa")
    parser.add_argument('--outputOthers', help="Also output stretches of sequence between annotated domains")
    args   = parser.parse_args()
    
    # ----------
    # Init
    
    yInc  = 40
    xScale = 0.1
    yScale = 0.5
    geneHeight = yInc * 0.4
    sliverCount = 0
    xMargin = 500
    sliverProteins = {}
    coords         = {}
     
    #-----------
    if "HMM_DIR" in os.environ:
        HMMDir = os.environ["HMM_DIR"]
    if args.hmmdir:
        HMMDir = args.hmmdir
    else:
        if not "HMM_DIR" in os.environ:
            sys.exit("ERROR: HMM_DIR not specified on command line or environment var")
    proteins    = fasta.getProteinSeqs(args.fastaFile) 
    proteinList = fasta.getSequenceOrder(args.fastaFile)

    annotations = annotateGenes.organismClass(args.fastaFile, HMMDir, annotateAll=False)
    
    annotations.getStraightforwardsProtein("X", [args.xdomain], HMMDir=HMMDir)
    annotations.annotateFlagellins()
         
    # ----------------- 
    SVGFile = SVGTools.SVGClass(imageHeight=1500)
     
    y = yInc

    for protein in proteinList:
        slivers = [[1,len(proteins[protein])]]
        Y = (y ) * yScale
        SVGFile.addItem(SVGTools.lineClass(x1=xMargin*xScale,y1=Y,x2=(xMargin+len(proteins[protein]))*xScale,y2=Y))
        SVGFile.addItem(SVGTools.textClass(protein,x=0,y=(geneHeight+y)*yScale),layer=1)
        y += yInc
        if protein in annotations.families["FliC"]:
            annotation = annotations.families["FliC"][protein]
            for FliCn in annotation.notes["FliCn"]:
                slivers = recalcSlivers(slivers,FliCn["from"],FliCn["to"])
                x1 = (xMargin+FliCn["from"])*xScale
                x2 = (xMargin+FliCn["to"])*xScale
                SVGFile.addItem(SVGTools.rectClass(x=x1, y=Y- (geneHeight * 0.5),width=x2-x1,height=geneHeight,colour=SVGTools.colourClass(r=255)))
            for FliCc in annotation.notes["FliCc"]:
                slivers = recalcSlivers(slivers,FliCc["from"],FliCc["to"])
                x1 = (xMargin+FliCc["from"])*xScale
                x2 = (xMargin+FliCc["to"])*xScale
                SVGFile.addItem(SVGTools.rectClass(x=x1, y=Y- (geneHeight * 0.5),width=x2-x1,height=geneHeight,colour=SVGTools.colourClass(b=255)))
        if protein in annotations.families["X"]:
            annotation = annotations.families["X"][protein]
            green = 255
            for domainHit in annotation.notes["notes"]:
                greenColour = SVGTools.colourClass(g=green)
                
                x1 = (xMargin+domainHit["from"])*xScale
                x2 = (xMargin+domainHit["to"])*xScale
                
                if "iEvalue" in domainHit:
                    if (domainHit["iEvalue"] < 0.001):
                        slivers = recalcSlivers(slivers,domainHit["from"],domainHit["to"])
                        SVGFile.addItem(SVGTools.rectClass(x=x1, y=Y- (geneHeight * 0.5),width=x2-x1,height=geneHeight,colour=greenColour))
                        #SVGFile.addItem(SVGTools.textClass(str(domainHit["iEvalue"]),fontsize=8, x=((x1+x2)/2),y=(geneHeight+y-yInc)*yScale),layer=1)
        
                #green = green * 0.75
        if args.extractSlivers:
            for sliver in slivers:
                if sliver[1]-sliver[0] >= int(args.extractSlivers):
                    _x = (xMargin+(sliver[0]+sliver[1])/2) * xScale
                    _y = (geneHeight+y-yInc)*yScale
                    coords[str(sliverCount)] = [_x,_y]
                    SVGFile.addItem(SVGTools.textClass(str(str(sliverCount)),fontsize=8, x=_x,y=_y),layer=1)
                    
                    sliverProteins[str(sliverCount)] = proteins[protein][sliver[0]:sliver[1]]
                    sliverCount+=1

    f = open("slivers.xy", 'w')
    for sliver in sliverProteins:
        f.write(sliver+"\t"+str(int(coords[sliver][0]))+","+str(int(coords[sliver][1]))+"\n")
    fasta.writeFASTA("slivers.faa", sliverProteins)
    SVGout = open(args.SVGout, 'w')
    SVGout.write(str(SVGFile)+"\n")

if __name__ == '__main__':
    main()
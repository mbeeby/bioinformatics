#!/usr/bin/python
'''
Takes an SVGTools tree output from seaview and 
adds circles at the tips designating metadata

'''

import argparse
import re
import sys
import random
import fasta
import numpy

class SVGlineClass:
    def __init__(self,m):
        self.x1 = float(m.group(1))        
        self.y1 = float(m.group(2))
        self.x2 = float(m.group(3))
        self.y2 = float(m.group(4))
        self.colour= m.group(5)
        self.strokewidth = float(m.group(6))
        #self.strokelinecap = m.group(7)
        self.strokelinecap = "round"
        
    def __str__(self):
        return '<line x1="'+str(self.x1)+'" y1="'+str(self.y1)+'" x2="'+str(self.x2)+'" y2="'+str(self.y2)+'" style="stroke:rgb('+self.colour+');stroke-width:'+ str(self.strokewidth)+';stroke-linecap:'+self.strokelinecap+'" />'
    
class SVGleafClass:
    def __init__(self,m):
        self.translatex = float(m.group(1))
        self.translatey = float(m.group(2))
        self.rotate     = float(m.group(3))
        self.x          = float(m.group(4))
        self.y          = float(m.group(5))
        self.fontfamily = m.group(6)
        self.fontsize   = float(m.group(7))
        self.fill       = m.group(8)
        self.text       = m.group(9)
        self.id         = None
        
    def __str__(self):
        output = ''
        if self.text:
            output  = '<g '
            output += 'id="'+self.id+'_text" '
            output += 'transform="translate(0,0) '
            #output += 'transform="translate('+str(self.translatex)+','+str(self.translatey)+') '
            #output += 'rotate('+str(self.rotate)+')'
            output += '">'
            output += '<text x="'+str(self.x)+'" y="'+str(self.y)+'" font-family="'+self.fontfamily+'" font-size="'+str(self.fontsize)+'" xml:space="preserve" fill="rgb('+self.fill+')">'+self.text+'</text></g>'
        return output
        
def getColourDependingOnLength(aa_len):
    if aa_len < 400:
        colour  = "#0000ff"
    elif aa_len < 600:
        colour  = "#00ff00"
    elif aa_len < 800:
        colour = "#ffff00"
    elif aa_len < 900:
        colour = "#ff9900"
    else:
        colour = "#ff0000"
    return colour

class labelClass:
    def __init__(self,x,y,id,leaf):
        self.x = x
        self.y = y
        self.leafX = x
        self.leafY = y
        self.start  = id+"_circle"
        self.end    = id+"_text"
        self.leaf   = leaf
        
    def radiateFromCentre(self,x0,y0,offset):
        h = (self.x-x0)
        v = (self.y-y0)
        mag = numpy.sqrt((h*h)+(v*v))
        
        h = h / mag
        v = v / mag
        
        self.x = self.x + offset * h
        self.y = self.y + offset * v
        self.leaf.x = self.x
        self.leaf.y = self.y
        
    def __str__(self):
        connector = "<path "
        connector += 'style="fill:none;fill-rule:evenodd;stroke:#c4c4c4;stroke-width:0.5;stroke-linecap:butt;stroke-linejoin:miter;strike-miterlimit:4;stroke-dasharray:1,1;stroke-dashoffset:0;stroke-opacity:1" '
        connector += 'd="m '+str(self.leafX)+','+str(self.leafY)+' '+str(self.x)+','+str(self.y)+'" '
        connector += 'inkscape:connector-type="polyline" '
        connector += 'inkscape:connector-curvature="0" '
        connector += 'inkscape:connection-start="#'+self.start+'" '
        connector += 'inkscape:connection-end="#'+self.end+'" />\n'
        connector += str(self.leaf)
        return connector

def processSVG(SVGFile, insertSeqs,sizeScaling, circleScaling,fontScaling, rename={}, alwaysLabel={}, labelThresh=1,labelOffset=50):
    line_re = re.compile('<line x1="([0-9\.\-]+)" y1="([0-9\.\-]+)" x2="([0-9\.\-]+)" y2="([0-9\.\-]+)" style="stroke:rgb\((.+)\);stroke-width:([0-9+]);stroke-linecap:([a-zA-Z]+)" />')
    leaf_re = re.compile('<g transform="translate\(([0-9\.\-]+),([0-9\.\-]+)\) rotate\(([0-9\.\-]+)\)"><text x="([0-9\.\-]+)" y="([0-9\.\-]+)" font-family="([A-Za-z]+)" font-size="([0-9\.\-]+)" xml:space="preserve"\s+fill="rgb\((.+)\)">\s*([^\s]+)\s*</text>')
    
    SVG       = open(SVGFile,'r')
    lines     = SVG.readlines()
    labels    = []
    Xs        = []
    Ys        = []

    lineCount = 5
    
    sys.stdout.write('<?xml version="1.0" encoding="UTF-8" standalone="no"?>\n')
    sys.stdout.write('<svg\n')
    sys.stdout.write('   xmlns:dc="http://purl.org/dc/elements/1.1/"\n')
    sys.stdout.write('   xmlns:cc="http://creativecommons.org/ns#"\n')
    sys.stdout.write('   xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#"\n')
    sys.stdout.write('   xmlns:svg="http://www.w3.org/2000/svg"\n')
    sys.stdout.write('   xmlns="http://www.w3.org/2000/svg"\n')
    sys.stdout.write('   xmlns:sodipodi="http://sodipodi.sourceforge.net/DTD/sodipodi-0.dtd"\n')
    sys.stdout.write('   xmlns:inkscape="http://www.inkscape.org/namespaces/inkscape"\n')
    sys.stdout.write('   width="595px"\n')
    sys.stdout.write('   height="842px"\n')
    sys.stdout.write('   viewBox="0 0 595 842"\n')
    sys.stdout.write('   version="1.1"\n')
    sys.stdout.write('   id="svg2"\n')
    sys.stdout.write('   inkscape:version="0.91 r"\n')
    sys.stdout.write('   sodipodi:docname="ends_inkscape.svg">\n')
    
    while lineCount < len(lines)-1:
        line = lines[lineCount]
        lineCount += 1
        # ------------------
        m = line_re.search(line)
        if m:
            curLine = SVGlineClass(m)
            curLine.strokewidth = curLine.strokewidth * sizeScaling
            line = str(curLine)
            prevLine = curLine
            
        m = leaf_re.search(line)
        if m:
            lineCount+= 1
            curLeaf = SVGleafClass(m)
            curLeaf.fontsize = curLeaf.fontsize * sizeScaling * fontScaling
            curLeaf.translatex = prevLine.x2
            curLeaf.translatey = prevLine.y2
            aa_len = len(insertSeqs[curLeaf.text])
            size   = aa_len * sizeScaling * circleScaling
            if curLeaf.text in rename:
                curLeaf.id   = curLeaf.text
                curLeaf.text = rename[curLeaf.text]

            if (aa_len < labelThresh) and (not (curLeaf.id in alwaysLabel)):
                curLeaf.text = None
                
            colour = getColourDependingOnLength(aa_len)
            opacity = 0.5 
            sys.stdout.write('<circle ')
            if curLeaf.text:
                sys.stdout.write('id="'+curLeaf.id+'_circle" ')
                labels.append(labelClass(prevLine.x2,prevLine.y2,curLeaf.id,curLeaf)) 
            sys.stdout.write('cx="'+str(prevLine.x2)+'" cy="'+str(prevLine.y2)+'" r="'+str(size)+'"  opacity="'+str(opacity)+'" fill="'+colour+'" />\n')
            Xs.append(prevLine.x2)
            Ys.append(prevLine.y2)

        else:
            sys.stdout.write(line+"\n")
    xCentre = numpy.mean(Xs)
    yCentre = numpy.mean(Ys)
    for label in labels:
        label.radiateFromCentre(xCentre,yCentre,labelOffset)
        sys.stdout.write(str(label))
    sys.stdout.write("</svg>")
        
def loadRenameFile(filename):
    '''
    Rename file: four columns:
    1. Gene name
    2. Gene rename (i.e., consecutive numbers)
    3. organism name
    4. organism taxID
    '''
    renames = {}
    orgs    = {}
    f = open(filename, 'r')
    lines = f.readlines()
    for line in lines:
        line.rstrip()
        (accession,leafName, org, taxID) = line.rstrip().split("\t")
        # Strip out quotes
        org = org[1:-1]
        renames[leafName] = accession
        orgs[leafName]    = org
    return (renames, orgs)

def loadAlwaysLabelFile(filename):
    '''
    Always label file: one column:
    1. Gene name
    '''
    
    accession = {}
    f = open(filename, 'r')
    lines = f.readlines()
    for line in lines:
        alwaysLabelID = re.split("\s*#",line.rstrip())[0]
        accession[alwaysLabelID] = ""
        sys.stderr.write(alwaysLabelID+"\n")
        
    return accession

def main():
    parser = argparse.ArgumentParser(description = "Annotate SVGTools tree")
    parser.add_argument('SVGFile', metavar = 'SVGTools file', type=str, \
        help = "The input SVGTools tree")
    parser.add_argument('insertFile', metavar = 'Insert file', type=str, \
        help = "The FASTA file of extracted insert sequences")
    parser.add_argument('--sizeScaling', nargs='?', const=1, type=float, \
        help = "Global size scaling factor ", default=0.1)
    parser.add_argument('--circleScaling', nargs='?', const=1, type=float, \
        help = "Circle scaling factor ", default=0.05)
    parser.add_argument('--fontScaling', nargs='?', const=1, type=float, \
        help = "Extra font scaling factor ", default=1)
    parser.add_argument('--renameFile', type=str, \
        help = "File for renaming leaves")
    parser.add_argument('--alwaysLabel', type=str, \
        help = "Always label leaves with this accession")
    parser.add_argument('--labelThresh', default=1, type=int, \
        help = "Minimum size of insert to label")
    
    
    args   = parser.parse_args()
    seqs = fasta.getProteinSeqs(args.insertFile)
    rename      = {}
    alwaysLabel = {}
    if args.renameFile:
        (rename,orgs) = loadRenameFile(args.renameFile)
    if args.alwaysLabel:
        alwaysLabel = loadAlwaysLabelFile(args.alwaysLabel)
    
    processSVG(args.SVGFile, seqs, args.sizeScaling, args.circleScaling, args.fontScaling, rename=orgs, alwaysLabel=alwaysLabel, labelThresh=args.labelThresh)

if __name__ == '__main__':
    main()

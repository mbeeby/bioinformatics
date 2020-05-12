#!/usr/bin/python
'''
Created on 20 Apr 2018

@author: mbeeby
'''

import sys

class colourClass:
    
    def __init__(self,r=0,g=0,b=0):
        self.r = r
        self.g = g
        self.b = b
        
    def __str__(self):
        return "rgb("+str(self.r)+","+str(self.g)+","+str(self.b)+")"
    
class lineClass:
    def __init__(self,x1=0,y1=0,x2=0,y2=0,colour=colourClass(),strokewidth=1,strokelinecap="round"):
        self.x1 = x1        
        self.y1 = y1
        self.x2 = x2
        self.y2 = y2
        self.colour= colour
        self.strokewidth = strokewidth
        self.strokelinecap = strokelinecap
        
    def __str__(self):
        return '<line x1="'+str(self.x1)+'" y1="'+str(self.y1)+'" x2="'+str(self.x2)+ \
               '" y2="'+str(self.y2)+'" style="stroke:'+str(self.colour)+ \
               ';stroke-width:'+ str(self.strokewidth)+';stroke-linecap:'+self.strokelinecap+'" />'
               
class rectClass:
    def __init__(self, x=0,y=0, width=0, height=0, opacity=1, colour=colourClass()):
        self.x = x
        self.y = y
        self.width = width
        self.height = height
        self.opacity = opacity
        self.colour = colour
        
    def __str__(self):
        return "<rect x='"+str(self.x)+"' y='"+str(self.y)+"' \
               width='"+str(self.width)+"' height='"+str(self.height)+"' \
               style='fill-opacity:"+str(self.opacity)+"; \
               fill:"+str(self.colour)+"' />\n"
                   
class textClass:
    
    def __init__(self, text, x=0,y=0, fontfamily="Verdana", fontsize="12", fill=colourClass()):
        self.text = text
        self.x = x
        self.y = y
        self.fontsize = fontsize
        self.fontfamily = fontfamily
        self.fill = fill
        
    def __str__(self):
        return "<text x='"+str(self.x)+"' y='"+str(self.y)+"' \
               font-family='"+str(self.fontfamily)+"' font-size='"+str(self.fontsize)+"' \
               fill='"+str(self.fill)+"'>"+self.text+"</text>\n"
 
class SVGClass:
    
    def __init__(self,imageWidth=1000,imageHeight=1000):
        self.imageWidth  = imageWidth
        self.imageHeight = imageHeight
        self.layer       = []    # An array of array of strings representing layers        
        
    def __str__(self):
        output = ""
        output += "<svg version='1.1' \n"
        output += "baseProfile='full' \n"
        output += "     width='"+str(self.imageWidth)+"' height='"+str(self.imageHeight)+"'\n"
        output += "     xmlns='http://www.w3.org/2000/svg'"
        output += "     xmlns:xlink='http://www.w3.org/1999/xlink'>\n"
        
        for layer in self.layer:
            for item in layer:
                output += str(item)
        
        output += "</svg>\n"
        return output
    
    def addItem(self, item, layer = 0):
        while len(self.layer) <= layer:
            self.layer.append([])
        self.layer[layer].append(item)

if __name__ == '__main__':
    sys.stdout.write("Hello, world!")
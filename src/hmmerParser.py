#!/usr/bin/python
#
# Python library for parsing hmmer outputs

import sys
import os
import subprocess
import random
import re

alignmentLine        = "\s+([0-9]+)\s+([!?])\s+([0-9e\-\+\.]+)\s+([0-9e\-\+\.]+)\s+([0-9e\-\+\.]+)\s+([0-9e\-\+\.]+)\s+([0-9e\-\+\.]+)\s+([0-9e\-\+\.]+)\s+([\]\[\.]{2})\s+([0-9e\-\+\.]+)\s+([0-9e\-\+\.]+)\s+\.{2}\s+([0-9e\-\+\.]+)\s+([0-9e\-\+\.]+)\s+\.{2}\s+([0-9e\-\+\.]+)"

class domainHitClass:
    def __init__(self,line):    
        #                                                                            --- full sequence --- -------------- this domain -------------   hmm coord   ali coord   env coord
        # target name        accession   tlen query name           accession   qlen   E-value  score  bias   #  of  c-Evalue  i-Evalue  score  bias  from    to  from    to  from    to  acc description of target
        #------------------- ---------- ----- -------------------- ---------- ----- --------- ------ ----- --- --- --------- --------- ------ ----- ----- ----- ----- ----- ----- ----- ---- ---------------------
        # WP_010911285.1       -            549 YscJ_FliF            PF01514.16   194   7.9e-61  203.8   0.1   1   1   5.5e-64   1.9e-60  202.5   0.1     2   193    27   215    26   216 0.96 flagellar M-ring protein FliF [Mesorhizobium loti]
        # WP_044549206.1       -            287 YscJ_FliF            PF01514.16   194   1.3e-27   95.4   0.0   1   1   5.4e-31   1.9e-27   94.9   0.0    21   185    33   198    17   202 0.87 EscJ/YscJ/HrcJ family type III secretion inner membrane ring protein [Mesorhizobium loti]

        data = re.split("\s+",line)

        self.targetname            = data[0]
        self.accession             = data[1]
        self.tlen                  = int(data[2])
        self.queryname             = data[3]
        self.accession             = data[4]
        self.qlen                  = int(data[5])
        self.Evalue                = float(data[6])
        self.score                 = float(data[7])
        self.bias                  = float(data[8])
        self.num                   = int(data[9])
        self.of                    = int(data[10])
        self.cEvalue               = float(data[11])
        self.iEvalue               = float(data[12])
        self.score                 = float(data[13])
        self.bias                  = float(data[14])
        self.hmm_from              = int(data[15])
        self.hmm_to                = int(data[16])
        self.ali_from              = int(data[17])
        self.ali_to                = int(data[18])
        self.env_from              = int(data[19])
        self.env_to                = int(data[20])
        self.acc                   = float(data[21])
        self.description_of_target = data[22]

class hmmhitClass:
   
    def __init__(self, line):    
        # Takes as input the text output from hmmsearch 
        # and populates class variables relevant to the 
        # search
        
        # tblout # #                                                               --- full sequence ---- --- best 1 domain ---- --- domain number estimation ----
        # target name        accession  query name           accession    E-value  score  bias   E-value  score  bias   exp reg clu  ov env dom rep inc description of target
        #------------------- ---------- -------------------- ---------- --------- ------ ----- --------- ------ -----   --- --- --- --- --- --- --- --- ---------------------
        # WP_010911285.1       -          YscJ_FliF            PF01514.16   7.9e-61  203.8   0.1   1.9e-60  202.5   0.1   1.7   1   0   0   1   1   1   1 flagellar M-ring protein FliF [Mesorhizobium loti]
        # WP_044549206.1       -          YscJ_FliF            PF01514.16   1.3e-27   95.4   0.0   1.9e-27   94.9   0.0   1.2   1   0   0   1   1   1   1 EscJ/YscJ/HrcJ family type III secretion inner membrane ring protein [Mesorhizobium loti]

        data = re.split("\s+", line)
        
        self.targetname            = data[0]
        self.accession             = data[1]
        
        self.query_name            = data[2]
        self.query_accession       = data[3]
        
        self.full_Evalue           = float(data[4])
        self.full_score            = float(data[5])
        self.full_bias             = float(data[6])
        
        self.best_Evalue           = float(data[7])
        self.best_score            = float(data[8])
        self.best_bias             = float(data[9])
        
        self.exp                   = float(data[10])
        
        self.reg                   = int(data[11])
        self.clu                   = int(data[12])
        self.ov                    = int(data[13])
        self.env                   = int(data[14])
        self.dom                   = int(data[15])
        self.rep                   = int(data[16])
        self.inc                   = int(data[17])
        self.description_of_target = data[18]
        self.domainHits            = []

class HMMSearchClass:

    def __init__(self,genomeFile,HMM,genomeFileFormat="FAA"):
        # genomeFileFormat can be GPFF or FAA

        #
        # Run hmmsearch and work through header
        #
        self.hits     = {}
        randomSuffix  = str(int(random.random()*1000000))
        tbloutFile    = ".tbloutFile"+randomSuffix
        domtbloutFile = ".domtbloutFile"+randomSuffix
        
        commands = ["hmmsearch", "--tblout", tbloutFile, "--domtblout", domtbloutFile, HMM, genomeFile]

        try:
            process  = subprocess.Popen(commands, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        except:
            sys.stderr.write("Error executing "+" ".join(commands)+"\n\n")
            sys.exit(1)
        
        stdout, stderr = process.communicate()
        if stderr != '':
            sys.stderr.write("Command: "+" ".join(commands)+"\n")
            sys.exit("ERRORS:\n"+stderr+"\n")
        
        #
        # Parse hits to different proteins into hits{}
        #
        tblout = open(tbloutFile, 'r')
        lines = tblout.readlines()
        lineCount= 0
        while lineCount < len(lines):
            line = lines[lineCount]
            if line[0] != "#":
                if genomeFileFormat == "GPFF":
                    # Deal with messed-up output format of GFF which inserts
                    # carriage-returns after the second and last columns
                    line = line.strip()
                    line = line + lines[lineCount+1]
                    line = line.strip()
                    lineCount += 2
                hit = hmmhitClass(line)
                self.hits[hit.targetname] = hit
            lineCount += 1
        tblout.close()

        #
        # Now parse more detailed domain hit information and feed more into hits{}
        #
        
        domtblout = open(domtbloutFile, 'r')
        lines = domtblout.readlines()
        
        
                
                
        lineCount= 0
        while lineCount < len(lines):
            line = lines[lineCount]
            if line[0] != "#":
                if genomeFileFormat == "GPFF":
                    # Deal with messed-up output format of GFF which inserts
                    # carriage-returns after the second and last columns
                    line = line.strip()
                    line = line + lines[lineCount+1]
                    line = line.strip()
                    lineCount += 2
                domhit = domainHitClass(line)
                self.hits[domhit.targetname].domainHits.append(domhit)
                self.hits[domhit.targetname].domainHits.sort(key=lambda x: x.iEvalue)
            lineCount += 1                 
        domtblout.close()

        # Delete temporary hmmsearch output files    
        os.system("rm "+tbloutFile)
        os.system("rm "+domtbloutFile)
        
    def inHmmSearch(self, candidate):
        '''
        returns True if the candidate is in the search
         
        'candidate' is of type hmmhitClass
               
        '''
    
        for result in self.hits:
            if candidate.targetname == result:
                return result
        return False
    
    def nameInHmmSearch(self, candidate):
        '''
        returns True if the candidate is in the search
         
        'candidate' is a gene name
        
        Returns relevant hmmhitClass, or False if not found
        
        '''
    
        for result in self.hits:
            if candidate == result:
                return self.hits[result]
        return False
                
if __name__ == "__main__":
    print "Hello, world!"
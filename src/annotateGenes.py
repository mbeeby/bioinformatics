#!/usr/bin/python
# Takes a two-column file (ignoring empty lines)
# with desired name and PFxxxxx number and (1) works
# out current PFxxxxx version, then outputs as
# <desired name>.hmm

import numpy
import random
import sys
from hmmerParser import HMMSearchClass

# --------------------------------------------

class organismClass:
    
    def __init__(self,genomeFile,HMMDir, annotateAll = True, organism=None, taxID=None): 
        # Get a list of families

        self.organism   = None
        self.taxID      = None
        self.families   = {}
        self.genomeFile = genomeFile
        self.HMMDir     = HMMDir
        
        if organism:
            self.organism = organism
        if taxID:
            self.taxID    = taxID
        if annotateAll:
            self.annotateAll()
                       
    def annotateAll(self):
        self.annotateCoreFlagellarGenes()
        self.annotateFlagellins()
        self.annotateRibosomalProteins()
        self.getFliY()
        self.getFliN()
        self.getFlgE()
        
    def annotateFlagellins(self): 
        self.getFliCFlgL()
        
    def annotateRibosomalProteins(self):
        straightforwardsRibosomalProteins = \
            ["L1", "L10", "L13", "L14", "L15", \
             "L16", "L17", "L18", "L2", "L3", "L4", "L5"]

        for ribosomalProtein in straightforwardsRibosomalProteins:
            self.getStraightforwardsProtein(ribosomalProtein, [ribosomalProtein+".hmm"])
        
    def annotateCoreFlagellarGenes(self):
        # Calculate the estimated number of flagellar systems by taking the
        # average number of core flagellar proteins
        
        self.getFliF()
        self.getFliG()
        self.getFliM()
        coreFlagellarProteins = ["FliF", "FliG", "FliM"]
        total = []
        for coreFlagellarProtein in coreFlagellarProteins:
            total.append(len(self.families[coreFlagellarProtein]))
        self.number_flagellar_systems = numpy.mean(total)
          
    def annotateSRP(self):
        self.getFlhFFtsYFfh()
    
    def getLpoB(self):
        search          = {}
        finalCandidates = {}
    
        # Perform search(es) and get list of all candidates
        search["LpoB"]  = HMMSearchClass(self.genomeFile, self.HMMDir+"/LpoB.hmm")
        candidates      = getCandidates(search.values())
    
        # Establish criteria match and not matched by candidate (as a dictionary);
        # Then use criteria to judge which are actual final candidates
    
        for candidate in candidates:
            evidence = {}
            notes    = {}
            evidence["Hits LpoB"] = search["LpoB"].inHmmSearch(candidate)
    
            if "Hits LpoB" in evidence:
                finalCandidates[candidate.targetname] = \
                                     candidateClass(evidence=evidence,notes=notes)
                
        self.families["LpoB"] = finalCandidates
    
    def getTolA(self):
    
        search = {}
        # Perform search
        search["TolAc"]  = HMMSearchClass(self.genomeFile, self.HMMDir+"/TolAc.hmm")
    
        # Get list of possible candidates
        candidates = getCandidates(search.values())
    
        # Establish criteria match and not matched by candidate (as a bit array):
    
        for candidate in candidates:
            crit={}
    
            critNum = 0
            crit[critNum] = "Hits TolAc"
            if search["TolAc"].inHmmSearch(candidate):
                candidates[candidate] = candidates[candidate] + 2 ** critNum
    
        # Finally use criteria to judge which are actual final candidates
        finalCandidates = {}
        for candidate in candidates:
            score = candidates[candidate]
            if critMatch(score, 0):
                finalCandidates[candidate.targetname] = ""
        self.families["TolA"] = finalCandidates
    
    def getRcsF(self):
    
        search = {}
        # Perform search
        search["RcsF"]  = HMMSearchClass(self.genomeFile, self.HMMDir+"/RcsF.hmm")
    
        # Get list of possible candidates
        candidates = getCandidates(search.values())
    
        # Establish criteria match and not matched by candidate (as a bit array):
    
        for candidate in candidates:
            crit={}
    
            critNum = 0
            crit[critNum] = "Hits RcsF"
            if search["RcsF"].inHmmSearch(candidate):
                candidates[candidate] = candidates[candidate] + 2 ** critNum
    
        # Finally use criteria to judge which are actual final candidates
        finalCandidates = {}
        for candidate in candidates:
            score = candidates[candidate]
            if critMatch(score, 0):
                finalCandidates[candidate.targetname] = ""
        self.families["RcsF"] = finalCandidates
    
    def getIgaA(self):
    
        # Perform search
        search  = HMMSearchClass(self.genomeFile, self.HMMDir+"/IgaA.hmm")
    
        # Get list of possible candidates
        candidates = getCandidates([search])
    
        # Establish criteria match and not matched by candidate (as a bit array):
    
        for candidate in candidates:
            crit={}
    
            critNum = 0
            crit[critNum] = "Hits IgaA"
            if search.inHmmSearch(candidate):
                candidates[candidate] = candidates[candidate] + 2 ** critNum
    
        # Finally use criteria to judge which are actual final candidates
        finalCandidates = {}
        for candidate in candidates:
            score = candidates[candidate]
            if critMatch(score, 0):
                finalCandidates[candidate.targetname] = ""
        self.families["IgaA"] = finalCandidates
        
    def getLpp(self):
    
        # Perform search
        Lpp_search  = HMMSearchClass(self.genomeFile, self.HMMDir+"/Lpp.hmm")
    
        # Get list of possible candidates
        candidates = getCandidates([Lpp_search])
    
        # Establish criteria match and not matched by candidate (as a bit array):
    
        for candidate in candidates:
            crit={}
    
            critNum = 0
            crit[critNum] = "Hits Lpp"
            if Lpp_search.inHmmSearch(candidate):
                candidates[candidate] = candidates[candidate] + 2 ** critNum
    
        # Finally use criteria to judge which are actual final candidates
        finalCandidates = {}
        for candidate in candidates:
            score = candidates[candidate]
            if critMatch(score, 0):
                finalCandidates[candidate.targetname] = ""
        self.families["Lpp"] = finalCandidates
    
    def getTonB(self):
    
        # Perform search
        TonBn_search  = HMMSearchClass(self.genomeFile, self.HMMDir+"/TonBn.hmm")
        TonBc_search  = HMMSearchClass(self.genomeFile, self.HMMDir+"/TonBc.hmm")
        TonB_2_search = HMMSearchClass(self.genomeFile, self.HMMDir+"/TonB_2.hmm")
    
        # Get list of possible candidates
        candidates = getCandidates([TonBn_search, TonBc_search, TonB_2_search])
    
        # Establish criteria match and not matched by candidate (as a bit array):
    
        for candidate in candidates:
            crit={}
    
            critNum = 0
            crit[critNum] = "Hits TonBn"
            if TonBn_search.inHmmSearch(candidate):
                candidates[candidate] = candidates[candidate] + 2 ** critNum
    
            critNum = 1
            crit[critNum] = "Hits TonBc"
            if TonBc_search.inHmmSearch(candidate):
                candidates[candidate] = candidates[candidate] + 2 ** critNum
    
            critNum = 2
            crit[critNum] = "Hits TonB_2 (alternative to TonBc)"
            if TonBc_search.inHmmSearch(candidate):
                candidates[candidate] = candidates[candidate] + 2 ** critNum
    
        # Finally use criteria to judge which are actual final candidates
        finalCandidates = {}
        for candidate in candidates:
            score = candidates[candidate]
            if critMatch(score, 0) and (critMatch(score, 1) or critMatch(score,2)):
                finalCandidates[candidate.targetname] = ""
        self.families["TonB"] = finalCandidates
    
    def getStraightforwardsProtein(self, proteinName, HMMFiles, HMMDir=""):
        # Returns all candidates who hit all HMMs in [HMMFiles]
    
        # Inititialize
        search     = {}
        candidates = {}
    
        if HMMDir == "":
            HMMDir = self.HMMDir
        
        # Perform search
        for HMMFile in HMMFiles:
            search[HMMFile]   = HMMSearchClass(self.genomeFile, HMMDir+"/"+HMMFile)
            # Get list of possible candidates
            for hit in search[HMMFile].hits:
                candidates[search[HMMFile].hits[hit]] = 0
    
        # Establish criteria match and not matched by candidate (as a bit array):
    
        for candidate in candidates:
            crit={}
            critNum = 0
    
            for HMMFile in HMMFiles:
                crit[critNum] = "Hits "+HMMFile
                if search[HMMFile].inHmmSearch(candidate):
                    candidates[candidate] = candidates[candidate] + 2 ** critNum
                critNum += 1
    
        # Finally use criteria to judge which are actual final candidates
        finalCandidates = {}
        for candidate in candidates:
            trueCandidate = True
            
            # Make sure that the candidate was hit by all of the HMMs:
            for critCount in range(critNum-1): 
                if not critMatch(candidates[candidate], critCount):
                    trueCandidate = False
                    
            if trueCandidate:
                finalCandidates[candidate.targetname] = candidateClass()
                finalCandidates[candidate.targetname].notes["notes"] = []
                                
                for domainHit in search[HMMFile].hits[candidate.targetname].domainHits:
                                        
                    finalCandidates[candidate.targetname].notes["notes"].append({"from":domainHit.ali_from,"to":domainHit.ali_to, "iEvalue":domainHit.iEvalue})
                
                    
                #ali_from = search[HMMFile].hits[candidate.targetname].domainHits[0].ali_from
                #ali_to   = search[HMMFile].hits[candidate.targetname].domainHits[0].ali_to
            
                #finalCandidates[candidate.targetname].notes["notes"] = {"from":ali_from,"to":ali_to}
                
        self.families[proteinName] = finalCandidates
                 
    
    def getFliY(self):
    
        # Perform search
        FliM_search   = HMMSearchClass(self.genomeFile, self.HMMDir+"/FliM.hmm")
        FliMNc_search = HMMSearchClass(self.genomeFile, self.HMMDir+"/FliMNc.hmm")
        # HMMDir+"/FliNn.hmm"  # I don't use FliNn because it doesn't match Campy FliN (or Campy FliY for that matter)
    
        # Get list of possible FliY candidates, and also get FliM candidates
        # to rule out from FliY list.
        FliY_candidates = {}
        FliM_candidates = {}
        for FliM in FliM_search.hits:
            FliM_candidates[FliM_search.hits[FliM]] = 0
        for FliMNc in FliMNc_search.hits:
            FliY_candidates[FliMNc_search.hits[FliMNc]] = 0
    
        # Establish criteria match and not matched by candidate (as a bit array):
    
        for candidate in FliY_candidates:
            crit={}
    
            critNum = 0
            crit[critNum] = "DOESN'T hit FliM"
            if not FliM_search.inHmmSearch(candidate):
                FliY_candidates[candidate] = FliY_candidates[candidate] + 2 ** critNum
    
            critNum = 1
            crit[critNum] = "Hits FliMNc"
            if FliMNc_search.inHmmSearch(candidate):
                FliY_candidates[candidate] = FliY_candidates[candidate] + 2 ** critNum
    
            critNum = 2
            crit[critNum] = "HMM alignment start > 100 aa into sequence"
            if candidate.domainHits[0].ali_from > 100:
                FliY_candidates[candidate] = FliY_candidates[candidate] + 2 ** critNum
    
        # Finally use criteria to judge which are actual final candidates
        finalCandidates = {}
        for candidate in FliY_candidates:
            score = FliY_candidates[candidate]
            if critMatch(score, 0) and critMatch(score,1) and critMatch(score,2):
                finalCandidates[candidate.targetname] = ""
    
        self.families["FliY"] = finalCandidates
    
    def getFliN(self):
    
        # Perform search
        FliM_search   = HMMSearchClass(self.genomeFile, self.HMMDir+"/FliM.hmm")
        FliMNc_search = HMMSearchClass(self.genomeFile, self.HMMDir+"/FliMNc.hmm")
        # HMMDir+"/FliNn.hmm"  # I don't use FliNn because it doesn't match Campy FliN (or Campy FliY for that matter)
    
        # Get list of possible FliN candidates, and also get FliM candidates
        # to rule out from FliN list.
        FliN_candidates = {}
        FliM_candidates = {}
        for FliM in FliM_search.hits:
            FliM_candidates[FliM_search.hits[FliM]] = 0
        for FliMNc in FliMNc_search.hits:
            FliN_candidates[FliMNc_search.hits[FliMNc]] = 0
    
        # Establish criteria match and not matched by candidate (as a bit array):
    
        for candidate in FliN_candidates:
            crit={}
    
            critNum = 0
            crit[critNum] = "DOESN'T hit FliM"
            if not FliM_search.inHmmSearch(candidate):
                FliN_candidates[candidate] = FliN_candidates[candidate] + 2 ** critNum
    
            critNum = 1
            crit[critNum] = "Hits FliMNc"
            if FliMNc_search.inHmmSearch(candidate):
                FliN_candidates[candidate] = FliN_candidates[candidate] + 2 ** critNum
    
            critNum = 2
            crit[critNum] = "HMM alignment start < 100 aa into sequence"
            if candidate.domainHits[0].ali_from < 100:
                FliN_candidates[candidate] = FliN_candidates[candidate] + 2 ** critNum
    
        # Finally use criteria to judge which are actual final candidates
        finalCandidates = {}
        for candidate in FliN_candidates:
            score = FliN_candidates[candidate]
            if critMatch(score, 0) and critMatch(score,1) and critMatch(score,2):
                finalCandidates[candidate.targetname] = ""
        self.families["FliN"] = finalCandidates
    
    def getFliF(self):
    
        # Perform search
        FliF_search  = HMMSearchClass(self.genomeFile, self.HMMDir+"/FliF.hmm")
        FliFc_search = HMMSearchClass(self.genomeFile, self.HMMDir+"/FliFc.hmm")
    
        # Get list of possible candidates
        FliF_candidates = {}
        for FliF in FliF_search.hits:
            FliF_candidates[FliF_search.hits[FliF]] = 0
        for FliFc in FliFc_search.hits:
            FliF_candidates[FliFc_search.hits[FliFc]] = 0
    
        # Establish criteria match and not matched by candidate (as a bit array):
    
        for candidate in FliF_candidates:
    
            # criterion number 0: hits FliF
            critNum = 0
            if FliF_search.inHmmSearch(candidate):
                FliF_candidates[candidate] = FliF_candidates[candidate] + 2 ** critNum
    
            # criterion number 1: hits FliFc  
            critNum = 1
            if FliFc_search.inHmmSearch(candidate):
                FliF_candidates[candidate] = FliF_candidates[candidate] + 2 ** critNum
    
        # Finally use criteria to judge which are actual final candidates
        finalCandidates = {}
        for candidate in FliF_candidates:
            score = FliF_candidates[candidate]
            if critMatch(score, 0) and critMatch(score,1):
                finalCandidates[candidate.targetname] = ""
        self.families["FliF"] = finalCandidates
    
    def getFliM(self):
    
        # Perform search
        FliM_search   = HMMSearchClass(self.genomeFile, self.HMMDir+"/FliM.hmm")
        FliMNc_search = HMMSearchClass(self.genomeFile, self.HMMDir+"/FliMNc.hmm")
    
        # Get list of possible candidates
        FliM_candidates = {}
        for FliM in FliM_search.hits:
            FliM_candidates[FliM_search.hits[FliM]] = 0
        for FliMNc in FliMNc_search.hits:
            FliM_candidates[FliMNc_search.hits[FliMNc]] = 0
    
        # Establish criteria match and not matched by candidate (as a bit array):
    
        for candidate in FliM_candidates:
    
            # criterion number 0: hits FliM
            critNum = 0
            if FliM_search.inHmmSearch(candidate):
                FliM_candidates[candidate] = FliM_candidates[candidate] + 2 ** critNum
    
            # criterion number 1: hits FliMNc  
            critNum = 1
            if FliMNc_search.inHmmSearch(candidate):
                FliM_candidates[candidate] = FliM_candidates[candidate] + 2 ** critNum
    
        # Finally use criteria to judge which are actual final candidates
        finalCandidates = {}
        for candidate in FliM_candidates:
            score = FliM_candidates[candidate]
            if critMatch(score, 0) and critMatch(score,1):
                finalCandidates[candidate.targetname] = ""
        self.families["FliM"] = finalCandidates
    
    def getFliG(self):
    
        # Perform search
        FliGn_search = HMMSearchClass(self.genomeFile, self.HMMDir+"/FliGn.hmm")
        FliGm_search = HMMSearchClass(self.genomeFile, self.HMMDir+"/FliGm.hmm")
        FliGc_search = HMMSearchClass(self.genomeFile, self.HMMDir+"/FliGc.hmm")
    
        # Get list of possible candidates
        FliG_candidates = {}
        for FliGn in FliGn_search.hits:
            FliG_candidates[FliGn_search.hits[FliGn]] = 0
        for FliGm in FliGm_search.hits:
            FliG_candidates[FliGm_search.hits[FliGm]] = 0
        for FliGc in FliGc_search.hits:
            FliG_candidates[FliGc_search.hits[FliGc]] = 0
    
        # Establish criteria match and not matched by candidate (as a bit array):
    
        for candidate in FliG_candidates:
    
            # criterion number 0: hits FliGn   
            critNum = 0
            if FliGn_search.inHmmSearch(candidate):
                FliG_candidates[candidate] = FliG_candidates[candidate] + 2 ** critNum
    
            # criterion number 1: hits FliGm   
            critNum = 1
            if FliGm_search.inHmmSearch(candidate):
                FliG_candidates[candidate] = FliG_candidates[candidate] + 2 ** critNum
    
            # criterion number 2: hits FliGc   
            critNum = 2
            if FliGc_search.inHmmSearch(candidate):
                FliG_candidates[candidate] = FliG_candidates[candidate] + 2 ** critNum
    
        # Finally use criteria to judge which are actual final candidates
        finalCandidates = {}
        for candidate in FliG_candidates:
            score = FliG_candidates[candidate]
            if critMatch(score, 0) and critMatch(score,1) and critMatch(score,2):
                finalCandidates[candidate.targetname] = ""
        self.families["FliG"] = finalCandidates
    
    def getFlhFFtsYFfh(self):
    
        # Perform search
        SRP54n_search = HMMSearchClass(self.genomeFile, self.HMMDir+"/SRP54n.hmm")
        SRP54_search  = HMMSearchClass(self.genomeFile, self.HMMDir+"/SRP54.hmm")
        SPB_search    = HMMSearchClass(self.genomeFile, self.HMMDir+"/SPB.hmm")
        
        CustomFlhF_search = HMMSearchClass(self.genomeFile, self.HMMDir+"/custom_FlhF.hmm")
        CustomFtsY_search = HMMSearchClass(self.genomeFile, self.HMMDir+"/custom_FtsY.hmm")
    
        # Get list of possible candidates. Keys are gene names, values hmmhitClass
        candidates = getCandidates([SRP54n_search, SRP54_search, SPB_search, CustomFlhF_search])
        
        # Establish criteria match and not matched by candidate (as a bit array):
    
        for candidate in candidates:
            if SRP54n_search.nameInHmmSearch(candidate):
                candidates[candidate].evidence["Hits SRP54n"] = True 
            if SRP54_search.nameInHmmSearch(candidate):
                candidates[candidate].evidence["Hits SRP54"] = True
            if SPB_search.nameInHmmSearch(candidate):
                candidates[candidate].evidence["Hits SPB"] = True
            if CustomFlhF_search.nameInHmmSearch(candidate):
                hit = CustomFlhF_search.nameInHmmSearch(candidate)
                if hit.full_Evalue < 1e-20:
                    candidates[candidate].evidence["Hits custom FlhF"] = True
        
        # Finally use criteria to judge which are actual final candidates
        FlhF_candidates = {}
        FtsY_candidates = {}
        Ffh_candidates  = {}
        
        
        for candidate in candidates:
            curCand = candidateClass()
            curCand.evidence = candidates[candidate].evidence
            if "Hits SRP54n" in candidates[candidate].evidence:
                SRP54n_first_domHit = SRP54n_search.hits[candidate].domainHits[0]
                curCand.notes['length'] = SRP54n_first_domHit.tlen
                curCand.notes['SRP54n'] = []
                curCand.notes['SRP54n'].append({"from":SRP54n_first_domHit.ali_from,"to":SRP54n_first_domHit.ali_to})
            if "Hits SRP54" in candidates[candidate].evidence:
                SRP54_first_domHit = SRP54_search.hits[candidate].domainHits[0]
                curCand.notes['length'] = SRP54_first_domHit.tlen
                curCand.notes['SRP54']  = []                                   
                curCand.notes['SRP54'].append({"from":SRP54_first_domHit.ali_from,"to":SRP54_first_domHit.ali_to})
            if "Hits SPB" in candidates[candidate].evidence:
                SPB_first_domHit = SPB_search.hits[candidate].domainHits[0]
                curCand.notes['length'] = SPB_first_domHit.tlen
                curCand.notes['SPB']  = []                                   
                curCand.notes['SPB'].append({"from":SPB_first_domHit.ali_from,"to":SPB_first_domHit.ali_to}) 
    
            if ("Hits SRP54" in candidates[candidate].evidence):
                if ("Hits SRP54n" in candidates[candidate].evidence) and \
                   ("Hits SPB" in candidates[candidate].evidence):
                    # This protein is therefore Ffh
                    Ffh_candidates[candidate] = curCand
                else:
                    if ("Hits SRP54n" in candidates[candidate].evidence) or \
                       ("Hits custom FlhF" in candidates[candidate].evidence):
                         
                        # Get e-values for both and use this to differentiate   
                        FlhF_eValue = "Not detected"                
                        FlhF_hit    = CustomFlhF_search.nameInHmmSearch(candidate)
                        if FlhF_hit:
                            FlhF_eValue = FlhF_hit.full_Evalue
                            
                        FtsY_eValue = "Not detected"                
                        FtsY_hit    = CustomFtsY_search.nameInHmmSearch(candidate)
                        if FtsY_hit:
                            FtsY_eValue = FtsY_hit.full_Evalue
                         
                        # We now have enough information to differentiate FliC from FlgL     
                        if FlhF_eValue < FtsY_eValue:
                            # This protein is therefore FlhF
                            curCand.evidence["FlhF e-value ("+str(FlhF_eValue)+") < FtsY e-value ("+str(FtsY_eValue)+")"] = True
                            curCand.notes['SRP54_start'] = \
                                SRP54_first_domHit.ali_from + 1 
                            curCand.notes['SRP54_end'] = \
                                SRP54_first_domHit.ali_to - 1
                            FlhF_candidates[candidate] = curCand
                        else:
                            # This protein is therefore FtsY        
                            curCand.evidence["FtsY e-value ("+str(FtsY_eValue)+") < FlhF e-value ("+str(FlhF_eValue)+")"] = True
                            curCand.notes['SRP54_start'] = \
                                SRP54_first_domHit.ali_from + 1 
                            curCand.notes['SRP54_end'] = \
                                SRP54_first_domHit.ali_to - 1
                            if "SRP54n" in curCand.notes:
                                curCand.notes['A-domain length'] = \
                                    SRP54n_first_domHit.ali_from - 1
                            FtsY_candidates[candidate] = curCand
                    
        self.families["FlhF"] = FlhF_candidates
        self.families["FtsY"] = FtsY_candidates
        self.families["FfH"]  = Ffh_candidates
    
    def getFliCFlgL(self):
    
        # Perform search
        FliCn_search       = HMMSearchClass(self.genomeFile, self.HMMDir+"/Flagellin_N.hmm")
        FliCc_search       = HMMSearchClass(self.genomeFile, self.HMMDir+"/Flagellin_C.hmm")
        FlorianFliC_search = HMMSearchClass(self.genomeFile, self.HMMDir+"/FlorianFliC.hmm")
        FlorianFlgL_search = HMMSearchClass(self.genomeFile, self.HMMDir+"/FlorianFlgL.hmm")
        
    
        # Get list of possible candidates. Keys are gene names, values hmmhitClass
        candidates = {}
        for FliCn in FliCn_search.hits:
            candidates[FliCn] = candidateClass()
        for FliCc in FliCc_search.hits:
            candidates[FliCc] = candidateClass()
    
        # Establish criteria match and not matched by candidate (as a bit array):
    
        for candidate in candidates:
            if FliCn_search.nameInHmmSearch(candidate):
                candidates[candidate].evidence["Hits FliCn"] = True 
            if FliCc_search.nameInHmmSearch(candidate):
                candidates[candidate].evidence["Hits FliCc"] = True
        
        # Finally use criteria to judge which are actual final candidates
        FliC_candidates = {}
        FlgL_candidates = {}
        for candidate in candidates:
            
            if ("Hits FliCc" in candidates[candidate].evidence) and \
               ("Hits FliCn" in candidates[candidate].evidence):
                
                FliC_eValue = "Not detected"
                FliC_hit    = FlorianFliC_search.nameInHmmSearch(candidate)
                if FliC_hit:
                    FliC_eValue = FliC_hit.full_Evalue
                    
                FlgL_eValue = "Not detected"                
                FlgL_hit    = FlorianFlgL_search.nameInHmmSearch(candidate)
                if FlgL_hit:
                    FlgL_eValue = FlgL_hit.full_Evalue
                 
                # ----------------------------------
                # We now have enough information to differentiate FliC from FlgL     
                if FliC_eValue < FlgL_eValue:
                    # This protein is therefore FliC
                    curCand = candidateClass()
                    FliCn_first_domHit = FliCn_search.hits[candidate].domainHits[0]
                    FliCc_first_domHit = FliCc_search.hits[candidate].domainHits[0]
                    
                    curCand.notes['length'] = FliCn_first_domHit.tlen
                    curCand.notes['FliCn']  = []
                    curCand.notes['FliCc']  = []
    
                    curCand.notes['FliCn'].append({"from":FliCn_first_domHit.ali_from,"to":FliCn_first_domHit.ali_to})
                    curCand.notes['FliCc'].append({"from":FliCc_first_domHit.ali_from,"to":FliCc_first_domHit.ali_to})
                                   
                    curCand.evidence["Hits to both N- and C-terminus of flagellin HMM"] = True
                    curCand.evidence["FliC e-value ("+str(FliC_eValue)+") < FlgL e-value ("+str(FlgL_eValue)+")"] = True
                    
                    curCand.notes['insert_length'] = \
                        FliCc_first_domHit.ali_from - \
                        FliCn_first_domHit.ali_to -2
                    curCand.notes['insert_start'] = \
                        FliCn_first_domHit.ali_to + 1 
                    curCand.notes['insert_end'] = \
                        FliCc_first_domHit.ali_from - 1
                    FliC_candidates[candidate] = curCand
                else:
                    # This protein is therefore FlgL
                    curCand = candidateClass()
                    curCand.evidence["Hits to both N- and C-terminus of flagellin HMM"] = True
                    curCand.evidence["FlgL e-value ("+str(FlgL_eValue)+") < FliC e-value ("+str(FliC_eValue)+")"] = True
                    FlgL_candidates[candidate] = curCand
        
        self.families["FliC"] = FliC_candidates
        self.families["FlgL"] = FlgL_candidates
    
    def getFlgE(self):
        search          = {}
        finalCandidates = {}
    
        # Perform search(es) and get list of all candidates
        search["FlgE"]  = HMMSearchClass(self.genomeFile, self.HMMDir+"/FlaE.hmm")
        candidates      = getCandidates(search.values())
    
        # Establish criteria match and not matched by candidate (as a dictionary);
        # Then use criteria to judge which are actual final candidates
    
        for candidate in candidates:
            if search["FlgE"].nameInHmmSearch(candidate):
                candidates[candidate].evidence["Hits FlaE HMM"] = True 
          
        # Finally use criteria to judge which are actual final candidates
        
        for candidate in candidates:      
            if "Hits FlaE HMM" in candidates[candidate].evidence:
                finalCandidates[candidate] = candidates[candidate]
            
        self.families["FlgE"] = finalCandidates
    

# --------------------------------------------

def critMatch(score,criterion):
    return (score & (2**criterion)) == (2**criterion)

def getCandidates(searches):
    ''' 
    Takes multiple HMMSeachClass and outputs a dictionary of candidates 
    with gene names as keys, and empty candidateClass's as values. 
    '''
    
    candidates = {}
    for search in searches:
        for hit in search.hits:
            candidates[search.hits[hit].targetname] = candidateClass()
    return candidates

class candidateClass:
    
    def __init__(self):
        self.evidence = {}
        self.warnings = {}
        self.notes    = {}

#--------------------------------------------


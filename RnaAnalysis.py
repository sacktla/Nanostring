##RNA analysis
import pandas as pd
from RccProcessing import *
import scipy.stats.mstats
import numpy as np
import math
class RnaAnalysis:


    def __init__(self, rccProcessing,houseKeepingGenes):#, rlfDict):

        ##Combine all data into one
        self.samplesInfo = rccProcessing.getSamplesInfo()
        self.lanesInfo = rccProcessing.getLanesInfo()
        self.data = rccProcessing.getSamples()
        self.rccKeys = self.data.keys()
        self.hkg = houseKeepingGenes
        
        #Change all 0's to 1's
        for i in range(len(self.rccKeys)):
            self.data[self.rccKeys[i]] = self.data[self.rccKeys[i]].applymap(lambda x: 1 if x == '0' else self.parseInt(x) )##makes sure it is a number for the export to come out correctly
            
        
            
        
        
    ##function for catching parse int
    def parseInt(self,valueToParse):
        try:
            return int(valueToParse, 10)
        except ValueError:
            return valueToParse
        

    def calcPosGeoMean(self):##GOOD
        geomeanAll = []
        for i in range(len(self.rccKeys)):
            geomeanAll.extend(self.data[self.rccKeys[i]].loc['Positive']['Count'].tolist())
        geomeanAll = map(float,geomeanAll)
        
        self.posGeomean = scipy.stats.mstats.gmean(geomeanAll)
        


    def calcGeoMeanPerSample(self):##GOOD

        self.samplesGeomeans = {}

        for i in range(len(self.rccKeys)):

            positiveValues = self.data[self.rccKeys[i]].loc['Positive']['Count'].tolist()
            positiveValues = map(float,positiveValues)
            geoMean = scipy.stats.mstats.gmean(positiveValues)
            self.samplesGeomeans[self.rccKeys[i]] = geoMean
        
        



    def calcPcn(self):##GOOD
        self.pcn = {}
        for i in range(len(self.rccKeys)):
            self.pcn[self.rccKeys[i]] = self.posGeomean/self.samplesGeomeans[self.rccKeys[i]]

    def pcnNormalization(self):#GOOD
        self.pcnNormalized = {}
        for i in range(len(self.rccKeys)):
###i am going to remove endogenous as we want everything
            self.pcnNormalized[self.rccKeys[i]] = self.data[self.rccKeys[i]][['Name','Count']].set_index('Name').applymap(float)\
                .applymap(lambda x:x*self.pcn[self.rccKeys[i]])


    def calcLog2(self):#GOOD
        self.log2 = {}
        for i in range(len(self.rccKeys)):
            self.log2[self.rccKeys[i]] = self.data[self.rccKeys[i]]['Count'].to_frame()\
                .applymap(float)\
                .applymap(np.log2).set_index(self.data[self.rccKeys[i]]['Name'])
        
                
    def calcLog2Endogenous(self):#good
        self.log2Endogenous = {}
        for i in range(len(self.rccKeys)):
            self.log2Endogenous[self.rccKeys[i]] = self.data[self.rccKeys[i]].loc['Endogenous'][['Name','Count']].set_index('Name')\
                .applymap(float)\
                .applymap(np.log2)
        




    #We will use the raw counts
    def calcHkgCount(self):#GOOD

        self.hkgCount = {}
        
        for i in range(len(self.rccKeys)):
            hkgValues = []
            for j in range(len(self.hkg)):
                hkgValues.extend(self.data[self.rccKeys[i]].loc[self.data[self.rccKeys[i]]['Name'] == self.hkg[j]]['Count'])
            hkgValues = map(float,hkgValues)
            self.hkgCount[self.rccKeys[i]] = scipy.stats.mstats.gmean(hkgValues)
       


    ##Calc HKG geomean using log2

    def calcHkgCount2(self):##GOOD
        self.hkgCount2 = {}
        #print self.hkg

        for i in range(len(self.rccKeys)):
            hkgValues = []
            for j in range(len(self.hkg)):
                hkgValues.append(float(self.data[self.rccKeys[i]].loc[self.data[self.rccKeys[i]]['Name'] == self.hkg[j]]['Count']))
            #print len(hkgValues)
            #print hkgValues##Do geomean first then log2
            self.hkgCount2[self.rccKeys[i]] = math.log(scipy.stats.mstats.gmean(hkgValues),2)
        #print self.hkgCount2
        

    def calcDeltaCounts(self):##GOOD
        self.deltaCounts = {}

        for i in range(len(self.rccKeys)):
            self.deltaCounts[self.rccKeys[i]]= self.log2Endogenous[self.rccKeys[i]].applymap(lambda x: x - self.hkgCount2[self.rccKeys[i]])
            
    def calcDeltaCounts2(self):##GOOD ##Calculates and becomes NEN if a count is less than 30
        self.deltaCounts2 = {}

        for i in range(len(self.rccKeys)):
            self.deltaCounts2[self.rccKeys[i]]= self.log2Endogenous[self.rccKeys[i]].applymap(lambda x: x - self.hkgCount2[self.rccKeys[i]] if 2**x > 30 else "NEN")
            ##need to change log2Endogenous to log2PCN


    """def calcCcnFactor(self):#GOOD
        self.ccnFactor = {}
        tempCcn = []
        for i in range(len(self.rccKeys)):
            for j in range(len(self.hkg)):
                tempCcn.extend(self.data[self.rccKeys[i]].loc[self.data[self.rccKeys[i]]['Name'] == self.hkg[j]]['Count'])
            tempCcn = map(float,tempCcn)
        geomeanHKGAll = scipy.stats.mstats.gmean(tempCcn)##This works
       

        for i in range(len(self.rccKeys)):
            self.ccnFactor[self.rccKeys[i]] = geomeanHKGAll/self.hkgCount[self.rccKeys[i]]"""
            
        



    def calcCcnPcnArray(self):##GOOD

        self.ccnNormalized = {}

        for i in range(len(self.rccKeys)):
            self.ccnNormalized[self.rccKeys[i]] = self.data[self.rccKeys[i]].loc['Endogenous'][['Name','Count']].set_index('Name')\
                                                    .applymap(float)\
                                                    .applymap(lambda x:x*self.ccnFactor[self.rccKeys[i]])
                                                    


    def getPCNNormalized(self):
        return self.pcnNormalized
        
    def getDeltaCounts(self):
        return self.deltaCounts
        
    def getCCNNormalized(self):
        return self.ccnNormalized
        
    def getPCNFactor(self):
        return self.pcn
        
    def getKeys(self):
        return self.rccKeys
        
    def getCCNFactor(self):
        return self.ccnFactor
        
    def getRawData(self):
        return self.data
        
    def getDeltaCounts2(self):
        return self.deltaCounts2




'''
test = RccProcessing()
test.addRcc('20151208_151208 PBMC2 Titration_F1289 Frozen 1 week_04.RCC')

test.addRcc('20151208_151208 PBMC2 Titration_CMK86 Frozen_01.RCC')
test.addRcc('20151208_151208 PBMC2 Titration_F1289 Fresh_02.RCC')

files = open('genes.csv','r')
genes = files.readline().rstrip()
genes = genes.split(',')
testRNA = RnaAnalysis(test, genes)
testRNA.calcPosGeoMean()
testRNA.calcGeoMeanPerSample()
testRNA.calcPcn()
testRNA.pcnNormalization()
testRNA.calcLog2()
testRNA.calcLog2Endogenous()
testRNA.calcHkgCount()
testRNA.calcHkgCount2()
testRNA.calcDeltaCounts()
testRNA.calcCcnFactor()
testRNA.calcCcnPcnArray()
'''

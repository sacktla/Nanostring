from RccProcessing import *
import math
import scipy.stats.mstats
class ProteinAnalysis:

    def __init__(self, rccProcessing):
        ##Combine all data into one
        self.samplesInfo = rccProcessing.getSamplesInfo()
        self.lanesInfo = rccProcessing.getLanesInfo()
        self.data = rccProcessing.getSamples()

        self.rccKeys = self.data.keys()
        
        self.protein = {}
        self.proteinNeg = {}
    
        
        for key in self.rccKeys:
            
            self.protein[key] = self.data[key].loc['Protein']
            self.proteinNeg[key] = self.data[key].loc['Protein_NEG']
            
    
        
    def calcCTP(self):
        self.proteinLog2 = {}
        for key in self.rccKeys:
            self.proteinLog2[key] = self.protein[key][['Name','Count']].set_index('Name')\
                                                    .applymap(float)\
                                                    .applymap(lambda x:math.log(x,2))
                                                    
    def calcCTHKG(self, hkgList):
        self.hkgLog2 = {}
        for key in self.rccKeys:
            values = []
            for element in hkgList:
                values.append(self.protein[key].loc[self.protein[key]['Name'] == element]['Count'])
            values = map(float,values)
                
            geomeanHKG = scipy.stats.mstats.gmean(values)
            
            self.hkgLog2[key] = math.log(geomeanHKG,2)
            
        
            
    def calcDeltaCt(self):
        self.deltaCt = {}
        for key in self.rccKeys:
            self.deltaCt[key] = self.proteinLog2[key].applymap(lambda x: x - self.hkgLog2[key])
        
       
    def calcNegGmean(self):
        #Negative - Positive
        self.negGeoMean = {}
        for key in self.rccKeys:
            values = self.proteinNeg[key]['Count'].tolist()
            values = map(float,values)
            geomean = scipy.stats.mstats.gmean(values)
            self.negGeoMean[key] = math.log(geomean,2)
            
    def calcPosNegValue(self):
        self.posNegValue = {}
        for key in self.rccKeys:
            self.posNegValue[key] = self.hkgLog2[key] - self.negGeoMean[key] 
        
        
    ##Now we compare
    def compareValues(self):
        self.compared = {}
        for key in self.rccKeys:
            self.compared[key] = self.deltaCt[key].applymap(lambda x:"NEN" if x<(-1*self.posNegValue[key]+ math.log(4,2)) else x)
        
    
    def compareValuesSpecific(self,keysOfGenes):
        listOfGenes = list(keysOfGenes.keys())
        
        for key in self.rccKeys:
            for gene in listOfGenes:
                for gene in listOfGenes:
                    if self.deltaCt[key].loc[gene]['Count'] < (self.posNegValue[key]*-1 + math.log(float(keysOfGenes[gene]),2)):
                        self.compared[key].loc[gene]['Count'] = "NEN"
                    else:
                         self.compared[key].loc[gene]['Count'] = self.deltaCt[key].loc[gene]['Count']
                         
    def getComparedValues(self):
        return self.compared
        
    def getDeltaValues(self):
        return self.deltaCt
            
    def getKeys(self):
        return self.rccKeys
     
'''        
test = RccProcessing()


test.addRcc('20151202_151201 PBMC Titration_F1289 2mil cells A_09.RCC')



proteinTest = ProteinAnalysis(test)
proteinTest.calcCTP()
proteinTest.calcCTHKG(['CD45(HI30)|CD45|PTPRC|HI30|0','CD9(HI9a)|CD9|CD9|HI9a|0'])
proteinTest.calcDeltaCt()
proteinTest.calcNegGmean()
proteinTest.calcPosNegValue()
proteinTest.compareValues()
proteinTest.compareValuesSpecific({'OX40(Ber-ACT35)|CD134|TNFRSF4|Ber-ACT35|0':.5})

'''            
      
            
            
            
    
            
    
        
    
    
        
        
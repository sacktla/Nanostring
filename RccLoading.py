import os.path
from os import walk
class RccLoading:

    def __init__(self, pathname):

        self.totalFiles = []

        for(dirpaths, dirnames, filenames) in walk(pathname):
            for f in filenames:
                fname = os.path.join(dirpaths, f)
                assert(os.path.exists(fname))
                self.totalFiles.append(fname)


    def getTotalFiles(self):
        return self.totalFiles
      
    
    def getControls(self):
        controls = []
        names = []
        for i in range(len(self.totalFiles)):
            if self.totalFiles[i].endswith(".RCC"):
                tempRunSplit = self.totalFiles[i].split("_")
                tempRun = tempRunSplit[len(tempRunSplit)-2]
                tempSplit = tempRun.split()[0] 
                if not tempSplit in controls:
                    names.append([])
                    controls.append(tempSplit)
        for i in range(len(controls)):
            names[i].append(controls[i])
            names[i].append(0)
        
        for i in range(len(self.totalFiles)):
            if self.totalFiles[i].endswith(".RCC"):
                tempRunSplit = self.totalFiles[i].split("_")
                tempRun = tempRunSplit[len(tempRunSplit)-2]
                tempSplit = tempRun.split()[0]
                for j in range(len(names)):
                    control, quantity = names[j]
                    if control == tempSplit:
                        names[j][1] += 1
                        break
        
        ##FInd the max
        numbers = []
        for i in range(len(names)):
            numbers.append(names[i][1])
        
        maxNumber = max(numbers)
        oneRun = False
        if maxNumber == 1:
            oneRun = True
        
        maxNames = []
        possibleNames = []
        for i in range(len(names)):
            if names[i][1] == maxNumber:
                maxNames.append(names[i][0])
            elif names[i][1] > 1:
                possibleNames.append(names[i][0])
        
        return maxNames, possibleNames, oneRun
        
#Test high quantity
'''        
test = RccLoading('LEE Runs')
testList = test.getTotalFiles()

testProcessing = RccProcessing()

for element in testList:
    if element.endswith(".RCC"):
        testProcessing.addRcc(element)


testRNA = RnaAnalysis(testProcessing, ['GUSB','HPRT1','PGK1','GAPDH','TUBB','CLTC'])
testRNA.calcPosGeoMean()
testRNA.calcGeoMeanPerSample()
testRNA.calcPcn()
testRNA.pcnNormalization()
testRNA.calcLog2()
testRNA.calcHkgCount()
testRNA.calcHkgCount2()
testRNA.calcDeltaCounts()
testRNA.calcCcnFactor()
testRNA.calcCcnPcnArray()



proteinTest = ProteinAnalysis(testProcessing)
proteinTest.calcCTP()
proteinTest.calcCTHKG(['CD45(HI30)|CD45|PTPRC|HI30|0','CD9(HI9a)|CD9|CD9|HI9a|0'])
proteinTest.calcDeltaCt()
proteinTest.calcNegGmean()
proteinTest.calcPosNegValue()
proteinTest.compareValues()
proteinTest.compareValuesSpecific({'OX40(Ber-ACT35)|CD134|TNFRSF4|Ber-ACT35|0':.5})
'''

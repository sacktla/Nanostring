import pandas as pd
from Rcc import *
from Rlf import *
from scipy import stats
import numpy as st
class QualityCheck:

    def __init__(self, fileName, rccObject, rlfObject):

        self.fileName = fileName

        self.rccSample = rccObject.getSampleInfo()
        self.rccLane = rccObject.getLaneInfo()
        self.rccExperiment = rccObject.getExperiment()

        self.rlfName = rlfObject.getRlfName()
        self.rlfClass = rlfObject.getClassInformation()
        self.rlfContent = rlfObject.getRlfContent()

    #CHECKED
    def rlfChecker(self):
        return str(self.rlfName).lower().endswith(str(self.rccSample.loc['GeneRLF']['Results']).lower() + ".rlf")


    def geneCount(self):
        endogenousId = self.rlfClass[self.rlfClass['Class Name']=='Endogenous'].index.tolist()[0]
        return len(self.rccExperiment.loc['Endogenous'].index) == len(self.rlfContent.loc[endogenousId].index)##Error is in the last part

    def reqIdChecker(self):
        check = True
           
        if "omit" in self.fileName.lower():
            check = False
        return check

    def hkgCountsChecker(self,minCount, percentCutoff):
        hkCountList = self.rccExperiment.loc['Housekeeping']['Count'].tolist()##THis needs to be added the flexibility of the hkg
        #hkCountList = map(float,hkgCountList)
        if min([float(i) for i in hkCountList])<minCount:
            temp = map(float,hkCountList)
            new = [i for i in temp if i<minCount]
            if (float(len(new))/float(len(temp)) )*100<percentCutoff:#This logic doesn't make sense
                return False
            else:
                return True
        else:
            return True



    def fovRatioChecker(self, minFovRatio):
        fovCounted = float(self.rccLane.loc['FovCounted']['Results'])
        fovCount = float(self.rccLane.loc['FovCount']['Results'])

        fovRatio = (fovCounted/fovCount) * 100.0

        if fovRatio < minFovRatio:
            return False
        else: return True

    def bindingDensityChecker(self, lowBound, upBound):

        bindingDensity = float(self.rccLane.loc['BindingDensity']['Results'])
        if bindingDensity <= lowBound or bindingDensity >= upBound:
            return False
        else: return True

    def positiveControlLinearityChecker(self, rSquared):
        positiveControlConcentrations = [128,32,8,2,0.5,0.125]
        positiveDF = pd.DataFrame.sort(self.rccExperiment.loc['Positive'].set_index('Name'))['Count']
        listPositive = [float(i) for i in positiveDF.tolist()]

        _,_,r_value,_,_ = stats.linregress(listPositive, positiveControlConcentrations)
        if r_value**2 < rSquared:
            return False
        else: return True

    ##Confirm with Nate about this one!
    def positiveControlLimitDetectionChecker(self, sdFactor):

        positiveDF = pd.DataFrame.sort(self.rccExperiment.loc['Positive'].set_index('Name'))['Count']
        negativeDF = pd.DataFrame.sort(self.rccExperiment.loc['Negative'].set_index('Name'))['Count']

        listPositive = [float(i) for i in positiveDF.tolist()]
        listNegative = [float(i) for i in negativeDF.tolist()]

        fmPositiveValue = listPositive[4]
        negativeValuesMean = st.mean(listNegative)

        negativeValuesSD = st.std(listNegative)

        errorWindow = negativeValuesMean + (sdFactor*negativeValuesSD)

        if fmPositiveValue < errorWindow:
            return False
        else: return True

'''
testRCC = Rcc('20151202_151201 PBMC Titration_F1289 1mil cells B_08omit.RCC')
testRLF = Rlf('NS_CancerImmune_RNAProtein_1.1.rlf')
testQC = QualityCheck('20151202_151201 PBMC Titration_F1289 1mil cells B_08.RCC',testRCC, testRLF)
print testQC.positiveControlLimitDetectionChecker(0.05)

'''



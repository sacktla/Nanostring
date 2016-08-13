'''This class is going to store the elements of an RCC'''

import pandas as pd
from bs4 import BeautifulSoup
class Rcc:

    def __init__(self,fileName):

        fileRead = open(fileName, 'r')
        listDate = fileName.split('\\')
        sampleName = listDate[len(listDate)-1]
        date = sampleName.split('_')[0]
        #date = listDate[len(listDate)-1]##SPlit for last
        content = fileRead.read()
        fileRead.close()
        soup = BeautifulSoup(content, "lxml")

        ##Get all the attributes in the file
        sampleAttributes = soup.find_all('sample_attributes')[0].contents[0]
        laneAttributes= soup.find_all('lane_attributes')[0].contents[0]
        summary = soup.find_all('code_summary')[0].contents[0]

        sampleAttr = str(sampleAttributes).splitlines()
        laneAttr = str(laneAttributes).splitlines()
        summ = str(summary).splitlines()

        sample = []
        lane = []
        summary = []

        for info in sampleAttr:
            if len(info)>0:
                sample.append(tuple(info.split(',')))

        for info in laneAttr:
            if len(info)>0:
                lane.append(tuple(info.split(',')))

        for info in summ:
            if len(info)>0:
                summary.append(tuple(info.split(',')))

        self.sampleDF = pd.DataFrame(data=sample, columns=['Attributes','Results'])
        self.sampleDF = self.sampleDF.set_index('Attributes')
        self.sampleDF.loc['Date']['Results'] = date
        self.laneDF = pd.DataFrame(data=lane, columns=['Attributes', 'Results'])
        self.laneDF = self.laneDF.applymap(lambda x: 'LANE' if x=='ID' else x)
        self.laneDF = self.laneDF.set_index('Attributes')
        self.summaryDF = pd.DataFrame(data=summary)
        self.summaryDF.columns = self.summaryDF.iloc[0]
        self.summaryDF = self.summaryDF[1:]
        self.summaryDF = self.summaryDF.sort(columns=['CodeClass','Name'], axis=0)
        self.summaryDF = self.summaryDF.set_index('CodeClass')
        #self.summaryDF = self.summaryDP.loc['Endogenous'].sort(columns='Name',axis=)
        

    ##Returners
    def getSampleInfo(self):
        return self.sampleDF

    def getLaneInfo(self):
        return self.laneDF

    def getExperiment(self):
        return self.summaryDF





#testSampleInfo.set_index(attributes)
#print testSampleInfo
#print pd.DataFrame.sort(testExperiment.loc['Positive'].set_index('Name'))
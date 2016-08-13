from Rcc import *
#This might not be needed
class RccProcessing:

    ##This one should grab one by one and create
    ##a static variable that gets updated everytime
    ##you add an Rcc
    ##I need to add a reference dictionary
    def __init__(self):

        ##start the array of RCC arrays
        self.samplesInfo = {}
        self.lanesInfo = {}
        self.samples= {}


    def addRcc(self, rcc):
        ##use scannerID to match all three elements
        """
        @type rcc:Rcc
        """
        rcc = Rcc(rcc)
        sampleInfo = rcc.getSampleInfo()
        laneInfo = rcc.getLaneInfo()
        experiment = rcc.getExperiment()

        scannerId = str(rcc.getSampleInfo().loc['ID']['Results'])+"_"+str(rcc.getLaneInfo().loc['LANE']['Results'])#Modified the identifier
       

        self.samplesInfo[scannerId] = sampleInfo
        self.lanesInfo[scannerId] = laneInfo
        self.samples[scannerId] = experiment

    def getSamplesInfo(self):
        return self.samplesInfo

    def getLanesInfo(self):
        return self.lanesInfo

    def getSamples(self):
        return self.samples


'''
test = RccProcessing()
test.addRcc('/Users/jesuszaragoza/Desktop/Mixed RCC Files/Raw Data/20140402_ WNT Pathway Exp2_CAPAN2- 260ng_04.RCC')

test.addRcc('/Users/jesuszaragoza/Desktop/Mixed RCC Files/Raw Data/20140402_ WNT Pathway Exp2_CAPAN2+ 260ng_03.RCC')
test.addRcc('/Users/jesuszaragoza/Desktop/Mixed RCC Files/Raw Data/20140402_ WNT Pathway Exp2_HPAFII- 260ng_02.RCC')
test.addRcc('/Users/jesuszaragoza/Desktop/Mixed RCC Files/Raw Data/20140402_ WNT Pathway Exp2_XG- 260ng_06.RCC')

testCon =  test.getSamples()
keys= testCon.keys()

#result = pd.concat(testCon, axis= 1,join_axes=[testCon[0].index])
#testCon2 = result[['Count','Name']].loc['Housekeeping']
testA = testCon[keys[0]][['Name','Count']]
testB = testCon[keys[1]][['Name','Count']]
testC= testCon[keys[2]][['Name','Count']]
testD=testCon[keys[3]][['Name','Count']]
result = pd.merge(testA,testB,on='Name')
result2 = pd.merge(result,testC,on='Name')
result3=pd.merge(result2,testD,on='Name')
print result3


#testCon2 = pd.merge(left=testA,right=testB,left_index=True,left_on='Name',right_on='Name',suffixes=('_Chill','_Nice'))
#print testCon2.loc['Housekeeping']




##I think I need to do a merging in the class code or the
##test code.
'''
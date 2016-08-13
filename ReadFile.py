#!usr/bin/env python

class ReadFile:
	##I need to find a better way to determine file
    def printFile(self, file, finalLocation):
        cartridgeId = "x,CartridgeID missing"
        sampleNameInputNg = "x,Sample ID missing"
        fovCounted = "x,0"
        fovCount = "x,0"
        bindingDensity = "x,0"
        
        dataDict = {}
        for key in valueDictionary:
		dataDict[valueDictionary[key]] = []
	
        readFile = open(file, 'r')
        contents = readFile.read().splitlines()
        readFile.close()
        sampleAttributes = False
      
        for line in contents:
            if line.startswith("<Sample_Attributes>"):
                sampleAttributes = True
            if sampleAttributes:
                if line.startswith("ID"):
                    sampleNameInputNg = line
                    sampleAttributes = False
            if line.startswith("CartridgeID"):
                cartridgeId = line
            elif line.startswith("FovCounted"):
                fovCounted = line
            elif line.startswith("FovCount,"):
                fovCount = line
            elif line.startswith("BindingDensity"):
                bindingDensity = line
            
        for line in contents:
            temp = line.split(",")
            if temp[0] in dataDict:
                dataDict[temp[0]].append(temp[1] + "," + temp[3])
                
        for key in dataDict:
            dataDict[key].sort()
         
        headers = [cartridgeId, sampleNameInputNg, fovCount, fovCounted, bindingDensity]

        tempSplitFile = file.split("\\")
      
        headersParsed = [["x" , tempSplitFile[len(tempSplitFile) - 1]]]
        for i in range(0, len(headers)):
            headersParsed.append(headers[i].split(','))
				

        printFile = open(finalLocation, 'a')
        
        for i in range(0, len(headersParsed)):
            printFile.write(headersParsed[i][1] + "/")
            if i == len(headersParsed)- 1:
                printFile.write(',')
		
        keyList = dataDict.keys()
        keyList.sort()
        for key in keyList:
            for item in dataDict[key]:
                printFile.write(item.split(",")[1] + ",")
		
        printFile.write('\n')
        printFile.close()
		
	#we need to print the headers function
	
    def printHeader(self, nameOfRlf, finalLocation):
        
        openFile = open(nameOfRlf, 'r')
        contents = openFile.read().splitlines()

        #Obtain class types
        classTypes = []
        for line in contents:
            if line.startswith("ClassName"):
                classTypes.append(line.split("=")[1] + ',' + line.split("=")[0][9])
            if line.startswith("[Content]"):
                break

        ##Delete unecessary class types
        bindingIndex = 0
        purificationIndex = 0
        reservedIndex = 0
        for i in range(len(classTypes)):
            if classTypes[i].startswith("Binding"):
                bindingIndex = i
            elif classTypes[i].startswith("Purification"):
                purificationIndex = i
            elif classTypes[i].startswith("Reserved"):
                reservedIndex = i

        
        clearedClassTypes = []
        for i in range(len(classTypes)):
            if i == bindingIndex or i == purificationIndex or i == reservedIndex:
                continue
            clearedClassTypes.append(classTypes[i])
       
       ##NOw that we have the class types we need to create the headers
        numTypes = len(clearedClassTypes)
        data = {}
        global valueDictionary
        valueDictionary = {}
        for i in range(numTypes):
            temp = clearedClassTypes[i].split(",")
            valueDictionary[temp[1]]=temp[0]
            data[temp[1]] = []
           ##Up to here we have a working dictionary.
           ##Next thing to do is to append it and print it
    
		
        for line in contents:
            if line == "":
               continue
            temp = line.split(',')
            recordValue = temp[0]
            if recordValue[len(recordValue)-1] in valueDictionary and recordValue.startswith('Record') and not recordValue.startswith('RecordCount'):
                data[recordValue[len(recordValue)-1]].append(temp[3])
        #Sorted data alphabetically
            for key in data:
                data[key].sort()

        printFile = open(finalLocation, 'w')

        printFile.write('Gene Name,') 
        for key, value in sorted(valueDictionary.iteritems(), key=lambda(k,v):(v,k)):
            for item in data[key]:
                printFile.write(item + ",")

        printFile.write("\n")

        printFile.close()
        
# Using pandas as an object creator for an Rlf class


import pandas as pd


class Rlf:
    def __init__(self, fileName):

        readFile = open(fileName, 'r')
        contents = readFile.readlines()
        #readFile.close()

        # First we need to create the Rlf information attributes
        # Name of rlf
        tempName = fileName.split('\\')
        self.rlfName = tempName[len(tempName) - 1]

        header = False
        headerInformation = []

        # Go through file
        # Get all the Header information in an array.
        for line in contents:
            # if line starts with content we got to the table.
            if line.startswith('[Header]'):
                header = True
                continue
            if line.startswith('[Content]'):
                break
            if header:
                headerInformation.append(line.rstrip())

        content = False
        contentInformation = []

        # Get all the Content information in an array.
        for line in contents:

            if line.startswith('[Content]'):
                content = True
                continue
            if content:
                contentInformation.append(line.rstrip())

        classKey = {}
        for i in range(len(headerInformation)):
            if headerInformation[i].startswith('ClassCount'):
                classCount = headerInformation[i].split('=')[1]
            if headerInformation[i].startswith('ClassKey'):
                classKey[headerInformation[i].split('=')[1]] = -1
            if headerInformation[i].startswith('ClassName'):
                temp = headerInformation[i].split('=')
                classKey[temp[0][len(temp[0]) - 1]] = temp[1]

        ##Now lets create a relational table we can refer to when comparing to the content
        self.classDF = pd.DataFrame(classKey.items(),columns=['Class Key', 'Class Name'])
        self.classDF = self.classDF.set_index('Class Key')

        contentKey = []

        for i in range(len(contentInformation)):
            if contentInformation[i].startswith('ColumnCount'):
                columnCount = contentInformation[i].split('=')[1]
            if contentInformation[i].startswith('RecordCount'):
                recordCount = contentInformation[i].split('=')[1]
                continue
            if contentInformation[i].startswith('Columns'):
                columnHeader = contentInformation[i].split('=')[1]
            if contentInformation[i].startswith('Record'):
                contentKey.append(contentInformation[i].split('=')[1])

        newList = []
        for item in contentKey:
            tempList = item.split(',')
            newItem = tuple(tempList)
            newList.append(newItem)


        columnHeader = columnHeader.split(',')
        self.contentDF = pd.DataFrame.from_records(data=newList,columns=columnHeader)
        self.contentDF = self.contentDF.set_index('Classification')
    #Returns a pandas Data Frame with the class information
    def getClassInformation(self):
        return self.classDF

    #Returns a pandas Data Frame with the content information
    def getRlfContent(self):
        return self.contentDF

    def getRlfName(self):
        return self.rlfName

'''

test = Rlf('test.rlf')

#print test.getClassInformation()
#print test.getRlfContent()
index =  test.getClassInformation()
proteinIndex = index.loc[index['Class Name'] == "Protein"].index.tolist()
print len(proteinIndex)'''
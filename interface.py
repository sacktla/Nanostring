# -*- coding: utf-8 -*-
"""
Created on Tue Mar 01 12:40:29 2016

@author: jzaragoza
Interface for the Nanostring Normalization version 3.0
"""
#.loc raises a KeyError that needs to be addressed and caught
from PySide.QtGui import *
from PySide.QtCore import *
import sys
from nanostringNormalization2 import *
from Rlf import *
from ProteinAnalysis import *
from QualityCheck import *
from RnaAnalysis import *
from RccProcessing import *
from RccLoading import *
import ReadFile as rf
import os
import datetime

class MyApplication(QtGui.QMainWindow, Ui_MainWindow):
        def __init__(self, parent=None):
                super(MyApplication, self).__init__(parent)
                self.setupUi(self)
                
                ###Initial disability of Buttons
                self.confirmSelectionBtn.setEnabled(False)
                self.confirmNanostringQC.setEnabled(False)
                self.confirmGenoptixQC.setEnabled(False)
                self.confirmControls.setEnabled(False)
                self.actionRUN_EXPERIMENT.setEnabled(False)
                self.tab_5.setEnabled(False)
                
                ##Page 1 of interface
                ##Selection of Files
                self.selectRCCBtn.clicked.connect(self.findRCC)
                self.RLFSelectionBtn.clicked.connect(self.findRLF)
                self.outputSelectionBtn.clicked.connect(self.findOutput)
                self.confirmSelectionBtn.clicked.connect(self.confirm)
                
                ##Page 2 of interface
                ##QC parameters determined
                self.confirmNanostringQC.clicked.connect(self.changeQCSettings)
                
                ##Page 3 of interface
                self.confirmGenoptixQC.clicked.connect(self.changeGenQC)
                
                
                ##Page 5 of interface
                self.confirmControls.clicked.connect(self.controlQC)
                
                ##UP and down buttons
                self.upRCCBtn.clicked.connect(self.upRCC)
                self.downRCCBtn.clicked.connect(self.downRCC)
                self.upHkgBtn.clicked.connect(self.upHkg)
                self.downHkgBtn.clicked.connect(self.downHkg)
                
                ##Left and right buttons for controls
                self.pushButton_2.clicked.connect(self.moveToNotControl)##Towards not control
                self.pushButton.clicked.connect(self.moveToControl)##Towards control
                
                self.upProteinBtn_2.clicked.connect(self.toHKPProtein)
                self.downProteinBtn_2.clicked.connect(self.toEndogenousProtein)
                
                ##Run experiment
                self.actionRUN_EXPERIMENT.setShortcut('Ctrl+R')
                self.actionRUN_EXPERIMENT.triggered.connect(self.runExperiment)
                
                
                
                
                
                
        def findRCC(self):
            global rccPath
            rccPath = QtGui.QFileDialog.getExistingDirectory(self, "Select RCC Directory")
            self.RCCText.setText(rccPath)
            if not self.RCCText.text() == '' and not self.rlfTxt.text() == '' and not self.outputTxt.text() == '':
                self.confirmSelectionBtn.setEnabled(True)
                  
        def findRLF(self):
            global rlfPath
            rlfPath, _ = QFileDialog.getOpenFileName(self, "Find RLF")
            self.rlfTxt.setText(rlfPath)
            if not self.RCCText.text() == '' and not self.rlfTxt.text() == '' and not self.outputTxt.text() == '':
                self.confirmSelectionBtn.setEnabled(True)
           
        def findOutput(self):
            global outputPath
            outputPath = QtGui.QFileDialog.getExistingDirectory(self, "Select Output Directory")
            self.outputTxt.setText(outputPath)
            if not self.RCCText.text() == '' and not self.rlfTxt.text() == '' and not self.outputTxt.text() == '':
                self.confirmSelectionBtn.setEnabled(True)
            
        def confirm(self):
            confirmation = QtGui.QMessageBox.Cancel
            if rlfPath.endswith(".rlf") or rlfPath.endswith(".RLF"):
                if self.isRCC(rccPath):
                    confirmation = QtGui.QMessageBox.question(self,"Confirm", "YOU HAVE THE FOLLOWING SELECTIONS:\nLOCATION OF FILES:\n %(filelocation)s\nRLF LOCATION:\n %(rlflocation)s\nOUTPUT DATA LOCATION:\n%(datalocation)s" %\
                                {"filelocation":rccPath, "rlflocation":rlfPath,  "datalocation":outputPath}, QtGui.QMessageBox.Ok | QtGui.QMessageBox.Cancel, QtGui.QMessageBox.Cancel)
                else:
                    noRCCMessage = QMessageBox()
                    noRCCMessage.setText("No RCC files exist in directory selected!")
                    noRCCMessage.exec_()                
            else:
                noRlfMessage = QMessageBox()
                noRlfMessage.setText("No RLF file was selected!")
                noRlfMessage.exec_()
                
        
                
            if confirmation == QtGui.QMessageBox.Ok:
                
                self.confirmNanostringQC.setEnabled(True)
                
                global outputPath2 ##This where the outputFiles will go
                global errorLogLocation ##This is where the error log will go
                global rnaLocation 
                #global deltaProteinLocationPlot
                global controlsFile
                global qcFailed
                qcFailed = {}#This one will get filled as we go with the errors and manipulated in case they need to be addressed
		
        	     #The time will be associated per process
                time = datetime.datetime.now().strftime("%Y%B%d_%I%M")
                
                ##FILES being created ahead of time
                os.mkdir(outputPath + "/OutputFiles_" + time)
                outputPath2 = outputPath + "/OutputFiles_" + time
                errorLogLocation = outputPath2 + "/" + time + "_ErrorLog.csv"
                rnaLocation = outputPath2 + "/" + time + "_Data_Analysis.xlsx"
                #deltaProteinLocationPlot = outputPath2 + "/" + time + "_PROTEIN.xlsx"
                controlsFile = outputPath2 + "/" + time
                
                ##Populate housekeeping genes and endogenous genes
                rlf = Rlf(rlfPath)
                self.populateGeneBoxes(rlf)
                
        def changeQCSettings(self):
            global qcParameters
            
            self.confirmGenoptixQC.setEnabled(True)
            self.processRCCList.clear()
            self.nonProcessedRCCList.clear()
            
            qcParameters = {}
            qcParameters["imagingQC"] = self.fovRegistrationCheckBox.checkState()
            if qcParameters["imagingQC"] == QtCore.Qt.CheckState.Checked:
                qcParameters["imagingQCValue"] = self.fovRegistrationSB.value()
            qcParameters["bindingDensityQC"] = self.bindingDensityCheckBox.checkState()
            if qcParameters["bindingDensityQC"] == QtCore.Qt.CheckState.Checked:
                qcParameters["bindingDensityQCUpperValue"] = self.upperBoundSpinBox.value()
                qcParameters["bindingDensityQCLowerValue"] = self.lowerBoundSpinBox.value()
            qcParameters["positiveLinearityQC"] = self.positiveControlCheckBox.checkState()
            if qcParameters["positiveLinearityQC"] == QtCore.Qt.CheckState.Checked:
                qcParameters["positiveLinearityQCValue"] = self.posControlSpinBox.value()
            qcParameters["positiveDetectionQC"] = self.negControlCheck.checkState()
            if qcParameters["positiveDetectionQC"] == QtCore.Qt.CheckState.Checked:
                qcParameters["positiveDetectionQCValue"] = self.negControlSpinBox.value()
            qcParameters["pcnFactorQC"] = self.pcnRange.checkState()
            if qcParameters["pcnFactorQC"] == QtCore.Qt.CheckState.Checked:
                qcParameters["pcnFactorQCLowerValue"] = self.lowBoundPCN.value()
                qcParameters["pcnFactorQCUpperValue"] = self.upBoundPCN.value()
            
        def changeGenQC(self):
            
            self.listWidget.clear()
            self.listWidget_2.clear()
            self.processRCCList.clear()
            self.nonProcessedRCCList.clear() 
            
            qcParameters["backgroundLimitQC"] = self.backgroundTresholdCheck.checkState()
            if qcParameters["backgroundLimitQC"] == QtCore.Qt.CheckState.Checked:
                
                qcParameters["backgroundLimitQCValue"] = self.backgroundThresSpinBox.value()
                qcParameters["backgroundPercentQCValue"] = self.percentCutoffSpinBox.value()
            qcParameters["omittedQC"] = self.reqIdCheck.checkState()
            if qcParameters["omittedQC"] == QtCore.Qt.CheckState.Checked:
                qcParameters["omittedQCValue"] = True
                
            if proteinsAvailable:
                global proteins
                proteins = []
                for i in range(self.hkpList.count()):
                    proteins.append(self.hkpList.item(i).text())
                    
            global hkgList
            hkgList = []
            for i in range(self.defaultHkgList.count()):
                hkgList.append(self.defaultHkgList.item(i).text())
                
#i need to determine how to populate the controls            
            global totalRCC
            rcc = RccLoading(rccPath)
            totalRCC = rcc.getTotalFiles()
            
            global oneRun
            controls,possibleControls,oneRun = rcc.getControls()
            self.populateControlList(controls,possibleControls,oneRun)
            
            global noProteinTest
            if proteinsAvailable:
                if len(proteins) == 0:
                    #self.confirmControls.setEnabled(False)
                    completeBox = QMessageBox()
                    completeBox.setText("No Housekeeping Proteins have been selected!\nNo Protein analysis will be completed!")
                    completeBox.exec_()
                    noProteinTest = True
                else:
                    noProteinTest = False
            
            self.confirmControls.setEnabled(True) #Nate does  not care if no proteins are available for HKP
                
        def controlQC(self):
            #empty the confirm area just in case they re do 
            self.processRCCList.clear()
            self.nonProcessedRCCList.clear()
            self.actionRUN_EXPERIMENT.setEnabled(True)
            
            global processedFiles
            global nonProcessedFiles
            
            global backgroundFailures
            backgroundFailures= {}
				
            processedFiles = []
            nonProcessedFiles = []
            
            ##Prepare RccProcession
            rlfObject = Rlf(rlfPath)
            
            if oneRun==True:
                for i in range(len(totalRCC)):
                    ##QCCHECK
                    if(totalRCC[i].lower().endswith(".rcc")):
                        tempRCC = Rcc(totalRCC[i])
                        tempQC = QualityCheck(totalRCC[i],tempRCC, rlfObject)
                        failed = []
                        ##qcFailed[totalRCC[i]] =
                        if "omittedQCValue" in qcParameters.keys() and not tempQC.reqIdChecker():#This is good
                            failed.append("Omitted added to name")
                        if not tempQC.geneCount():#This is good
                            failed.append("Endogenous genes does not match RLF")
                        if not tempQC.rlfChecker():#This is good
                            failed.append("RCC does not match RLF name")
                        if "imagingQCValue" in qcParameters.keys() and not tempQC.fovRatioChecker(qcParameters["imagingQCValue"]):#This is good
                            failed.append("FOV ratio failed")
                        if "bindingDensityQCLowerValue" in qcParameters.keys() and not tempQC.bindingDensityChecker(qcParameters["bindingDensityQCLowerValue"], qcParameters["bindingDensityQCUpperValue"]):
                            failed.append("Binding density failed")
                        if "positiveLinearityQCValue" in qcParameters.keys() and not tempQC.positiveControlLinearityChecker(qcParameters["positiveLinearityQCValue"]):
                            failed.append("Positive Linearity failed")
                        if "positiveDetectionQCValue" in qcParameters.keys() and not tempQC.positiveControlLimitDetectionChecker(qcParameters["positiveDetectionQCValue"]):
                            failed.append("Positive detection failed")
                            
                        ##For the pcnFactor range we will take care of it after putting all the files together. ask nate
                        if "backgroundLimitQCValue" in qcParameters.keys() and not tempQC.hkgCountsChecker(qcParameters["backgroundLimitQCValue"],qcParameters["backgroundPercentQCValue"]):
                            ##Nate just wants a statement in the file
                            backgroundFailures[str(tempRCC.getSampleInfo().loc['ID']['Results'])] = "Failed background limit QC"
                        
                        if len(failed)==0:
                            processedFiles.append(totalRCC[i])
                        else:
                            nonProcessedFiles.append(totalRCC[i])
                            qcFailed[str(tempRCC.getSampleInfo().loc['ID']['Results'])] = failed
                            
            elif oneRun == False:
                global controls
                controls = []
                acceptedControl = self.controlLists() ##NEed to complete functin to match
                for i in range(len(totalRCC)):
                    if(totalRCC[i].endswith(".RCC")):
                        tempRCC = Rcc(totalRCC[i])
                        tempQC = QualityCheck(totalRCC[i],tempRCC, rlfObject)
                        failed = []
                        ##THere is a bug with RNA due to the the control AHRR not been complete
                        passit = True                       
                        for controlNames in acceptedControl:
                            if str(tempRCC.getSampleInfo().loc['ID']['Results']).startswith(controlNames):
                                passit = False
                        if passit:#str(tempRCC.getSampleInfo().loc['ID']['Results']) not in acceptedControl:
                            if "omittedQCValue" in qcParameters.keys() and not tempQC.reqIdChecker():
                                failed.append("Omitted added to name")
                            if not tempQC.geneCount():
                                failed.append("Endogenous genes does not match RLF")
                            if not tempQC.rlfChecker():
                                failed.append("RCC does not match RLF name")
                            if "imagingQCValue" in qcParameters.keys() and not tempQC.fovRatioChecker(qcParameters["imagingQCValue"]):
                                failed.append("FOV ratio failed")
                            if "bindingDensityQCLowerValue" in qcParameters.keys() and not tempQC.bindingDensityChecker(qcParameters["bindingDensityQCLowerValue"], qcParameters["bindingDensityQCUpperValue"]):
                                failed.append("Binding density failed")
                            if "positiveLinearityQCValue" in qcParameters.keys() and not tempQC.positiveControlLinearityChecker(qcParameters["positiveLinearityQCValue"]):
                                failed.append("Positive Linearity failed")
                            if "positiveDetectionQCValue" in qcParameters.keys() and not tempQC.positiveControlLimitDetectionChecker(qcParameters["positiveDetectionQCValue"]):
                                failed.append("Positive detection failed")
                                
                            ##For the pcnFactor range we will take care of it after putting all the files together. ask nate
                            if "backgroundLimitQCValue" in qcParameters.keys() and not tempQC.hkgCountsChecker(qcParameters["backgroundLimitQCValue"],qcParameters["backgroundPercentQCValue"]):
                                ##Nate just wants a statement in the file
                                #print "Going into the background Limit"
                                backgroundFailures[str(tempRCC.getSampleInfo().loc['ID']['Results'])] = "Failed background limit QC"
                            
                            if len(failed)==0:
                                processedFiles.append(totalRCC[i])
                            else:
                                nonProcessedFiles.append(totalRCC[i])
                                qcFailed[str(tempRCC.getSampleInfo().loc['ID']['Results'])] = failed
                        else:
                            #nonProcessedFiles.append(totalRCC[i])
                            controls.append(totalRCC[i])
                            
                            
            self.populateList()
            
            ##Need to add a Log
            writer = open(errorLogLocation, 'w')
            writer.write('ID,Error\n')
            ##I just need a file
            keysError = qcFailed.keys()
            for key in keysError:
                for error in qcFailed[key]:
                    writer.write(key + ',' + error + '\n')
            
            writer.close()
            
            
            for selectedControl in acceptedControl:
                fileToWrite = controlsFile + "_" + selectedControl + ".csv"
                writer = rf.ReadFile()
                writer.printHeader(rlfPath,fileToWrite)
                for files in controls:
                    nameControlList = files.split("_")
                    nameControl = nameControlList[len(nameControlList)-2]
                    if nameControl.startswith(selectedControl):
                        writer.printFile(files,fileToWrite)
                        
        def populateList(self):
            for i in range (0,len(processedFiles)):
                self.processRCCList.addItem(processedFiles[i])
			
            for i in range (0,len(nonProcessedFiles)):
                self.nonProcessedRCCList.addItem(nonProcessedFiles[i])
                
            self.processRCCList.sortItems(Qt.AscendingOrder)
            self.nonProcessedRCCList.sortItems(Qt.AscendingOrder)
        
        def controlLists(self):
            acceptedControl = []
            for i in range(self.listWidget_2.count()):
                acceptedControl.append(self.listWidget_2.item(i).text())
            
            return acceptedControl
                            
                    
                
        ##This function gets used to check whether there is an RCC files in the RCC path given
        def isRCC(self,path):
            import os.path
            from os import walk
            rccExists = False
            
            for(dirpaths, dirnames, filenames) in walk(path):
                for f in filenames:
                    fname = os.path.join(dirpaths, f)
                    assert(os.path.exists(fname))
                    if fname.endswith(".rcc") or fname.endswith(".RCC"):
                        rccExists = True
                        return rccExists
            return rccExists
            
        def populateGeneBoxes(self, rlfObject):
            index = rlfObject.getClassInformation()
            
            housekeepingIndex = index.loc[index['Class Name'] == "Housekeeping"].index.tolist()
            endogenousIndex = index.loc[index['Class Name'] == "Endogenous"].index.tolist()
            proteinIndex = index.loc[index['Class Name'] == "Protein"].index.tolist()
            ##Need to add to contain protein as well.
            proteinHkpIndex = index.loc[index['Class Name'] == "Protein_CELL_NORM"].index.tolist()
            
            
            
            
            rlfContent = rlfObject.getRlfContent()
            
            housekeepingGenes = rlfContent.loc[housekeepingIndex[0]]['GeneName'].tolist()
            endogenousGenes = rlfContent.loc[endogenousIndex[0]]['GeneName'].tolist()
            global proteinsAvailable
            proteinsAvailable = False
            if len(proteinIndex) > 0:
                proteinsAvailable = True
                proteins2 = rlfContent.loc[proteinIndex[0]]['GeneName'].tolist() #Returns a series if plural but a string if singular?
                #print proteins2
                proteinsHKP = rlfContent.loc[proteinHkpIndex[0]]['GeneName']#.tolist()
                #print type(proteinsHKP)
                if isinstance(proteinsHKP, basestring):
                    proteinsHKP = proteinsHKP.split("*-*")#.tolist() #I cannot get an object since its only one? #NEED to do a considilation of the code
                else:
                    proteinsHKP = proteinsHKP.tolist()
                #print proteinsHKP
            for gene in housekeepingGenes:
                self.defaultHkgList.addItem(gene)
            self.defaultHkgList.sortItems()
            
            for gene in endogenousGenes:
                self.genesList.addItem(gene)
            self.genesList.sortItems()
                
            if proteinsAvailable:
                self.tab_5.setEnabled(True)
                self.proteinTable.setEnabled(True)
                self.proteinTable.setRowCount(len(proteins2))
                self.proteinTable.setHorizontalHeaderItem(0,QTableWidgetItem("Protein"))
                self.proteinTable.setHorizontalHeaderItem(1,QTableWidgetItem("Fold Change"))
                
                i = 0
                for protein in proteins2:
                    self.proteinList.addItem(protein)
                    self.proteinTable.setItem(i,0,QTableWidgetItem(protein))
                    self.proteinTable.setItem(i,1,QTableWidgetItem("4"))
                    i += 1
                self.proteinList.sortItems()
                
                for protein in proteinsHKP:
                    self.hkpList.addItem(protein)##This can be reduced to steps if i think of it a little harder
                    self.proteinTable.setItem(i,0,QTableWidgetItem(protein))
                    self.proteinTable.setItem(i,1,QTableWidgetItem("4"))
                    i += 1
                self.hkpList.sortItems()
                
        def populateControlList(self, controls, possibleControls, oneRun):
            if not oneRun:
                for i in range (0,len(controls)):
                    self.listWidget_2.addItem(controls[i])
			
                for i in range (0,len(possibleControls)):
                    self.listWidget.addItem(possibleControls[i])
            else:
                self.listWidget_2.addItem("One run registered")
                
        def upHkg(self):
            widgetItem = self.genesList.selectedItems()
            if not widgetItem: return
            for item in widgetItem:
                self.genesList.takeItem(self.genesList.row(item))
            for i in range(len(widgetItem)):
                self.defaultHkgList.addItem(widgetItem[i].text())
            self.defaultHkgList.sortItems(Qt.AscendingOrder)
                
        def downHkg(self):
            widgetItem = self.defaultHkgList.selectedItems()
            if not widgetItem: return
            for item in widgetItem:
                self.defaultHkgList.takeItem(self.defaultHkgList.row(item))
            for i in range(len(widgetItem)):
                self.genesList.addItem(widgetItem[i].text())
            self.genesList.sortItems(Qt.AscendingOrder)
            
        def toEndogenousProtein(self):
            widgetItem = self.hkpList.selectedItems()
            if not widgetItem: return
            for item in widgetItem:
                self.hkpList.takeItem(self.hkpList.row(item))
            for i in range(len(widgetItem)):
                self.proteinList.addItem(widgetItem[i].text())
            self.proteinList.sortItems(Qt.AscendingOrder)
        
        def toHKPProtein(self):
            widgetItem = self.proteinList.selectedItems()
            if not widgetItem:return
            for item in widgetItem:
                self.proteinList.takeItem(self.proteinList.row(item))
            for i in range(len(widgetItem)):
                self.hkpList.addItem(widgetItem[i].text())
            self.hkpList.sortItems(Qt.AscendingOrder)
            
        def upRCC(self):
            ###Need to add a function to make sure that the selected Items that
            ###Try to go to the process list are of the same RLF
            ##rightRLF(self, fileToMove, rlfUsed)
            widgetItem = self.nonProcessedRCCList.selectedItems()
            if not widgetItem: return
            for item in widgetItem:
                self.nonProcessedRCCList.takeItem(self.nonProcessedRCCList.row(item))
            for i in range(len(widgetItem)):
                if self.rightRLF(widgetItem[i].text(), rlfPath):
                    self.processRCCList.addItem(widgetItem[i].text())
                else:
                    noRCCMessage = QMessageBox()
                    noRCCMessage.setText(widgetItem[i].text() + " does not match RLF being tested!")
                    noRCCMessage.exec_() 
            self.processRCCList.sortItems(Qt.AscendingOrder)
            
        def rightRLF(self, fileToMove, rlfUsed):
            ##Take an rlfName that should be decide at the time way early
            ##so that you dont have to call it everytime.
            readFile = open(fileToMove, 'r')
            contents = readFile.read().splitlines()
            readFile.close()
            check = False
            for line in contents:
                if line.startswith("GeneRLF"):
                    tempSplit = line.split(",")
                    break
            rlfName = tempSplit[1] + ".rlf"
        
            if rlfUsed.lower().endswith(rlfName.lower()):
                check = True
                
            return check
        
        def controlList(self):
            acceptedControl = []
            for i in range(self.listWidget_2.count()):
                acceptedControl.append(self.listWidget_2.item(i).text())
            return acceptedControl
            
        def downRCC(self):
            widgetItem = self.processRCCList.selectedItems()
            if not widgetItem: return
            for item in widgetItem:
                self.processRCCList.takeItem(self.processRCCList.row(item))
                
            for i in range(len(widgetItem)):
                self.nonProcessedRCCList.addItem(widgetItem[i].text())
            
            self.nonProcessedRCCList.sortItems(Qt.AscendingOrder)
    		
        def moveToNotControl(self):
            widgetItem = self.listWidget_2.selectedItems()
            if not widgetItem: return
            for item in widgetItem:
                self.listWidget_2.takeItem(self.listWidget_2.row(item))
                
            for i in range(len(widgetItem)):
                self.listWidget.addItem(widgetItem[i].text())
                
            self.listWidget.sortItems(Qt.AscendingOrder)
        
        def moveToControl(self):
            widgetItem = self.listWidget.selectedItems()
            if not widgetItem: return
            for item in widgetItem:
                self.listWidget.takeItem(self.listWidget.row(item))
                
            for i in range(len(widgetItem)):
                self.listWidget_2.addItem(widgetItem[i].text())
                
            self.listWidget_2.sortItems(Qt.AscendingOrder)
            
        def runExperiment(self):
            ##print qcParameters
            rcc = RccProcessing()
            
            
            
            #&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
            
            for i in range(self.processRCCList.count()):
                rcc.addRcc(self.processRCCList.item(i).text())
            
            
            testRNA = RnaAnalysis(rcc, hkgList)
            rawData = testRNA.getRawData()
            testRNA.calcPosGeoMean()
            testRNA.calcGeoMeanPerSample()
            testRNA.calcPcn()
            testRNA.pcnNormalization()
            testRNA.calcLog2()
            testRNA.calcLog2Endogenous()
            testRNA.calcHkgCount()
            testRNA.calcHkgCount2()
            testRNA.calcDeltaCounts()
            testRNA.calcDeltaCounts2()
            #testRNA.calcCcnFactor()
            #testRNA.calcCcnPcnArray()
            
            
            pcnFactors = testRNA.getPCNFactor()
            pcnNormalizedData = testRNA.getPCNNormalized()
            #ccnNormalizedData = testRNA.getCCNNormalized()
            rnaDeltaAnalysis = testRNA.getDeltaCounts2() ##this is the extra one
            rnaDeltaCounts = testRNA.getDeltaCounts()
            keys = testRNA.getKeys()
            
            
            ##############################HEADER#########################
            
            countHeader = ["Count" for i in range(len(keys))]            
            sampleInfo = rcc.getSamplesInfo()
            laneInfo = rcc.getLanesInfo()
            
            #print pcnFactors
            ##Add the PCN parameter to the laneInfo Dataframe
            if qcParameters.has_key("pcnFactorQCLowerValue"):
                        lowBound = float(qcParameters["pcnFactorQCLowerValue"])
                        #print "Low bound: ",lowBound
                        upBound = float(qcParameters["pcnFactorQCUpperValue"])
                        #print "Up bound: ", upBound
                        for i in range(len(keys)):
                            if float(pcnFactors[keys[i]]) > lowBound and float(pcnFactors[keys[i]]) < upBound: 
                                laneInfo[keys[i]].loc['PCN FACTOR VALID'] = "TRUE"
                            else:
                                laneInfo[keys[i]].loc['PCN FACTOR VALID'] = "FALSE"
                                
            ##Add the background failure message
            #print backgroundFailures
            if qcParameters.has_key("backgroundLimitQCValue"):##I still need to fix
                for key in keys:    
                    if backgroundFailures.has_key(key):
                        laneInfo[key].loc['BACKGROUND PASS'] = "FALSE"
                    else:
                        laneInfo[key].loc['BACKGROUND PASS'] = "TRUE"
            
            ##We need to create a row with the date and batch reference
            ##The first part is messed up. I need to do the date ordering somehow.
            dateOrder = []
            dateContainer = {}
            counter = 0
            for key in keys:
                dateOrder.append(sampleInfo[key].loc['Date']['Results'])
            dateOrder.sort()
            for i in range(len(dateOrder)):
                if not dateContainer.has_key(dateOrder[i]):
                    counter += 1
                    dateContainer[dateOrder[i]] = counter
            
            for key in keys:
                sampleInfo[key].loc["Batch"] = dateContainer[sampleInfo[key].loc["Date"]["Results"]]
                
            ##we need to get batch I
            import re
            for key in keys:
                string = str(laneInfo[key].loc['CartridgeID']['Results']).lower()
                batchList = re.findall('batch([0-9 ]*)',string)#str(re.findall('batch([0-9 ]*)',string)[0].strip())
                
                if len(batchList)==0:
                    batch = string#"000"
                else:
                    batch = str(batchList[0]).strip()
                    if len(batch) == 1:
                        batch = "00" + batch
                    elif len(batch) == 2:
                        batch = "0" + batch
                    elif len(batch) == 0:##If there is no batch number
                        batch = "000"##safe mechanism
                    
                date = str(sampleInfo[key].loc["Batch"]["Results"])
                if len(date) == 1:
                    date = "00" + date
                elif len(date) == 2:
                    date = "0" + date
                
                lane = str(laneInfo[key].loc["LANE"]["Results"])
                if len(lane) == 1:
                    lane = "00" + lane
                elif len(lane) == 2:
                    lane = "0" + lane
                    
                sampleInfo[key].loc["Batch"] = date + "-" +  batch + "-" + lane
                
                
            
            ##Make Sample info merger
            if len(keys) > 1:
                sampleInfoMerge = pd.concat([sampleInfo[keys[0]],sampleInfo[keys[1]]], axis = 1)
            
                for i in range(2,len(keys)):
                    sampleInfoMerge = pd.concat([sampleInfoMerge,sampleInfo[keys[i]]],axis=1)
            
            sampleInfoMerge.columns = countHeader
            #sampleInfoMerge.to_string()
            
            ##Make Lane info merger
            if len(keys) > 1:
                laneInfoMerge = pd.concat([laneInfo[keys[0]],laneInfo[keys[1]]], axis = 1)
            
                for i in range(2,len(keys)):
                    laneInfoMerge = pd.concat([laneInfoMerge,laneInfo[keys[i]]],axis=1)
            
            laneInfoMerge.columns = countHeader
            
            ##Concatenate headers
            title = pd.concat([sampleInfoMerge,laneInfoMerge])
            
            
            ################RNA###################################
            
            ##RAW DATA
            cleanRaw = {}
            for key in keys:
                cleanRaw[key] = rawData[key][['Name','Count']].set_index('Name')
            
            if len(keys) >1:
                rawMerger = pd.concat([cleanRaw[keys[0]], cleanRaw[keys[1]]], axis = 1)
                for i in range(2,len(keys)):
                    rawMerger = pd.concat([rawMerger, cleanRaw[keys[i]]], axis = 1)
            #print rawMerger.to_string()
            
            ##PCN
            if len(keys) > 1:##I need to add a warning message if PCN factor is outside the range
                    pcnMerger = pd.concat([pcnNormalizedData[keys[0]],pcnNormalizedData[keys[1]]], axis = 1)
                    for i in range(2,len(keys)):
                        pcnMerger = pd.concat([pcnMerger,pcnNormalizedData[keys[i]]], axis=1)
                        
            '''#CCN
            if len(keys) > 1:
                    ccnMerger = pd.concat([ccnNormalizedData[keys[0]],ccnNormalizedData[keys[1]]], axis = 1)
                    for i in range(2,len(keys)):
                        ccnMerger = pd.concat([ccnMerger,ccnNormalizedData[keys[i]]], axis=1)'''
                        
            ##Delta count
            if len(keys) > 1:
                    deltaMerger = pd.concat([rnaDeltaCounts[keys[0]],rnaDeltaCounts[keys[1]]], axis = 1)
                    for i in range(2,len(keys)):
                        deltaMerger = pd.concat([deltaMerger,rnaDeltaCounts[keys[i]]], axis=1)
                        
            ##Delta analysis
            if len(keys) > 1:
                deltaAnalysisMerger = pd.concat([rnaDeltaAnalysis[keys[0]],rnaDeltaAnalysis[keys[1]]], axis = 1)
                for i in range(2,len(keys)):
                    deltaAnalysisMerger = pd.concat([deltaAnalysisMerger,rnaDeltaAnalysis[keys[i]]], axis=1)
                        
            ##Print
            writer = pd.ExcelWriter(rnaLocation, engine='xlsxwriter')
            
                
            rawMerger = pd.concat([title,rawMerger])
            rawMerger.columns = rawMerger.loc["Batch"]
            rawMerger = rawMerger.sort(axis=1)
            
            rawMerger.to_excel(writer,sheet_name='RAW_DATA')
            
            pcnMerger = pd.concat([title,pcnMerger])
            pcnMerger.columns = pcnMerger.loc["Batch"]
            pcnMerger = pcnMerger.sort(axis=1)
            pcnMerger.to_excel(writer,sheet_name='PCN')   
           
            deltaMerger = pd.concat([title,deltaMerger])
            deltaMerger.columns = deltaMerger.loc["Batch"]
            deltaMerger = deltaMerger.sort(axis=1)
            deltaMerger.to_excel(writer,sheet_name='RNA Delta Counts')
            
            deltaAnalysis = pd.concat([title,deltaAnalysisMerger])
            deltaAnalysis.columns = deltaAnalysis.loc["Batch"]
            deltaAnalysis = deltaAnalysis.sort(axis=1)
            deltaAnalysis.to_excel(writer, sheet_name='RNA Delta Analysis')
            
                
            #writer.close()
            
            
            ################PROTEIN###############################
            if proteinsAvailable:
                
                if not noProteinTest:
                    proteinTest = ProteinAnalysis(rcc)
                    proteinTest.calcCTP()
                    proteinTest.calcCTHKG(proteins)##This is not right
                    proteinTest.calcDeltaCt()
                    proteinTest.calcNegGmean()
                    proteinTest.calcPosNegValue()
                    proteinTest.compareValues()
                    ##Here is where we need to do the integration of the table
                    if self.proteinTable.rowCount() > 0:
                        proteinKey = {}
                        for i in range(self.proteinTable.rowCount()):
                            if float(self.proteinTable.item(i,1).text()) != 4:
                                proteinKey[str(self.proteinTable.item(i,0).text())] = float(self.proteinTable.item(i,1).text())
                        
                        proteinTest.compareValuesSpecific(proteinKey)
                        
                    proteinDeltaCounts = proteinTest.getDeltaValues()
                    proteinCompare = proteinTest.getComparedValues()
    
    
    
                    
                        
                    #Make deltacount protein merger
                    if len(keys) > 1:
                        deltaProteinmerger = pd.concat([proteinDeltaCounts[keys[0]],proteinDeltaCounts[keys[1]]], axis = 1)
                        for i in range(2,len(keys)):
                            deltaProteinmerger = pd.concat([deltaProteinmerger,proteinDeltaCounts[keys[i]]], axis=1)
                    deltaProteinmerger = pd.concat([title,deltaProteinmerger])
                    deltaProteinmerger.columns = deltaProteinmerger.loc["Batch"]
                    deltaProteinmerger = deltaProteinmerger.sort(axis=1)
                    deltaProteinmerger.to_excel(writer,sheet_name='Protein Delta Count')
                    
                    #Make compare protein merger
                    if len(keys) >1:
                        compareMerged = pd.concat([proteinCompare[keys[0]],proteinCompare[keys[1]]], axis=1)
                        for i in range(2,len(keys)):
                            compareMerged = pd.concat([compareMerged,proteinCompare[keys[i]]], axis=1)
                    compareMerged = pd.concat([title,compareMerged])
                    compareMerged.columns = compareMerged.loc["Batch"]
                    compareMerged = compareMerged.sort(axis=1)
                    compareMerged.to_excel(writer,sheet_name='Protein Delta Analysis')
                
                
                
                
                
            writer.close()
                
            
            ##RESTARTING 
            completeBox = QMessageBox()
            completeBox.setText("The process has been completed!")
            completeBox.exec_()
            
            self.processRCCList.clear()
            self.nonProcessedRCCList.clear() 
            self.listWidget.clear()
            self.listWidget_2.clear()
            self.genesList.clear()
            self.defaultHkgList.clear()
            self.RCCText.clear()
            self.rlfTxt.clear()
            self.outputTxt.clear()
            self.proteinList.clear()
            self.hkpList.clear()
            self.confirmSelectionBtn.setEnabled(False)
            self.confirmNanostringQC.setEnabled(False)
            self.confirmGenoptixQC.setEnabled(False)
            self.confirmControls.setEnabled(False)
            self.actionRUN_EXPERIMENT.setEnabled(False)
            self.tab_5.setEnabled(False)
            self.proteinTable.setEnabled(False)
            self.proteinTable.clear()
                
                
if __name__ == "__main__":
        app = QtGui.QApplication(sys.argv)
        window = MyApplication()
        window.show()
        sys.exit(app.exec_())		
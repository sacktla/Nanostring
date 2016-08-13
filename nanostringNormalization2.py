# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'test.ui'
#
# Created: Tue Apr 26 16:48:25 2016
#      by: pyside-uic 0.2.15 running on PySide 1.2.2
#
# WARNING! All changes made in this file will be lost!

from PySide import QtCore, QtGui

class Ui_MainWindow(object):
    def setupUi(self, MainWindow):
        MainWindow.setObjectName("MainWindow")
        MainWindow.resize(794, 504)
        self.centralwidget = QtGui.QWidget(MainWindow)
        self.centralwidget.setObjectName("centralwidget")
        self.horizontalLayout = QtGui.QHBoxLayout(self.centralwidget)
        self.horizontalLayout.setObjectName("horizontalLayout")
        self.proteinTab = QtGui.QTabWidget(self.centralwidget)
        self.proteinTab.setEnabled(True)
        self.proteinTab.setObjectName("proteinTab")
        self.widget = QtGui.QWidget()
        self.widget.setObjectName("widget")
        self.verticalLayout = QtGui.QVBoxLayout(self.widget)
        self.verticalLayout.setObjectName("verticalLayout")
        self.selectRCCBtn = QtGui.QPushButton(self.widget)
        self.selectRCCBtn.setObjectName("selectRCCBtn")
        self.verticalLayout.addWidget(self.selectRCCBtn)
        self.RCCText = QtGui.QLineEdit(self.widget)
        self.RCCText.setReadOnly(True)
        self.RCCText.setObjectName("RCCText")
        self.verticalLayout.addWidget(self.RCCText)
        self.RLFSelectionBtn = QtGui.QPushButton(self.widget)
        self.RLFSelectionBtn.setObjectName("RLFSelectionBtn")
        self.verticalLayout.addWidget(self.RLFSelectionBtn)
        self.rlfTxt = QtGui.QLineEdit(self.widget)
        self.rlfTxt.setReadOnly(True)
        self.rlfTxt.setObjectName("rlfTxt")
        self.verticalLayout.addWidget(self.rlfTxt)
        self.outputSelectionBtn = QtGui.QPushButton(self.widget)
        self.outputSelectionBtn.setObjectName("outputSelectionBtn")
        self.verticalLayout.addWidget(self.outputSelectionBtn)
        self.outputTxt = QtGui.QLineEdit(self.widget)
        self.outputTxt.setText("")
        self.outputTxt.setReadOnly(True)
        self.outputTxt.setObjectName("outputTxt")
        self.verticalLayout.addWidget(self.outputTxt)
        self.confirmSelectionBtn = QtGui.QPushButton(self.widget)
        self.confirmSelectionBtn.setObjectName("confirmSelectionBtn")
        self.verticalLayout.addWidget(self.confirmSelectionBtn)
        self.proteinTab.addTab(self.widget, "")
        self.tab_2 = QtGui.QWidget()
        self.tab_2.setObjectName("tab_2")
        self.fovRegistrationCheckBox = QtGui.QCheckBox(self.tab_2)
        self.fovRegistrationCheckBox.setGeometry(QtCore.QRect(9, 48, 384, 17))
        self.fovRegistrationCheckBox.setChecked(True)
        self.fovRegistrationCheckBox.setTristate(False)
        self.fovRegistrationCheckBox.setObjectName("fovRegistrationCheckBox")
        self.fovRegistrationSB = QtGui.QSpinBox(self.tab_2)
        self.fovRegistrationSB.setGeometry(QtCore.QRect(410, 50, 41, 21))
        self.fovRegistrationSB.setProperty("value", 75)
        self.fovRegistrationSB.setObjectName("fovRegistrationSB")
        self.bindingDensityCheckBox = QtGui.QCheckBox(self.tab_2)
        self.bindingDensityCheckBox.setGeometry(QtCore.QRect(10, 90, 477, 17))
        self.bindingDensityCheckBox.setChecked(True)
        self.bindingDensityCheckBox.setObjectName("bindingDensityCheckBox")
        self.upperBoundLabel = QtGui.QLabel(self.tab_2)
        self.upperBoundLabel.setGeometry(QtCore.QRect(160, 130, 66, 13))
        self.upperBoundLabel.setMaximumSize(QtCore.QSize(16777215, 13))
        self.upperBoundLabel.setObjectName("upperBoundLabel")
        self.upperBoundSpinBox = QtGui.QDoubleSpinBox(self.tab_2)
        self.upperBoundSpinBox.setGeometry(QtCore.QRect(230, 130, 49, 20))
        self.upperBoundSpinBox.setProperty("value", 2.25)
        self.upperBoundSpinBox.setObjectName("upperBoundSpinBox")
        self.lowerBoundLabel = QtGui.QLabel(self.tab_2)
        self.lowerBoundLabel.setGeometry(QtCore.QRect(20, 130, 66, 13))
        self.lowerBoundLabel.setMaximumSize(QtCore.QSize(16777215, 13))
        self.lowerBoundLabel.setObjectName("lowerBoundLabel")
        self.lowerBoundSpinBox = QtGui.QDoubleSpinBox(self.tab_2)
        self.lowerBoundSpinBox.setGeometry(QtCore.QRect(90, 130, 49, 20))
        self.lowerBoundSpinBox.setProperty("value", 0.05)
        self.lowerBoundSpinBox.setObjectName("lowerBoundSpinBox")
        self.positiveControlCheckBox = QtGui.QCheckBox(self.tab_2)
        self.positiveControlCheckBox.setGeometry(QtCore.QRect(10, 170, 474, 17))
        self.positiveControlCheckBox.setChecked(True)
        self.positiveControlCheckBox.setObjectName("positiveControlCheckBox")
        self.posControlSpinBox = QtGui.QDoubleSpinBox(self.tab_2)
        self.posControlSpinBox.setGeometry(QtCore.QRect(500, 170, 49, 20))
        self.posControlSpinBox.setProperty("value", 0.95)
        self.posControlSpinBox.setObjectName("posControlSpinBox")
        self.negControlCheck = QtGui.QCheckBox(self.tab_2)
        self.negControlCheck.setGeometry(QtCore.QRect(10, 210, 752, 17))
        self.negControlCheck.setChecked(True)
        self.negControlCheck.setObjectName("negControlCheck")
        self.negControlSpinBox = QtGui.QSpinBox(self.tab_2)
        self.negControlSpinBox.setGeometry(QtCore.QRect(10, 240, 41, 20))
        self.negControlSpinBox.setProperty("value", 2)
        self.negControlSpinBox.setObjectName("negControlSpinBox")
        self.confirmNanostringQC = QtGui.QPushButton(self.tab_2)
        self.confirmNanostringQC.setGeometry(QtCore.QRect(330, 360, 75, 23))
        self.confirmNanostringQC.setObjectName("confirmNanostringQC")
        self.pcnRange = QtGui.QCheckBox(self.tab_2)
        self.pcnRange.setGeometry(QtCore.QRect(10, 280, 181, 17))
        self.pcnRange.setChecked(True)
        self.pcnRange.setObjectName("pcnRange")
        self.upperBoundLabel_2 = QtGui.QLabel(self.tab_2)
        self.upperBoundLabel_2.setGeometry(QtCore.QRect(160, 310, 66, 13))
        self.upperBoundLabel_2.setMaximumSize(QtCore.QSize(16777215, 13))
        self.upperBoundLabel_2.setObjectName("upperBoundLabel_2")
        self.lowerBoundLabel_2 = QtGui.QLabel(self.tab_2)
        self.lowerBoundLabel_2.setGeometry(QtCore.QRect(20, 310, 66, 13))
        self.lowerBoundLabel_2.setMaximumSize(QtCore.QSize(16777215, 13))
        self.lowerBoundLabel_2.setObjectName("lowerBoundLabel_2")
        self.lowBoundPCN = QtGui.QDoubleSpinBox(self.tab_2)
        self.lowBoundPCN.setGeometry(QtCore.QRect(90, 310, 51, 22))
        self.lowBoundPCN.setProperty("value", 0.5)
        self.lowBoundPCN.setObjectName("lowBoundPCN")
        self.upBoundPCN = QtGui.QDoubleSpinBox(self.tab_2)
        self.upBoundPCN.setGeometry(QtCore.QRect(230, 310, 51, 22))
        self.upBoundPCN.setProperty("value", 2.0)
        self.upBoundPCN.setObjectName("upBoundPCN")
        self.proteinTab.addTab(self.tab_2, "")
        self.tab = QtGui.QWidget()
        self.tab.setObjectName("tab")
        self.confirmGenoptixQC = QtGui.QPushButton(self.tab)
        self.confirmGenoptixQC.setGeometry(QtCore.QRect(340, 390, 75, 23))
        self.confirmGenoptixQC.setObjectName("confirmGenoptixQC")
        self.backgroundTresholdCheck = QtGui.QCheckBox(self.tab)
        self.backgroundTresholdCheck.setGeometry(QtCore.QRect(9, 9, 154, 17))
        self.backgroundTresholdCheck.setChecked(True)
        self.backgroundTresholdCheck.setObjectName("backgroundTresholdCheck")
        self.downHkgBtn = QtGui.QPushButton(self.tab)
        self.downHkgBtn.setGeometry(QtCore.QRect(273, 160, 221, 23))
        self.downHkgBtn.setObjectName("downHkgBtn")
        self.backgroundThresSpinBox = QtGui.QSpinBox(self.tab)
        self.backgroundThresSpinBox.setGeometry(QtCore.QRect(270, 10, 51, 21))
        self.backgroundThresSpinBox.setProperty("value", 65)
        self.backgroundThresSpinBox.setMaximum(1000)
        self.backgroundThresSpinBox.setObjectName("backgroundThresSpinBox")
        self.reqIdCheck = QtGui.QCheckBox(self.tab)
        self.reqIdCheck.setGeometry(QtCore.QRect(510, 10, 171, 17))
        self.reqIdCheck.setChecked(True)
        self.reqIdCheck.setObjectName("reqIdCheck")
        self.genesList = QtGui.QListWidget(self.tab)
        self.genesList.setGeometry(QtCore.QRect(10, 110, 241, 99))
        self.genesList.setSelectionMode(QtGui.QAbstractItemView.MultiSelection)
        self.genesList.setObjectName("genesList")
        self.endogenousGenesLabel = QtGui.QLabel(self.tab)
        self.endogenousGenesLabel.setGeometry(QtCore.QRect(80, 90, 116, 16))
        self.endogenousGenesLabel.setObjectName("endogenousGenesLabel")
        self.defaultHkgList = QtGui.QListWidget(self.tab)
        self.defaultHkgList.setGeometry(QtCore.QRect(520, 110, 241, 99))
        self.defaultHkgList.setSelectionMode(QtGui.QAbstractItemView.MultiSelection)
        self.defaultHkgList.setObjectName("defaultHkgList")
        self.upHkgBtn = QtGui.QPushButton(self.tab)
        self.upHkgBtn.setGeometry(QtCore.QRect(273, 120, 221, 23))
        self.upHkgBtn.setObjectName("upHkgBtn")
        self.hkgLabel = QtGui.QLabel(self.tab)
        self.hkgLabel.setGeometry(QtCore.QRect(590, 90, 111, 16))
        self.hkgLabel.setObjectName("hkgLabel")
        self.minCountLbl = QtGui.QLabel(self.tab)
        self.minCountLbl.setGeometry(QtCore.QRect(210, 10, 71, 16))
        self.minCountLbl.setObjectName("minCountLbl")
        self.percentCutoffLbl = QtGui.QLabel(self.tab)
        self.percentCutoffLbl.setGeometry(QtCore.QRect(210, 40, 51, 16))
        self.percentCutoffLbl.setObjectName("percentCutoffLbl")
        self.percentCutoffSpinBox = QtGui.QSpinBox(self.tab)
        self.percentCutoffSpinBox.setGeometry(QtCore.QRect(270, 40, 51, 21))
        self.percentCutoffSpinBox.setMaximum(100)
        self.percentCutoffSpinBox.setProperty("value", 50)
        self.percentCutoffSpinBox.setObjectName("percentCutoffSpinBox")
        self.hkpLbl = QtGui.QLabel(self.tab)
        self.hkpLbl.setGeometry(QtCore.QRect(580, 240, 131, 16))
        self.hkpLbl.setObjectName("hkpLbl")
        self.hkpList = QtGui.QListWidget(self.tab)
        self.hkpList.setGeometry(QtCore.QRect(520, 260, 241, 99))
        self.hkpList.setSelectionMode(QtGui.QAbstractItemView.MultiSelection)
        self.hkpList.setObjectName("hkpList")
        self.upProteinBtn_2 = QtGui.QPushButton(self.tab)
        self.upProteinBtn_2.setGeometry(QtCore.QRect(260, 270, 251, 23))
        self.upProteinBtn_2.setObjectName("upProteinBtn_2")
        self.proteinList = QtGui.QListWidget(self.tab)
        self.proteinList.setGeometry(QtCore.QRect(10, 260, 241, 99))
        self.proteinList.setSelectionMode(QtGui.QAbstractItemView.MultiSelection)
        self.proteinList.setObjectName("proteinList")
        self.downProteinBtn_2 = QtGui.QPushButton(self.tab)
        self.downProteinBtn_2.setGeometry(QtCore.QRect(260, 310, 251, 23))
        self.downProteinBtn_2.setObjectName("downProteinBtn_2")
        self.endogenousGenesLabel_2 = QtGui.QLabel(self.tab)
        self.endogenousGenesLabel_2.setGeometry(QtCore.QRect(70, 240, 116, 16))
        self.endogenousGenesLabel_2.setObjectName("endogenousGenesLabel_2")
        self.label_6 = QtGui.QLabel(self.tab)
        self.label_6.setGeometry(QtCore.QRect(30, 40, 141, 21))
        self.label_6.setObjectName("label_6")
        self.proteinTab.addTab(self.tab, "")
        self.tab_5 = QtGui.QWidget()
        self.tab_5.setObjectName("tab_5")
        self.proteinTable = QtGui.QTableWidget(self.tab_5)
        self.proteinTable.setEnabled(False)
        self.proteinTable.setGeometry(QtCore.QRect(140, 60, 431, 331))
        self.proteinTable.setRowCount(1000)
        self.proteinTable.setColumnCount(2)
        self.proteinTable.setObjectName("proteinTable")
        self.proteinTable.setColumnCount(2)
        self.proteinTable.setRowCount(1000)
        self.proteinTable.horizontalHeader().setCascadingSectionResizes(True)
        self.proteinTable.horizontalHeader().setDefaultSectionSize(200)
        self.proteinTable.horizontalHeader().setStretchLastSection(True)
        self.proteinLabl = QtGui.QLabel(self.tab_5)
        self.proteinLabl.setGeometry(QtCore.QRect(100, 20, 531, 20))
        self.proteinLabl.setObjectName("proteinLabl")
        self.proteinTab.addTab(self.tab_5, "")
        self.tab_4 = QtGui.QWidget()
        self.tab_4.setObjectName("tab_4")
        self.listWidget = QtGui.QListWidget(self.tab_4)
        self.listWidget.setGeometry(QtCore.QRect(50, 110, 231, 192))
        self.listWidget.setSelectionMode(QtGui.QAbstractItemView.MultiSelection)
        self.listWidget.setObjectName("listWidget")
        self.listWidget_2 = QtGui.QListWidget(self.tab_4)
        self.listWidget_2.setGeometry(QtCore.QRect(490, 110, 231, 192))
        self.listWidget_2.setSelectionMode(QtGui.QAbstractItemView.MultiSelection)
        self.listWidget_2.setObjectName("listWidget_2")
        self.pushButton = QtGui.QPushButton(self.tab_4)
        self.pushButton.setGeometry(QtCore.QRect(310, 140, 151, 23))
        self.pushButton.setObjectName("pushButton")
        self.pushButton_2 = QtGui.QPushButton(self.tab_4)
        self.pushButton_2.setGeometry(QtCore.QRect(310, 220, 151, 23))
        self.pushButton_2.setObjectName("pushButton_2")
        self.label = QtGui.QLabel(self.tab_4)
        self.label.setGeometry(QtCore.QRect(100, 40, 571, 20))
        self.label.setObjectName("label")
        self.confirmControls = QtGui.QPushButton(self.tab_4)
        self.confirmControls.setGeometry(QtCore.QRect(290, 370, 191, 23))
        self.confirmControls.setObjectName("confirmControls")
        self.label_2 = QtGui.QLabel(self.tab_4)
        self.label_2.setGeometry(QtCore.QRect(130, 80, 46, 13))
        self.label_2.setText("")
        self.label_2.setObjectName("label_2")
        self.label_3 = QtGui.QLabel(self.tab_4)
        self.label_3.setGeometry(QtCore.QRect(110, 80, 121, 20))
        self.label_3.setObjectName("label_3")
        self.label_4 = QtGui.QLabel(self.tab_4)
        self.label_4.setGeometry(QtCore.QRect(550, 80, 121, 20))
        self.label_4.setObjectName("label_4")
        self.proteinTab.addTab(self.tab_4, "")
        self.tab_3 = QtGui.QWidget()
        self.tab_3.setObjectName("tab_3")
        self.verticalLayout_2 = QtGui.QVBoxLayout(self.tab_3)
        self.verticalLayout_2.setObjectName("verticalLayout_2")
        self.analyzedLbl = QtGui.QLabel(self.tab_3)
        self.analyzedLbl.setObjectName("analyzedLbl")
        self.verticalLayout_2.addWidget(self.analyzedLbl)
        self.processRCCList = QtGui.QListWidget(self.tab_3)
        self.processRCCList.setSelectionMode(QtGui.QAbstractItemView.MultiSelection)
        self.processRCCList.setObjectName("processRCCList")
        self.verticalLayout_2.addWidget(self.processRCCList)
        self.upRCCBtn = QtGui.QPushButton(self.tab_3)
        self.upRCCBtn.setObjectName("upRCCBtn")
        self.verticalLayout_2.addWidget(self.upRCCBtn)
        self.downRCCBtn = QtGui.QPushButton(self.tab_3)
        self.downRCCBtn.setObjectName("downRCCBtn")
        self.verticalLayout_2.addWidget(self.downRCCBtn)
        self.label_5 = QtGui.QLabel(self.tab_3)
        self.label_5.setObjectName("label_5")
        self.verticalLayout_2.addWidget(self.label_5)
        self.nonProcessedRCCList = QtGui.QListWidget(self.tab_3)
        self.nonProcessedRCCList.setSelectionMode(QtGui.QAbstractItemView.MultiSelection)
        self.nonProcessedRCCList.setObjectName("nonProcessedRCCList")
        self.verticalLayout_2.addWidget(self.nonProcessedRCCList)
        self.proteinTab.addTab(self.tab_3, "")
        self.horizontalLayout.addWidget(self.proteinTab)
        MainWindow.setCentralWidget(self.centralwidget)
        self.menubar = QtGui.QMenuBar(MainWindow)
        self.menubar.setGeometry(QtCore.QRect(0, 0, 794, 21))
        self.menubar.setObjectName("menubar")
        self.menuMENU = QtGui.QMenu(self.menubar)
        self.menuMENU.setObjectName("menuMENU")
        self.menuABOUT_ME = QtGui.QMenu(self.menubar)
        self.menuABOUT_ME.setObjectName("menuABOUT_ME")
        MainWindow.setMenuBar(self.menubar)
        self.statusbar = QtGui.QStatusBar(MainWindow)
        self.statusbar.setObjectName("statusbar")
        MainWindow.setStatusBar(self.statusbar)
        self.actionABOUT_ME = QtGui.QAction(MainWindow)
        self.actionABOUT_ME.setObjectName("actionABOUT_ME")
        self.actionNANOSTRING_QC = QtGui.QAction(MainWindow)
        self.actionNANOSTRING_QC.setObjectName("actionNANOSTRING_QC")
        self.actionFILE_SELECTION = QtGui.QAction(MainWindow)
        self.actionFILE_SELECTION.setObjectName("actionFILE_SELECTION")
        self.actionCONFIRM = QtGui.QAction(MainWindow)
        self.actionCONFIRM.setObjectName("actionCONFIRM")
        self.actionGENOPTIX_PARAMETERS = QtGui.QAction(MainWindow)
        self.actionGENOPTIX_PARAMETERS.setObjectName("actionGENOPTIX_PARAMETERS")
        self.actionNANOSTRING_PARAMETERS = QtGui.QAction(MainWindow)
        self.actionNANOSTRING_PARAMETERS.setObjectName("actionNANOSTRING_PARAMETERS")
        self.actionQUIT = QtGui.QAction(MainWindow)
        self.actionQUIT.setObjectName("actionQUIT")
        self.actionRUN_EXPERIMENT = QtGui.QAction(MainWindow)
        self.actionRUN_EXPERIMENT.setObjectName("actionRUN_EXPERIMENT")
        self.menuMENU.addAction(self.actionRUN_EXPERIMENT)
        self.menuABOUT_ME.addAction(self.actionABOUT_ME)
        self.menubar.addAction(self.menuMENU.menuAction())
        self.menubar.addAction(self.menuABOUT_ME.menuAction())

        self.retranslateUi(MainWindow)
        self.proteinTab.setCurrentIndex(0)
        QtCore.QMetaObject.connectSlotsByName(MainWindow)

    def retranslateUi(self, MainWindow):
        MainWindow.setWindowTitle(QtGui.QApplication.translate("MainWindow", "JZ SOLVER", None, QtGui.QApplication.UnicodeUTF8))
        self.selectRCCBtn.setText(QtGui.QApplication.translate("MainWindow", "SELECT RCC", None, QtGui.QApplication.UnicodeUTF8))
        self.RLFSelectionBtn.setText(QtGui.QApplication.translate("MainWindow", "SELECT RLF", None, QtGui.QApplication.UnicodeUTF8))
        self.outputSelectionBtn.setText(QtGui.QApplication.translate("MainWindow", "SELECT OUTPUT LOCATION", None, QtGui.QApplication.UnicodeUTF8))
        self.confirmSelectionBtn.setText(QtGui.QApplication.translate("MainWindow", "CONFIRM", None, QtGui.QApplication.UnicodeUTF8))
        self.proteinTab.setTabText(self.proteinTab.indexOf(self.widget), QtGui.QApplication.translate("MainWindow", "FILE_SELECTION", None, QtGui.QApplication.UnicodeUTF8))
        self.fovRegistrationCheckBox.setText(QtGui.QApplication.translate("MainWindow", "Imaging QC: Do not include RCC when percent FOV registration is less than", None, QtGui.QApplication.UnicodeUTF8))
        self.bindingDensityCheckBox.setText(QtGui.QApplication.translate("MainWindow", "Binding Density QC: Do not include RCC when binding density is outside of the following range:", None, QtGui.QApplication.UnicodeUTF8))
        self.upperBoundLabel.setText(QtGui.QApplication.translate("MainWindow", "Upper Bound:", None, QtGui.QApplication.UnicodeUTF8))
        self.lowerBoundLabel.setText(QtGui.QApplication.translate("MainWindow", "Lower Bound:", None, QtGui.QApplication.UnicodeUTF8))
        self.positiveControlCheckBox.setText(QtGui.QApplication.translate("MainWindow", "Positive Control Linearity QC: Do not include RCC when Positive Control R2 value is less than: ", None, QtGui.QApplication.UnicodeUTF8))
        self.negControlCheck.setText(QtGui.QApplication.translate("MainWindow", "Positive Control Limit of Detection QC: Do not include RCC if 0.5fM positive control is less than the following standard deviations of the negative control:", None, QtGui.QApplication.UnicodeUTF8))
        self.confirmNanostringQC.setText(QtGui.QApplication.translate("MainWindow", "Confirm QC", None, QtGui.QApplication.UnicodeUTF8))
        self.pcnRange.setText(QtGui.QApplication.translate("MainWindow", "Determine PCN factor range", None, QtGui.QApplication.UnicodeUTF8))
        self.upperBoundLabel_2.setText(QtGui.QApplication.translate("MainWindow", "Upper Bound:", None, QtGui.QApplication.UnicodeUTF8))
        self.lowerBoundLabel_2.setText(QtGui.QApplication.translate("MainWindow", "Lower Bound:", None, QtGui.QApplication.UnicodeUTF8))
        self.proteinTab.setTabText(self.proteinTab.indexOf(self.tab_2), QtGui.QApplication.translate("MainWindow", "NANOSTRING QC METRICS", None, QtGui.QApplication.UnicodeUTF8))
        self.confirmGenoptixQC.setText(QtGui.QApplication.translate("MainWindow", "Confirm QC", None, QtGui.QApplication.UnicodeUTF8))
        self.backgroundTresholdCheck.setText(QtGui.QApplication.translate("MainWindow", "BACKGROUND THRESHOLD", None, QtGui.QApplication.UnicodeUTF8))
        self.downHkgBtn.setText(QtGui.QApplication.translate("MainWindow", "REMOVE GENES AS HOUSEKEEPING GENES", None, QtGui.QApplication.UnicodeUTF8))
        self.reqIdCheck.setText(QtGui.QApplication.translate("MainWindow", "EXCLUDE OMITTED SAMPLES", None, QtGui.QApplication.UnicodeUTF8))
        self.endogenousGenesLabel.setText(QtGui.QApplication.translate("MainWindow", "ENDOGENOUS GENES", None, QtGui.QApplication.UnicodeUTF8))
        self.upHkgBtn.setText(QtGui.QApplication.translate("MainWindow", "ADD GENES AS HOUSEKEEPING GENES", None, QtGui.QApplication.UnicodeUTF8))
        self.hkgLabel.setText(QtGui.QApplication.translate("MainWindow", "HOUSEKEEPING GENES", None, QtGui.QApplication.UnicodeUTF8))
        self.minCountLbl.setText(QtGui.QApplication.translate("MainWindow", "Min Count", None, QtGui.QApplication.UnicodeUTF8))
        self.percentCutoffLbl.setText(QtGui.QApplication.translate("MainWindow", "% Cutoff", None, QtGui.QApplication.UnicodeUTF8))
        self.hkpLbl.setText(QtGui.QApplication.translate("MainWindow", "HOUSEKEEPING PROTEINS", None, QtGui.QApplication.UnicodeUTF8))
        self.upProteinBtn_2.setText(QtGui.QApplication.translate("MainWindow", "ADD PROTEINS AS HOUSEKEEPING PROTEINS", None, QtGui.QApplication.UnicodeUTF8))
        self.downProteinBtn_2.setText(QtGui.QApplication.translate("MainWindow", "REMOVE PROTEINS AS HOUSEKEEPING PROTEINS", None, QtGui.QApplication.UnicodeUTF8))
        self.endogenousGenesLabel_2.setText(QtGui.QApplication.translate("MainWindow", "ENDOGENOUS PROTEINS", None, QtGui.QApplication.UnicodeUTF8))
        self.label_6.setText(QtGui.QApplication.translate("MainWindow", "DETECTION RNA HKGS", None, QtGui.QApplication.UnicodeUTF8))
        self.proteinTab.setTabText(self.proteinTab.indexOf(self.tab), QtGui.QApplication.translate("MainWindow", "GENOPTIX QC", None, QtGui.QApplication.UnicodeUTF8))
        self.proteinLabl.setText(QtGui.QApplication.translate("MainWindow", "IF YOU WOULD LIKE TO MODIFY THE DEFAULT BACKGROUND LIMIT PLEASE USE THE TABLE BELOW:", None, QtGui.QApplication.UnicodeUTF8))
        self.proteinTab.setTabText(self.proteinTab.indexOf(self.tab_5), QtGui.QApplication.translate("MainWindow", "PROTEIN SNR QC", None, QtGui.QApplication.UnicodeUTF8))
        self.pushButton.setText(QtGui.QApplication.translate("MainWindow", "ADD AS CONTROL", None, QtGui.QApplication.UnicodeUTF8))
        self.pushButton_2.setText(QtGui.QApplication.translate("MainWindow", "REMOVE AS CONTROL", None, QtGui.QApplication.UnicodeUTF8))
        self.label.setText(QtGui.QApplication.translate("MainWindow", "WE HAVE DETERMINED BASED ON QUANTITY WHAT YOUR CONTROLS SHOULD BE. PLEASE MODIFY IF NECESSARY", None, QtGui.QApplication.UnicodeUTF8))
        self.confirmControls.setText(QtGui.QApplication.translate("MainWindow", "CONFIRM", None, QtGui.QApplication.UnicodeUTF8))
        self.label_3.setText(QtGui.QApplication.translate("MainWindow", "POTENTIAL CONTROLS", None, QtGui.QApplication.UnicodeUTF8))
        self.label_4.setText(QtGui.QApplication.translate("MainWindow", "SELECTED CONTROLS", None, QtGui.QApplication.UnicodeUTF8))
        self.proteinTab.setTabText(self.proteinTab.indexOf(self.tab_4), QtGui.QApplication.translate("MainWindow", "RUN CONTROLS", None, QtGui.QApplication.UnicodeUTF8))
        self.analyzedLbl.setText(QtGui.QApplication.translate("MainWindow", "SAMPLES TO BE ANALYZED:", None, QtGui.QApplication.UnicodeUTF8))
        self.upRCCBtn.setText(QtGui.QApplication.translate("MainWindow", "ADD TO ANALYSIS", None, QtGui.QApplication.UnicodeUTF8))
        self.downRCCBtn.setText(QtGui.QApplication.translate("MainWindow", "REMOVE FROM ANALYSIS", None, QtGui.QApplication.UnicodeUTF8))
        self.label_5.setText(QtGui.QApplication.translate("MainWindow", "SAMPLES NOT INCLUDED IN ANALYSIS:", None, QtGui.QApplication.UnicodeUTF8))
        self.proteinTab.setTabText(self.proteinTab.indexOf(self.tab_3), QtGui.QApplication.translate("MainWindow", "CONFIRM RCC", None, QtGui.QApplication.UnicodeUTF8))
        self.menuMENU.setTitle(QtGui.QApplication.translate("MainWindow", "MENU", None, QtGui.QApplication.UnicodeUTF8))
        self.menuABOUT_ME.setTitle(QtGui.QApplication.translate("MainWindow", "HELP", None, QtGui.QApplication.UnicodeUTF8))
        self.actionABOUT_ME.setText(QtGui.QApplication.translate("MainWindow", "ABOUT ME", None, QtGui.QApplication.UnicodeUTF8))
        self.actionNANOSTRING_QC.setText(QtGui.QApplication.translate("MainWindow", "NANOSTRING QC", None, QtGui.QApplication.UnicodeUTF8))
        self.actionFILE_SELECTION.setText(QtGui.QApplication.translate("MainWindow", "FILE SELECTION", None, QtGui.QApplication.UnicodeUTF8))
        self.actionCONFIRM.setText(QtGui.QApplication.translate("MainWindow", "CONFIRM", None, QtGui.QApplication.UnicodeUTF8))
        self.actionGENOPTIX_PARAMETERS.setText(QtGui.QApplication.translate("MainWindow", "GENOPTIX PARAMETERS", None, QtGui.QApplication.UnicodeUTF8))
        self.actionNANOSTRING_PARAMETERS.setText(QtGui.QApplication.translate("MainWindow", "NANOSTRING PARAMETERS", None, QtGui.QApplication.UnicodeUTF8))
        self.actionQUIT.setText(QtGui.QApplication.translate("MainWindow", "QUIT", None, QtGui.QApplication.UnicodeUTF8))
        self.actionRUN_EXPERIMENT.setText(QtGui.QApplication.translate("MainWindow", "RUN EXPERIMENT", None, QtGui.QApplication.UnicodeUTF8))

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 15 14:25:06 2019

@author: earmla
"""
from PyQt5 import QtCore, QtWidgets, QtGui
import sys, os
from fnmatch import fnmatch
#import subprocess as subp
import pickle
import copy
import requests
import numpy as np
from lxml import html
import wget

class var():
    outDir = None

class ClickLabel(QtWidgets.QLabel):
    clicked = QtCore.pyqtSignal()
    imagepath = None
    def mousePressEvent(self, event):
        self.clicked.emit()
        QtWidgets.QLabel.mousePressEvent(self, event)

class Ui_dialogChooseFrames(object):
    #this way I do a class variable, i.e. shared through all classes.. maybe should be private one?
    frameList = None
    frameListAll = None
        
    def runChecker(self):
        #print(str(self.listFrames.selectedItems()[0].text()))
        track = self.listTracks.selectedItems()[0].text()
        frame = self.listFrames.selectedItems()[0].text()
        print('Downloading images, wait a bit')
        self.labeldown.setVisible(True)
        self.buttonBox.setDisabled(True)
        try:
            downloadIfgs(track, frame)
            downloaded = True
        except:
            downloaded = False
        if downloaded:
            self.DialogCheck = QtWidgets.QDialog()
            self.UiCheckIfg = Ui_DialogCheckIfg()
            self.UiCheckIfg.setupUi(self.DialogCheck)
            self.UiCheckIfg.selectedFrame = frame
            self.UiCheckIfg.selectedTrack = track
            self.UiCheckIfg.ifgs = getIfgList(frame)
            self.UiCheckIfg.ifgerrors, self.UiCheckIfg.posIfg = getIfgErrors(frame)
            #self.UiCheckIfg.posIfg = 0
            self.UiCheckIfg.retranslateUi(self.DialogCheck)
            self.UiCheckIfg.drawIfgs(self.DialogCheck)
            #DialogCheck.selectedFrame=selectedFrame
            self.checkOnly.setChecked(False)
            self.DialogCheck.show()
            self.DialogCheck.raise_()
        else:
            error='some error happened (maybe the frame has no data ready?)'
            error_dialog = QtWidgets.QErrorMessage()
            error_dialog.showMessage(str(error))
            error_dialog.exec_()
        self.buttonBox.setEnabled(True)
        self.labeldown.setVisible(False)
    
    def trackSelected(self, item):
        #once user selects a track, he should get a list of available frames
        track = int(str(item.text()))
        framelist = self.frameList[track]
        #QtWidgets.QListView.clear .QListWidget.clear(self.listFrames)
        self.listFrames.clear()
        for item in framelist:
            self.listFrames.addItem(item)
        self.listFrames.takeItem(0)
        #self.listFrames.addItem()
    
    def checkonly(self, clicked):
        if self.checkOnly.isChecked():
            self.frameList = getUncheckedFrames(self.frameListAll)
        else:
            self.frameList = self.frameListAll
        if len(self.listTracks.selectedItems()) > 0:
            self.trackSelected(self.listTracks.selectedItems()[0])

    def setupUi(self, dialogChooseFrames):
        dialogChooseFrames.setObjectName("dialogChooseFrames")
        dialogChooseFrames.resize(400, 300)
        self.buttonBox = QtWidgets.QDialogButtonBox(dialogChooseFrames)
        self.buttonBox.setGeometry(QtCore.QRect(30, 240, 341, 32))
        self.buttonBox.setOrientation(QtCore.Qt.Horizontal)
        self.buttonBox.setStandardButtons(QtWidgets.QDialogButtonBox.Cancel|QtWidgets.QDialogButtonBox.Ok)
        self.buttonBox.setObjectName("buttonBox")
        self.labeldown = QtWidgets.QLabel(dialogChooseFrames)
        self.labeldown.setGeometry(QtCore.QRect(10, 250, 221, 16))
        self.labeldown.setObjectName("labelDown") 
        #self.labeldown.setVisible(False)
        self.labeldown.setVisible(True)
        self.label = QtWidgets.QLabel(dialogChooseFrames)
        self.label.setGeometry(QtCore.QRect(10, 20, 221, 16))
        self.label.setObjectName("label")
        self.checkOnly = QtWidgets.QCheckBox(dialogChooseFrames)
        self.checkOnly.setGeometry(QtCore.QRect(260, 20, 131, 21))
        self.checkOnly.setObjectName("checkOnly")
        self.checkOnly.setChecked(True)
        self.checkOnly.clicked.connect(self.checkonly)
        self.listTracks = QtWidgets.QListWidget(dialogChooseFrames)
        self.listTracks.setGeometry(QtCore.QRect(10, 40, 101, 192))
        self.listTracks.setObjectName("listTracks")
        for i in range(1,176):
            self.listTracks.addItem(str(i))
        self.listTracks.itemClicked.connect(self.trackSelected)

        self.listFrames = QtWidgets.QListWidget(dialogChooseFrames)
        self.listFrames.setGeometry(QtCore.QRect(130, 40, 256, 192))
        self.listFrames.setObjectName("listFrames")
        #self.listFrames.clear

        self.retranslateUi(dialogChooseFrames)
        #self.buttonBox.accepted.connect(dialogChooseFrames.accept)
        #self.buttonBox.accepted.connect(self.dummy)
        self.buttonBox.accepted.connect(self.runChecker)
        self.buttonBox.rejected.connect(dialogChooseFrames.reject)
        QtCore.QMetaObject.connectSlotsByName(dialogChooseFrames)

    def retranslateUi(self, dialogChooseFrames):
        _translate = QtCore.QCoreApplication.translate
        dialogChooseFrames.setWindowTitle(_translate("dialogChooseFrames", "LiCS data checker"))
        self.label.setText(_translate("dialogChooseFrames", "Please select track and frame"))
        #self.labeldown.setText(_translate("dialogChooseFrames", "Please wait, downloading..."))
        self.labeldown.setText(_translate("dialogChooseFrames", "wait for download after clicking"))
        self.checkOnly.setText(_translate("dialogChooseFrames", "exclude checked"))

class Ui_DialogCheckIfg(object):
    selectedFrame = None
    selectedTrack = None
    curIfg = ''
    posIfg = 0
    totalIfgs = 0
    curIfgErr = 0
    ifgs = ['']
    ifgerrors = ['']
    #def closeme(self):
    #    self.
    def setupUi(self, Dialog):
        Dialog.setObjectName("Dialog")
        #Dialog.resize(717, 590)
        Dialog.resize(790, 840)
        Dialog.setFocusPolicy(QtCore.Qt.StrongFocus)
        #sizePolicy = QSizePolicy(QSizePolicy.Preferred, QSizePolicy.Preferred)
        #sizePolicy.setHeightForWidth(True)
        Dialog.setSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed)
        self.graphicsViewIfg = ClickLabel(Dialog)
        #self.graphicsViewIfg = QtWidgets.QLabel(Dialog)
        #self.graphicsViewIfg.setGeometry(QtCore.QRect(10, 60, 331, 251))
        self.graphicsViewIfg.setGeometry(QtCore.QRect(10, 60, 460, 370))
        self.graphicsViewIfg.setObjectName("graphicsViewIfg")
        self.graphicsViewIfg.setScaledContents(True)
        self.graphicsViewIfg.clicked.connect(self.viewFullImageIfg)
        
        self.graphicsViewUnw = ClickLabel(Dialog)
        #self.graphicsViewUnw = QtWidgets.QLabel(Dialog)
        #self.graphicsViewUnw.setGeometry(QtCore.QRect(10, 330, 331, 251))
        self.graphicsViewUnw.setGeometry(QtCore.QRect(10, 460, 460, 370))
        self.graphicsViewUnw.setObjectName("graphicsViewUnw")
        self.graphicsViewUnw.setScaledContents(True)
        self.graphicsViewUnw.clicked.connect(self.viewFullImageUnw)
        
        self.label = QtWidgets.QLabel(Dialog)
        self.label.setGeometry(QtCore.QRect(10, 40, 81, 16))
        self.label.setObjectName("label")
        self.label_2 = QtWidgets.QLabel(Dialog)
        self.label_2.setGeometry(QtCore.QRect(10, 440, 101, 16))
        self.label_2.setObjectName("label_2")
        #self.btnCancel = QtWidgets.QPushButton(Dialog)
        #self.btnCancel.setGeometry(QtCore.QRect(610, 540, 80, 23))
        #self.btnCancel.setObjectName("btnCancel")
        #self.btnCancel.clicked.connect(self.close())
        #self.btnCancel.clicked.connect(self.reject())
        self.formLayoutWidget = QtWidgets.QWidget(Dialog)
        self.formLayoutWidget.setGeometry(QtCore.QRect(490, 40, 251, 271))
        self.formLayoutWidget.setObjectName("formLayoutWidget")
        self.formLayout = QtWidgets.QFormLayout(self.formLayoutWidget)
        self.formLayout.setContentsMargins(0, 0, 0, 0)
        self.formLayout.setObjectName("formLayout")
        self.btnError1 = QtWidgets.QRadioButton(self.formLayoutWidget)
        self.btnError1.setObjectName("btnError1")
        self.btnError1.toggled.connect(self.setIfgError)
        self.formLayout.setWidget(2, QtWidgets.QFormLayout.LabelRole, self.btnError1)
        self.btnError2 = QtWidgets.QRadioButton(self.formLayoutWidget)
        self.btnError2.setObjectName("btnError2")
        self.btnError2.toggled.connect(self.setIfgError)
        self.formLayout.setWidget(3, QtWidgets.QFormLayout.LabelRole, self.btnError2)
        self.label_3 = QtWidgets.QLabel(self.formLayoutWidget)
        self.label_3.setObjectName("label_3")
        self.formLayout.setWidget(0, QtWidgets.QFormLayout.LabelRole, self.label_3)
        self.btnError3 = QtWidgets.QRadioButton(self.formLayoutWidget)
        self.btnError3.setObjectName("btnError3")
        self.btnError3.toggled.connect(self.setIfgError)
        self.formLayout.setWidget(4, QtWidgets.QFormLayout.LabelRole, self.btnError3)
        self.btnError7 = QtWidgets.QRadioButton(self.formLayoutWidget)
        self.btnError7.setObjectName("btnError7")
        self.btnError7.toggled.connect(self.setIfgError)
        self.formLayout.setWidget(5, QtWidgets.QFormLayout.LabelRole, self.btnError7)
        self.btnError4 = QtWidgets.QRadioButton(self.formLayoutWidget)
        self.btnError4.setObjectName("btnError4")
        self.btnError4.toggled.connect(self.setIfgError)
        self.formLayout.setWidget(6, QtWidgets.QFormLayout.LabelRole, self.btnError4)
        self.btnError5 = QtWidgets.QRadioButton(self.formLayoutWidget)
        self.btnError5.setObjectName("btnError5")
        self.btnError5.toggled.connect(self.setIfgError)
        self.formLayout.setWidget(7, QtWidgets.QFormLayout.LabelRole, self.btnError5)
        self.btnError6 = QtWidgets.QRadioButton(self.formLayoutWidget)
        self.btnError6.setObjectName("btnError6")
        self.btnError6.toggled.connect(self.setIfgError)
        self.formLayout.setWidget(8, QtWidgets.QFormLayout.LabelRole, self.btnError6)
        self.rbNoError = QtWidgets.QRadioButton(self.formLayoutWidget)
        self.rbNoError.setObjectName("rbNoError")
        self.rbNoError.setChecked(True)
        self.rbNoError.toggled.connect(self.setIfgError)
        self.formLayout.setWidget(1, QtWidgets.QFormLayout.LabelRole, self.rbNoError)
        self.btnNextImage = QtWidgets.QCommandLinkButton(Dialog)
        self.btnNextImage.setGeometry(QtCore.QRect(490, 320, 178, 41))
        self.btnNextImage.setObjectName("btnNextImage")
        self.btnPrevImage = QtWidgets.QPushButton(Dialog)
        self.btnPrevImage.setGeometry(QtCore.QRect(690, 320, 80, 23))
        self.btnPrevImage.setObjectName("btnPrevImage")
        #self.btnNextImage.toggled.connect(self.showNextIfg)
        self.btnNextImage.clicked.connect(self.showNextIfg)
        self.btnPrevImage.clicked.connect(self.showPrevIfg)
        shortcutleft = QtWidgets.QShortcut(QtGui.QKeySequence.MoveToPreviousChar, self.btnPrevImage)
        shortcutright = QtWidgets.QShortcut(QtGui.QKeySequence.MoveToNextChar, self.btnNextImage)
        shortcutleft.activated.connect(self.showPrevIfg)
        shortcutright.activated.connect(self.showNextIfg)
        
        shortcut0 = QtWidgets.QShortcut(QtGui.QKeySequence('0'), self.rbNoError)
        shortcut0.activated.connect(self.rbNoError.toggle)
        shortcut1 = QtWidgets.QShortcut(QtGui.QKeySequence('1'), self.btnError1)
        shortcut1.activated.connect(self.btnError1.toggle)
        shortcut2 = QtWidgets.QShortcut(QtGui.QKeySequence('2'), self.btnError2)
        shortcut2.activated.connect(self.btnError2.toggle)
        shortcut3 = QtWidgets.QShortcut(QtGui.QKeySequence('3'), self.btnError3)
        shortcut3.activated.connect(self.btnError3.toggle)
        shortcut4 = QtWidgets.QShortcut(QtGui.QKeySequence('5'), self.btnError4)
        shortcut4.activated.connect(self.btnError4.toggle)
        shortcut5 = QtWidgets.QShortcut(QtGui.QKeySequence('6'), self.btnError5)
        shortcut5.activated.connect(self.btnError5.toggle)
        shortcut6 = QtWidgets.QShortcut(QtGui.QKeySequence('7'), self.btnError6)
        shortcut6.activated.connect(self.btnError6.toggle)
        shortcut7 = QtWidgets.QShortcut(QtGui.QKeySequence('4'), self.btnError7)
        shortcut7.activated.connect(self.btnError7.toggle)

        shortcutpgup = QtWidgets.QShortcut(QtGui.QKeySequence.MoveToPreviousPage, self.label)
        shortcutpgdwn = QtWidgets.QShortcut(QtGui.QKeySequence.MoveToNextPage, self.label)
        shortcutpgup.activated.connect(self.showNextIfgMore)
        shortcutpgdwn.activated.connect(self.showPrevIfgMore)

        shortcutshpgup = QtWidgets.QShortcut(QtGui.QKeySequence.SelectPreviousPage, self.label)
        shortcutshpgdwn = QtWidgets.QShortcut(QtGui.QKeySequence.SelectNextPage, self.label)
        shortcutshpgup.activated.connect(self.showNextIfgMore100)
        shortcutshpgdwn.activated.connect(self.showPrevIfgMore100)
        
        self.labelCheckingIfg = QtWidgets.QLabel(Dialog)
        self.labelCheckingIfg.setGeometry(QtCore.QRect(10, 10, 581, 16))
        self.labelCheckingIfg.setObjectName("labelCheckingIfg")
        self.labelIfgNo = QtWidgets.QLabel(Dialog)
        self.labelIfgNo.setGeometry(QtCore.QRect(490, 370, 201, 16))
        self.labelIfgNo.setObjectName("labelIfgNo")
        self.labelHelp = QtWidgets.QLabel(Dialog)
        self.labelHelp.setGeometry(QtCore.QRect(490, 395, 281, 16))
        self.labelHelp.setObjectName("labelHelp")
        self.totalIfgs = len(self.ifgs)
        self.curIfg = self.ifgs[self.posIfg]
        
        self.retranslateUi(Dialog)
        QtCore.QMetaObject.connectSlotsByName(Dialog)
    
    def setIfgError(self):
        errCode = 0
        if self.rbNoError.isChecked()==True: errCode=0
        if self.btnError1.isChecked()==True: errCode=1
        if self.btnError2.isChecked()==True: errCode=2
        if self.btnError3.isChecked()==True: errCode=3
        if self.btnError4.isChecked()==True: errCode=4
        if self.btnError5.isChecked()==True: errCode=5
        if self.btnError6.isChecked()==True: errCode=6
        if self.btnError7.isChecked()==True: errCode=7
        self.ifgerrors[self.posIfg] = errCode

    def retranslateUi(self, Dialog):
        _translate = QtCore.QCoreApplication.translate
        Dialog.setWindowTitle(_translate("Dialog", "Manual check of ifgs for frame "+str(self.selectedFrame)))
        self.label.setText(_translate("Dialog", "wrapped ifg:"))
        self.label_2.setText(_translate("Dialog", "unwrapped ifg:"))
        #self.btnCancel.setText(_translate("Dialog", "Cancel"))
        self.label_3.setText(_translate("Dialog", "Type of error identified"))
        self.rbNoError.setText(_translate("Dialog", "0: no errors"))
        self.btnError1.setText(_translate("Dialog", "1: swath discontinuity ('vertical' line)"))
        self.btnError2.setText(_translate("Dialog", "2: unwrapping error (tiles visible)"))
        self.btnError3.setText(_translate("Dialog", "3: coregistration error ('ghosting')"))
        self.btnError7.setText(_translate("Dialog", "4: spectral div. error (azimuth fringes)"))
        self.btnError4.setText(_translate("Dialog", "5: missing bursts (gray/blue areas)"))
        self.btnError5.setText(_translate("Dialog", "6: missing (only) unwrapped data"))
        self.btnError6.setText(_translate("Dialog", "7: other (or mixed) error"))
        self.btnNextImage.setText(_translate("Dialog", "Next image pair"))
        self.btnPrevImage.setText(_translate("Dialog", "Previous"))
        
        self.totalIfgs = len(self.ifgs)
        pairpath = str(self.selectedTrack)+'/'+str(self.selectedFrame)+"/"+str(self.ifgs[self.posIfg])
        self.labelCheckingIfg.setText(_translate("Dialog", "Checking combination "+pairpath))
        #self.labelIfgNo.setText(_translate("Dialog", "ifg no "))
        self.labelIfgNo.setText(_translate("Dialog", "ifg no. "+str(int(self.posIfg)+1)+" from "+str(len(self.ifgs))))
        self.labelHelp.setText(_translate("Dialog", "try arrows, (shift+)PgUp, numbers, click image"))

    def showNextIfgMore(self):
        self.showNextIfg(20)

    def showPrevIfgMore(self):
        self.showPrevIfg(20)

    def showNextIfgMore100(self):
        self.showNextIfg(100)

    def showPrevIfgMore100(self):
        self.showPrevIfg(100)

    def showNextIfg(self, num = 1):
        _translate = QtCore.QCoreApplication.translate
        if self.posIfg+num >= self.totalIfgs:
            print('You have reached the end of the dataset')
            print('Will save your current inputs - thank you')
            saveIfgErrors(self.selectedFrame,self.ifgs,self.ifgerrors)
            self.labelIfgNo.setText(_translate("Dialog", "results saved"))
        else:
            if (self.posIfg > 0 and (np.mod(self.posIfg, 25)) == 0):
                print('autosaving in position {0}'.format(str(self.posIfg)))
                saveIfgErrors(self.selectedFrame,self.ifgs,self.ifgerrors)
            self.posIfg = self.posIfg+num
            self.curIfg = self.ifgs[self.posIfg]
            self.curIfgErr = self.ifgerrors[self.posIfg]
            pairpath = str(self.selectedTrack)+'/'+str(self.selectedFrame)+"/"+str(self.ifgs[self.posIfg])
            self.labelCheckingIfg.setText(_translate("Dialog", "Checking combination "+pairpath))
            self.labelIfgNo.setText(_translate("Dialog", "ifg no. "+str(int(self.posIfg)+1)+" from "+str(len(self.ifgs))))
            self.toggleErrors(self.posIfg)
            if self.posIfg+1 == self.totalIfgs:
                self.btnNextImage.setText(_translate("Dialog", "Save results"))
            if self.posIfg == 1:
                self.btnPrevImage.setText(_translate("Dialog", "Previous"))
            #self.rbNoError.toggle()
            self.drawIfgs(self)

    def showPrevIfg(self, num = 1):
        _translate = QtCore.QCoreApplication.translate
        if self.posIfg-num < 0:
            print('You have reached the end of the dataset')
            print('Will save your current inputs - thank you')
            saveIfgErrors(self.selectedFrame,self.ifgs,self.ifgerrors)
            self.labelIfgNo.setText(_translate("Dialog", "results saved"))
        else:
            self.posIfg = self.posIfg-num
            self.curIfg = self.ifgs[self.posIfg]
            self.curIfgErr = self.ifgerrors[self.posIfg]
            pairpath = str(self.selectedTrack)+'/'+str(self.selectedFrame)+"/"+str(self.ifgs[self.posIfg])
            self.labelCheckingIfg.setText(_translate("Dialog", "Checking combination "+pairpath))
            self.labelIfgNo.setText(_translate("Dialog", "ifg no. "+str(int(self.posIfg)+1)+" from "+str(len(self.ifgs))))
            self.toggleErrors(self.posIfg)
            if self.posIfg+1 == self.totalIfgs-1:
                self.btnNextImage.setText(_translate("Dialog", "Next image pair"))
            if self.posIfg == 0:
                self.btnPrevImage.setText(_translate("Dialog", "Save"))
            #self.rbNoError.toggle()
            self.drawIfgs(self)   
    
    def toggleErrors(self, posIfg):
        errNo = self.ifgerrors[posIfg]
        options = {0: self.rbNoError.toggle,
                   1: self.btnError1.toggle,
                   2: self.btnError2.toggle,
                   3: self.btnError3.toggle,
                   4: self.btnError4.toggle,
                   5: self.btnError5.toggle,
                   6: self.btnError6.toggle,
                   7: self.btnError7.toggle,
        }
        options[errNo]()
        #self.rbNoError.toggle()

    def drawIfgs(self, Dialog):
        ifgfile = os.path.join(var.outDir,self.selectedFrame,self.ifgs[self.posIfg]+'.geo.diff.png')
        unwfile = os.path.join(var.outDir,self.selectedFrame,self.ifgs[self.posIfg]+'.geo.unw.png')
        if os.path.exists(ifgfile):
            ifgimage = QtGui.QImage(ifgfile)
            self.graphicsViewIfg.setPixmap(QtGui.QPixmap(ifgimage))
            #self.graphicsViewIfg.mousePressEvent = self.viewFullImage(ifgfile)
            self.graphicsViewIfg.imagepath = ifgfile
        else:
            self.graphicsViewIfg.clear()
            self.graphicsViewIfg.imagepath = None
        if os.path.exists(unwfile):
            unwimage = QtGui.QImage(unwfile)
            self.graphicsViewUnw.setPixmap(QtGui.QPixmap(unwimage))
            #self.graphicsViewUnw.mousePressEvent = self.viewFullImage(unwfile)
            self.graphicsViewUnw.imagepath = unwfile
        else:
            self.graphicsViewUnw.clear()
            self.graphicsViewUnw.imagepath = None

    class Ui_ImageView(object):
        imagepath = None #'/nfs/a1/insar/lics_check/003D_09757_111111/20170510_20170522.geo.unw.bmp'
        def setupUi(self, Dialog):
            Dialog.setObjectName("Dialog")
            Dialog.setFocusPolicy(QtCore.Qt.StrongFocus)
            self.imageView = QtWidgets.QLabel(Dialog)
            image = QtGui.QPixmap(self.imagepath)
            self.imageView.setPixmap(image)
            Dialog.resize(image.width(), image.height())
    
    def viewFullImageIfg(self):
        print(self.graphicsViewIfg.imagepath)
        imagepath = self.graphicsViewIfg.imagepath
        if imagepath is not None:
            self.DialogView = QtWidgets.QDialog()
            self.UiImView = self.Ui_ImageView()
            #imagepath = '/nfs/a1/insar/lics_check/003D_09757_111111/20170510_20170522.geo.unw.bmp'
            self.UiImView.imagepath = imagepath
            self.UiImView.setupUi(self.DialogView)
            self.DialogView.show()
            self.DialogView.raise_()

    def viewFullImageUnw(self):
        print(self.graphicsViewUnw.imagepath)
        imagepath = self.graphicsViewUnw.imagepath
        if imagepath is not None:
            self.DialogView = QtWidgets.QDialog()
            self.UiImView = self.Ui_ImageView()
            #imagepath = '/nfs/a1/insar/lics_check/003D_09757_111111/20170510_20170522.geo.unw.bmp'
            self.UiImView.imagepath = imagepath
            self.UiImView.setupUi(self.DialogView)
            self.DialogView.show()
            self.DialogView.raise_()

def get_framelist():
    import requests
    a=requests.get('http://gws-access.ceda.ac.uk/public/nceo_geohazards/LiCSAR_products/framelist.txt')
    #initiate framelist as dictionary d
    framelist = dict()
    for i in range(1,175+1):
        framelist.update({i:['']})
    for frame in a.text.split('\n'):
        if frame is not '':
            track = int(frame[0:3])
            framelist[track].append(frame)
    return framelist

def getUncheckedFrames(frameList):
    resFiles = [file for file in os.listdir(var.outDir) if fnmatch(file, '*.savedResults')]
    checkedFrames = [f.split('.')[0] for f in resFiles]
    framelist = frameList.copy()
    framelist = copy.deepcopy(frameList)
    for track in framelist:
        frames = framelist[track]
        for item in frames:
            if item in checkedFrames: frames.remove(item)
        #framelist.update({track:frames})
    return framelist

def getIfgList(frame):
    #unws=var.tempDir
    #ifgs=var.tempDir
    fullpath=os.path.join(var.outDir,frame)
    #unws =  [file for file in os.listdir(fullpath) if fnmatch(file, '*unw.bmp')]
    #ifgs =  [file for file in os.listdir(fullpath) if fnmatch(file, '*diff.bmp')]
    allbmps = [file for file in os.listdir(fullpath) if fnmatch(file, '*.png')]
    ifglist = []
    for item in allbmps:
        ifglist.append(item.split('.')[0])
    ifglist = list(set(ifglist))
    ifglist.sort()
    #print(len(ifglist))
    #print(ifglist[0])
    #os.listdir(os.path.join(var.tempDir,frame))
    return ifglist

def saveIfgErrors(frame,ifglist,ifgerrors):
    ifgset=[ifglist,ifgerrors]
    #fullpath=os.path.join(var.outDir,frame)
    pickleFile=frame+'.savedResults'
    pickleFilePath=os.path.join(var.outDir,pickleFile)
    outfile = open(pickleFilePath,'wb')
    pickle.dump(ifgset,outfile)
    outfile.close()
    #os.chmod(pickleFilePath,stat.S_IRWXG)
    try:
        os.chmod(pickleFilePath,0o777)
    except:
        print('warning, could not change savedResults permissions')
    print('Saved to '+pickleFilePath)
    #print('Please send this file to M.Lazecky@leeds.ac.uk')
    
def getIfgErrors(frame):
    #fullpath=os.path.join(var.outDir,frame)
    pickleFile=frame+'.savedResults'
    ifglist = getIfgList(frame)
    lastpos = 0
    if os.path.exists(os.path.join(var.outDir,pickleFile)):
        #print('yes - just load the pickle')
        infile = open(os.path.join(var.outDir,pickleFile),'rb')
        ifgset = pickle.load(infile)
        infile.close()
        ifgerrorslist = ifgset[1]
        #position of last erroneous ifg
        try:
            lastpos = max([i for i, e in enumerate(ifgerrorlist) if e != 0])
        except:
            lastpos = 0
        if len(ifgerrorslist)<len(ifglist):
            lastpos = len(ifgerrorslist) - 1
            bb = [0] * (len(ifglist)-len(ifgerrorslist))
            ifgerrorslist = ifgerrorslist + bb
    else:
        ifgerrorslist = [0] * len(ifglist)
    return ifgerrorslist, lastpos

def downloadIfgs(track,frame):
    outDir = os.path.join(var.outDir,str(frame))
    if not os.path.exists(outDir):
        os.mkdir(outDir)
    rc = os.system('chmod -R 777 '+outDir+' 2>/dev/null')
    #get list of ifgs:
    webaddr = 'http://gws-access.ceda.ac.uk/public/nceo_geohazards/LiCSAR_products/'+str(track)+'/'+str(frame)+'/interferograms'
    print('getting list of interferograms in LiCSAR')
    a = requests.get(webaddr)
    b = html.fromstring(a.content)
    ifgs = b.xpath('//tr/td/a/text()')
    ifgs.remove(ifgs[0])
    ifgs = [w.split('/')[0] for w in ifgs]
    for ifg in ifgs:
        print('checking/downloading '+ifg)
        for ext in ['unw.png', 'diff.png']:
            webpath = webaddr+'/'+ifg+'/'+ifg+'.geo.'+ext
            outpath = os.path.join(outDir,ifg+'.geo.'+ext)
            ##cmd = ['wget','-O',outpath,'-nc',webpath]
            #cmd = 'wget -O {0} -nc {1} 2>/dev/null'.format(outpath,webpath)
            #rc = os.system(cmd)
            #wget.download(url, '/Users/scott/Downloads/cat4.jpg')
            wget.download(webpath, outpath)
    #wgetcmd = 'wget -r -nd -c -e robots=off -A unw.png,diff.png -P '+outDir
    #cmd = wgetcmd+' '+webadd
    #cmd = [wgetcmd,webadd]
    #cmd = ['wget','-r','-nd','-c','-np','-e','robots=off','-A','unw.png,diff.png','-P',outDir,webaddr]
    #logfilename = os.path.join(outDir,'wget.log')
    #with open(logfilename,'a') as f:
    #    try:
    #        subp.check_call(cmd,stdout=f)
    #        #print('jojojo')
    #    except:
    #        print('some error happened downloading the products')
    #        #messageBox('some error happened (maybe the frame has no data ready?)')
    #print('rc je')
    #print(rc)
    #if rc != 0:
    #    return False
    #else:
    #    return True
    return True

if __name__ == "__main__":
    outDir = '/nfs/a1/insar/lics_check'
    if not os.path.exists(outDir):
        print('warning, directory {} does not exist - perhaps you do not run it at Leeds Uni. No problem though'.format(outDir))
        #try:
            #os.mkdir(outDir)
        #except:
        outDir=os.path.expanduser('licscheck')
        print('we will use directory {} instead'.format(outDir))
        if not os.path.exists(outDir): os.mkdir(outDir)
    var.outDir = outDir
    #app = QtWidgets.QApplication(sys.argv)
    app = QtWidgets.QApplication([])
    #first of all - load all frames that exist within LiCS database
    Dialog = QtWidgets.QDialog()
    uiChooseFrames = Ui_dialogChooseFrames()
    uiChooseFrames.setupUi(Dialog)
    print('Getting list of existing frames')
    uiChooseFrames.frameListAll = get_framelist()
    uiChooseFrames.frameList = getUncheckedFrames(uiChooseFrames.frameListAll)
    #print(len(uiChooseFrames.frameListAll))
    #print(len(uiChooseFrames.frameList))
    print('done')
    #setting tempdir
    #user = os.environ['USER']
    #tempDir = '/tmp/licscheck-'+user
    #ui.listFrames.addI
            #modality not working in Qt5??
    #self.setWindowModality(QtCore.Qt.ApplicationModal)
    #self.setWindowModality(QtCore.Qt.WindowModal)
    Dialog.show()
    sys.exit(app.exec_())

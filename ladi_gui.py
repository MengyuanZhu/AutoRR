#!/usr/bin/python
import sys
from PyQt4.QtCore import *
from PyQt4.QtGui import *
from rdkit import Chem
from rdkit.Chem import AllChem
import openbabel
import threading

# create our window
app = QApplication(sys.argv)
w = QWidget()
w.setWindowTitle('LADI v1.0')
w.setGeometry(50,50,230,440)

#open and save buttons
btnOpen = QPushButton('Open', w)
btnOpen.setGeometry(10,10,100,30)
btnSave = QPushButton('Save', w)
btnSave.setGeometry(120,10,100,30)

#library selection
labelFragment=QLabel("Fragment libraries: ",w )
labelFragment.move(5,45)
chkbox1=QCheckBox('hydrophobic',w)
chkbox1.move(20,145)
chkbox2=QCheckBox('aromatic',w)
chkbox2.move(20,125)
chkbox3=QCheckBox('polar_negative',w)
chkbox3.move(20,105)
chkbox4=QCheckBox('polar_positive',w)
chkbox4.move(20,85)
chkbox5=QCheckBox('polar_uncharged',w)
chkbox5.move(20,65)

#optimization
chkboxOpt=QCheckBox('Optimize 3D structure',w)
chkboxOpt.move(5,175)
radioUFF=QRadioButton('Universal Force Field',w)
radioUFF.move(20,195)
radioMMFF=QRadioButton('MMFF Force Field',w)
radioMMFF.move(20,215)
radioUFF.setEnabled(False)
radioMMFF.setEnabled(False)

#protonation
chkboxProt=QCheckBox('Fixpka \nBased on Evan\'s pKa table',w)
chkboxProt.move(5,245)

#format
lOutput=QLabel("Output format: ",w)
lOutput.move(5,305)
cbOutput=QComboBox(w)
cbOutput.addItem("mol2")
cbOutput.addItem("sdf")
cbOutput.addItem("mol")
cbOutput.addItem("smi")
cbOutput.addItem("pdbqt")
cbOutput.move(120,300)

#progressbar
lProgress=QLabel("Progress:",w)
lProgress.move(5,335)
pbLibrary=QProgressBar(w)
pbLibrary.setGeometry(5,5,220,20)
pbLibrary.setRange(0,100)
pbLibrary.setValue(0)
pbLibrary.move(5,355)

#start
btnStart = QPushButton('Start', w)
btnStart.setGeometry(90,90,90,50)
btnStart.move(75,380)

#global variables
openfilename=""
savefilename=""
librarySize=0
libraryDone=0
obConversion=openbabel.OBConversion()
obmol=openbabel.OBMol()
mutexlock=threading.Lock()
# Create the actions
@pyqtSlot()
def onClickOpen():
	global openfilename
	openfilename=QFileDialog.getOpenFileName()

@pyqtSlot()
def onClickSave():
	global savefilename
	savefilename= QFileDialog.getSaveFileName()

@pyqtSlot()
def onClickOpt():
	if (chkboxOpt.isChecked()):
		radioUFF.setEnabled(True)
		radioMMFF.setEnabled(True)
		radioUFF.setChecked(True)
	else:
		radioUFF.setEnabled(False)
		radioMMFF.setEnabled(False)


@pyqtSlot()
def updateProgress(val):
	pbLibrary.setValue(val)


@pyqtSlot()
def onClickStart():
	global openfilename
	global savefilename
	global obConversion
	global librarySize

	#check input file
	if (openfilename==""):
		msg = QMessageBox()
		msg.setIcon(QMessageBox.Warning)
		msg.setText("Lead is not defined.")
   		msg.setWindowTitle("Warning")
  		msg.setDetailedText("Click \"Open\" button to select a lead compound.")
   		msg.setStandardButtons(QMessageBox.Ok)
   		msg.exec_()
		return

	#check output file
	if (savefilename==""):
		msg = QMessageBox()
		msg.setIcon(QMessageBox.Warning)
		msg.setText("Save file is not defined")
   		msg.setWindowTitle("Warning")
  		msg.setDetailedText("Click \"save\" button to select the save location.")
   		msg.setStandardButtons(QMessageBox.Ok)
   		msg.exec_()
		return

	#check format
	fileextension=openfilename.split(".")[-1]
	supported=["mol2","sdf","mol","pdb"]
	if str(fileextension).lower() not in supported:
		msg = QMessageBox()
		msg.setIcon(QMessageBox.Warning)
		msg.setText("Input format not supported")
   		msg.setWindowTitle("Warning")
  		msg.setDetailedText("We support formats mol2, pdb, sdf and mol.")
   		msg.setStandardButtons(QMessageBox.Ok)
   		msg.exec_()

	lProgress.setText("Progress")
	obConversion.SetInAndOutFormats("pdb",str(cbOutput.currentText()))

	#library files
	file1="hydrophobic.lib"
	file2="aromatic.lib"
	file3="polar_negative.lib"
	file4="polar_positive.lib"
	file5="polar_uncharged.lib"

	#process aromatic
	arolink=file("aromatic-linker.lib","r")
	aroend=file("aromatic-end.lib","r")
	arocom=file("aromatic-complex.lib","r")
	arrlink=[]
	arrend=[]
	line=arolink.readline()
	while True:
		if not line:
			break
		arrlink.append(line)
		line=arolink.readline()
	line=aroend.readline()
	while True:
		if not line:
			break
		arrend.append(line)
		line=aroend.readline()
	aro=file("aromatic.lib","w")
	for link in arrlink:
		for end in arrend:
			aro.write(link[:len(link)-1]+end)
	line=arocom.readline()
	while True:
		if not line:
			break
		aro.write(line)
		line=arocom.readline()
	#end of process aromatic

	#calculate library sizes
	sizeHydrophobic = sum(1 for line in open(file1))
	sizeAromatic = sum(1 for line in open(file2))
	sizeNegative = sum(1 for line in open(file3))
	sizePositive = sum(1 for line in open(file4))
	sizeUncharged = sum(1 for line in open(file5))


	if str(fileextension)=="mol2":
		m=Chem.MolFromMol2File(str(openfilename))  #readfile in mol2 format
	if str(fileextension)=="mol":
		m=Chem.MolFromMolFile(str(openfilename))   #readfile in mol format
	if str(fileextension)=="pdb":
		m=Chem.MolFromPDBFile(str(openfilename))   #readfile in pdb format
	if str(fileextension)=="sdf":
		suppl = SDMolSupplier(str(openfilename))   #readfile in sdf format
		m= suppl.next()

	patt=Chem.MolFromSmarts("[Au]")
	if (str(savefilename).split(".")[-1]!=str(cbOutput.currentText())):
		savefilename=str(savefilename)+"."+str(cbOutput.currentText())
	writer=open(str(savefilename),"w")
	librarySize=0
	if chkbox1.checkState()!=0:#hydrophobic
		librarySize=librarySize+sizeHydrophobic
	if chkbox2.checkState()!=0:#aromatic
		librarySize=librarySize+sizeAromatic
	if chkbox3.checkState()!=0:#negative
		librarySize=librarySize+sizeNegative
	if chkbox4.checkState()!=0:#positivie
		librarySize=librarySize+sizePositive
	if chkbox5.checkState()!=0:#uncharged
		librarySize=librarySize+sizeUncharged


	if chkbox1.checkState()!=0:
		thread1=threading.Thread(target=replace, args=(m,patt,writer,file1,"hydrophobic"))
		thread1.start()
		
	if chkbox2.checkState()!=0:
		thread2=threading.Thread(target=replace, args=(m,patt,writer,file2,"aromatic"))
		thread2.start()
		
	if chkbox3.checkState()!=0:
		thread3=threading.Thread(target=replace, args=(m,patt,writer, file3,"negative"))
		thread3.start()
	
	if chkbox4.checkState()!=0:
		thread4=threading.Thread(target=replace, args=(m,patt,writer, file4,"positive"))
		thread4.start()
		
	if chkbox5.checkState()!=0:
		thread5=threading.Thread(target=replace, args=(m,patt,writer,file5,"uncharged"))
		thread5.start()
		

#molecules generation
def replace(m,patt,writer,library,fragmentType):
	f=open(library,"r")
	line=f.readline()
	global libraryDone
	global mutexlock
	i=1
	while True:
		if not line:
			break
		print line
		repl=Chem.MolFromSmiles(line)

		rms=AllChem.ReplaceSubstructs(m,patt,repl)
		#rms[0].SetProp("_Name","ligand_"+"fragmentType"+str(i))
		Chem.SanitizeMol(rms[0])

		smiles=Chem.MolToSmiles(rms[0])
		mol=Chem.MolFromSmiles(smiles)
		mol=Chem.AddHs(mol)
		Chem.rdPartialCharges.ComputeGasteigerCharges(mol)
		AllChem.EmbedMolecule(mol)
		if (radioUFF.isChecked()):
			AllChem.UFFOptimizeMolecule(mol)
		if (radioMMFF.isChecked()):
			AllChem.MMFFOptimizeMolecule(mol)
		obConversion.ReadString(obmol, Chem.MolToPDBBlock(mol))
		obmol.SetTitle("ligand_"+fragmentType+str(i))
		writer.write(obConversion.WriteString(obmol))
		mutexlock.acquire()  #mutex lock for the critical section
		libraryDone=libraryDone+1
		mutexlock.release()
		val=int(libraryDone*100/librarySize)
		if (val%2==0):
			w.emit(SIGNAL("test"),val)
		i=i+1
		line=f.readline()

# connect the signals to the slots
btnOpen.clicked.connect(onClickOpen)
btnSave.clicked.connect(onClickSave)
btnStart.clicked.connect(onClickStart)
chkboxOpt.clicked.connect(onClickOpt)
QObject.connect(w,SIGNAL("test"),updateProgress)

w.show()

app.exec_()

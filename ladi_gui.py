#!/usr/bin/python
import sys
from PyQt4.QtCore import pyqtSlot
from PyQt4.QtGui import *
from rdkit import Chem
from rdkit.Chem import AllChem
import openbabel

obConversion=openbabel.OBConversion()

obmol=openbabel.OBMol()
 
# create our window
app = QApplication(sys.argv)
w = QWidget()
w.setWindowTitle('LADI v1.0')
w.setGeometry(50,50,230,450)

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
radioMMFF=QRadioButton('MMFF Field',w)
radioMMFF.move(20,215)
radioUFF.setEnabled(False)
radioMMFF.setEnabled(False)

#partical charge
chkboxParChar=QCheckBox('Computer partial charge \nGasteigerCharge',w)
chkboxParChar.move(5,245)

#protonation
chkboxProt=QCheckBox('Fixpka \nBased on Evan\'s pKa table',w)
chkboxProt.move(5,295)

#format
lOutput=QLabel("Output format: ",w)
lOutput.move(5,355)
cbOutput=QComboBox(w)
cbOutput.addItem("mol2")
cbOutput.addItem("sdf")
cbOutput.addItem("mol")
cbOutput.addItem("smi")
cbOutput.addItem("pdbqt")
cbOutput.move(120,350)

#start
btnStart = QPushButton('Start', w)
btnStart.setGeometry(90,90,90,50)
btnStart.move(75,390)

#global variables
openfilename=""
savefilename=""

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
def onClickStart():
	global openfilename
	global savefilename
	global obConversion
	

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


   	print str(cbOutput.currentText())

	obConversion.SetInAndOutFormats("pdb",str(cbOutput.currentText()))
	file1="hydrophobic.lib"
	file2="hydrophobic_aromatic.lib"
	file3="polar_negative.lib"
	file4="polar_positive.lib"
	file5="polar_uncharged.lib"

	if str(fileextension)=="mol2":
		m=Chem.MolFromMol2File(str(openfilename))  #readfile in mol2 format
	if str(fileextension)=="mol":
		m=Chem.MolFromMolFile(str(openfilename))   #readfile in mol format
	if str(fileextension)=="pdb":
		m=Chem.MolFromPDBFile(str(openfilename))   #readfile in pdb format
	if str(fileextension)=="sdf":
		suppl = SDMolSupplier(str(openfilename))
		m= suppl.next()   						   #readfile in sdf format



	patt=Chem.MolFromSmarts("[Au]")
	

	writer=open(str(savefilename)+"."+str(cbOutput.currentText()),"w")
	
	i=1

	if chkbox1.checkState()!=0:
		i=replace(m,patt,writer,i, file1)
	if chkbox2.checkState()!=0:
		i=replace(m,patt,writer,i, file2)
	if chkbox3.checkState()!=0:
		i=replace(m,patt,writer,i, file3)
	if chkbox4.checkState()!=0:
		i=replace(m,patt,writer,i, file4)
	if chkbox5.checkState()!=0:
		i=replace(m,patt,writer,i, file5)


#molecules generation	
def replace(m,patt,writer,i, library):
	f=open(library,"r")
	line=f.readline()
	while True:
		if not line:
			break
		print line	
		repl=Chem.MolFromSmiles(line)
	
		rms=AllChem.ReplaceSubstructs(m,patt,repl)
		rms[0].SetProp("_Name","ligand"+str(i))
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
		obmol.SetTitle("ligand"+str(i))	
		writer.write(obConversion.WriteString(obmol))

		i=i+1
		line=f.readline()
	return i
	
# connect the signals to the slots
btnOpen.clicked.connect(onClickOpen)
btnSave.clicked.connect(onClickSave)
btnStart.clicked.connect(onClickStart)
chkboxOpt.clicked.connect(onClickOpt)

w.show()

app.exec_()

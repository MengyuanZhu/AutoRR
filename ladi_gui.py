import sys
from PyQt4.QtCore import pyqtSlot
from PyQt4.QtGui import *
from rdkit import Chem
from rdkit.Chem import AllChem

 
# create our window
app = QApplication(sys.argv)
w = QWidget()
w.setWindowTitle('LADI v1.0')
w.setGeometry(50,50,200,300)
# Create a button in the window
btn = QPushButton('Open a File', w)
btn.move(50,2)
label1=QLabel("Fragment libraries: ",w )
label1.move(0,35)
chkbox1=QCheckBox('hydrophobic',w)
chkbox1.move(0,135)
chkbox2=QCheckBox('aromatic',w)
chkbox2.move(0,115)
chkbox3=QCheckBox('polar_negative',w)
chkbox3.move(0,95)
chkbox4=QCheckBox('polar_positive',w)
chkbox4.move(0,75)
chkbox5=QCheckBox('polar_uncharged',w)
chkbox5.move(0,55)

btn1 = QPushButton('Save file location', w)

btn1.move(35,170)

btn2 = QPushButton('Start', w)
btn2.setGeometry(90,90,90,90)
btn2.move(55,200)

openfilename=""
savefilename=""




# Create the actions 
@pyqtSlot()
def on_click():
	global openfilename
	openfilename=QFileDialog.getOpenFileName()
 
@pyqtSlot()
def on_click_save():
	global savefilename
	savefilename= QFileDialog.getSaveFileName()

@pyqtSlot()
def on_click_start():
	global openfilename
	global savefilename
	file1="hydrophobic.lib"
	file2="hydrophobic_aromatic.lib"
	file3="polar_negative.lib"
	file4="polar_positive.lib"
	file5="polar_uncharged.lib"

	print openfilename
	print savefilename
	m=Chem.MolFromMol2File(str(openfilename))
	patt=Chem.MolFromSmarts("[Au]")
	
	writer=Chem.SDWriter(str(savefilename))
	
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
	
		writer.write(rms[0])
		i=i+1
		line=f.readline()
	return i

	
 
# connect the signals to the slots
btn.clicked.connect(on_click)
btn1.clicked.connect(on_click_save)
btn2.clicked.connect(on_click_start)


w.show()

app.exec_()

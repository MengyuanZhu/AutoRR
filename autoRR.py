#!/usr/bin/python
from rdkit import Chem
from rdkit.Chem import AllChem
from openeye.oechem import *
from openeye.oeomega import *
import sys

pkaLib={}

libraryFile=["library/aromatic.lib","library/hydrophobic.lib","library/polar_negative.lib","library/polar_positive.lib","library/polar_uncharged.lib"]
library=["aromatic","hydrophobic","polar_negative","polar_positive","polar_uncharged"]

def preparepKa():
	pkafile=open("pka.lib","r")
	line=pkafile.readline()
	while True:
		if not line:
			break
		pkaItem=line.split(" ")
		pkaLib[pkaItem[0]]=pkaItem[1]
		line=pkafile.readline()

def prepareAromatic():
	#Generate aromatic library
	#It will read linker, end and complex to generate aromatic.lib
	arolink=file("library/aromatic-linker.lib","r")
	aroend=file("library/aromatic-end.lib","r")
	arocom=file("library/aromatic-complex.lib","r")
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
	aro=file("library/aromatic.lib","w")
	for link in arrlink:
		for end in arrend:
			aro.write(link[:len(link)-1]+end)
	line=arocom.readline()
	while True:
		if not line:
			break
		aro.write(line)
		line=arocom.readline()

def main(argv=[__name__]):
	
	if len(argv) != 3:
		print("%s <infile> <outfile>" % argv[0])
		exit()

	ofs = oemolostream()
	if not ofs.open(argv[2]):
		OEThrow.Fatal("Unable to open %s for writing" % argv[2])
	
	fragmentIndex=input("Library 0: Aromatic, 1: Hydrophobic, 2: Polar Negative, 3: Polar Positive, 4: Polar Uncharged\nLibrary:")

	numConformation=input("Maximum number of conformation:")
	m=Chem.MolFromMol2File(argv[1])
	patt=Chem.MolFromSmarts("[Au]")

	prepareAromatic()

	numFragments = sum(1 for line in open(libraryFile[fragmentIndex]))

	fragmentFile=open(libraryFile[fragmentIndex],"r")
	line=fragmentFile.readline()
	
	i=1

	while True:
		if not line:
			break

		repl=Chem.MolFromSmiles(line)
	
		rms=AllChem.ReplaceSubstructs(m,patt,repl)
		smiles=Chem.MolToSmiles(rms[0])

		oemol=OEMol()
		omega = OEOmega()
		omega.SetMaxConfs(numConformation)
		omega.SetStrictStereo(False)	
		if OESmilesToMol(oemol,smiles):
			print str(round(100.0*i/numFragments,2))+"%"
		oemol.SetTitle("molecule_"+library[fragmentIndex]+"_"+str(i))
		if omega(oemol):
			OEWriteMolecule(ofs, oemol)
		i=i+1
		line=fragmentFile.readline()
	
if __name__ == "__main__":
    sys.exit(main(sys.argv))








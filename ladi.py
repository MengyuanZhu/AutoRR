#!/usr/bin/python
from rdkit import Chem
from rdkit.Chem import AllChem
from openeye.oechem import *
from openeye.oeomega import *
import sys

def main(argv=[__name__]):

	if len(argv) != 3:
		raise ValueError("%s <infile> <outfile>" % argv[0])

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
	#end of process aromatic



	m=Chem.MolFromMol2File("example.mol2")
	patt=Chem.MolFromSmarts("[Au]")
	filename="result"
	file1="library/aromatic.lib"
	file2="library/hydrophobic.lib"
	file3="library/polar_negative.lib"
	file4="library/polar_positive.lib"
	file5="library/polar_uncharged.lib"


	f=open(file2,"r")

	pkafile=open("pka.lib","r")

	line=pkafile.readline()
	pkaLib={}

	while True:
		if not line:
			break
		pkaItem=line.split(" ")
		pkaLib[pkaItem[0]]=pkaItem[1]
		line=pkafile.readline()

	line=f.readline()
	i=1

	ofs = oemolostream()
	if not ofs.open(argv[1]):
		OEThrow.Fatal("Unable to open %s for writing" % argv[1])

	while True:
		if not line:
			break

		repl=Chem.MolFromSmiles(line)
	
		rms=AllChem.ReplaceSubstructs(m,patt,repl)
		smiles=Chem.MolToSmiles(rms[0])

		oemol=OEMol()
		omega = OEOmega()
		omega.SetMaxConfs(1)
		omega.SetStrictStereo(False)	
		if OESmilesToMol(oemol,smiles):
			print oemol

		if omega(oemol):
			OEWriteMolecule(ofs, oemol)
		i=i+1
		line=f.readline()
	
if __name__ == "__main__":
    sys.exit(main(sys.argv))


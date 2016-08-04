from rdkit import Chem
from rdkit.Chem import AllChem
import openbabel

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

obConversion=openbabel.OBConversion()
obConversion.SetInAndOutFormats("pdb","mol2")
obmol=openbabel.OBMol()
ff = openbabel.OBForceField.FindForceField("mmff94")
obbuilder=openbabel.OBBuilder()

m=Chem.MolFromMol2File("ligand_au.mol2")
patt=Chem.MolFromSmarts("[Au]")
filename="result"
file1="aromatic-linker.lib"
file2="hydrophobic.lib"
file3="polar_negative.lib"
file4="polar_positive.lib"
file5="polar_uncharged.lib"
file6="aromatic-complex.lib"

f=open(file5,"r")

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

while True:
	if not line:
		break
	#print line	
	repl=Chem.MolFromSmiles(line)
	
	rms=AllChem.ReplaceSubstructs(m,patt,repl)
	smiles=Chem.MolToSmiles(rms[0])
	
	
	#for key, value in pkaLib.iteritems():
	#	print
	rms=Chem.MolFromSmiles(smiles)
	rms=AllChem.ReplaceSubstructs(rms,Chem.MolFromSmarts("C(=O)[OH]"),Chem.MolFromSmiles("C([O-])=O"), True)
	smiles=Chem.MolToSmiles(rms[0])

	rms=Chem.MolFromSmiles(smiles)
	rms=AllChem.ReplaceSubstructs(rms,Chem.MolFromSmarts("[N;H2;!$(Nc)]"),Chem.MolFromSmiles("[NH3+]"), True)
	smiles=Chem.MolToSmiles(rms[0])

	rms=Chem.MolFromSmiles(smiles)
	rms=AllChem.ReplaceSubstructs(rms,Chem.MolFromSmarts("[N;H1;!$(NC=O)]"),Chem.MolFromSmiles("[NH2+]"), True)
	smiles=Chem.MolToSmiles(rms[0])
	#Chem.SanitizeMol(rms[0])
	
	mol=Chem.MolFromSmiles(smiles)
	mol=Chem.AddHs(mol)
	AllChem.EmbedMolecule(mol)
	AllChem.MMFFOptimizeMolecule(mol)
	obConversion.ReadString(obmol, Chem.MolToPDBBlock(mol))
	print obConversion.WriteString(obmol)
	i=i+1
	line=f.readline()
	



from rdkit import Chem
from rdkit.Chem import AllChem
import openbabel

obConversion=openbabel.OBConversion()
obConversion.SetInAndOutFormats("pdb","mol2")
obmol=openbabel.OBMol()
ff = openbabel.OBForceField.FindForceField("mmff94")
obbuilder=openbabel.OBBuilder()

m=Chem.MolFromMol2File("ligand_au.mol2")
patt=Chem.MolFromSmarts("[Au]")
filename="result"
file1="hydrophobic_aromatic.lib"
file2="hydrophobic.lib"
file3="polar_negative.lib"
file4="polar_positive.lib"
file5="polar_uncharged.lib"

f=open(file4,"r")

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
	



for key, value in pkaLib.items():
	print "SB"+key+"SB"+value+"SB"
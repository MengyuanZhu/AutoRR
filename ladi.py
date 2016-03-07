from rdkit import Chem
from rdkit.Chem import AllChem


m=Chem.MolFromMol2File("ligand_au.mol2")
patt=Chem.MolFromSmarts("[Au]")
filename="result"
f=open("hydrophobic_aromatic.lib","r")


w=Chem.SDWriter("results.sdf")
line=f.readline()
i=1

while True:
	if not line:
		break
	print line	
	repl=Chem.MolFromSmiles(line)
	
	rms=AllChem.ReplaceSubstructs(m,patt,repl)
	rms[0].SetProp("_Name","ligand"+str(i))
	
	w.write(rms[0])
	i=i+1
	line=f.readline()
	




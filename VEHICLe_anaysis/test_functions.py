import pandas as pd
from VEHICLe_read import read_csv_to_dataframe
from VEHICLe_read import read_list_of_smiles
from VEHICLe_read import read_mol_from_smiles
import rdkit
from mol_translator.structure.structure_write import write_mol_toxyz
from openbabel import pybel
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdmolfiles
VEHICLe=read_csv_to_dataframe('VEHICLe.csv')
head=VEHICLe.head()

sample=VEHICLe.sample()
smiles=sample['Smiles']
for j in smiles:
    print(j)
    y=read_mol_from_smiles(j)

    AllChem.Compute2DCoords(y)
    for c in y.GetConformers():
        p=c.GetPositions()

        ### prints numpy array with coords for each atom but how do we get them out of that? ###


#j=['1','2','3','4','5', '6', '7', '8', '9', '10']
#for k in j:
#    x=VEHICLe.sample()
#    y=read_list_of_smiles(x)#

 #   for i in y:
 #       print(i)
 #       mol=read_mol_from_smiles(i)
 #       molHs=Chem.AddHs(mol)
 #       xyz=Chem.MolToXYZFile(molHs)
 #       print(xyz)
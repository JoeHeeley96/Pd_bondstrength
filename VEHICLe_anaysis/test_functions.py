import pandas as pd
from VEHICLe_read import read_csv_to_dataframe
from VEHICLe_read import read_list_of_smiles
from VEHICLe_read import read_mol_from_smiles
import numpy as np
import rdkit
import sys
#from mol_translator.structure.structure_write import write_mol_toxyz
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
    x=Chem.AddHs(y)

    z=AllChem.Compute2DCoords(x)
    type_array = np.zeros(x.GetNumAtoms(), dtype=np.int32)


    for j, atoms in enumerate(x.GetAtoms()):
        type_array[j] = atoms.GetAtomicNum()

    for c in x.GetConformers():
        xyz=c.GetPositions()
        xyz_list=xyz.tolist()
        print(xyz_list)

    match_coords=list(zip(type_array, xyz_list))
    print(match_coords)
    with open('galoo.txt', 'w') as f:
        for i in match_coords:
            print(*i, sep=', ', file=f)




            #str=np.savetxt(sys.stdout, j, newline='', fmt='%.3f')




           # with open('test.txt', 'w') as f:
            #    print(str, file=f)

    #    for j in c.GetPositions():
    #        for i in type_array:
    #            with open('oogy.xyz', 'w') as f:
    #               print(i, j, file=f)


    #for c in y.GetConformers():
     #   with open('oogy.xyz', 'w') as f:
      #      print(c.GetPositions(), file=f)

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
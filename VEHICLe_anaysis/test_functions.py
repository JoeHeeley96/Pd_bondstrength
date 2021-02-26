import pandas as pd
import glob
from VEHICLe_read import read_csv_to_dataframe
from VEHICLe_read import read_list_of_smiles
from VEHICLe_read import read_mol_from_smiles
import numpy as np
import rdkit
import sys
import cirpy
from datetime import date
from file_generation import VEHICLe_string_to_com
#from mol_translator.structure.structure_write import write_mol_toxyz
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdmolfiles
VEHICLe=read_csv_to_dataframe('VEHICLe.csv')
head=VEHICLe.head(30)



#VEHICLe_string_to_com(head)

with open('com_files/S1_2021-02-26_wb97xd_631gd_opt.com') as f:
    H_num=[]
    H_coord=[]
    for line in f:
        if 'H' in line:
            H_coord.append(line)

    for i in range(len(H_coord)):
        H_num.append(i + 1)

    zip=zip(H_num, H_coord)


    with open('test_anions/S1_'  + 'anion_2021-02-26_wb97xd_631gd_opt.com', 'w') as p:
        for i in zip:
            print(i[1])
            lines=f.readlines()
            f.seek(0)
            for j in lines:
                print(j)
                if j != i[1]:
                    p.write(line)



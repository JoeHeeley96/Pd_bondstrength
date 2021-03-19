import glob
import sys
import cirpy
from datetime import date
import rdkit
import pandas as pd
import glob
from VEHICLe_read import read_csv_to_dataframe
from VEHICLe_read import read_list_of_smiles
from VEHICLe_read import read_mol_from_smiles
import numpy as np
import rdkit
from com_files_from_smiles import VEHICLe_string_to_com
#from mol_translator.structure.structure_write import write_mol_toxyz
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdmolfiles
from anion_generation import anion_from_com

comfilenames=glob.glob('neutral_comfiles/*')

for file in comfilenames:
    anion_from_com(file)


#VEHICLe_string_to_com(head)

#with open('neutral_comfiles/S1_2021-02-26_wb97xd_631gd_opt.com') as f:
 #   H_num=[]
  #  H_coord=[]
   # for line in f:
    #    if 'H' in line:
     #       H_coord.append(line)
#
 #   for i in range(len(H_coord)):
  #      H_num.append(i + 1)


   # zip=zip(H_num, H_coord)

    #for j, k in zip:
#
 #       with open('anion_comfiles/S1_' + str(j) + 'anion_2021-02-26_wb97xd_631gd_opt.com', 'w') as p:
  #              with open('neutral_comfiles/S1_2021-02-26_wb97xd_631gd_opt.com') as f:
   #                     for line in f:
    #                        if line == '0 1\n':
     #                           print('-1 1', file =p)
      #                      elif line != k:
       #                         print(line.strip('\n'), file=p)





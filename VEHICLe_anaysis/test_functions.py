import pandas as pd
import glob
from datetime import date
from bromine_generation import bromines_from_com
from com_files_from_smiles import VEHICLe_string_to_com
from xyzfile_generation import xyz_from_com
from xyz2mol import read_xyz_file
from xyz2mol import xyz2mol
from rdkit.Chem import Draw
from rdkit import Chem
from rdkit.Chem import rdchem
from rdkit.Chem import rdmolfiles
from anion_generation import anion_from_com
from rdkit.Chem import AllChem
from com_files_from_smiles import write_indexfile

VEHICLe=pd.read_csv('VEHICLe.csv')

VEHICLe_1=VEHICLe[:1]

#write_indexfile(VEHICLe_1)

smiles=VEHICLe_1['Smiles']
regid = VEHICLe_1['Regid']
C_index_list=[]

for l in regid:
    comfiles = glob.glob('neutral_comfiles/' + str(l) + '_*')
    index_file=glob.glob('txt_files/' + str(l) + '_indexing.txt')

    for i in smiles:
        print(i)
        mol=Chem.MolFromSmiles(i)
        molH=Chem.AddHs(mol)
        molc=AllChem.Compute2DCoords(molH)
        a=molH.GetAtoms()

        #this bit isn't needed
        #for c in molH.GetConformers():
         #   index=[]
          #  atomicnum = []
           # coords = list(c.GetPositions())

        for j, atoms in enumerate(a):
            b=atoms.GetAtomicNum()
            d = atoms.GetIdx()
            if b == 1:
                e = atoms.GetNeighbors()
                for j, k in enumerate(e):
                    if k.GetAtomicNum() == 6:
                        if str(k.GetHybridization()) == 'SP2':
                            C_index = k.GetIdx()
                            C_index_list.append(C_index)

    C_numbering=[]
    for x in range(len(C_index_list)):
        C_numbering.append(x + 1)

    match_Cnum_Cindex = zip(C_index_list, C_numbering)
    ### Use same numbering systen in anion generation ###

    for file in comfiles:
        todays_date = str(date.today())
        split = file.split('\\')
        name = split[1].split('_')
        #print(name)
        for m in index_file:
            for n, o in match_Cnum_Cindex:
                    with open(m) as p:
                            for line in p:
                                comindex=line[2]
                                if str(n) == comindex:
                                    bromine_comfilename='bromine_comfiles/' + name[0] + '_' + str(o) + 'bromine_' + todays_date + '_' + name[2] + '-' + name[3] + '_opt.com'
                                    bromine_chkfilename= name[0] + '_' + str(o) + 'bromine_' + todays_date + '_' + name[2] + '-' + name[3] + '_opt.chk'
                                    split = line.split(' ')
                                    Br_xcoord = float(split[2])
                                    Br_ycoord = float(split[3])
                                    Br_zcoord = float(split[4])
                                    Br_coords=(str(Br_xcoord) +  ' ' + str(Br_ycoord) + ' ' + str(Br_zcoord))
                                    print(Br_coords)

                                with open(bromine_comfilename, 'w') as q:
                                    with open(file) as f:
                                        for comline in f:

                                            if comline.startswith('%chk='):
                                                print('%chk=' + bromine_chkfilename, file=q)

                                            elif comline == '0 1\n':
                                                print('1 1', file=q)

                                            elif Br_coords in comline:
                                                print(comline.strip('\n'), file=q)
                                                print('Br', Br_xcoord, Br_ycoord, (Br_zcoord + 2), file=q)


                                            else:
                                                print(comline.strip('\n'), file=q)




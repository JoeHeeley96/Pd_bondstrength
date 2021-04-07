from VEHICLe_read import read_mol_from_smiles
from xyzfile_generation import xyzcoords_from_type_array
from rdkit import Chem
from rdkit.Chem import AllChem
from datetime import date
import numpy as np


def write_gaussian_command_with_regid(Regid, chkfilename):
    with open('Gaussian_command.txt', 'w') as p:
        p.write('%chk=' + chkfilename)
        p.write('\n')
        p.write('# opt wb97xd/6-31g(d)')
        p.write('\n')
        p.write('\n')
        p.write(Regid)
        p.write('\n')
        p.write('\n')
        p.write('0 1')

def com_from_xyz_coords(xyz_coord_filename, comfilename):
    with open('Gaussian_command.txt') as fp:
        data = fp.read()

    with open(xyz_coord_filename) as fd:
        data2 = fd.read()

        data += '\n'
        data += data2

    with open(comfilename, 'w') as fp:
        fp.write(data)
        fp.write('\n')
        fp.write('\n')
        fp.write('\n')


def VEHICLe_string_to_com(dataframe):
    smiles = dataframe['Smiles']
    for a in smiles:
        #print(a)
        row=dataframe[dataframe.Smiles == a]
        Regid=row.Regid.item()
        todays_date=date.today()
        xyz_coord_filename ='txt_files/' + str(Regid) + '_xyzcoords.txt'
        comfilename ='neutral_comfiles/' + str(Regid) + '_' + str(todays_date) + '_wb97xd_631gd_opt.com'
        chkfilename= str(Regid) + '_' + str(todays_date) + '_wb97xd_631gd_opt.chk'
        print(comfilename)

        x = read_mol_from_smiles(a)
        y = Chem.AddHs(x)

        z = AllChem.Compute2DCoords(y)
        type_array = np.zeros(y.GetNumAtoms(), dtype=np.int32)

        for j, atoms in enumerate(y.GetAtoms()):
            type_array[j] = atoms.GetAtomicNum()

        for c in y.GetConformers():
            xyz = c.GetPositions()
            xyz_list = xyz.tolist()

        match_coords = zip(type_array, xyz_list)
        xyzcoords_from_type_array(xyz_coord_filename, match_coords)
        write_gaussian_command_with_regid(Regid, chkfilename)
        com_from_xyz_coords(xyz_coord_filename, comfilename)


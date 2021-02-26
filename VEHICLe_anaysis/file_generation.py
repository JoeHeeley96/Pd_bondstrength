from VEHICLe_read import read_mol_from_smiles
from rdkit import Chem
from rdkit.Chem import AllChem
from datetime import date
import numpy as np

def write_xyz_from_type_array(filename, zipped_type_array):
    Atoms = {1: 'H', 6: 'C', 7: 'N', 8: 'O', 16: 'S'}
    with open(filename, 'w') as f:
        for i in zipped_type_array:
            for atomicnum, atom in Atoms.items():
                if i[0] == atomicnum:
                    print(atom, *i[1], sep=' ', file=f)

def write_gaussian_command_with_regid(Regid, filename):
    with open('Gaussian_command.txt', 'w') as p:
        p.write('%chk=' + filename)
        p.write('\n')
        p.write('# opt wb97xd/6-31g(d)')
        p.write('\n')
        p.write('\n')
        p.write(Regid)
        p.write('\n')
        p.write('\n')
        p.write('0 1')

def VEHICLe_string_to_com(dataframe):
    smiles = dataframe['Smiles']
    for a in smiles:
        #print(a)
        row=dataframe[dataframe.Smiles == a]
        Regid=row.Regid.item()
        todays_date=date.today()
        xyzfilename ='xyz_files/' + str(Regid) + '.xyz'
        comfilename ='com_files/' + str(Regid) + '_' + str(todays_date) + '_wb97xd_631gd_opt.com'
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
        write_xyz_from_type_array(xyzfilename, match_coords)
        write_gaussian_command_with_regid(Regid, chkfilename)

        with open('Gaussian_command.txt') as fp:
            data=fp.read()

        with open(xyzfilename) as fd:
            data2=fd.read()

            data += '\n'
            data += data2

        with open(comfilename, 'w') as fp:
           fp.write(data)
           fp.write('\n')
           fp.write('\n')
           fp.write('\n')


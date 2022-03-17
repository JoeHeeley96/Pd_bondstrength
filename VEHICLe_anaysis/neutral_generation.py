from VEHICLe_read import read_mol_from_smiles
from xyzfile_generation import xyzcoords_from_type_array
from rdkit import Chem
from rdkit.Chem import AllChem
from datetime import date
import numpy as np
from tqdm import tqdm


def type_dict(warning=True):

    if warning:
        print('Type Dict includes: H, C, N, O, S ONLY')

    type_dict = {1 : 'H', 6: 'C', 7: 'N', 8: 'O', 16: 'S'}

    return type_dict

def write_gaussian_command_with_regid(Regid, chkfilename, modredundant=False):
    with open('Gaussian_command.txt', 'w') as p:
        p.write('%chk=' + chkfilename)
        p.write('\n')

        if modredundant:
            p.write('# opt=modredundant wb97xd/6-31g(d)')
        else:
            p.write('# opt wb97xd/6-31g(d)')

        p.write('\n')
        p.write('\n')
        p.write(Regid)
        p.write('\n')
        p.write('\n')
        p.write('0 1')

def com_from_xyz_coords(indexfilename, comfilename):

    with open(comfilename, 'w') as fn:
        with open('Gaussian_command.txt') as fp:
            with open(indexfilename) as fd:

                for gaussline in fp:
                    print(gaussline.strip('\n'), file=fn)

                for indexline in fd:
                    print(indexline[0].strip('\n'), indexline[4:].strip('\n').lstrip(), file=fn)


def write_indexfile(dataframe):
    smiles = dataframe['Smiles']
    for a in smiles:
        row = dataframe[dataframe.Smiles == a]
        Regid = row.Regid.item()
        indexfilename = 'txt_files/' + str(Regid) + '_indexing.txt'

        mol = read_mol_from_smiles(a)
        molH = Chem.AddHs(mol)
        AllChem.EmbedMolecule(molH)
        AllChem.MMFFOptimizeMolecule(molH)
        conf = molH.GetConformer()
        coordinates = conf.GetPositions()
        xyz_list = coordinates.tolist()

        coords_2d = AllChem.Compute2DCoords(molH)
        type_array = np.zeros(molH.GetNumAtoms(), dtype=np.int32)
        index_list = []

        for j, atoms in enumerate(molH.GetAtoms()):
            type_array[j] = atoms.GetAtomicNum()
            index = atoms.GetIdx()
            index_list.append(index)

        match_coords = zip(type_array, index_list, xyz_list)

        xyzcoords_from_type_array(indexfilename, match_coords)

def basehet_from_smiles(dataframe):

    print('PLEASE NOTE: If you are using neutral_generation.basehet_from_smiles outside of '
          'workflow.comfile_generation_workflow you need to append blank lines to all generated comfiles')

    regid = dataframe['Regid']
    todays_date = date.today()
    write_indexfile(dataframe)

    for i in tqdm(regid):

        indexfilename = 'txt_files/' + str(i) + '_indexing.txt'
        comfilename = 'comfiles/' + str(i) + '_' + str(todays_date) + '_wb97xd_631gd_opt.com'
        chkfilename = str(i) + '_' + str(todays_date) + '_wb97xd-631gd_opt.chk'

        write_gaussian_command_with_regid(i, chkfilename)
        com_from_xyz_coords(indexfilename, comfilename)





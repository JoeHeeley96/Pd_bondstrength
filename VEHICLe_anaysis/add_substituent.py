from neutral_generation import type_dict
from datetime import date
from random import randint
import numpy as np
from mol_translator.structure.structure_write import write_mol_tosdf
import pandas as pd
from neutral_generation import write_gaussian_command_with_regid
from xyzfile_generation import xyz_from_com
from mol_translator.aemol import aemol
from mol_translator.structure import structure_write as structwrt
from neutral_generation import write_indexfile
import rdkit.Chem as rd
from anion_generation import anions_from_smiles
from bromine_generation import bromines_from_smiles

def add_methyls(amol):

    types = type_dict(warning=False)
    mol_types = amol.structure['types']
    xyz = amol.structure['xyz']
    conn = amol.structure['conn']

    h_list = [x for x in mol_types if x == 1]

    h_to_change = randint((len(mol_types) - len(h_list) + 1), len(mol_types))
    mol_types[h_to_change - 1] = 6
    c_coords = xyz[h_to_change - 1]

    if c_coords[0] < 0:
        sum = -1

    else:
        sum = 1

    methyl_h1_coords = [c_coords[0], c_coords[1], c_coords[2] - sum]
    methyl_h2_coords = [c_coords[0], c_coords[1], c_coords[2] + sum]
    methyl_h3_coords = [c_coords[0], c_coords[1] + sum, c_coords[2]]

    extended_types = np.append(mol_types, [1, 1, 1])
    h_conn = np.zeros(len(extended_types), dtype=np.int)

    np.put(h_conn, [h_to_change - 1], [1])

    h_conn_list = h_conn.tolist()

    extended_conn = []
    extended_xyz = []

    for i, j in zip(conn.tolist(), xyz.tolist()):
        i.extend([0 ,0 ,0])
        extended_conn.append(i)

        extended_xyz.append(j)

    extended_conn.extend([h_conn_list, h_conn_list, h_conn_list])
    extended_xyz.extend([methyl_h1_coords, methyl_h2_coords, methyl_h3_coords])

    amol.structure.update({'size': len(extended_types), 'xyz': np.array(extended_xyz), 'types': extended_types,
                         'conn': np.array(extended_conn)})


    structwrt.write_mol_tosdf(amol, 'tmp.sdf')

    molblock = open('tmp.sdf', 'r').read()

    rdmol = rd.MolBlockTo(molblock)

    print(rdmol)

    #df = pd.DataFrame({'Regid': amol.info['molid'] + '-methylated', 'Smiles': rd.MolToSmiles(rdmol)})



''' x = amol.info['filepath'].split('\\')
y = x[2].split('_')
comfile = x[0] + '/' + 'comfiles/' + x[2].split('.')[0] + '.com'

methylated_comfile = x[0] + '\\' + 'comfiles\\' + y[0] + '-methylated_' + str(date.today()) + y[2] + y[3] + '.com'
methylated_chkfile = y[0] + '-methylated_' + str(date.today()) + y[2] + y[3] + '.chk'

write_gaussian_command_with_regid(Regid=y[0], chkfilename=methylated_chkfile, modredundant=False)

with open(methylated_comfile, 'w') as f:
    with open('Gaussian_command.txt') as fp:
            for gaussline in fp:
                print(gaussline.strip('\n'), file=f)

            for i, j in zip(mol_types, xyz):
                print(types[i], *j, file=f)

            print('H', *methyl_h1_coords, file=f)
            print('H', *methyl_h2_coords, file=f)
            print('H', *methyl_h3_coords, file=f)
            print('\n', file=f)
            print('\n', file=f)
            print('\n', file=f)

methylated_xyzfile = xyz_from_com(comfilename=methylated_comfile, save_loc='xyzfiles')

methylated_amol = aemol(amol.info['molid'] + '-methylated', filepath=methylated_xyzfile)
methylated_amol.from_file(file=methylated_xyzfile, ftype='xyz')
rdmol = methylated_amol.to_rdkit()
print(rdmol)

data = {'Regid': methylated_amol.info['molid'], 'Smiles': rd.MolToSmiles(rdmol)}
df = pd.DataFrame(data)

write_indexfile(dataframe=df)
anions_from_smiles(dataframe=df)
bromines_from_smiles(dataframe=df)'''
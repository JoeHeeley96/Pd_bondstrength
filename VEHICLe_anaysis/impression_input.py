import numpy as np
from mol_translator.imp_converter import dataframe_write as dfw
from mol_translator.aemol import aemol
from mol_translator.structure.structure_write import write_mol_tosdf
from tqdm import tqdm
import os

def logfile_to_aemol(logfiles, write=False, ftype='log'):

    array = np.unique(np.array(logfiles))
    amols = []

    for file in tqdm(list(array)):
        x = file.split('\\')
        p = x[1].split('_')[0]

        outfile = 'aemols/' + str(p) + '_AEMOL.sdf'

        amol = aemol(p)
        amol.from_file(file, ftype=ftype)
        amol.prop_fromfile(file, 'g16', 'scf')

        assert amol.mol_properties['energy'] < 1000000, print(amol.mol_properties)

        amol.get_bonds()
        amol.get_path_lengths()

        if amol not in amols:
            amols.append(amol)
        else:
            print('amol error', amol.info['molid'])
        if write:
            write_mol_tosdf(amol, outfile)

    return amols

def xyzfile_to_aemol(xyzfiles, write=False):

    amols = []
    for file in tqdm(xyzfiles):
        x = file.split('\\')
        p = x[1].split('_')[0]
        try:
            outfile = 'aemols/' + str(p) + '_AEMOL.sdf'
            if os.path.isfile(outfile):

                amol = aemol(p)
                amol.from_file(file, ftype='xyz')

                assert amol.mol_properties['energy'] < 1000000, print(amol.mol_properties)

                amol.get_bonds()
                amol.get_path_lengths()
        finally:
            amols.append(amol)
            if write:
                write_mol_tosdf(outfile, amol)
    return amols

def write_imp_input(aemols, fulldata_df, outname, write=True):
    atom_df = dfw.make_atom_df(aemols)
    pair_df = dfw.make_pair_df(aemols)

    regid = list(set(atom_df.molecule_name.unique()))

    props=np.zeros(len(atom_df))

    for j in regid:
        mol_prop = fulldata_df[fulldata_df['Regid'].str.fullmatch(j)]
        acidities = mol_prop.filter(regex='acidity').values
        elec_affs = mol_prop.filter(regex='electrophile_affinity').values
        mol_df = atom_df[atom_df['molecule_name'].str.fullmatch(j)]
        catom_index = mol_df[mol_df['typestr'].str.fullmatch('C')].atom_index
        hatom_index = mol_df[mol_df['typestr'].str.fullmatch('H')].atom_index
        # because of my stupid numbering I had to identify the non-quaternary carbons using the connectivity array
        # Did this by finding the positions of the H's and seeing if the atom in question was attached
        carbons = []
        hydrogens = []

        for ind in catom_index.values:
            c_conn=mol_df.iloc[ind].conn
            for i in hatom_index.values:
                c = str(c_conn).split(' ')
                if '1' in c[i]:
                    carbons.append(ind)
                    hydrogens.append(i)

        for a in acidities:
                for u,v, in zip(hydrogens, a):
                    hdf_index = mol_df.loc[mol_df['atom_index'] == u].index.values
                    np.put(props, [hdf_index], [v])

        for e in elec_affs:
                for x,y in zip(carbons, e):
                    cdf_index = mol_df.loc[mol_df['atom_index'] == x].index.values
                    np.put(props, [cdf_index], [y])

    atom_df['shift'] = props
    final_atom_df = atom_df.fillna(0)

    if write:

        with open(outname + '_atom_df.csv', 'w') as p:
            print(final_atom_df.to_csv(sep=','), file=p)

        with open(outname + '_pair_df.csv', 'w') as f:
            print(pair_df.to_csv(sep=','), file=f)

    return final_atom_df, pair_df

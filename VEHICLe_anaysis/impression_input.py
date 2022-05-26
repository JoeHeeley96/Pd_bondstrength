import numpy as np
import pandas as pd
from mol_translator.imp_converter import dataframe_write as dfw
from mol_translator.aemol import aemol
from mol_translator.structure.structure_write import write_mol_tosdf
from tqdm import tqdm
import os

def file_to_aemol(logfiles, write=False, ftype='log'):

    array = np.unique(np.array(logfiles))
    amols = []

    for file in tqdm(list(array)):
        x = file.split('\\')
        p = x[1].split('_')[0]

        outfile = 'aemols/' + str(p) + '_AEMOL.sdf'

        amol = aemol(p, filepath=file)
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
            connected_hs = []
            c_conn = mol_df.iloc[ind].conn
            for i in hatom_index.values:
                c = str(c_conn).split(' ')
                if '1' in c[i]:
                    connected_hs.append(i)

            if len(connected_hs) == 1:
                carbons.append(ind)
                hydrogens.extend(connected_hs)

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

def classifier_input(vehicle_dataframe, calculate_properties_dataframe, outname, djr_lim=1.2, write=True):

    '''
    Writes a dataframe for the input of a classification machine. Will set acidic structures to 1 and nucleophilic
    structures to 0 depending on the djr limit given

    :param vehicle_dataframe: The VEHICLe dataframe - must contain same structures as calculate_properties_dataframe
    :param calculate_properties_dataframe: Output from calculate_properties
    :param outname: Name to save the output dataframe under
    :param write: Bool: saves if True
    :param djr_lim: limit above which acidic activation is seen (experimentally determined)
    :return:
    '''

    regid = calculate_properties_dataframe['Regid']

    class_dataframe = pd.DataFrame(columns=['Regid', 'Smiles', 'djr', 'Activation'])

    for i in regid:

        acidities = {}
        elec_affs = {}
        relative_acidities = []
        relative_electro_affs = []

        df = calculate_properties_dataframe.loc[calculate_properties_dataframe['Regid'] == i]
        for j in df.columns:
            if 'anion' in j:
                acidities[j] = float(df[j].values)
            elif 'bromine' in j:
                elec_affs[j] = float(df[j].values)

        for m, c in acidities.items():
            if c == (min(acidities.values())):
                relative_acidities.append(
                    24.044391262487466 / c) if 24.044391262487466 / c not in relative_acidities else relative_acidities

                if len(relative_acidities) == 1:
                    break

        for f, l in elec_affs.items():
            if l == max(elec_affs.values()):
                relative_electro_affs.append(
                    l / 183.64724482984454) if l / 183.64724482984454 not in relative_electro_affs else relative_electro_affs

        if len(relative_acidities) != len(relative_electro_affs):
            print('There is a problem with: ', i, 'please check!')

        djr = round(relative_acidities[0] / relative_electro_affs[0], 2)

        if djr > djr_lim:
            activation = 1

        else:
            activation = 0

        data = {'Regid': i, 'Smiles': vehicle_dataframe[vehicle_dataframe['Regid'] == i].Smiles.values[0], 'djr': djr,
                'Activation': activation}

        class_dataframe = class_dataframe.append(data, ignore_index=True)

    with open(outname, 'w') as f:
        print(class_dataframe.to_csv(sep=','), file=f)

def exclude_structures(structures_to_exclude, atom_df, pair_df, outname, write=True):

    excluded_atom = atom_df.drop(atom_df.index[atom_df['molecule_name'].isin(structures_to_exclude)])
    excluded_pair = pair_df.drop(pair_df.index[pair_df['molecule_name'].isin(structures_to_exclude)])

    print('Removed:', len([x for x in atom_df['molecule_name'].unique() if x not in excluded_atom['molecule_name'].unique()]), 'structures')

    if write:
        with open(outname + '_atom_df.csv', 'w') as f:
            print(excluded_atom.to_csv(), file=f)

        with open(outname + '_pair_df.csv', 'w') as p:
            print(excluded_pair.to_csv(), file=p)

    return excluded_atom, excluded_pair

#def merge_datasets(dataset1_atom, dataset1_pair, dataset2_atom, dataset2_pair):


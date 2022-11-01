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

def classifier_input(vehicle_dataframe, relative_properties_dataframe, outname,
                     djr_lim=1.21, acid_lim=0.67, nuc_lim=0.84, write=True):

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

    regid = relative_properties_dataframe['Regid']

    class_dataframe = pd.DataFrame(columns=['Regid', 'Smiles', 'djr', 'Activation'], index=[0])

    for i in regid:

        acidities = {}
        elec_affs = {}

        df = relative_properties_dataframe.loc[relative_properties_dataframe['Regid'] == i]

        for j in df.columns:
            if 'anion' in j:
                if df[j].values != 0:
                    acidities[j] = float(df[j].values)

            elif 'bromine' in j:
                elec_affs[j] = float(df[j].values)

        most_acidic = round(max(acidities.values()), 2)
        most_nuc = round(max(elec_affs.values()), 2)

        djr = round(most_acidic / most_nuc, 2)

        activation = -404.404

        if most_acidic < acid_lim and most_nuc < nuc_lim:
            activation = 0

        elif most_acidic > acid_lim and most_nuc < nuc_lim:
            activation = 1

        elif most_acidic > acid_lim and most_nuc > nuc_lim:
            if djr > djr_lim:
                activation = 1

            else:
                activation = 2

        else:
            activation = 2

        if activation == -404.404:
            print(i, most_acidic, most_nuc, djr)

        data = pd.DataFrame({'Regid': i, 'Smiles': vehicle_dataframe[vehicle_dataframe['Regid'] == i].Smiles.values[0], 'djr': djr,
                'Activation': activation}, index=[0])

        class_dataframe = pd.concat([class_dataframe, data])

    with open(outname +'.csv', 'w') as f:
        print(class_dataframe.to_csv(sep=','), file=f)

def exclude_structures_imp_dfs(structures_to_exclude, atom_df, pair_df, outname, write=True):

    excluded_atom = atom_df.drop(atom_df.index[atom_df['molecule_name'].isin(structures_to_exclude)])
    excluded_pair = pair_df.drop(pair_df.index[pair_df['molecule_name'].isin(structures_to_exclude)])

    print('Removed:', len([x for x in atom_df['molecule_name'].unique() if x not in excluded_atom['molecule_name'].unique()]), 'structures')

    if write:
        with open(outname + '_atom_df.csv', 'w') as f:
            print(excluded_atom.to_csv(), file=f)

        with open(outname + '_pair_df.csv', 'w') as p:
            print(excluded_pair.to_csv(), file=p)

    return excluded_atom, excluded_pair

def exclude_structures_prop_dfs(structures_to_exclude, df1, df2, write=True):

    excluded_df1 = df1.drop(df1.index[df1['Regid'].isin(structures_to_exclude)])
    excluded_df2 = df2.drop(df2.index[df2['Regid'].isin(structures_to_exclude)])

    print('Removed:',
          len([x for x in df1['Regid'].unique() if x not in excluded_df1['Regid'].unique()]),
          'structures from df1')

    print('Removed:',
          len([x for x in df2['Regid'].unique() if x not in excluded_df2['Regid'].unique()]),
          'structures from df2')

    if write:
        print('you need to write the code to save the excluded dfs ya dingus')

    return excluded_df1, excluded_df2

def merge_datasets(dataset1_atom, dataset1_pair, dataset2_atom, dataset2_pair, outname, write=True):

    #use the exclude_structures function to remove any overlapping structures between the two input dataframes
    # before concatenation
    structures_to_exclude = dataset1_atom['molecule_name'].unique()

    unique_dfs = exclude_structures_imp_dfs(structures_to_exclude, dataset2_atom, dataset2_pair, outname=None, write=False)

    join_atom = pd.concat([dataset1_atom, unique_dfs[0]], axis=0)
    join_pair = pd.concat([dataset1_pair, unique_dfs[1]], axis=0)

    if write:
        with open(outname + '_atom_df.csv', 'w') as f:
            print(join_atom.to_csv(), file=f)

        with open(outname + '_pair_df.csv', 'w') as p:
            print(join_pair.to_csv(), file=p)

    return join_atom, join_pair

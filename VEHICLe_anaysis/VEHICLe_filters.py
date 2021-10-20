import numpy as np
from rdkit import Chem
import pandas as pd
from VEHICLe_read import read_list_of_mol_objects
from print_file import print_to_csv
from rdkit.Chem import Draw

def nonetype_filter(dataframe):
    for i in read_list_of_mol_objects(dataframe):
        print(i.GetNumAtoms())

def synthesis_filter(dataframe, export):
    synthesis_filtered=dataframe[dataframe.Pgood!=0]

    if export == 'yes':
        print_to_csv(synthesis_filtered, 'Synthesis_filter')
    return synthesis_filtered

def beilstein_filter(dataframe, export):
    beilstein_filtered =dataframe[dataframe.Beilstein_hits_June_08!=0]
    if export == 'yes':
        print_to_csv(beilstein_filtered, 'Beilstein_filter')

    return beilstein_filtered

def dj1_filter(dataframe, export):
    '''
    apply this to the VEHICLe dataframe to get the dj1 dataset
    '''

    filtered_data = pd.DataFrame({})

    for i in dataframe:
        filtered_data[i] = i

    for j in range(len(dataframe)):
        row = dataframe.iloc[j]
        smiles = row['Smiles']
        mol = Chem.MolFromSmiles(smiles)
        type_array = np.zeros(mol.GetNumAtoms(), dtype=np.int32)

        for x, atoms in enumerate(mol.GetAtoms()):
            type_array[x] = atoms.GetAtomicNum()

        if np.count_nonzero(type_array == 6, axis=0) > 0:
            if np.count_nonzero(type_array == 7, axis=0) <= 4:
                if np.count_nonzero(type_array == 8, axis=0) == 0:
                    if np.count_nonzero(type_array == 16, axis=0) == 0:
                        filtered_data = filtered_data.append(row, ignore_index=True)

    if export == 'yes':
        print_to_csv(filtered_data, 'Nitrogen_only')

def dj2_filter(dataframe, write=True):

    filtered_data = pd.DataFrame({})

    for i in dataframe:
        filtered_data[i] = i

    for j in range(len(dataframe)):
        row = dataframe.iloc[j]
        smiles = row['Smiles']
        mol = Chem.MolFromSmiles(smiles)
        Hmol = Chem.AddHs(mol)
        type_array = np.zeros(Hmol.GetNumAtoms(), dtype=np.int32)

        for j, atoms in enumerate(Hmol.GetAtoms()):
            type_array[j] = atoms.GetAtomicNum()

        if np.count_nonzero(type_array == 6, axis=0) >0:
            if np.count_nonzero(type_array == 1, axis=0) >0:
                if np.count_nonzero(type_array == 7, axis=0) <=4:
                    if np.count_nonzero(type_array == 8, axis=0) <=4:
                        if np.count_nonzero(type_array == 16, axis=0) <=4:
                            filtered_data = filtered_data.append(row, ignore_index=True)

    if write:
        with open('dj2.csv', 'w') as f:
            print(filtered_data.to_csv(sep=','), file=f)

    return filtered_data



def NCN_filter(dataframe, export):
    ### This filter extracts molecules with C-C cores.
    filtered_data = pd.DataFrame({})

    for i in dataframe:
        filtered_data[i] = i

    for j in range(len(dataframe)):
        row = dataframe.iloc[j]
        smiles = row['Smiles']
        mol = Chem.MolFromSmiles(smiles)
        type_array = np.zeros(mol.GetNumAtoms(), dtype=np.int32)

        for j, atoms in enumerate(mol.GetAtoms()):
            if atoms.GetAtomicNum() == 6:
                if len(atoms.GetNeighbors()) == 3:
                    for i in atoms.GetNeighbors():
                        if i.GetAtomicNum() == 6: ### Change condiction to == 7 for CN cores ###
                            if len(i.GetNeighbors()) == 3:
                                type_array = np.zeros(len(atoms.GetNeighbors()), dtype=np.int32)

                                for x, y in enumerate(atoms.GetNeighbors()):
                                    type_array[x] = y.GetAtomicNum()

        if np.count_nonzero(type_array == 7, axis=0) >= 2: ### Change condition to >= 3 for C-N core ###
            filtered_data = filtered_data.append(row, ignore_index=True)


    if export == 'yes':
        print_to_csv(filtered_data, 'NCN')

def fulldata_filter(VEHICLe_dataframe, fulldata_dataframe, outname):
    regid = list(set(fulldata_dataframe.Regid))
    structures=[]
    for i in regid:
        acidities = {}
        elec_affs = {}
        relative_acidities = []
        relative_electro_affs = []

        df = fulldata_dataframe.loc[fulldata_dataframe['Regid'] == i]
        for j in df.columns:
            if 'anion' in j:
                acidities[j] = float(df[j].values)
            elif 'bromine' in j:
                elec_affs[j] = float(df[j].values)

        for m, c in acidities.items():
            if c == (min(acidities.values())):
                relative_acidities.append(
                    24.038503 / c) if 24.038503 / c not in relative_acidities else relative_acidities
                acidic_position = list(m.split('_')[0])

                if len(relative_acidities) == 1:
                    break

        for f, l in elec_affs.items():
            if l == max(elec_affs.values()):
                relative_electro_affs.append(
                    l / 183.651714) if l / 183.651714 not in relative_electro_affs else relative_electro_affs
                nucleophilic_position = list(f.split('_')[0])

        if len(relative_acidities) != len(relative_electro_affs):
            print('There is a problem with: ', i, 'please check!')

        #if (relative_acidities[0] < 1.0) and 0.5 <= relative_electro_affs[0] <= 0.7:
        print(i, 'acidic position is:', acidic_position[0], 'with relative acidity:', relative_acidities,
                                      'nucleophilic position is:', nucleophilic_position[0],
                                            'with relative electrophile affinity:', relative_electro_affs)
        structures.append(i)

    mol_list = []
    for i in structures:
        row = VEHICLe_dataframe[VEHICLe_dataframe['Regid'].str.fullmatch(i)]
        for j in row.Smiles.values:
            mol = Chem.MolFromSmiles(j)
            mol_list.append(mol)
        img = Draw.MolsToGridImage(mol_list, molsPerRow=5, legends=[x for x in structures])
        img.save(outname + '.png')




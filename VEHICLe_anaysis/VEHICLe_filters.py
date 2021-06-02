import numpy as np
from rdkit import Chem
import pandas as pd
from VEHICLe_read import read_list_of_mol_objects
from print_file import print_to_csv

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

def atom_filter(dataframe, export):
    filtered_data=pd.DataFrame({})

    for i in dataframe:
        filtered_data[i] = i

    for j in range(len(dataframe)):
        row = dataframe.iloc[j]
        smiles = row['Smiles']
        mol = Chem.MolFromSmiles(smiles)
        type_array = np.zeros(mol.GetNumAtoms(), dtype=np.int32)

        for j, atoms in enumerate(mol.GetAtoms()):
            type_array[j] = atoms.GetAtomicNum()

        if np.count_nonzero(type_array == 6, axis=0) >0:
            if np.count_nonzero(type_array == 7, axis=0) <=4:
                if np.count_nonzero(type_array == 8, axis=0) <=4:
                    if np.count_nonzero(type_array == 16, axis=0) <=4:
                        filtered_data = filtered_data.append(row, ignore_index=True)

    if export == 'yes':
        print_to_csv(filtered_data, 'Atom_filter')

    return filtered_data

def nitrogen_only_filter(dataframe, export):
    filtered_data = pd.DataFrame({})

    for i in dataframe:
        filtered_data[i] = i

    for j in range(len(dataframe)):
        row = dataframe.iloc[j]
        smiles = row['Smiles']
        mol = Chem.MolFromSmiles(smiles)
        type_array = np.zeros(mol.GetNumAtoms(), dtype=np.int32)

        for j, atoms in enumerate(mol.GetAtoms()):
            type_array[j] = atoms.GetAtomicNum()

        if np.count_nonzero(type_array == 6, axis=0) > 0:
            if np.count_nonzero(type_array == 7, axis=0) <= 4:
                if np.count_nonzero(type_array == 8, axis=0) == 0:
                    if np.count_nonzero(type_array == 16, axis=0) == 0:
                        filtered_data = filtered_data.append(row, ignore_index=True)

    if export == 'yes':
        print_to_csv(filtered_data, 'Nitrogen_only')


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

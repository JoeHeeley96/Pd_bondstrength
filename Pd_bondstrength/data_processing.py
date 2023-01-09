import pandas as pd

def assign_directed_activation_to_atom_df(atom_df, directed_activation_df):

    assigned_atom_df = pd.DataFrame()

    for id in atom_df['molecule_name'].unique():

        idx = []
        ncs = []

        mol_atom = atom_df[atom_df['molecule_name'] == id]
        mol_directed = directed_activation_df[directed_activation_df['Structure'].str.endswith(id)]

        for atom in mol_directed.Structure.values:
            idx.append(int(atom.split('-')[0]) - 34)
            ncs.append(mol_directed[mol_directed['Structure'] == atom].iloc[0]['ExchG_kcalpermol'])

        for a, b in zip(idx, ncs):
            mol_atom.loc[mol_atom['atom_index'] == a, 'shift'] = b

        assigned_atom_df = pd.concat([assigned_atom_df, mol_atom])

    with open('NCS_atom_df.csv', 'w') as f:
        print(assigned_atom_df.to_csv(), file=f)

    return assigned_atom_df



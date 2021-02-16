from rdkit import Chem
from rdkit.Chem import Draw

def molecule_write_image(VEHICLe_as_dataframe, desired_filename_png):
    list_of_smiles = list(VEHICLe_as_dataframe['Smiles'])
    list_of_mol_objects = []

    for i in list_of_smiles:
        mol_object = Chem.MolFromSmiles(i)
        list_of_mol_objects.append(mol_object)

    image=Draw.MolsToGridImage(list_of_mol_objects)
    image.save(desired_filename_png)
### NEED TO FIX THE WRITING SECTION OF THIS###
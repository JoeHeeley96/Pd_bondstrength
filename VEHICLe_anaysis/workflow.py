from neutral_generation import VEHICLe_string_to_com
from anion_generation import anions_from_smiles
from bromine_generation import bromines_from_smiles

def comfile_workflow(dataframe):
    VEHICLe_string_to_com(dataframe)
    anions_from_smiles(datframe)
    bromines_from_smiles(datframe)
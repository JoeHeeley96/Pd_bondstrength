from bromine_generation import bromines_from_smiles
from anion_generation import anions_from_smiles
from com_files_from_smiles import VEHICLe_string_to_com
from com_files_from_smiles import write_indexfile


def vehicle_workflow(dataframe):

    VEHICLe_string_to_com(dataframe)
    write_indexfile(dataframe)
    anions_from_smiles(dataframe)
    bromines_from_smiles(dataframe)

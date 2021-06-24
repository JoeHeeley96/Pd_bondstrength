import glob
import pandas as pd
from neutral_generation import VEHICLe_string_to_com
from anion_generation import anions_from_smiles
from bromine_generation import bromines_from_smiles
from logfile_read import energy_readdata
from calculation_check import opt_check
from property_calculate import calculate_properties
from plots import plot_activation_map


def comfile_generation_workflow(dataframe):
    VEHICLe_string_to_com(dataframe)
    anions_from_smiles(datframe)
    bromines_from_smiles(datframe)

def property_calculate_workflow(filelocations, calc_check, plot_map):
    if calc_check != 'no':
        files = glob.glob(filelocations + '/*')
        if opt_check(i, 'SCF Done'):
            pass
        else:
            print('Error! SCF not found in file: ', i)

    energy_readdata(filelocations, 'Nonly_rawdata.csv')
    calculate_properties(pd.read_csv('Nonly_rawdata.csv'), 'VEHICLe_NOnly_fulldata.csv')

    if plot_map != 'no':
        plot_activation_map(pd.read_csv('VEHICLe_Nonly_fulldata.csv'), 'Nonly_activation_map.png')
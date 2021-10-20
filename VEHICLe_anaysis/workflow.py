import glob
import pandas as pd
from neutral_generation import VEHICLe_string_to_com
from anion_generation import anions_from_smiles
from bromine_generation import bromines_from_smiles
from logfile_read import energy_readdata
from calculation_check import opt_check
from property_calculate import calculate_properties
from plots import plot_activation_map
from impression_input import logfile_to_aemol
from impression_input import write_imp_input


def comfile_generation_workflow(dataframe):
    VEHICLe_string_to_com(dataframe)
    print('do you need to write the index files here?')
    anions_from_smiles(dataframe)
    bromines_from_smiles(dataframe)

def property_calculate_workflow(filelocations, outname, calc_check=True, plot_map=False, write=True):

    if calc_check:
        files = glob.glob(filelocations + '/*')
        for i in files:
            if opt_check(i, 'SCF Done'):
                pass
            else:
                print('Error! SCF not found in file: ', i)

    energy_readdata(filelocations, outname=(outname + '_rawdata.csv'), write=write)
    calculate_properties(pd.read_csv(outname + '_rawdata.csv'), outname=(outname + '_fulldata.csv'), write=write)

    if plot_map:
        plot_activation_map(pd.read_csv(outname + '_fulldata.csv'), outname + '_activation_map.png')


def impression_input_workflow(logfiles, fulldata_df, outname, write=True):

    print('aemols are written by default but not saved! To change please edit the '
          'impression_input_workflow function in workflow.py')

    aemols = logfile_to_aemol(logfiles, write=False)
    imp_input = write_imp_input(aemols=aemols, fulldata_df=fulldata_df, outname=outname, write=write)

    return imp_input[0], imp_input[1]


import glob
import pandas as pd
from neutral_generation import basehet_from_smiles
from anion_generation import anions_from_smiles
from bromine_generation import bromines_from_smiles
from logfile_read import energy_readdata
from calculation_check import opt_check, output_structure_check
from property_calculate import calculate_properties
from property_calculate import calculate_relative_properties
from plots import plot_activation_map
from impression_input import file_to_aemol, write_imp_input
from sampling import get_tm_df, get_fps_ids
from calculation_check import comfile_check
from tqdm import tqdm

def comfile_generation_workflow(dataframe):
    basehet_from_smiles(dataframe, workflow=True)
    anions_from_smiles(dataframe, workflow=True)
    bromines_from_smiles(dataframe, workflow=True)

    comfiles = glob.glob('comfiles/*')
    regids = dataframe['Regid']

    for i in comfiles:
        with open(i, 'a') as f:
            f.write('\n\n')
            f.close()

    comfile_check(comfiles=comfiles, sampled_structures=regids)

def logfile_analysis_workflow(logfilelocation, xyzfilelocation, outname, calc_check=True, calc_data=True, calc_rel_data=True,
                              plot_map=False, write=True):

    if calc_check:
        failed_runs = []
        files = glob.glob(logfilelocation + '/*')
        xyzfiles = glob.glob(xyzfilelocation + '/*')
        for i in tqdm(files):
            if opt_check(i, 'Error'):
                failed_runs.append(i)

        bad_structures = output_structure_check(xyzfiles=xyzfiles, logfiles=files)
        print('Found', len(failed_runs), 'failed structures:', failed_runs)

        with open(outname + '_failed_structures.txt', 'w') as f:
            for i in failed_runs:
                print(i, file=f)

        with open(outname + '_bad_structures.txt', 'w') as p:
            for j in [x for x in bad_structures if x not in failed_runs]:
                print(j, file=p)

    if calc_data:
        read_data = energy_readdata(logfilelocation, outname=(outname + '_rawdata.csv'), write=write)
        calculate_properties(pd.read_csv(outname + '_rawdata.csv'), outname=(outname + '_fulldata.csv'), write=write)

        if calc_rel_data:
            calculate_properties_df = pd.read_csv(outname + '_fulldata.csv')
            calculate_relative_properties(calculate_properties_dataframe=calculate_properties_df,
                                          outname= outname + '_fulldata_relative.csv', write=True)

    if plot_map:
        plot_activation_map(pd.read_csv(outname + '_fulldata.csv'), outname + '_activation_map.png')


def impression_input_workflow(neutral_logfiles, fulldata_df, outname, write=True):
    '''

    :param logfiles: should be the neutral logfiles for the structures you want to use
    :param fulldata_df: the output of calculate_properties or logfile_analysis_workflow
    :param outname: What you want to name the output atom and pair dataframes
    :param write: Bool, will write the dataframes to .csv files if True
    :return: atom_df, pair_df for structures you inputted at the start
    '''

    aemols = file_to_aemol(neutral_logfiles, write=False)
    imp_input = write_imp_input(aemols=aemols, fulldata_df=fulldata_df, outname=outname, write=write)

    return imp_input[0], imp_input[1]


def logfile_fingerprint_sampling_workflow(neutral_logfiles, training_structures, vehicle_df, outname,
                                  fulldata_df, num=100, imp_input=True, write=True):
    '''
    :param neutral_logfiles: list of logfiles that you want to sample
    :param training_structures: vehicle regid's of the structures you want to sample from
    :param vehicle_df: the full vehicledatabase as a dataframe
    :param outname: what you want to save your files as
    :param fulldata_df: the output of the properties_calculate.calculate_properties function
    :param num: the number of structures you want to sample
    :param imp_input: boolean, writes the input atom and pari dataframes for IMPRESSION Gen I if True
    :param write: boolean, saves all generated files if True

    :return: a list of regids that correspond to the sampled structures
    '''

    tm_df_array = get_tm_df(training_structures=training_structures, vehicle_df=vehicle_df)
    fpsample = get_fps_ids(tm_df=tm_df_array[0], tm_array=tm_df_array[1], num=num, outname=outname, write=write)

    sample_structures = []
    for i in neutral_logfiles:
        name = i.split('\\')[1].split('_')[0]
        if name in fpsample:
            sample_structures.append(i)

    if imp_input:
        impression_input_workflow(logfiles=sample_structures, fulldata_df=fulldata_df, outname=outname + str(num),
                                    write=write)

    return fpsample


def vehicle_fingerprint_sampling_workflow(vehicle_dataframe, selectdf, outname, num=100, write=True):

    sampled_dataframe = pd.DataFrame([])
    regids = vehicle_dataframe['Regid']
    tm_df_array = get_tm_df(training_structures=regids, vehicle_df=vehicle_dataframe)
    fpsample = get_fps_ids(tm_df=tm_df_array[0], tm_array=tm_df_array[1], num=num, outname=outname, write=write)

    for i in fpsample:
        row = vehicle_dataframe[vehicle_dataframe.Regid == i]
        sampled_dataframe = sampled_dataframe.append(row)

    if write:
        with open(outname + str(num) +'.csv', 'w') as f:
            print(sampled_dataframe.to_csv(sep=','), file=f)

    return sampled_dataframe

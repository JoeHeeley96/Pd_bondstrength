import pandas as pd
import glob
from neutral_generation import write_gaussian_command_with_regid
from calculation_check import opt_check
from file_generation import type_dict
from mol_translator import aemol

def g16_scfread(file):
    energy = -404
    with open(file, 'r') as f:
        for line in f:
            if 'SCF Done' in line:
                items = line.split()
                energy = float(items[4])

    return energy

def orca_scfread(file):
    energy = -404
    with open(file, 'r') as f:
        for line in f:
            if 'FINAL SINGLE POINT ENERGY' in line:
                items = line.split()
                energy = float(items[-1])

    return energy

def orcaxyz_to_comfile(file):

    structure = file.split('.')[0]

    comfilename = structure + '.gjf'
    write_gaussian_command_with_regid(structure, 'tmp.chk')

    with open(comfilename, 'w') as fn:
        with open(file, 'r') as f:
            lines = f.readlines()[2:]
            with open('Gaussian_command.txt') as fp:

                    for gaussline in fp:
                        print(gaussline.strip('\n'), file=fn)

                    for xyzline in lines:
                        print(xyzline.strip('\n'), file=fn)


def energy_readdata(file_locations, outname, write=True):

        energies = []
        date=[]
        structure=[]
        RegId = []
        theory=[]
        calc=[]

        list_of_files = glob.glob(file_locations + '/*')
        for file in list_of_files:
            energy = g16_scfread(file)
            if energy != -404:
                energies.append(energy)
                x = file.split('\\')
                z = x[1].split('_')
                RegId.append(z[0])
                theory.append(z[3])
                calc.append(z[-1])

                if 'anion' in z[1]:
                    structure.append(z[1])
                    date.append(z[2])
                elif 'bromine' in z[1]:
                    structure.append(z[1])
                    date.append(z[2])

                else:
                    structure.append('Base_het')
                    date.append(z[1])

            else:
                print(file, 'Energy not found')

        dict = {'Regid': RegId, 'Structure': structure, 'Energy': energies, 'Date': date, 'Functional/basisset': theory, 'Calculation': calc}
        data = pd.DataFrame(dict)

        if write:
            with open(outname, 'w') as f:
                print(data.to_csv(sep=','), file=f)

        return data

def read_freq_calcs(freq_logfile_location):

    files = glob.glob(freq_logfile_location + '/*.log')
    failed_runs = []
    imaginary_freqs= []
    imaginary_files = []

    for i in files:
        freq = []
        if opt_check(i, 'Error'):
            failed_runs.append(i)

        with open(i, 'r') as f:
            for line in f:
                if 'Frequencies' in line:
                    for item in [-1, -2, -3]:
                        freq.append(float(line.split()[item]))

        negatives = [x for x in freq if x < 0]

        if len(negatives) != 0:
            imaginary_freqs.append(negatives)
            print(i)
            imaginary_files.append(i)

    print('Found', len(failed_runs), 'Failed runs')
    print('Found', len(imaginary_freqs), 'Imaginary Frequencies')

def comfile_to_aemol(comfile):

    file = open(comfile)
    typ_dict = type_dict()

    lines = file.readlines()

    id = lines[3]

    types = []
    xyz = []

    xyz = lines[6:]

    #for i in xyz:


import numpy as np
from mol_translator import aemol
from mol_translator.properties.structure.bond_angle import get_bond_angles
from tqdm import tqdm
import glob
import os
from logfile_read import comfile_to_aemol
from logfile_read import opt_check

def comfile_check(comfiles, sampled_structures):

    print('NEED TO FIX COMFILE CHECK')
    for i in tqdm(sampled_structures):
        molcom_list = [w for w in comfiles if i + '_' in w]

        anion = [x for x in molcom_list if 'anion' in x]
        bromine = [z for z in molcom_list if 'bromine' in z]
        neutral = [y for y in bromine if 'anion' not in y]

        with open(neutral[0], 'r') as f:
            nlines = [l for l in (line.strip() for line in f) if l]
            ncoords = nlines[4:]


#                    if len(alines) != (len(neutralxyz) - 1):
#                        print(alines)

#                        if len(blines) != (len(neutralxyz) + 1):
#                            print(blines)

def get_dist_array(aemol):
    num_atoms = len(aemol.structure['types'])
    dist_array = np.zeros(num_atoms, num_atoms, dtype=np.float64)
    xyz_array = aemol.structure['xyz']
    for i in range(num_atoms):
        for j in range(num_atoms):
            dist_array[i][j] = np.absolute(np.linalg.norm(xyz_array[i] - xyz_array[j]))

    return dist_array

def output_structure_check(comfiles, logfiles):
    '''
    This function checks the MAE between bond angles
    for the input and output structures.

    Typically, successful anions range between 0 - 12 (approx)
    and successful bromines range between 0 - 30 (approx)
    so any structures with MAE >= 12 or 30 are thrown up

    PLEASE CHECK THE STRUCTURES BEFORE DELETING THEM

    If this doesnt work then just change the condition:
    if np.mean(np.absolute(m-c)) > *set new number here*:
    '''

    badfiles = []
    for i in tqdm(logfiles):
        out_file = i.split('\\')
        out_split = out_file[2].split('_')

        for j in comfiles:
            in_file = j.split('\\')
            in_split = in_file[2].split('_')

            if in_split[0] + '_' + in_split[1] == out_split[0] + '_' + out_split[1]:
                input_aemol = comfile_to_aemol(j)

                output_aemol = aemol(in_file[1].split('_')[0])
                output_aemol.from_file(i, ftype='log')

                input_structure = get_bond_angles(input_aemol)
                output_structure = get_bond_angles(output_aemol)

                err = []
                for x, y in zip(input_structure, output_structure):
                    for m, c in zip(x, y):
                        mnan = np.asarray([0 if x != x else x for x in m])
                        cnan = np.asarray([0 if x != x else x for x in c])

                        if 'anion' in i:
                            if np.mean(np.absolute(mnan-cnan)) > 12.0:
                                err.append(np.mean(np.absolute(mnan-cnan)))
                        elif 'bromine' in i:
                            if np.mean(np.absolute(mnan-cnan)) > 32.0:
                                err.append(np.mean(np.absolute(mnan-cnan)))
                        else:
                            if np.mean(np.absolute(mnan - cnan)) > 10.0:
                                err.append(np.mean(np.absolute(mnan - cnan)))

                if len(err) > 0:
                    badfiles.append(i)

    print('Large variances detected in ', len(badfiles), ' output structures:', badfiles)

    return badfiles

def quick_dj1_sub_check():

    with open('quick_dj1.2_check.txt', 'w') as f:
        for i in glob.glob('comfiles/*'):

            if os.path.isdir:
                logs = glob.glob(i + '/*.log')
                coms = glob.glob(i + '/*.com')

                failed_runs = []

                for i in logs:
                    if opt_check(i, 'Error'):
                        failed_runs.append(i)

                bad_files = output_structure_check(coms, logs)
                print('Found:', len(failed_runs), 'Failed Runs')

                print(i.split('\\')[1], file=f)
                print('FAILED STRUCTURES', file=f)
                for i in failed_runs:
                    print(i.split('\\')[-1], file=f)
                print('LARGE VARIENCES', file=f)
                for i in bad_files:
                    print(i.split('\\')[-1], file=f)
                print('\n', file=f)
                print('\n', file=f)

                    #bad_structures = output_structure_check(xyzfiles=xyzfiles, logfiles=files)
                    #print('Found', len(failed_runs), 'failed structures:', failed_runs)

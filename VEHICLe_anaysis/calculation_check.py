import numpy as np
from mol_translator import aemol
from mol_translator.properties.structure.bond_angle import get_bond_angles
from tqdm import tqdm

def comfile_check(comfiles, sampled_structures):

    print('need to check this against broken files to see if it picks it up')
    for i in tqdm(sampled_structures):
        molcom_list = [w for w in comfiles if i + '_' in w]

        neutral = [y for y in molcom_list if 'neutral' in y]
        anion = [x for x in molcom_list if 'anion' in x]
        bromine = [z for z in molcom_list if 'bromine' in z]

        with open(neutral[0], 'r') as f:
            neutralxyz = f.readlines()[6:]
            for j, k in zip(anion, bromine):
                with open(j, 'r') as p, open(k, 'r') as q:
                    alines = p.readlines()[6:]
                    blines = q.readlines()[6:]


                anionxyz = set(alines).intersection(neutralxyz)
                brominexyz = set(blines).intersection(neutralxyz)

                if len(anionxyz) != len(neutralxyz) - 1:
                    print(i)
                    print(neutral[0])
                    print(j)

                    if len(brominexyz) != len(neutralxyz) + 1:
                        pass
                        #print('CHECK FILE:', k)


def opt_check(file, string_to_search):
    with open(file, 'r') as f:
        for line in f:
            if string_to_search in line:
                return True
    return False

def get_dist_array(aemol):
    num_atoms = len(aemol.structure['types'])
    dist_array = np.zeros(num_atoms, num_atoms, dtype=np.float64)
    xyz_array = aemol.structure['xyz']
    for i in range(num_atoms):
        for j in range(num_atoms):
            dist_array[i][j] = np.absolute(np.linalg.norm(xyz_array[i] - xyz_array[j]))

    return dist_array

def output_structure_check(xyzfiles, logfiles, failed_runs=[]):
    '''
    This function checks the MAE between bond angles
    for the input and output structures.

    Typically, successful anions range between 0 - 12 (approx)
    and successful bromines range between 0 - 30 (approx)
    so any structures with MAE >= 30 are thrown up

    PLEASE CHECK THE STRUCTURES BEFORE DELETING THEM

    If this doesnt work then just change the condition:
    if np.mean(np.absolute(m-c)) > *set new number here*:
    '''

    badfiles = []
    for i in tqdm(xyzfiles):
        in_file = i.split('\\')
        in_split = in_file[1].split('-')

        for j in logfiles:

            out_file = j.split('\\')
            out_split = out_file[1].split('-')

            if in_split[0] == out_split[0]:
                input_aemol = aemol(in_file[1].split('_')[0])
                input_aemol.from_file(i, ftype = 'xyz')

                output_aemol = aemol(in_file[1].split('_')[0])
                output_aemol.from_file(j, ftype='log')

                input_structure = get_bond_angles(input_aemol)
                output_structure = get_bond_angles(output_aemol)

                err = []
                for x,y in zip(input_structure, output_structure):
                    for m, c in zip(x,y):
                        mnan = np.asarray([0 if x != x else x for x in m])
                        cnan = np.asarray([0 if x != x else x for x in c])

                        if 'anion' in i:
                            if np.mean(np.absolute(mnan-cnan)) > 12.0:
                                err.append(np.mean(np.absolute(mnan-cnan)))
                        elif 'bromine' in i:
                            if np.mean(np.absolute(mnan-cnan)) > 22.0:
                                err.append(np.mean(np.absolute(mnan-cnan)))

                if len(err) > 0:
                    badfiles.append(j)

    print('Large variances detected in ', len(badfiles), ' output structures:', badfiles)
    return badfiles
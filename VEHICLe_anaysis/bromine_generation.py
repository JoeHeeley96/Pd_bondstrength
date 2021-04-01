from xyzfile_generation import xyz_from_com
from xyz2mol import xyz2mol
from datetime import date

def bromine_check(bromine_comfilename):
    xyz_from_com(bromine_comfilename)



def bromines_from_com(comfilename):
    ### Probelm with this approach is that you end up appending a Br to each carbon, despite the nature of the atom in question (eg/ quaternary, carbonyl etc)
    split = comfilename.split('\\')
    name = split[1].split('_')

    with open(comfilename) as f:

        C_num = []
        C_coord = []
        Atoms=[]
        for line in f:
            if 'C' in line:
                C_coord.append(line)

        for i in range(len(C_coord)):
            C_num.append(i + 1)

        num_C_coords = zip(C_num, C_coord)

        for j, k in num_C_coords:
            bromine_xyzfilename = name[0] + str(j) + 'bromine.xyz'
            bromine_comfilename = name[0] + '_' + str(j) + 'bromine_' + name[1] + '_' + name[2] + '_' + name[3] + '_' + name[4]
            bromine_chkfilename = name[0] + '_' + str(j) + 'bromine_' + name[1] + '_' + name[2] + '_' + name[3] + '_opt.chk'


            with open('bromine_comfiles/' + bromine_comfilename, 'w') as p:
                with open(comfilename) as f:
                    for line in f:
                        if line.startswith('%chk='):
                            print('%chk=' + bromine_chkfilename, file=p)

                        elif line == '0 1\n':
                            print('1 1', file=p)

                        elif line == k:
                            exp = line.split(' ')
                            Br_xcoord = float(exp[1])
                            Br_ycoord = float(exp[2])
                            Br_zcoord = float(exp[3])
                            print(line.strip('\n'), file=p)
                            print('Br', (Br_xcoord), (Br_ycoord), (Br_zcoord + 2.0), file=p)


                        else:
                            print(line.strip('\n'), file=p)


def bromine_comfiles_from_xyz(xyzfilename):
    todays_date = str(date.today())
    split = xyzfilename.split('\\')
    name = split[1].split('_')

    with open(xyzfilename) as f:
            C_num = []
            C_coord = []
            for line in f:
                if 'C' in line:
                    C_coord.append(line)

            for i in range(len(C_coord)):
                C_num.append(i + 1)

            num_C_coords = zip(C_num, C_coord)

            for j, k in num_C_coords:
                bromine_xyzfilename = name[0] + '_' + str(j) + 'bromine_' + todays_date + '.xyz'
                bromine_comfilename = name[0] + '_' + str(j) + 'bromine_' + todays_date + '_wb97xd_631gd_opt.com'
                bromine_chkfilename = name[0] + '_' + str(j) + 'bromine_' + todays_date + '_wb97xd_631gd_opt.chk'

                with open('xyz_files/' + bromine_xyzfilename, 'w') as p:
                    with open(xyzfilename) as f:
                        for line in f:
                            ### How do we change the charge from here? ##


                            if line == k:
                                exp = line.split(' ')
                                Br_zcoord = float(exp[3])
                                print(line.strip('\n'), file=p)
                                print('Br', exp[1], exp[2], (Br_zcoord + 2.0), file=p)

                            else:
                                print(line.strip('\n'), file=p)
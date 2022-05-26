from datetime import date

def type_dict(warning=True):

    if warning:
        print('Type Dict includes: H, C, N, O, S and F ONLY! Please append it if you are including other atoms!'
              'See: neutral_generation.type_dict)')

    type_dict = {1: 'H', 6: 'C', 7: 'N', 8: 'O', 9: 'F', 16: 'S', 35: 'Br'}

    return type_dict

def xyzcoords_from_type_array(filename, zipped_type_array, warning=True):

    Atoms = type_dict(warning=warning)
    with open(filename, 'w') as f:
        for i in zipped_type_array:
            for atomicnum, atom in Atoms.items():
                if i[0] == atomicnum:
                    print(atom, i[1], *i[2], file=f)

def xyz_from_com(comfilename, save_loc, warning=True):

    split = comfilename.split('\\')
    name = split[-1].split('.')
    atoms = type_dict(warning=warning).values()
    with open(comfilename) as p:
        xyz = []
        read = p.readlines()[5:]
        remove_command=read[1:] ### Why did I do this? ###
        charge_multip= read[0]
        charge = int(charge_multip[0] + charge_multip[1])

        for line in remove_command:
            if not line.strip():
                continue

            elif 'F' in line and 'B' in line:
                pass

            else:
                for i in atoms:
                    if i in line:
                        xyz.append(line.strip('\n'))
        numatoms= len(xyz)
        xyzfilename = save_loc + '/' + name[0] + '.xyz'

        with open(xyzfilename, 'w') as g:
            print(numatoms,'\n', 'charge=' + str(charge) + '=', file=g, sep='')

            for x in xyz:
                line=str(x)
                print(line.strip('\n'), file=g)

    return xyzfilename


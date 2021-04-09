

def anion_from_com(comfilename):
    ### If we rewrite this we can make sure that only H's bound to carbons are stripped! ###
    split = comfilename.split('\\')
    name = split[1].split('_')
    with open(comfilename) as f:
        H_num=[]
        H_coord=[]
        for line in f:
            if 'H' in line:
                H_coord.append(line)

        for i in range(len(H_coord)):
            H_num.append(i + 1)

        num_H_coords=zip(H_num, H_coord)

        for j, k in num_H_coords:
            anion_comfilename = name[0] + '_' + str(j) + 'anion_' + name[1] + '_' + name[2] + '_' + name[3] + '_' + name[4]
            anion_chkfilename = name[0] + '_' + str(j) + 'anion_' + name[1] + '_' + name[2] + '_' + name[3] + '_opt.chk'
            with open('anion_comfiles/' + anion_comfilename, 'w') as p:
                    with open(comfilename) as f:
                            for line in f:
                                if line.startswith('%chk='):
                                    print('%chk=' + anion_chkfilename, file=p)
                                elif line == '0 1\n':
                                    print('-1 1', file =p)
                                elif line != k:
                                    print(line.strip('\n'), file=p)

def anions_from_smiles(dataframe, neutral_comfilename):
    print('not written yet. Use getneighbors() to find Hs only attached to Cs and the use rdgid to find comfile and to strip out those coords')
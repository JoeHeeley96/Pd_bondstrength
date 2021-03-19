

def anion_from_com(comfilename):
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
            with open('anion_comfiles/' + anion_comfilename, 'w') as p:
                    with open(comfilename) as f:
                            for line in f:
                                if line == '0 1\n':
                                    print('-1 1', file =p)
                                elif line != k:
                                    print(line.strip('\n'), file=p)
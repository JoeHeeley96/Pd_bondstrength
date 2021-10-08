import pandas as pd


def calculate_properties(dataset, outname):
    regId = list(set(dataset.Regid))
    property_data=pd.DataFrame(columns=['Regid'])


    for id in regId:
        regid_dict= {'Regid': id}
        regid_data = dataset[dataset['Regid'].str.fullmatch(id)]
        alist = list(regid_data.itertuples(index=False))
        basehet = (regid_data[regid_data['Structure'].str.contains('Base_het')]).Energy

        for j in basehet:
            for i in alist:
                if 'anion' in i.Structure:
                    acidity_au =(i[3] + -229.010318493) - (j + -228.425492536)
                    acidity_kcalpermol = acidity_au * 627.5
                    regid_dict[i[2] + '_acidity_kcalpermol'] = acidity_kcalpermol

                elif 'bromine' in i.Structure:
                    elec_affin_au = (j + -2571.17518775) - i[3]
                    elec_affin_kcalpermol = elec_affin_au * 627.5
                    regid_dict[i[2] + '_electrophile_affinity_kcalpermol'] = elec_affin_kcalpermol

        property_data=property_data.append(regid_dict, ignore_index=True)


    with open(outname, 'w') as f:
        print(property_data.to_csv(sep=','), file=f)


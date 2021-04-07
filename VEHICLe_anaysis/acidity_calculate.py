import pandas as pd


def calculate_properties(dataset, reference_data):

    basehets = dataset[dataset['Structure'].str.contains('(?!$)base_het')]
    anions = dataset[dataset['Structure'].str.contains('(?!$)anion')]
    bromines = dataset[dataset['Structure'].str.contains('(?!$)bromine')]

    RegIds = list(set(dataset.Regid))

    VEHICLe_AcidityData = pd.DataFrame([])
    acidities_au = []
    acidities_kcal_permol = []

    VEHICLe_EaData = pd.DataFrame([])
    electrophile_affinities_au = []
    electrophile_affinities_kcal_permol = []

    for Id in RegIds:
        hetanions_byId = anions[anions['Id'].str.contains(Id)]
        basehets_byId = basehets[basehets['Id'].str.contains(Id)]
        hetbromines_byId = bromines[bromines['Id'].str.contains(Id)]

        VEHICLe_AcidityData = VEHICLe_AcidityData.append([hetanions_byId])
        VEHICLe_EaData = VEHICLe_EaData.append([hetbromines_byId])

        print('WARNING: CHECK LABELS FOR REFERENCES AGAINST GENERATED REFERENCE DF')

        aceticacid = reference_data[reference_data['ExpNo'].str.contains('aceticacid')]
        acetate = reference_data[reference_data['ExpNo'].str.contains('acetate')]
        bromonium = reference_data[reference_data['ExpNo'].str.contains('bromonium')]

        anion_energy = hetanions_byId['Energy']
        aceticacid_energy = aceticacid['Energy']
        basehet_energy = basehets_byId['Energy']
        acetate_energy = acetate['Energy']
        bromine_energy = hetbromines_byId['Energy']
        bromonium_energy = bromonium['Energy']

        for acetateE in acetate_energy.items():
            for basehetE in basehet_energy.items():
                for aceticacidE in aceticacid_energy.items():
                    for anionE in anion_energy.items():
                        w = acetateE + basehetE + aceticacidE + anionE
                        x = list(w)
                        y = (x[-1] + x[-3]) - (x[1] + x[3])
                        z = y * 627.5

                        acidities_au.append(y)
                        acidities_kcal_permol.append(z)

        for bromoniumE in bromonium_energy.items():
            for basehetE_Ea in basehet_energy.items():
                for bromineE in bromine_energy.items():
                    j = bromoniumE + basehetE_Ea + bromineE
                    k = list(j)
                    l = (k[3] + k[1]) - k[-1]
                    m = l * 627.5
                    electrophile_affinities_au.append(l)
                    electrophile_affinities_kcal_permol.append(m)

    VEHICLe_AcidityData['Acidity_au'] = acidities_au
    VEHICLe_AcidityData['Acidity_kcalpermol'] = acidities_kcal_permol

    VEHICLe_EaData['Electrophile_Affinity_au'] = electrophile_affinities_au
    VEHICLe_EaData['Electrophile_Affinity_kcalpermol'] = electrophile_affinities_kcal_permol

    with open('VEHICLe_AcidityData.csv', 'w') as f:
        print(VEHICLe_AcidityData.to_csv(sep=','), file=f)

    with open('VEHICLe_EaDAta.csv', 'w') as p:
        print(VEHICLe_EaData.to_csv(sep=','), file=p)

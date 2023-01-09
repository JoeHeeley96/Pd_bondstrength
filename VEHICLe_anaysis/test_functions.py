import pandas as pd
import glob
import numpy as np

dj3_overlapped_structures = ['DR001', 'DR002', 'DR003', 'DR005', 'DR006', 'DR007', 'DR010', 'DR011', 'DR013', 'DR015', 'DR020', 'DR021',
                  'DR028', 'DR032', 'DR034', 'DR040', 'DR041', 'DR045', 'DR046', 'DR056', 'DR057', 'DR058', 'DR059',
                  'DR063', 'DR064']

succ_act = ['S1080', 'S1875', 'S1863', 'S5857', 'S6018', 'S5854', 'S5703', 'S1869', 'S5853', 'S1072']

unsucc_act = ['S5704', 'S2741']

targeted = ['S291', 'S5706', 'S5707', 'S2741', 'S7959', 'S54', 'S8063', 'S1873',
            'S5706', 'S5855', 'S2746', 'S5704', 'S49', 'S50']

dj3_exc = ['DR014', 'DR042', 'DR043', 'DR044', 'DR079', 'DR080', 'DR081', 'DR082', 'DR083', 'DR084', 'DR085',
           'DR086', 'DR087', 'DR088', 'DR001', 'DR002', 'DR003', 'DR005', 'DR006', 'DR007', 'DR010', 'DR011', 'DR013',
           'DR015', 'DR020', 'DR021', 'DR028', 'DR032', 'DR034', 'DR040', 'DR041', 'DR045', 'DR046', 'DR056', 'DR057',
           'DR058', 'DR059', 'DR063', 'DR064', 'DR009', ]


from calculation_check import quick_dj1_sub_check
#quick_dj1_sub_check()

from chembl_filters import dj1_point2_filter

dj1_point2_filter('chembl_31_chemreps.txt', outname='dj1_point2_chembl')





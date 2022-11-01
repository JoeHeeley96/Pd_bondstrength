from data_extraction import pd_readdata
import glob
import pandas as pd
from data_extraction import calculate_directed_activation
from data_extraction import plot_directed_activation
from data_extraction import calc_rel_exchange_data
from data_extraction import plot_3D_map_with_exchange_E

freqfiles = glob.glob(r'C:/Users/yv19151/OneDrive - University of Bristol/Documents/Orca_files/directed_activation/*_freq.out')
outfiles = glob.glob(r'C:/Users/yv19151/OneDrive - University of Bristol/Documents/Orca_files/directed_activation/*.out')
het_files = [x for x in outfiles if 'Pd' not in x]
opt_files = [x for x in outfiles if '_freq' not in x]
exchange_energies = pd.read_csv('directed_activation_Exchange_Energies.csv')
Pd_G = pd.read_csv('directed_activation_deltaG_Energies.csv')
#relE = pd.read_csv('directed_activation_Relative_Energies.csv')
rel_props = pd.read_csv('dj1_fulldata_relative.csv')

#df = calculate_directed_activation(fileloc, outname='directed_activation')
#plot_directed_activation(df, outname='directed_activation_results')
#calc_rel_exchange_data(df, 'directed_activation')

from data_extraction import read_deltaG
from data_extraction import plot_G_and_E
from data_extraction import calculate_deltaGexch
freqs = freqfiles + het_files

df1 = read_deltaG(freqs, outname='directed_activation')
df2 = calculate_deltaGexch(df1, 'directed_activation')

plot_G_and_E(df2, exchange_energies, outname='directed_activation_E_vs_G_comparison')


#plot_3D_map_with_exchange_E(df2, rel_props, outname='rel_props3D', animate=False)

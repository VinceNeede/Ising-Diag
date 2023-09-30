import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

plt.rcParams.update({'font.size': 22, 'figure.figsize':[11, 8]})
data=pd.read_csv("Magnetization and gap/chain_2.out")
figs=[None for _ in range(3)]
axs=[None for _ in range(3)]
for i in range(3):
    figs[i],axs[i]=plt.subplots()
for ispin in range(3,21):
    filter_data=data[data['L']==ispin]
    filter_data=filter_data.sort_values(by='kg')
    axs[0].plot(filter_data['kg'], abs(filter_data['broken_mag'])*ispin**(0.125),label=f"L= {ispin}",linewidth=1.5)
    axs[1].plot(filter_data['g'], abs(filter_data['broken_mag']),label=f"L= {ispin}",linewidth=1.5)
    axs[2].plot(filter_data['g'],abs(filter_data['E_m']-filter_data['E_p']),label=f"L= {ispin}",linewidth=1.5)
axs[0].set_ylabel(r'$m L^{1/8}$')
axs[1].set_ylabel(r'$m')
axs[2].set_ylabel(r'$\Delta$')
axs[0].set_xlabel(r'$(g-1)L$')
axs[1].set_xlabel(r'$g$')
axs[2].set_xlabel(r'$g$')
axs[0].set_xlim((-1.,1.))

figs[0].savefig("Magnetization and gap/universal_fun.png")
figs[1].savefig("Magnetization and gap/magn.png")
figs[2].savefig("Magnetization and gap/gap.png")
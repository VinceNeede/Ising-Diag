import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

data=pd.read_csv("Magnetization and gap/chain.out")
fig, axs = plt.subplots(2)

for ispin in range(3,20):
    axs[0].plot(data['g'].where(data['L']==ispin), np.abs(data['Mx'].where(data['L']==ispin)),label=f"L= {ispin}")
    axs[1].plot(data['g'].where(data['L']==ispin),abs((data['E_m0']-data['E_p0']).where(data['L']==ispin)))
plt.rcParams['text.usetex'] = True
axs[0].set_ylabel(r'$M_x$')
axs[0].set_ylim((0.0,1.05))
axs[1].set_ylabel(r'$\Delta$')
axs[1].set_xlabel('g/J')
plt.show()
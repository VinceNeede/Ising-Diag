import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

data=pd.read_csv("Magnetization and gap/chain.out")
fig, axs = plt.subplots(2)
for ispin in range(3,21):
    filter_data=data[data['L']==ispin]
    axs[0].plot(filter_data['kg'], abs(filter_data['broken_mag'])*ispin**0.125,label=f"L= {ispin}")
    axs[1].plot(filter_data['kg']/ispin+1.,abs(filter_data['E_m']-filter_data['E_p'])*ispin)
fig.legend()
axs[0].set_ylabel(r'$M L^{1/8}$')
#axs[0].set_xlim((-1.,1.))
axs[1].set_ylabel(r'$\Delta L$')
axs[1].set_xlabel(r'$k_g=(g-g_c)L$')
plt.show()
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

data=pd.read_csv("fixed_g_1.out")
fig, axs = plt.subplots(2)

axs[0].plot(data['L'], np.abs(data['Mx']))
axs[1].plot(data['L'],abs((data['E_m0']-data['E_p0'])))
for ax in axs:
    ax.set_xticks(data['L'])
axs[0].set_ylabel(r'$M_x$')
axs[1].set_ylabel(r'$\Delta$')
axs[1].set_xlabel('L')
plt.show()
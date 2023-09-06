import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

from scipy.optimize import curve_fit
data=pd.read_csv("Magnetization and gap/fixed_g_1.out")
fig, axs = plt.subplots(2)

def mag_scaling(x,a,b):
    return a*x**(-b)


popt, pcov = curve_fit(mag_scaling,data['L'],np.abs(data['Mx']),p0=[1.0, 1./8.])
print('beta= ', popt[1],'+-', np.sqrt(pcov[1,1]))
axs[0].plot(data['L'], np.abs(data['Mx']),'.',markersize=9)
x=np.linspace(3,23,100)
axs[0].plot(x, mag_scaling(x,*(popt)))

axs[1].plot(data['L'],abs((data['E_m0']-data['E_p0'])))
for ax in axs:
    ax.set_xticks(data['L'])
axs[0].set_ylabel(r'$M_x$')
axs[1].set_ylabel(r'$\Delta$')
axs[1].set_xlabel('L')
plt.plot
plt.show()
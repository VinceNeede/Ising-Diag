import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
plt.rcParams.update({'font.size': 22, 'figure.figsize':[11, 8]})

data=pd.read_csv('suscettività/susc.dat')
elles=[4,6,8,10,12,14,16,18,20] 
fig,axs=plt.subplots(1,2)
for ispin in elles:
    filtered_data = data[data['L'] == ispin]
    filtered_data=filtered_data.sort_values(by='g')
    axs[0].plot(filtered_data['g'], filtered_data['chi']*ispin**(-7./4.), label=f"L={ispin}",linewidth=2.5) 
    axs[1].plot((filtered_data['g']-1.)*ispin, filtered_data['chi']*ispin**(-7./4.), label=f"L={ispin}",linewidth=2.5)
axs[0].set_ylabel(r"$\chi L^{-7/4}$")
axs[0].set_xlabel(r"$g$")
axs[1].set_xlabel(r"$(g-1)L$")   
axs[1].legend([f"L={ell}" for ell in elles],loc="upper right",framealpha=0.)
#plt.legend()
#plt.xlim(0.9,1.1)
fig.savefig("suscettività/susc.png")

# from scipy.optimize import curve_fit
# data=pd.read_csv('suscettività/susc_fixed_g.dat')

# def mag_scaling(x,a,b):
#     return a*x**(b)

# popt, pcov = curve_fit(mag_scaling,data['L'], np.abs(data['chi']),p0=[1.0, 7./4.])
# print(popt)
# x=np.linspace(3,20,100)
# plt.plot(data['L'], np.abs(data['chi']))
# plt.plot(x,mag_scaling(x,*(popt)))
# plt.show()
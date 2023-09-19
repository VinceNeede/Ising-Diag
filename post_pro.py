import pandas as pd
import matplotlib.pyplot as plt
data=pd.read_csv("quench.dat")
elles=[10,12,15]#range(7,17)
#plt.rcParams['text.usetex'] = True
for ell in elles:
    filtered_data=data[(data['ell']==ell)]
    plt.plot(filtered_data['theta'], filtered_data['broken_mag'],label=f'L={ell}')
plt.legend()
plt.xlabel(r'$\Theta$')
plt.ylabel(r'$m L^{1/8}$')
#plt.xscale("log")
plt.show()
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
data=pd.read_csv("quench4.dat")
elles=[3,4,5,6]

for ell in elles:
    filtered_data=data[(data['ell']==ell)]
    #print(filtered_data)
    # print(filtered_data['theta'],filtered_data['magz'])
    plt.plot(filtered_data['steps'], np.abs((filtered_data['magz']-filtered_data.iloc[-1]['magz'])/filtered_data.iloc[-1]['magz'])*ell**(0.125),label=f'L={ell}')
plt.plot([1000,100000],[1.e-5,1.e-5])
plt.legend()
plt.xlabel(r'$\Theta$')
plt.ylabel(r'$m L^{1/8}$')
#plt.xlim(0.,39.)
#plt.ylim(-1.,0.6)
plt.yscale("log")
plt.show()
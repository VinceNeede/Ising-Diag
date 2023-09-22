import pandas as pd
import matplotlib.pyplot as plt
data=pd.read_csv("quench.dat")
elles=[5,8,10,13,15]
for ell in elles:
    filtered_data=data[(data['ell']==ell)&(data['theta']>=0.)]
    plt.plot(filtered_data['theta'], filtered_data['magz']*ell**(0.125),label=f'L={ell}')

plt.legend()
plt.xlabel(r'$\Theta$')
plt.ylabel(r'$m L^{1/8}$')
plt.xlim(34.,39.)
plt.ylim(-1.,0.6)
#plt.xscale("log")
plt.show()
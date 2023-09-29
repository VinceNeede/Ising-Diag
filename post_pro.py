import pandas as pd
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 22, 'figure.figsize':[11, 8]})
fig, ax=plt.subplots()

elles=[5,9,13,16,17]
for ell in elles:
    data=pd.read_csv(f"quench_{ell}.dat")
    filtered_data=data[(data['ell']==ell)&(data['theta']>=0.)]
    ax.plot(filtered_data['theta'], filtered_data['magz']*ell**(0.125),label=f'L={ell}',linewidth=2.5)

ax.legend(loc="upper left",framealpha=0.)
ax.set_xlabel(r'$\Theta$')
ax.set_ylabel(r'$m L^{1/8}$')
ax.set_xlim(34.,39.)
ax.set_ylim(-1.,0.6)
#plt.xscale("log")
fig.savefig("mag_vs_theta.png")
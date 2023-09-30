import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
plt.rcParams.update({'font.size': 22, 'figure.figsize':[11, 8]})

data=pd.read_csv('transverse susceptance/susc2.dat')
figs=[None for _ in range(2)]
axs=[None for _ in range(2)]
for i in range(2):
    figs[i],axs[i]=plt.subplots()
elles=[2*i+2 for i in range(1,10)]
for ispin in elles:
    filtered_data = data[data['L'] == ispin]
    filtered_data=filtered_data.sort_values(by='g')
    axs[0].plot(filtered_data['g'], filtered_data['chi']*filtered_data['g'],label=f"L={ispin}",linewidth=2.5)
axs[0].set_xlabel(r'$g$')
axs[0].set_ylabel(r'$g\chi_x$')
axs[0].legend(loc="upper right",framealpha=0.,ncol=1)
axs[0].set_xticks(np.arange(0., 2.5, step=0.5))
figs[0].savefig("transverse susceptance/specific_heat.png")

filtered_data = data[data['g'] == 1.]
filtered_data=filtered_data.sort_values(by='L')
axs[1].scatter(filtered_data['L'], filtered_data['chi'],color='red',marker='*',label="Lanczos")
axs[1].set_xlabel(r'$L$')
axs[1].set_ylabel(r'$g\chi_x(g=1)$')
x=np.linspace(1,100,10,endpoint=True)
def fun (x):
    return 0.32*np.log(x)+0.16
axs[1].semilogx()
axs[1].set_xticks([1,10,100],[r'$1$',r'$10$',r'$100$'])
axs[1].plot(x,fun(x),label=r'$0.32\log (L) + 0.16$')
axs[1].legend(loc="upper left",framealpha=0.,ncol=1)

figs[1].savefig("transverse susceptance/specific_heat_scaling.png")
plt.show()

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# data=pd.read_csv('suscettività/susc.dat')

# for ispin in range(3,20):
#     filtered_data = data[data['L'] == ispin]
#     plt.plot(filtered_data['g'], filtered_data['chi'])

# plt.show()

from scipy.optimize import curve_fit
data=pd.read_csv('suscettività/susc_fixed_g.dat')

def mag_scaling(x,a,b):
    return a*x**(-b)

popt, pcov = curve_fit(mag_scaling,data['L'], np.abs(data['chi']),p0=[1.0, 7./4.])
print(popt)
x=np.linspace(3,20,100)
plt.plot(data['L'], np.abs(data['chi']))
#plt.plot(x,mag_scaling(x,*(popt)))
plt.show()
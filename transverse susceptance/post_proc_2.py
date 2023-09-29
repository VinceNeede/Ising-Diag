import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

data=pd.read_csv('transverse susceptance/susc2.dat')

for ispin in range(3,21):
    filtered_data = data[data['L'] == ispin]
    filtered_data=filtered_data.sort_values(by='g')
    plt.plot(filtered_data['g'], filtered_data['chi']*filtered_data['g'])

plt.show()
filtered_data = data[data['g'] == 1.]
filtered_data=filtered_data.sort_values(by='L')
plt.scatter(filtered_data['L'], filtered_data['chi'])
x=np.linspace(3,20,10,endpoint=True)
def fun (x):
    return 0.32*np.log(x)+0.16
plt.semilogx(x,fun(x))
plt.show()

# from scipy.optimize import curve_fit
# #data=pd.read_csv('transverse susceptance/susc_fixed_g.dat')
# filtered_data = data[data['g'] == 1.]

# def mag_scaling(x,a,b):
#     return b*np.log(x)+a

# popt, pcov = curve_fit(mag_scaling,filtered_data['L'], filtered_data['chi'],p0=[1.0, 7./4.])
# print(popt)
# x=np.linspace(3,20,100)
# plt.semilogx()
# plt.plot(filtered_data['L'], filtered_data['chi'])
# plt.plot(x,mag_scaling(x,*(popt)))
# #plt.plot(x,mag_scaling(x,*(popt)))
# plt.show()
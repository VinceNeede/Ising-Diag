import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
data=pd.read_csv("quench_error_step.dat")
elles=[3,4,5,6,7,8,9,10,11,12]

def bubbleSort(x,y):
    len_array = len(x)
    for j in range(0, len_array-1):
      flag = False
      for i in range(0, len_array-1):
        if x[i]>x[i+1]:
            x[i], x[i+1] = x[i+1], x[i]
            y[i], y[i+1] = y[i+1], y[i]
            flag = True
            pos = i + 1
      if not(flag):
        break
      else: 
        len_array = pos
        

for ell in elles:
    filtered_data=data[(data['ell']==ell)]
    filtered_data=filtered_data.sort_values('steps')
    #print(filtered_data)
    # print(filtered_data['theta'],filtered_data['magz'])
    x=filtered_data['steps']
    y=filtered_data['magz']
    plt.plot(x, np.abs((y-y.iloc[-1])/y.iloc[-1]),label=f'L={ell}')
plt.plot([1000,300000],[1.e-5,1.e-5])
plt.legend()
plt.xlabel(r'$\Theta$')
plt.ylabel(r'$m L^{1/8}$')
plt.xlim(0.,3.e5)
#plt.ylim(-1.,0.6)
plt.yscale("log")
plt.savefig("step_test.png")
plt.show()
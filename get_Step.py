import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
data=pd.read_csv("quench4.dat")
elles=[3,4,5,6,7,8,9,10,11,12]
 
steps=[]
for ell in elles:
    filtered_data=data[(data['ell']==ell)&(data['steps']>=1000)]
    filtered_data=filtered_data.sort_values('steps')

    x=filtered_data['steps']
    y=filtered_data['magz']
    z=np.abs((y-y.iloc[-1])/y.iloc[-1])
    z=z-1.e-5
    for i,el in enumerate(z):
        if el<0:
            steps.append(x.iloc[i])
            #print(ell,x.iloc[i],y.iloc[i], z.iloc[i])
            break
        
from scipy.optimize import curve_fit
def model(x,a,b,c):
    return a*np.exp(b*x)+c
pops,pcov = curve_fit(model,elles,steps,p0=[5000,0.3,2000])
print(pops)
print(elles,steps)
plt.plot(elles,steps)
xp=np.linspace(3,12,100)
yp=model(xp,*pops)
plt.plot(xp,yp)
print("estrapolated:")
for i in range(11,18):
    print(i,[model(i,*pops)])
plt.savefig("prova.png")
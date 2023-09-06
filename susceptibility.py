from Ising import Ising
import numpy as np
import numdifftools as nd
import matplotlib.pyplot as plt

x=np.linspace(0.,2.,50)
x2=[]
y2=[]
def fun(g,ell):
    E_p, E_m, mag=Ising(ell,g,0.0,True,True)
    return mag
for ispin in range(3,10):
    f1=nd.Derivative(lambda x: fun(x,ispin), n=1, step=1.e-5)

    y= np.array([f1(el) for el in x])
    x2.append(ispin)
    y2.append(f1(1.))
    plt.plot(x,y)
plt.show()
plt.plot(x2,y2)
plt.show()
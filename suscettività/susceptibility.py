import numpy as np
import numdifftools as nd
import matplotlib.pyplot as plt
import sys
 
sys.path.append('./')
 
# importing
from Ising import Ising
x=np.linspace(0.,2.,50)
if 1. not in x:
    x=np.append(x,1.)
    x=np.sort(x)
x2=[]
y2=[]
def fun(g,ell):
    E_p, E_m, mag=Ising(ell,g,0.0,True,True)
    return mag
fout=open("suscettività/susc.dat","a")
#fout.write("L,g,chi\n")
from tqdm import tqdm
for ispin in range(16,20):
    f1=nd.Derivative(lambda x: fun(x,ispin), n=1, step=1.e-5)
    for el in tqdm(x):
        fout.write(str(ispin)+','+str(el)+','+str(f1(el)))
        fout.write('\n')
    fout.flush()
fout.close()
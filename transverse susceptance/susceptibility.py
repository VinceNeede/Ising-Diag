import numpy as np
import numdifftools as nd
import matplotlib.pyplot as plt
import sys
 
sys.path.append('./')
 
# importing
from Ising import Ising, Ising_Parity_Result
x=np.linspace(0.,2.,50)
if 1. not in x:
    x=np.append(x,1.)
    x=np.sort(x)
x2=[]
y2=[]
def fun(g,ell):
    res=Ising(ell,g,0.0,True,True)
    return res.tran_mag_p if res.E_p[0] < res.E_m[0] else res.tran_mag_m
fout=open("transverse susceptance/susc.dat","a")
fout.write("L,g,chi\n")
from tqdm import tqdm
for ispin in range(3,20):
    f1=nd.Derivative(lambda x: fun(x,ispin), n=1, step=1.e-5)
    for el in tqdm(x):
        fout.write(str(ispin)+','+str(el)+','+str(f1(el)))
        fout.write('\n')
    fout.flush()
fout.close()
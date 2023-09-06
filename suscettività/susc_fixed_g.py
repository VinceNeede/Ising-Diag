import numpy as np
import numdifftools as nd
import matplotlib.pyplot as plt
import sys
 
sys.path.append('./')
 
# importing
from Ising import Ising
import pandas as pd

data=pd.read_csv("critic g/g_max.dat")

def fun(h,ell,g=1.0):
    E, mag=Ising(ell,g,h,True,False)
    return mag
fout=open("suscettività/susc_fixed_g.dat","w")
fout.write("L,g,chi\n")
from tqdm import tqdm

def fun1(x,ell,h):
    return (fun(x+h,ell) - fun(x-h,ell))/(2*h)
from tqdm import tqdm

print(fun(0.5,12), fun(-0.5,12))

for ispin,gfield in tqdm(zip(data['L'],data['g_max'])):
    if ispin<18:
        f1=nd.Derivative(lambda x: fun(x,ispin,1.), n=1, step=1.e-5,method='forward')  
        fout.write(str(ispin)+','+str(gfield)+','+str(f1(0.))+'\n')
        fout.flush()
fout.close()
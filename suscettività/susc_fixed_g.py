import numpy as np
import numdifftools as nd
import matplotlib.pyplot as plt
import sys
 
sys.path.append('./')
 
# importing
from Ising import Ising
import pandas as pd

data=pd.read_csv("critic g/g_max.dat")

def fun(g,ell):
    E_p, E_m, mag=Ising(ell,g,0.0,True,True)
    return mag
fout=open("suscettività/susc_fixed_g.dat","w")
fout.write("L,g,chi\n")
from tqdm import tqdm
for ispin,gfield in zip(data['L'],data['g_max']):
    f1=nd.Derivative(lambda x: fun(x,ispin), n=1, step=1.e-5)
    fout.write(str(ispin)+','+str(gfield)+','+str(f1(gfield))+'\n')
    fout.flush()
fout.close()
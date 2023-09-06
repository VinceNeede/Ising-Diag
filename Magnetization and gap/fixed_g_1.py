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
fout=open("Magnetization and gap/fixed_g_1.dat","w")
fout.write("L,g,mag\n")
from tqdm import tqdm

for ispin in data['L']:
    gfield=1.0
    f1=nd.Derivative(lambda x: fun(x,ispin), n=1, step=1.e-5)
    fout.write(str(ispin)+','+str(fun(gfield, ispin))+'\n')
    fout.flush()
fout.close()
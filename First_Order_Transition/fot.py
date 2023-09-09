import numpy as np
import numdifftools as nd
import matplotlib.pyplot as plt
from tqdm import tqdm
import os
import sys

sys.path.append("./")

from Ising import Ising, Ising_Result


Ls=range(11,12)
gs=np.linspace(0.0,1.0,11)
def hs(L):

    hs_pos=np.linspace(start= 4.11e-5, stop=0.00019, num=11)
    hs_neg=-hs_pos
    hs=np.sort(np.concatenate((hs_pos,hs_neg)))
    return hs



def long_mag(L: int,g: np.double,h: np.double) -> np.double:
    res=Ising(L,g,h,True,False)
    return res._long_mag



fileout="First_Order_Transition/long_mag.dat"
if os.path.isfile(fileout):
    fout=open(fileout,"a")
else:
    fout=open(fileout,"w")
    fout.write("L,g,h,long_mag\n")


for L in tqdm(Ls, desc= "L", position=0):
    
    for g in tqdm(gs, desc="g", position=1, leave=False):
        for h in tqdm(hs(L), desc="h", position=2, leave=False):
            long_mag_result=long_mag(L,g,h)
            fout.write(str(L) +"," + str(g) + "," + str(h) + "," + str(long_mag_result) + "\n")
        fout.flush()
fout.close()
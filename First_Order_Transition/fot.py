import numpy as np
import numdifftools as nd
import matplotlib.pyplot as plt
from tqdm import tqdm
import os
import sys

sys.path.append("./")

from Ising import Ising, Ising_Result

gs=np.linspace(0.0,1.0,11)
hs_pos=np.logspace(start= np.log10(0.03), stop= np.log10(0.0002), num=16+1)
hs_neg=-hs_pos
hs=np.sort(np.concatenate((hs_pos,hs_neg)))
Ls=range(12,21)

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
    for g in tqdm(gs, desc="outer_loop", position=1, leave=False):
        for h in tqdm(hs, desc="inner_loop", position=2, leave=False):
            long_mag_result=long_mag(L,g,h)
            fout.write(str(L) +"," + str(g) + "," + str(h) + "," + str(long_mag_result) + "\n")
        fout.flush()
fout.close()
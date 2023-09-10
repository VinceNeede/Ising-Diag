import numpy as np
import numdifftools as nd
import matplotlib.pyplot as plt
from tqdm import tqdm
import os
import sys

sys.path.append("./")

from Ising import Ising, Ising_Result


Ls=range(17,18)
gs=np.linspace(0.0,1.0,11)
import pandas as pd


print("before hcrits")
hcrits=pd.read_csv("First_Order_Transition/critical_h/predicted_critical_hs.out", sep=",")
print(hcrits)


def h_critfunction(L):
    if 3 <= L <= 24:
        # Check if L is within the valid range
        row = hcrits[hcrits['L'] == L]  # Find the row where L matches
        if not row.empty:
            return row.iloc[0]['h_crit']  # Return the corresponding h_crit value
        else:
            return "L value not found in the DataFrame"
    else:
        return "L is outside the valid range"


def hs(L):
    hcrit=h_critfunction(L)
    hs_pos=np.logspace(start= np.log10(hcrit*0.1), stop=np.log10(hcrit*1000), num=40)
    hs=np.sort(hs_pos)
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
import numpy as np
import numdifftools as nd
import matplotlib.pyplot as plt
from tqdm import tqdm
import os
import sys

sys.path.append("./")

from Ising import Ising, Ising_Result


Ls=range(3,12)
gs=np.linspace(0.0,0.49,30)
hs=np.logspace(np.log10(1e-3),np.log10(1e-0),21)
hs=np.insert(hs,0,0)
PBC=True



fileout="FOTFSS/fot.dat"
if os.path.isfile(fileout):
    fout=open(fileout,"a")
else:
    fout=open(fileout,"w")
    fout.write("L,g,h,PBC,long_mag,E_0,E_1,E_2,tran_mag,broken_mag\n")


for L in tqdm(Ls, desc= "L", position=0):
    
    for g in tqdm(gs, desc="g", position=1, leave=False):
        for h in tqdm(hs, desc="h", position=2, leave=False):
            data=Ising(L,g,h,PBC=PBC,parity=False)
            string=f"{L},{g},{h},{PBC},{data.long_mag},{data.E[0]},{data.E[1]},{data.E[2]},{data.tran_mag},{data.broken_mag}\n"
            fout.write(string)
        fout.flush()
fout.close()
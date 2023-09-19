import numpy as np
import matplotlib.pyplot as plt
import sys
from tqdm import tqdm
sys.path.append('./')
from Ising import Ising, Ising_Parity_Result
 
elles=range(19,21)
kgs=np.linspace(-1.,1.,10)
if 0. not in kgs:
    kgs=np.append(kgs,0.)
    kgs=np.sort(kgs)
f=open("Magnetization and gap/chain.out","a")
#f.write("L,kg,E_m,E_p,broken_mag\n")
for ell in tqdm(elles):
    for kg in tqdm(kgs,leave=False):
        g=kg/ell+1.
        res=Ising(ell,g,0.,True,True)
        f.write(f"{ell},{kg},{res.E_m[0]},{res.E_p[0]},{res.broken_mag}\n")
    f.flush()
f.close()
        
        
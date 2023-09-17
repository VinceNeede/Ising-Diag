import numpy as np
import sys
from tqdm import tqdm
 
sys.path.append('./')
 
# importing
from Ising import quenching
fo=open("quenching.dat","w")
fo.write("L,mL,theta\n")
elles=[3,7,8,10]
thetas=np.linspace(0.,1.,500)
for ell in tqdm(elles):
    for theta in tqdm(thetas,leave=False):
        res=quenching(ell,0.8,0.9999,0.,theta,PBC=True)
        fo.write(f"{ell},{res.broken_mag*10**0.125},{theta}\n")
    fo.flush()
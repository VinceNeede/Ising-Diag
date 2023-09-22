import numpy as np
import sys
from tqdm import tqdm
 
sys.path.append('./')
 
# importing
from Quenching import quenching
elles=[3,4,5,6]
kg0=0.
kg=0.

kh0=-1.
kh=1.

steps=range(10000,100000,100)
f=open("quench4.dat","a")
#print("ell,steps,magz")
for ell in tqdm(elles):
    g0=kg0/float(ell)+1.
    g=kg/float(ell)+1.
    h0=kh0/(float(ell)**float(15./8.))
    h=kh/(float(ell)**float(15./8.))
    for i in tqdm(steps,leave=False):
        res = quenching(ell,g0,g,h0,h,39.,PBC=True,steps=i)
        f.write(f"{res.ell},{i},{res.magz}\n")
        #print(res.ell,i,res.magz,sep=',')
    f.flush()
f.close()
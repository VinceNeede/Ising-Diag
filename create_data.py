import numpy as np
import sys
from tqdm import tqdm
 
sys.path.append('./')
 
# importing
from Quenching import quenching
elles=[13]
kg0=-1.
kg=1.

for ell in tqdm(elles):
    g0=kg0/float(ell)+1.
    g=kg/float(ell)+1.
    quenching(ell,g0,g,0.,10.,PBC=True,steps=10000,OUTPUTFILE="quench.dat")

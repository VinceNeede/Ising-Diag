import numpy as np
import matplotlib.pyplot as plt
import sys
from tqdm import tqdm
sys.path.append('./')
from Ising import Ising, Ising_Parity_Result
import multiprocessing as mp

elles=range(3,5)
gs=np.linspace(0.,2.,10)
if 1. not in gs:
    gs=np.append(gs,1.)
    gs=np.sort(gs)
f=open("Magnetization and gap/chain_2.out","w")
f.write("L,kg,g,E_m,E_p,broken_mag\n")
f.flush()
points=np.ndarray(len(gs)*len(elles))
points=np.asarray([(ell, g) for ell in elles for g in gs])

def run(point):
    ell=int(point[0])
    g=np.double(point[1])
    kg=(g-1.)*ell
    res=Ising(ell,g,0.,True,True)
    f.write(f"{ell},{kg},{g},{res.E_m[0]},{res.E_p[0]},{res.broken_mag}\n")
    f.flush()
if __name__=='__main__':
    pool=mp.Pool(mp.cpu_count())
    for point in points:
        # run(point)
        pool.apply_async(run,args=(point,))
    pool.close()
    pool.join()

    f.close()
        
        
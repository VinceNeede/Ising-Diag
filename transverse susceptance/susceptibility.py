import numpy as np
import numdifftools as nd
import matplotlib.pyplot as plt
import sys
import multiprocessing as mp
sys.path.append('./')
 
# importing
from Ising import Ising, Ising_Parity_Result
x=np.linspace(0.9,1.1,10)
if 1. not in x:
    x=np.append(x,1.)
    x=np.sort(x)

def fun(g,ell):
    res=Ising(ell,g,0.0,True,True)
    return res.tran_mag_p if res.E_p[0] < res.E_m[0] else res.tran_mag_m
fout=open("transverse susceptance/susc2.dat","a")
#fout.write("L,g,chi\n")
#fout.flush()

elles=np.asarray([i for i in range(3,21)])
points=np.ndarray(len(x)*len(elles))
points=np.asarray([(ell, g) for ell in elles for g in x])


def run(point):
    ell=int(point[0])
    g=np.double(point[1])
    f1=nd.Derivative(lambda x: fun(x,ell), n=1, step=1.e-5)
    res=f1(g)
    fout.write(f"{ell},{g},{res}\n")
    fout.flush()

pool=mp.Pool(mp.cpu_count())
if __name__=='__main__':
    pool=mp.Pool(mp.cpu_count())
    for point in points:
        pool.apply_async(run,args=(point,))
    pool.close()
    pool.join()

    fout.close()

# for point in points:
#     f1=nd.Derivative(lambda x: fun(x,point[0]), n=1, step=1.e-5)
#     pool.apply_async()
    
# for ispin in range(3,21):
#     f1=nd.Derivative(lambda x: fun(x,ispin), n=1, step=1.e-5)
#     for el in x:
#         pool.apply_async(f1,)
#         fout.write(str(ispin)+','+str(el)+','+str(f1(el)))
#         fout.write('\n')
#     fout.flush()
# fout.close()
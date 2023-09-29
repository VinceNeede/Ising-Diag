import numpy as np
import numdifftools as nd
import matplotlib.pyplot as plt
import sys
import multiprocessing as mp
sys.path.append('./')
 
# importing
from Ising import Ising
x=np.linspace(0.9,1.1,50)

if 1. not in x:
    x=np.append(x,1.)
    x=np.sort(x)
x2=[]
y2=[]
def fun(g,h,ell):
    res=Ising(ell,g,h,True,False)
    # print(res.long_mag)
    return res.long_mag
fout=open("suscettività/susc.dat","w")
fout.write("L,g,chi\n")
fout.flush()
from tqdm import tqdm

elles=np.asarray([i for i in range(3,21)])

points=np.ndarray(len(x)*len(elles))
points=np.asarray([(ell, g) for ell in elles for g in x])
# yp=[] 
# ell=8
# for el in x:
#     dx=1.e-3
#     def f1(x):
#         return (fun(el,x+dx,ell) - fun(el,x-dx,ell))/(2.*dx)
#     yp.append(f1(0.)*ell**(-7./4.))
# plt.scatter(x,yp)
# plt.xlim(0.9,1.1)
# plt.show()

def run(point):
    dx=1.e-3
    ell=int(point[0]) 
    g=np.double(point[1])
    def f1(x):
        return (fun(g,x+dx,ell) - fun(g,x-dx,ell))/(2.*dx)
    res=f1(0.)
    fout.write(f"{ell},{g},{res}\n")
    fout.flush()
if __name__=='__main__':
    pool=mp.Pool(mp.cpu_count())
    for point in points:
        pool.apply_async(run,args=(point,))
    pool.close()
    pool.join()

    fout.close()
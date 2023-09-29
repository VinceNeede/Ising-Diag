import numpy as np
import sys
from tqdm import tqdm
import multiprocessing as mp
sys.path.append('./')
 
# importing
from Quenching import quenching
elles=[17]#=range(3,18)
steps=[410000]#=[10000, 20000, 40000, 50000, 80000, 100000, 130000, 150000, 180000, 200000, 240000, \
   # 270000, 330000, 360000, 410000]
kg0=0.
kg=0.

kh0=-1.
kh=1.

n_cpu=mp.cpu_count()


    
pool=mp.Pool(n_cpu)
for ell ,step in tqdm(zip(elles,steps)):
    g0=kg0/float(ell)+1.
    g=kg/float(ell)+1.
    h0=kh0/(float(ell)**float(15./8.))
    h=kh/(float(ell)**float(15./8.))
    quenching(ell,g0,g,h0,h,39., PBC=True, steps=step, OUTPUTFILE= f"quench_{ell}.dat")
   # pool.apply_async(quenching,args=(ell,g0,g,h0,h,39.),\
   #     kwds={'PBC': True,'steps': step, 'OUTPUTFILE': f"quench_{ell}.dat"})
pool.close()
pool.join()  

# for i in tqdm(steps,leave=False):
#     res = quenching(ell,g0,g,h0,h,39.,PBC=True,steps=i)
#     f.write(f"{res.ell},{i},{res.magz}\n")
#     #print(res.ell,i,res.magz,sep=',')

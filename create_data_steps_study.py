import numpy as np
import sys
from tqdm import tqdm
import multiprocessing as mp
sys.path.append('./')
 
# importing
from Quenching import quenching
elles=[13,14,15,16,17]
kg0=0.
kg=0.

kh0=-1.
kh=1.

n_cpu=mp.cpu_count()
#steps=[i for i in range(185000,190000,250)]#+[500000]
steps=[222286,259053,300000,350000,400000]
num_batches = len(steps) // n_cpu
if len(steps) % n_cpu != 0:
    num_batches += 1
# split the points
batched_points = np.array_split(np.asarray(steps), num_batches)

f=open("quench4.dat","a")
#f.write("ell,steps,magz\n")
#f.flush()
results=""
def coll_res(res):
    global results
    results+=f"{res.ell},{res.steps},{res.magz}\n"
for ell in tqdm(elles):
    g0=kg0/float(ell)+1.
    g=kg/float(ell)+1.
    h0=kh0/(float(ell)**float(15./8.))
    h=kh/(float(ell)**float(15./8.))

    for batch in tqdm(batched_points,leave=False):
        results=""
        pool=mp.Pool(n_cpu)
        for i in batch:
            pool.apply_async(quenching,args=(ell,g0,g,h0,h,39.),kwds={'PBC': True,'steps': i},callback=coll_res)
        pool.close()
        pool.join()   
        f.write(results)     
    # for i in tqdm(steps,leave=False):
    #     res = quenching(ell,g0,g,h0,h,39.,PBC=True,steps=i)
    #     f.write(f"{res.ell},{i},{res.magz}\n")
    #     #print(res.ell,i,res.magz,sep=',')
    f.flush()
f.close()
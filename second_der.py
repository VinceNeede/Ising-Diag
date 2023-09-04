import scipy as sc
import numdifftools as nd
import subprocess
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
for ispin in range(19,23):
    def fun(x0):
        f=open("chain.in", "w")
        #print(float(x0))
        f.write(
    f"{ispin}	    ! Number of spins in the system		(ell)\n\
    {float(x0)}	! Transverse magnetic field strength	(g)\n\
    0.d0	! Longitudinal magnetic field strength	(h)\n\
    .true.	! Type of boundary conditions  		(.true. -> PBC,        .false. -> OBC)\n\
    .true.	! Type of diagonalization      		(.true. -> Davidson,   .false. -> Lapack full diag) NOT IMPLEMENTED\n"\
            )
        f.close()
        p=subprocess.Popen("./_results/Ising_Parity", shell=True, stdout=subprocess.PIPE)
        out=p.communicate()[0].decode()
        outs=out.split(',')
        #print(outs)
        return abs(float(outs[8].replace('\n','')))

    f2=nd.Derivative(fun, n=2, step=1.e-5)
    # x=np.linspace(0., 2., 50)
    # y=np.array([f2(el) for el in x])
    # plt.plot(x,y)
    # plt.show()
    # print(f2(0.85),f2(1.05))
    sol=sc.optimize.root_scalar(f2,bracket=(0.88,1.03),method='bisect',rtol=1.e-7)
    print(ispin, sol.root, sep=',', flush=True)

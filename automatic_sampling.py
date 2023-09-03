import numpy as np
import matplotlib.pyplot as plt

def norm_pdf (x, mu=1.0, sigma=0.5):
    return 1/(sigma * np.sqrt(2 * np.pi))*np.exp( - (x - mu)**2 / (2 * sigma**2))
    
start=0.
end=2.
npoints=20

x=np.linspace(start, end, npoints)
fx=norm_pdf(x)
g=0
density=(1+g*(fx-1))/fx
x_density=np.cumsum(density)
x_density-=x_density.min()
x_density/=x_density.max()
x_density*=(end-start)+start
if 1. not in x_density:
    x_density=np.append(x_density,1.)
    x_density=np.sort(x_density)
# print(x)
# print(x_density)


# fx_density = norm_pdf(x_density)


# plt.plot(x,x_density,'ok',ms = 10,label = 'linear')
# # plt.plot(x_density,fx_density,'or',ms = 10,label = 'warped')
# # plt.legend(loc = 'upper right')
# plt.show()

import subprocess


for ispin in range(3,10):
    for gfield in x_density:
        f=open("chain.in", "w")
        f.write(
f"{ispin}	    ! Number of spins in the system		(ell)\n\
{gfield}d0	! Transverse magnetic field strength	(g)\n\
0.d0	! Longitudinal magnetic field strength	(h)\n\
.true.	! Type of boundary conditions  		(.true. -> PBC,        .false. -> OBC)\n\
.true.	! Type of diagonalization      		(.true. -> Davidson,   .false. -> Lapack full diag) NOT IMPLEMENTED\n"\
            )
        f.close()
        p=subprocess.Popen("./_results/Ising_Parity", shell=True, stdout=subprocess.PIPE)
        out=p.communicate()[0].decode()
        if (out!=""):
            print(ispin, gfield, out)
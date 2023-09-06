import subprocess
import numpy as np
import sys

if sys.platform=='linux':
    indir='./Linux_results/'
else:
    indir='win_results/'

np.set_printoptions(precision=15)

def Ising(ell: int, gfield: np.double, hfield:np.double, PBC: bool, parity: bool):
    if parity:
        if hfield!=0.:
            print("WARNING: hfield has to be zero in order to implement the parity simmetry")
        fileout='Ising_Parity'
    else:
        fileout='Ising'
        
    f=open("chain.in", "w")
    f.write(
f"{ell}	    ! Number of spins in the system		(ell)\n\
{gfield}	! Transverse magnetic field strength	(g)\n\
{hfield}	! Longitudinal magnetic field strength	(h)\n\
.{'true' if PBC else 'false'}.	! Type of boundary conditions  		(.true. -> PBC,        .false. -> OBC)\n\
.true.	! Type of diagonalization      		(.true. -> Davidson,   .false. -> Lapack full diag) NOT IMPLEMENTED\n"\
        )
    f.close()
    p=subprocess.Popen(indir+fileout, shell=True, stdout=subprocess.PIPE)
    out=p.communicate()[0].decode().replace('\n','').split(sep=',')
    out_ell=float(out.pop(0))
    out_gfield=float(out.pop(0))
    if parity:
        E_p=np.array([np.double(out.pop(0)) for _ in range(3)],dtype=np.double)
        E_m=np.array([np.double(out.pop(0)) for _ in range(3)],dtype=np.double)
        mag=abs(np.double(out.pop(0)))
        return E_p, E_m, mag        
    E=np.array([np.double(out.pop(0)) for _ in range(3)])
    mag=abs(np.double(out.pop(0)))
    return E, mag

""" def Ising_Parity(ell: int, gfield: float, PBC: bool):
    f=open("chain.in", "w")
    f.write(
f"{ell}	    ! Number of spins in the system		(ell)\n\
{gfield}d0	! Transverse magnetic field strength	(g)\n\
0.d0	! Longitudinal magnetic field strength	(h)\n\
.{'true' if PBC else 'false'}.	! Type of boundary conditions  		(.true. -> PBC,        .false. -> OBC)\n\
.true.	! Type of diagonalization      		(.true. -> Davidson,   .false. -> Lapack full diag) NOT IMPLEMENTED\n"\
        )
    f.close()
    p=subprocess.Popen(indir+"Ising_Parity", shell=True, stdout=subprocess.PIPE)
    out=p.communicate()[0].decode().replace('\n','').split(sep=',')
    out_ell=float(out.pop(0))
    out_gfield=float(out.pop(0))
    E_p=np.array([float(out.pop(0)) for _ in range(3)])
    E_m=np.array([float(out.pop(0)) for _ in range(3)])
    mag=abs(float(out.pop(0)))
 """

if __name__=='__main__':
    print(Ising(15,0.2, 0.5, False, parity=False))
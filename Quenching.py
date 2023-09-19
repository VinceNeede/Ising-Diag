import subprocess
import numpy as np
import sys
import pandas as pd
import os

if sys.platform=='linux':
    indir='./Linux_results/'
else:
    indir='.\\win_results\\'

np.set_printoptions(precision=15)

class QuenchResult:
    def __init__(self, ell, gfield0, gfield1, hfield, theta, PBC, 
            thetaf, broken_mag, magz, magx) -> None:
        self.ell=int(ell)
        self.gfield0=np.double(gfield0)
        self.gfield1=np.double(gfield1)
        self.hfield=np.double(hfield)
        self.theta=np.double(theta)
        self.PBC = True if PBC == ' T' else False
        self.tehtaf=np.double(thetaf)
        self.broken_mag=np.double(broken_mag)
        self.magz=np.double(magz)
        self.magx=np.double(magx)
        

def quenching(ell: int, gfield0: np.double, gfield1: np.double, hfield:np.double, theta: np.double, PBC: bool, 
            steps: int=1000, OUTPUTFILE:str=None):
    fileout="quenching"        
    f=open("chain.in", "w")
    f.write(
f"{ell}	    ! Number of spins in the system		(ell)\n\
{gfield0}	! Starting magnetic field strength	(g0)\n\
{gfield1}	! Ending magnetic field strength	(g1)\n\
{hfield}	! Longitudinal magnetic field strength	(h)\n\
{theta}	! Rescaled time	(theta)\n\
.{'true' if PBC else 'false'}.	! Type of boundary conditions  		(.true. -> PBC,        .false. -> OBC)\n\
.true.	! Type of diagonalization      		(.true. -> Davidson,   .false. -> Lapack full diag) NOT IMPLEMENTED\n\
{steps}    ! Number of steps for the runge-kutta\n\
.{'true' if OUTPUTFILE is not None else 'false'}.  ! .true. prints each step, .false. prints only the ending step"\
        )
    f.close()

    
    if OUTPUTFILE is not None:
        if os.path.isfile(OUTPUTFILE):
            f=open(OUTPUTFILE,"a")
        else:
            f=open(OUTPUTFILE,"w")
            f.write("ell,g0,g1,h,theta_fin,PBC,theta,broken_mag,magz,magX\n")
            f.flush()
        p=subprocess.Popen(indir+fileout, shell=True, stdout=f, text=True)
        p.wait()
        f.flush()
        f.close()
        return
    p=subprocess.Popen(indir+fileout, shell=True, stdout=subprocess.PIPE)
    return QuenchResult(*p.communicate()[0].decode().replace('\n','').split(sep=','))

if __name__=='__main__':
    res=quenching(10,0.8,0.9999,0.,1.,PBC=True, OUTPUTFILE="pippp.dat")
    print(res['ell'])

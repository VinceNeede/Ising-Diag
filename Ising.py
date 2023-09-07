import subprocess
import numpy as np
import sys

if sys.platform=='linux':
    indir='./Linux_results/'
else:
    indir='win_results/'

np.set_printoptions(precision=15)


class IRes:
    
    @property
    def ell(self)->int:
        return self._ell
    @ell.setter
    def ell(self, value: int):
        self._ell=value
        
    @property
    def long_field(self)->np.double:
        return self._long_field
    @long_field.setter
    def long_field(self, value: np.double):
        self._long_field=value
    
    @property
    def tran_field(self)->np.double:
        return self._tran_field
    @tran_field.setter
    def tran_field(self, value: np.double):
        self._tran_field=value
        
    @property
    def broken_mag(self):
        return self._broken_mag
    @broken_mag.setter
    def broken_mag(self, value):
        self._broken_mag=value
    
    @property
    def PBC(self)->bool:
        return self._PBC
    @PBC.setter
    def PBC(self, value: bool):
        self._PBC=value
        

class Ising_Parity_Result(IRes):
    def __init__(self) -> None:
        pass

    @property
    def E_p(self)->np.ndarray:
        return self._E_p
    @E_p.setter
    def E_p(self, value: np.ndarray):
        self._E_p=value
    
    @property
    def E_m(self)->np.ndarray:
        return self._E_m
    @E_m.setter
    def E_m(self, value: np.ndarray):
        self._E_m=value

    @property
    def tran_mag_p(self):
        return self._tran_mag_p
    @tran_mag_p.setter
    def tran_mag_p(self, value):
        self._tran_mag_p=value
    
    @property
    def tran_mag_m(self):
        return self._tran_mag_m
    @tran_mag_m.setter
    def tran_mag_m(self, value):
        self._tran_mag_m=value    
        

class Ising_Result(IRes):
    
    def __init__(self) -> None:
        pass
    
    @property
    def E(self)->np.ndarray:
        return self._E
    @E.setter
    def E(self, value: np.ndarray):
        self._E=value

    @property
    def long_mag(self):
        return self._long_mag
    @long_mag.setter
    def long_mag(self, value):
        self._long_mag=value
    
    @property
    def tran_mag(self):
        return self._tran_mag
    @tran_mag.setter
    def tran_mag(self, value):
        self._tran_mag=value    
        
        
        
        

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
        if PBC:
            broken_mag=np.double(out.pop(0))
            tran_mag_p=np.double(out.pop(0))
            tran_mag_m=np.double(out.pop(0))
        else:
            broken_mag=np.array([np.double(out.pop(0)) for _ in range(ell)])
            tran_mag_p=np.array([np.double(out.pop(0)) for _ in range(ell)])
            tran_mag_m=np.array([np.double(out.pop(0)) for _ in range(ell)])
        res=Ising_Parity_Result()
        res.ell=ell
        res.long_field=hfield
        res.tran_field=gfield
        res.PBC=PBC
        res.E_p=E_p
        res.E_m=E_m
        res.broken_mag=broken_mag
        res.tran_mag_p=tran_mag_p
        res.tran_mag_m=tran_mag_m
        return res    
    E=np.array([np.double(out.pop(0)) for _ in range(3)])
    broken_mag=np.abs(np.double(out.pop(0)))
    if PBC:
        long_mag=np.array([np.double(out.pop(0)) for _ in range(ell)])
        tran_mag=np.array([np.double(out.pop(0)) for _ in range(ell)])
    else:
        long_mag=np.double(out.pop(0))
        tran_mag=np.double(out.pop(0))
    res=Ising_Result()
    res.ell=ell
    res.long_field=hfield
    res.tran_field=gfield
    res.PBC=PBC
    res.E=E
    res.broken_mag=broken_mag
    res.long_mag=long_mag
    res.tran_mag=tran_mag
    return res

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
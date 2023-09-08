import subprocess
import numpy as np
import sys

if sys.platform=='linux':
    indir='./Linux_results/'
else:
    indir='.\\win_results\\'

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
    else:
        E=np.array([np.double(out.pop(0)) for _ in range(3)])
        broken_mag=np.abs(np.double(out.pop(0)))
        if PBC:
            long_mag=np.double(out.pop(0))
            tran_mag=np.double(out.pop(0))
        else:
            long_mag=np.array([np.double(out.pop(0)) for _ in range(ell)])
            tran_mag=np.array([np.double(out.pop(0)) for _ in range(ell)])
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







import os
def automatic_sampler (Ls, gfields, hfields, PBC: bool, parity:bool):
    from collections.abc import Iterable
    print("IT IS WIDELY SUGGESTED TO BUILD YOUR OWN SAMPLER AND USE THIS ONE ONLY AS AN EXAMPLE")
    if parity:
        print("It can raise some errors if parity=True given the difference in the outputs variables")
    #Idea: Ls, gfields, hfields can also be a single number, rather than an one-dimensional array. In that case, that variable is fixed.
    
    OUTPUTFILE="computed_stuff.out"
    if os.path.isfile(OUTPUTFILE):
        fout=open(OUTPUTFILE,"a")
    else:
        fout=open(OUTPUTFILE,"w")
        fout.write("L,g,h,E_ground,long_mag,trans_mag\n") #Qui ho deciso di mettere solo queste variabili e solo in questo ordine, poi giustamente se tu vuoi cambiarlo o se vuoi aggiungere altre energie oltre a quella di ground ne parliamo e possiamo modificare questa funzione di prova
    
    for L in Ls if isinstance(Ls,Iterable) else [Ls]:
        for gfield in gfields if isinstance(gfields,Iterable) else [gfields]:
            for hfield in hfields if isinstance(hfields,Iterable) else [hfields]:
                res=Ising(L, gfield, hfield, PBC, parity)
                E_ground=res.E[0]
                
                newrow= f"{res.ell},{res.tran_field},{res.long_field},{E_ground},{res.long_mag},{res.tran_mag}\n"
                fout.write(newrow)
                fout.flush()
    

    fout.close()
    return 0

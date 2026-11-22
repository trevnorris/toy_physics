"""Independent audit of closed-cavity readout saved before optional unscaled-Bessel assertion."""
import json, sys, hashlib
from pathlib import Path
import numpy as np
from scipy.special import hankel1e
sys.path.insert(0,'software/stage1_solver/tools')
from mouth_exchange_shear import OutgoingSlab
root=Path('_scratch/s11c/exchange_independent01')
for n in [32,64]:
    data=root/f'wave{n}'/'stdout'
    recs=[json.loads(line) for line in data.read_text().splitlines()]
    assert len(recs)==2 and {rec['kind'] for rec in recs}=={'interior','exchange'}
    for rec in recs:
        e=np.array(rec['clamped_modes'])
        assert e.shape==(3,2) and np.all(e[:,1]-e[:,0]>.77)
        assert abs((e[2,0]-e[0,0])/.01-rec['slope'])<1e-10
        assert abs((e[2,0]+e[0,0]-2*e[1,0])/.005**2-rec['curvature'])<1e-9
        print('closed',rec['kind'],n,rec['slope'],rec['curvature'],hashlib.sha256(data.read_bytes()).hexdigest())
r=3.
p=OutgoingSlab(np.c_[np.full(25,r),np.linspace(-1,0,25)],r)
omega=4.20802-3.9e-8j
k=np.sqrt(omega**2-p.cutoff**2)
k=np.where((p.cutoff>omega.real)&(k.imag<0),-k,k)
z=k*r;nu=4.5
exact=.5/r+.5*k*(hankel1e(nu-1,z)-hankel1e(nu+1,z))/hankel1e(nu,z)
err=np.abs(exact-p.traction(omega))
print('channels',len(err),'port max abs',float(max(err)))
assert max(err)<1e-11 and np.all(exact[:2].imag>0)

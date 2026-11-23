"""Guarded field and power check for the R2 curved duct and receiver."""
import argparse
import json
import numpy as np
from mouth_curved_sound import CurvedSound,mean_front,power_matrix
from mouth_curved_steady import reference, SteadyMouth

def main():
    p=argparse.ArgumentParser()
    p.add_argument('--n',type=int,required=True)
    p.add_argument('--box',type=float,default=6.)
    p.add_argument('--pml',type=float,default=3.)
    p.add_argument('--omega',type=float,default=3.06656)
    p.add_argument('--mean',action='store_true')
    a=p.parse_args()
    shape,chi,mesh,scale,_,rows,digest=reference()
    st=SteadyMouth(shape,chi,mesh)
    j=st.info['terminal_current']/(4*np.pi*shape.exit**3/3)
    sheet=mean_front(j) if a.mean else None
    op=CurvedSound(shape,scale,a.n,a.box,a.pml,sheet,st)
    res=op.evaluate(a.omega)
    ct=min([r for r in rows if r.get('kind')=='spectrum' and r['relative_exit']==0],key=lambda r:r['shape_step'])
    v=np.array(ct['mass_normalized_vectors'][0])
    wk=lambda m: (v@m@v).item()
    complex_pair=lambda z:[z.real,z.imag]
    print(json.dumps(dict(n=a.n,box=a.box,pml=a.pml,mean=a.mean,omega=a.omega,mesh=mesh,
         digest=digest,j=j,matrix_size=int(op.size),vector=v.tolist(),W=complex_pair(wk(res['W'])),
         work=wk(power_matrix(res['W'],a.omega)).real,feed=wk(res['feed_power']).real,
         receiver=wk(res['bulk_power']).real,feed_cut=wk(res['feed_cut']).real,receiver_cut=wk(res['bulk_cut']).real,
         radii=res['radii'].tolist(),area=res['disk_weights'].tolist(),
         currents=list(map(complex_pair,res['absolute_current']@v)),
         moving_current=list(map(complex_pair,res['relative_current']@v)),
         sheet_force=list(map(complex_pair,res['front_force']@v)),
         fed_pressure=list(map(complex_pair,res['pressure_minus']@v)),
         receiving_pressure=list(map(complex_pair,res['pressure_plus']@v)))))

if __name__=='__main__':main()

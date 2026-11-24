"""Regional R4 crystal energy, order derivative and signed displacement controls."""
import argparse
from pathlib import Path
from hashlib import sha256
import json
import sys

import numpy as np

root=Path.cwd()
parser=argparse.ArgumentParser(description='Check the three added R4 displacement forces by signed work on a smooth trial.')
parser.add_argument('output',type=Path)
area=parser.parse_args().output.resolve()
area.mkdir(parents=True,exist_ok=True)
sc=root/'_scratch/s11c'
sys.path[:0]=list(map(str,[sc/'r4_exact_support01/patch',sc/'r4_melted_survival01/pin/software/stage1_solver/tools',root/'software/stage1_solver/tools']))
from mouth_free_shape_solid import Solid
from mouth_r4_refined_solid import RefinedSolid
from mouth_r4_fold_order import axes_model
from mouth_sealed_shape_relax import Meridian

snapshot=sc/'r4_feasible_descent01/coarse/stalled_00011.npz'
assert sha256(snapshot.read_bytes()).hexdigest()=='c99706009efdaf3ce94f5c8a626a737927ff2ce4705cb030c843b5862d908cef'
prep=sc/'supercore_joined_front01/adopt309/state.npz'
with np.load(snapshot) as f:
    assert sha256(prep.read_bytes()).hexdigest()==str(f['state_sha'])
    shape=Meridian(4.,f['shape_params'],float(f['shape_extra']))
    g=axes_model(shape,json.loads((sc/'r4_whole_mouth01/coarse/endpoint.json').read_text())['depth'])
    body=Solid(prep,shape,float(f['x'][5]),8,24,6,g,degree=3,subdivisions=2)
    q=np.r_[f['q'],np.zeros(3)];order=f['field'][2]
body.buffers['base'].fill(0)
body.melted_labels=True
body.data.exact_melt_support=1
child=RefinedSolid(body,area/'build')
out={'seam':child.seam,'fraction':child.seam_fraction,'regions':{},'virtual_work':{}}
for label,reg in [('target',(5,18,0,0)),('star','star'),('outside',(7,18,0,0))]:
    sides=[child.evaluate_region(q,order,reg,refined=r,load=True) for r in (False,True)]
    old,new=sides
    delta=dict(energy=abs(old[0]-new[0]),
        g=float(np.max(abs(old[1][:-3]-new[1][:-3]))),
        deposition=float(np.max(abs(old[2]-new[2]))), extra=new[1][-3:].tolist())
    assert max(delta[k] for k in ('energy','g','deposition'))<1.e-14
    out['regions'][label]=delta

# Positive and spatially nonuniform order exercises the advected-order term
# and the moving cut quadrature. The nonzero state has positive determinants.
neutral=np.zeros_like(q)
neutral[-3:]=[.0012,-.0007,.002]
block=np.empty(body.ids.shape)
x,a,t=np.meshgrid(body.buffers['field_x'],body.buffers['field_a'],body.buffers['field_t'],indexing='ij')
block=.25+.55*x*x*(.8+.2*np.cos(a))*(.65+.35*t*t)
spatial=np.full_like(order,.25);spatial[body.ids]=block
e,grad,dep,info=child.evaluate_region(neutral,spatial,'star',load=True,diagnose=True)
assert info['min_J']>0.
out['neutral_min_J']=float(info['min_J'])
for c in range(3):
    direction=np.zeros_like(q); direction[-3+c]=1.
    def energy(h):
        return child.evaluate_region(neutral+h*direction,spatial,'star',gradient=False)[0]
    deltas=[]
    for h in [2e-6,1e-6]:
        slope=(energy(h)-energy(-h))/(2*h)
        deltas.append([slope,slope-grad[-3+c]])
    assert abs(deltas[-1][1])<1e-7*max(1.,abs(grad[-3+c]))
    out['virtual_work'][str(c)]=dict(analytic=grad[-3+c],signed=deltas)
vary=np.zeros_like(order); vary[body.ids]=.07*x*(1-t)*np.sin(a)
h=2e-6
eplus=child.evaluate_region(neutral,spatial+h*vary,'star',gradient=False)[0]
eminus=child.evaluate_region(neutral,spatial-h*vary,'star',gradient=False)[0]
slope=(eplus-eminus)/(2*h); adj=float(dep@vary)
assert abs(slope-adj)<1e-8
out['order_work_error']=slope-adj
(area/'result.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
print(json.dumps(out))

"""Plot the saved R4 order field at the deformed polar angle of the folding corner."""
import argparse
import hashlib
import json
from pathlib import Path
import sys

import numpy as np


def read(path):
    return json.loads(path.read_text())


def main(gate, output, figure):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.colors import PowerNorm
    from scipy.interpolate import RegularGridInterpolator
    d=read(gate)
    path=Path(d['checkpoint'])
    work=path.parents[1]
    if d['action']!='paused' or hashlib.sha256(path.read_bytes()).hexdigest()!=d['sha256']:
        raise ValueError('Expected the paused state identified by the stop gate')
    pin=work/'pin/software/stage1_solver/tools'
    sys.path.insert(0,str(pin))
    from mouth_continuous_depth_mesh import axes_for_depth
    from mouth_last_rings_fluid import GRIDS
    from mouth_sealed_shape_relax import Meridian
    from mouth_supercore_joined_front import Joined

    class Layout(Joined):
        def matrix(self,*unused):
            pass  # Use the producer's numbering, without assembling a new field solve.

    original=work.parent/'r4_whole_mouth01/coarse'
    reference_depth=read(original/'endpoint.json')['depth']
    state=Path(read(original/'invocation.json')['command'][5])
    with np.load(state) as f:
        cut=np.exp(-.5*np.exp(float(f['registration'].item())))
    with np.load(path) as f:
        shape=Meridian(4.,f['shape_params'],extra=f['shape_extra'].item())
        layout=Layout(shape, None, 1, GRIDS['double'],precondition=False)
        rr, zs=axes_for_depth(shape.length,reference_depth)
        order=f['field'][2,layout.ids[1,0]]
        solid=f['solid']
    if shape.neck!=4. or np.any(order<0.) or np.any(order>1.):
        raise ValueError('Invalid radius-four order input')
    n=(solid.shape[1]-1)//3
    sx=np.r_[np.linspace(0.,cut,2*n+1)[:-1],np.linspace(cut,1.,n+1)]
    st=np.linspace(0.,1.,solid.shape[0])
    sa=np.linspace(0.,np.pi/2,solid.shape[2])
    key=d['latest']['orientation']['extremizers']['min']
    j,i,k=key['cell_t_x_angle']
    corners=((0,0),(0,1),(1,0)) if key['triangle']==0 else ((1,1),(0,1),(1,0))
    dj,di=corners[key['vertex']]
    i,j=i+di,j+dj
    frac=key['angular_fraction']
    angle=(1-frac)*sa[k]+frac*sa[k+1]
    c=np.array([sx[i],st[j],angle])+(1-frac)*solid[j,i,k]+frac*solid[j,i,k+1]
    x,a,t=rr[0],layout.a,zs[1]/shape.length
    interp=RegularGridInterpolator((x,a,t),order)
    if abs(np.asarray(interp(c[[0,2,1]])).item())>1.+1.e-13:
        raise ValueError('Corner lies outside the order field')
    X,T=np.meshgrid(np.linspace(0.,1.,220),np.linspace(0.,1.,220))
    section=interp((X,np.full_like(X,c[2]),T))
    upper=order.max(axis=(0,1))
    lower=order.min(axis=(0,1))
    thresholds={}
    for val in (0.,.01,.1,.5,.9):
        bad=np.flatnonzero(upper>val)
        i0=0 if len(bad)==0 else bad[-1]+1
        thresholds[str(val)]={'first_knot_bounding_all_downstream':float(t[i0]) if i0<len(t) else None}
    table=dict(checkpoint=d['sha256'],folded_current_chart_x_t_theta=c.tolist(),depth=t.tolist(),
               section_max=upper.tolist(),section_min=lower.tolist(),bound_to_disk=thresholds,
               region='entire interior shaft; Q1 maximum over radius and polar angle at each depth knot')
    output.write_text(json.dumps(table,indent=2,allow_nan=False)+'\n')

    fig,(ax,side)=plt.subplots(1,2,figsize=(8.5,4.9),width_ratios=(1.25,.95),sharey=True,
                               layout='constrained')
    norm=PowerNorm(.6,vmin=0.,vmax=1.)
    mesh=ax.pcolormesh(X,T,section,cmap='cividis',norm=norm,shading='auto',rasterized=True)
    ax.scatter(c[0],c[1],s=83,marker='x',color='#f27342',linewidths=2.3,label='Folding corner')
    ax.set(xlim=(0,1),ylim=(1.015,-.015),xlabel='Radial mesh coordinate (wall at 1)',
           ylabel='Depth / mouth depth',title='Section through the corner')
    ax.legend(loc='lower left',frameon=True,fontsize=9)
    cb=fig.colorbar(mesh,ax=ax,orientation='horizontal',fraction=.07,pad=.09,ticks=[0,.1,.5,1])
    cb.set_label('Crystalline order χ')
    side.fill_betweenx(t,lower,upper,color='#d2e1ed',label='Range over radius and angle')
    side.plot(upper,t,color='#215476',lw=1.7,label='Maximum order')
    side.plot(lower,t,color='#718b9a',lw=1,ls='--',label='Minimum order')
    side.axhline(.97765,color='#a04b2c',ls=':',lw=1.3,label='Exit-rim sample')
    side.set(xlim=(-.016,1.015),xlabel='Crystalline order χ',title='Full shaft')
    side.grid(axis='x',alpha=.22)
    side.legend(loc='center right',fontsize=8,frameon=False)
    fig.savefig(figure,dpi=170)
    plt.close(fig)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for key in ('gate','output','figure'):
        p.add_argument(key,type=Path)
    a=p.parse_args()
    main(a.gate,a.output,a.figure)

"""Local work for the difference between the two interior extensions of the same moving R4 rim."""
import argparse
import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np

from mouth_closed_rim_surface_review import Reader, quadrature
from mouth_disk_pillbox import dot
from mouth_exchange_carrier import mapped_nodes
from mouth_rim_bank_work import Bank
from mouth_rim_disk_work_review import Disk
from mouth_rim_pinned_comparison import digest, validate_guard
from mouth_shared_rim_shift import SharedRimChart
from mouth_supercore_joined_front import ORDER_INERTIA, state_from


def gradient(disk, slab, side, coords, values):
    """Gradient of the registered Q1 order at a specified physical point."""
    r, theta, z = np.broadcast_arrays(*coords)
    b, slope = disk.sh.radius(z)
    width = 28-b if side else b
    x = (r-b)/width if side else r/width
    axes = (disk.model.rr[side], disk.model.a, disk.model.z[slab])
    loc = (x, theta, z)
    idx = [np.clip(np.searchsorted(a,v)-1, 0, len(a)-2) for a,v in zip(axes,loc)]
    gap = [a[i+1]-a[i] for a,i in zip(axes,idx)]
    frac = [(v-a[i])/d for a,i,d,v in zip(axes,idx,gap,loc)]
    g = np.zeros((r.size,3))
    for i,j,k in np.ndindex(2,2,2):
        bits=(i,j,k)
        v=values[idx[0]+i,idx[1]+j,idx[2]+k]
        w=[f if s else 1-f for s,f in zip(bits,frac)]
        for d in range(3):
            g[:,d] += v*(2*bits[d]-1)*w[(d+1)%3]*w[(d+2)%3]/gap[d]
    g[:,0] /= width
    g[:,1] /= r
    g[:,2] -= slope*(1-x if side else x)*g[:,0]
    return g


def measure_block(disk,chart,slab,side,nq):
    m=disk.model
    view = SimpleNamespace(model=m, sh=disk.sh, reader=SimpleNamespace(vortex=disk.reader.vortex,
                            fields={(1,side):disk.reader.fields[slab,side]}))
    bank=Bank(view,side,None,0.)
    # For W_4-W_8 the two extensions coincide outside this radial band.
    # Per-extension entries are therefore restricted to their differing support.
    edges=sorted([chart.edges[n][side] for n in (4,8)])
    xmin,xmax=(0.,edges[-1]) if side else (edges[0],1.)
    rx=m.rr[side]
    xs,wx=quadrature(np.unique(np.r_[xmin,rx[(rx>xmin)&(rx<xmax)],xmax]),nq)
    aa,wa=quadrature(m.a,nq)
    x,angle=np.meshgrid(xs,aa,indexing='ij')
    x,angle=x.ravel(),angle.ravel()
    weights=(wx[:,None]*wa).ravel()
    zmin,zmax=chart.tents[[0,2] if slab==1 else [2,4]]
    zgrid=m.z[slab]
    zz,wz=quadrature(np.unique(np.r_[zmin,zgrid[(zgrid>zmin)&(zgrid<zmax)],zmax]),nq)
    hvals=(2e-5,-2e-5,1e-5,-1e-5)
    aggregate=np.zeros((2,5,3))
    for z,dz in zip(zz,wz):
        b=disk.sh.radius(z)[0]
        width=28-b if side else b
        r=b+x*width if side else x*width
        volume=4*np.pi*weights*width*r*r*np.sin(angle)*dz
        flows=list(bank.physical(x,angle,z))
        register=chart.physical(side,r,z)
        for ix,n in enumerate((4,8)):
            eta,er,ez=np.subtract(chart.carrier(side,r,z,n),register)
            for sign,before,source in zip((1.,-1.),flows,disk.reader.fields[slab,side]):
                (_,p,c,rho,_,_),amp,a,F=before
                y=amp**2
                potential=F**10*(y-1)**2*(y**3+2*y**2+3*y+4)/20
                pc=dot(p,c)
                denominator=rho+ORDER_INERTIA*dot(c,c)
                projection=-.5*ORDER_INERTIA*rho*pc**2/denominator
                order=source.values[:,:,:,2]
                for it,h in enumerate(hvals):
                    jac=1+h*ez
                    c1=gradient(disk,slab,side,(r,angle,z+h*eta),order)
                    def pull(g):
                        out=g.copy()
                        out[:,2]/=jac
                        out[:,0]-=h*er*out[:,2]
                        return out
                    p1,a1=pull(p),pull(a)
                    projection1=-.5*ORDER_INERTIA*rho*dot(p1,c1)**2/(rho+ORDER_INERTIA*dot(c1,c1))
                    parts=[.5*F**2*(jac*dot(a1,a1)-dot(a,a))+(jac-1)*potential,
                           .5*rho*(jac*dot(p1,p1)-dot(p,p)),jac*projection1-projection]
                    aggregate[ix,it]-=sign*(np.array(parts)@volume)/h
                # Differentiate the integrand before integration. The smooth
                # spatial derivative of the material FE gradient uses a local
                # physical-coordinate centered probe within its own element.
                eps=1e-7
                dc=(gradient(disk,slab,side,(r,angle,z+eps),order)
                    -gradient(disk,slab,side,(r,angle,z-eps),order))*(eta/(2*eps))[:,None]
                dp=np.stack((-er*p[:,2],np.zeros_like(r),-ez*p[:,2]),axis=-1)
                da=np.stack((-er*a[:,2],np.zeros_like(r),-ez*a[:,2]),axis=-1)
                dpj=-ORDER_INERTIA*rho*pc/denominator*(dot(dp,c)+dot(p,dc))
                dpj-=projection*2*ORDER_INERTIA*dot(c,dc)/denominator
                derivative=[ez*(.5*F**2*dot(a,a)+potential)+F**2*dot(a,da),
                            .5*rho*ez*dot(p,p)+rho*dot(p,dp),ez*projection+dpj]
                aggregate[ix,4]-=sign*(np.array(derivative)@volume)
    return aggregate[0]-aggregate[1]


def run(a):
    validate_guard(a.job)
    validate_guard(a.faces.parent)
    previous=json.loads(a.manifest.read_text())
    c=next(c for c in previous['cases'] if c['inputs']['grid']==a.job.name)
    m=c['inputs']
    face=json.loads(a.faces.read_text())
    if (digest(a.job/'zero.npz') != c['zero_state_sha256'] or face['field']!=c['zero_state_sha256']
        or digest(a.original)!=m['field_sha256'] or digest(a.axes)!=m['reference_axes_sha256']
        or digest(a.wall)!=m['wall_sha256'] or digest(a.state)!=m['state_sha256']):
        raise ValueError('Unmatched input')
    disk=Disk(a.state,a.wall,a.original)
    with np.load(a.job/'zero.npz') as f:
        disk.reader=Reader(disk.model,state_from(a.state)[1],f['wound'],f['rest'])
    # Use the actual original physical taper and doubled cell extents from the independently checked chart.
    native_radial=mapped_nodes(13,disk.sh.length,28.,22.)[0]
    chart=object.__new__(SharedRimChart)
    chart.shape=disk.sh
    chart.outer_taper=native_radial[1][4:6]
    chart.tents=np.array(face['heights'])
    chart.window=np.array(disk.windows[0])
    chart.edges={n:(m['duct_bounds'][i],m['receiver_bounds'][i]) for i,n in enumerate((4,8))}
    rows=[]
    for read in face['readings']:
        nq=read['gauss']
        contributions=read['cases']
        swept=np.subtract(contributions['rim_four']['breakdown'],contributions['rim_eight']['breakdown']).sum(axis=1)
        cells=[]
        for slab,side in np.ndindex(2,2):
            w=measure_block(disk,chart,slab+1,side,nq)
            fneg,fpos=swept[:,slab,side]
            c=w[4].sum()+(fneg+fpos)/2
            cells.append(dict(slab=slab+1,outer=side,smooth_rows=w.tolist(),swept_negative=fneg,
                              swept_positive=fpos,centered_first_order=c,
                              true_one_sided=[w[4].sum()+fneg,w[4].sum()+fpos]))
        rows.append(dict(gauss=nq,regions=cells,
                         totals=np.sum([c['true_one_sided'] for c in cells],axis=0).tolist(),
                         centered=sum(c['centered_first_order'] for c in cells)))
    return dict(grid=a.job.name,zero_sha256=face['field'],face_sha=digest(a.faces),
                motion='Relative V4-S vs V8-S. Common material/disk motion cancels at first order.',
                row_order=['+2e-5','-2e-5','+1e-5','-1e-5','smooth first derivative'],
                sectors=['amplitude and native density','pure phase','phase-order projection'],checks=rows)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('state','wall','original','axes','job','manifest','faces'):
        p.add_argument(k,type=Path)
    print(json.dumps(run(p.parse_args()),allow_nan=False))

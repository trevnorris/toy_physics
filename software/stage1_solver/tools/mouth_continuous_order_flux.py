"""Two-sided weak and face-integrated order fluxes on a held absolute-Q1 R4 field."""
import argparse
import copy
import ctypes
from pathlib import Path
import sys

import numpy as np

from mouth_held_continuous_order import digest, read, write


def norm(v):
    return float(np.max(abs(v), initial=0.))


def split_load(t, q, j):
    """Duplicate the one material *test* trace at the former lower edge.

    Passing two equal axial knots makes its endpoint basis independent on
    either side. The native upper_bound lookup in the crystal kernel selects a
    positive-length cell on both sides, including a point exactly on the cut.
    The reference quadrature, placement and its physical integration partition
    are unchanged. Summing both work covectors must reproduce the original.
    """
    from mouth_free_shape_solid import Data, Grid, pointer
    b = copy.copy(t.b)
    b.data = Data.from_buffer_copy(t.b.data)
    old = t.model.z[1]/t.shape.length
    z = np.insert(old,j,old[j])
    b.data.field_t = Grid(len(z),pointer(z))
    count = t.model.size
    b.ids = np.concatenate((b.ids[:,:,:j+1],b.ids[:,:,j:]+count),axis=-1)
    b.zero_dc = np.zeros(b.ids.shape)
    _, _, load, _ = b.evaluate(q, gradient=False)
    loads = np.zeros(2*count)
    loads[:len(load)]=load
    return loads.reshape(2,count)


def terms(t, u):
    from mouth_breathing_mobile_patch import KAPPA, WELL, ORDER_INERTIA as I
    m=t.model
    c=m.V@u[2]
    g=np.stack([mat@u[2] for mat in m.G],axis=1)
    p=np.stack([mat@u[1] for mat in m.G],axis=1)
    n=(m.D[:,4]*(m.V@u[0]))**2
    ratio=I*np.einsum('ij,ij->i',g,p)/(n+I*np.einsum('ij,ij->i',g,g))
    pi=n*ratio
    velocity=p-ratio[:,None]*g
    return 2*WELL*c*(1-c)*(1-2*c),KAPPA*g-pi[:,None]*velocity


def weak(t, u, q, supplied):
    from mouth_breathing_mobile_patch import WELL, KAPPA
    m=t.model
    lower=int(np.argmin(abs(m.z[1]-t.shape.length+.12)))
    inner=split_load(t,q,lower)
    source_error=norm(inner.sum(axis=0)+t.ext-supplied)
    if source_error>2.e-11:
        raise ArithmeticError(f'Split material test load differs by {source_error}')
    well,flux=terms(t,u)
    pieces={}
    def block(slab,side,j,above):
        rows = m.blocks[slab,side]
        start = sum(8*len(nodes) for nodes,offset in m.blocks.values() if offset<rows[1])
        ids = np.arange(start,start+8*len(rows[0]))
        if j is not None:
            take = (m.D[ids,3]>m.z[slab][j])
            ids=ids[take if above else ~take]
        w=m.hq[ids]
        bulk=m.V[ids].T@(w*well[ids])
        for k,G in enumerate(m.G):
            bulk+=G[ids].T@(w*flux[ids,k])
        if (slab,side)==(1,0):
            bulk+=inner[int(above)]
        elif (slab,side)==(1,1):
            bulk+=t.ext
        return t.pull.T@bulk
    for slab,side in m.ids:
        if side==0 and slab in (1,2):
            j=lower if slab==1 else int(np.argmin(abs(m.z[2]-t.shape.length-.12)))
            for upper in (False,True):
                pieces[slab,side,int(upper)] = block(slab,side,j,upper)
        else:
            pieces[slab,side,0]=block(slab,side,None,False)
    c = u[2,t.first]
    total=sum(pieces.values())
    m.material_load=supplied
    grad=t.pull.T@m.evaluate(u)[1][2]
    difference=norm(total-grad)
    if difference>2.e-8:
        raise ArithmeticError(f'Independent weak assembly differs by {difference}')
    def compare(ix,left,right):
        ix=np.unique(t.number[ix])
        d0,d1=(sum(pieces[key] for key in group)[ix] for group in (left,right))
        # Exclude intersection lines; the two open sides exhaust each nodal star.
        omitted=norm(total[ix]-d0-d1)
        if omitted>2.e-10:
            raise ArithmeticError(f'Other domain enters the interface star: {omitted}')
        interior=(c[ix]>0)&(c[ix]<1)
        lo=c[ix]==0
        hi=c[ix]==1
        summed=d0+d1
        def row(mask):
            a,b=d0[mask],d1[mask]
            return dict(nodes=int(sum(mask)),lower_or_inner_outward_sum=float(a.sum()),
                        upper_or_outer_outward_sum=float(b.sum()),signed_mismatch=float((a+b).sum()),
                        nodal_mismatch_L1=float(abs(a+b).sum()),nodal_mismatch_max=norm(a+b),
                        both_sides_L1=float(abs(a).sum()+abs(b).sum()))
        return dict(interior=row(interior), lower_obstacle=row(lo),upper_obstacle=row(hi),
                    wrong_sign_lower_max=norm(np.minimum(summed[lo],0.)),
                    wrong_sign_upper_max=norm(np.maximum(summed[hi],0.)),
                    order_range=[float(min(c[ix])),float(max(c[ix]))])
    left_core=[(1,0,0),(1,0,1)]
    right_core=[(2,0,0),(2,0,1)]
    inside=m.ids[1,0]
    outer=m.ids[1,1]
    jtop=int(np.argmin(abs(m.z[2]-t.shape.length-.12)))
    tests={
        'former_lower_edge':(inside[:-1,:,lower],[(1,0,0)],[(1,0,1)]),
        'former_upper_edge':(m.ids[2,0][:-1,:,jtop],[(2,0,0)],[(2,0,1)]),
        'terminal_disk':(inside[:-1,:,-1],left_core,right_core),
        'entrance_open_disk':(inside[:-1,:,0],[(0,0,0)],left_core),
        'intact_slab':(outer[1:-1,:,0],[(0,1,0)],[(1,1,0)]),
        'lateral_two_banks':(inside[-1,:,1:-1],left_core,[(1,1,0)]),
        'post_disk_radial_join':(m.ids[2,0][-1,:,1:-1],right_core,[(2,1,0)]),
        'feed_radial_join':(m.ids[0,0][-1,:,:-1],[(0,0,0)],[(0,1,0)]),
        'outer_terminal_plane':(outer[1:-1,:,-1],[(1,1,0)],[(2,1,0)])}
    return dict(raw_gradient_reconstruction_max=difference,split_material_load_error=source_error,
                surfaces={name:compare(*args) for name,args in tests.items()})


def direct_face(t,u,slab,side,axis,plane,high_cell,degree=6):
    """Integrate physical Q·dS on one face, oriented in +w or +r.

    Values are direct Q1 one-sided traces; they need not coincide at finite
    resolution even when the common weak flux balances.
    """
    from mouth_joined_angular_work import bump
    from mouth_breathing_mobile_patch import KAPPA, ORDER_INERTIA as I
    m=t.model
    xs, angles, zs=m.rr[side], m.a, m.z[slab]
    ci=plane if high_cell else plane-1
    iv=np.array([ci]) if axis==0 else np.arange(len(xs)-1)
    jv=np.array([ci]) if axis==2 else np.arange(len(zs)-1)
    ix,ia,iz=np.meshgrid(iv,np.arange(len(angles)-1),jv,indexing='ij')
    dx, da, dz=np.diff(xs)[ix],np.diff(angles)[ia],np.diff(zs)[iz]
    points,weights=np.polynomial.legendre.leggauss(degree)
    qw=list(zip((1+points)/2,weights/2))
    fixed=[(0. if high_cell else 1.,1.)]
    sums=np.zeros_like(dx)
    for fx,wx in (fixed if axis==0 else qw):
        for fa,wa in qw:
            for fz,wz in (fixed if axis==2 else qw):
                x,a,z=xs[ix]+fx*dx,angles[ia]+fa*da,zs[iz]+fz*dz
                if axis == 2:
                    z = np.full_like(z, zs[plane])
                b,bz=t.shape.radius(z)
                # At w=0 and L the meridian returns the duct's one-sided
                # derivative. Both bath charts have constant radial scale.
                if slab != 1:
                    bz = 0.
                width=(28-b) if side else b
                r=(b+x*width) if side else x*b
                h,hz=bump(z,t.shape.length)
                P=(3*np.cos(a)**2-1)/2
                Pa=-3*np.sin(a)*np.cos(a)
                k=1-x if side else x**3
                kx=-1 if side else 3*x*x
                e=t.polar_shift
                rad=r+e*h*k*P
                rx=width+e*h*kx*P
                ra=e*h*k*Pa
                rz=bz*(1-x if side else x)+e*hz*k*P
                vals=np.zeros((3,*x.shape))
                raw=np.zeros((3,3,*x.shape))
                for i in (0,1):
                    for j in (0,1):
                        for k0 in (0,1):
                            bx=fx if i else 1-fx
                            ba=fa if j else 1-fa
                            bz0=fz if k0 else 1-fz
                            v=u[:,m.ids[slab,side][ix+i,ia+j,iz+k0]]
                            vals+=v*bx*ba*bz0
                            raw[:,0]+=v*(2*i-1)*ba*bz0/dx
                            raw[:,1]+=v*(2*j-1)*bx*bz0/da
                            raw[:,2]+=v*(2*k0-1)*bx*ba/dz
                raw[:,0]/=rx
                raw[:,1]=(raw[:,1]-ra*raw[:,0])/rad
                raw[:,2]-=rz*raw[:,0]
                n=(m._vortex(rad*np.sin(a))[0]*vals[0])**2
                g,p=raw[2],raw[1]
                ratio=I*(g*p).sum(axis=0)/(n+I*(g*g).sum(axis=0))
                flux=KAPPA*g-n*ratio*(p-ratio*g)
                if axis==2:
                    current=flux[2]*rx*dx
                else:
                    current=(flux[0]-ra/rad*flux[1]-rz*flux[2])*dz
                sums+=4*np.pi*rad**2*np.sin(a)*da*wx*wa*wz*current
    return sums.squeeze()


def manufactured_axial_traces(t):
    """Check each physical bank with chi=x+side and zero meridional phase.

    On the spherical disk faces the polar deformation vanishes. On either
    bath bank chi is independent of w at fixed physical radius. In the duct
    its physical w derivative follows from x=r/b or x=(r-b)/(28-b).
    Integrating these derivatives against 4*pi*r**2 dr gives the checks below.
    """
    from mouth_breathing_mobile_patch import KAPPA
    m = t.model
    u = np.zeros((3, m.size))
    u[0] = 1.
    for (_, side), ids in m.ids.items():
        u[2, ids] = side + m.rr[side][:, None, None]
    rows = {}
    for slab, high_cell, edge in ((0, False, "feed"), (1, True, "duct_entrance"),
                                  (1, False, "duct_exit"), (2, True, "receiver")):
        plane = 0 if high_cell else len(m.z[slab]) - 1
        b, db, _ = t.shape.profile(0. if edge in ("feed", "duct_entrance") else 1.)
        width = 28. - b
        for side in (0, 1):
            integral = (np.pi*b*b if side == 0 else
                        4*np.pi*(b*b/2 + b*width/3 + width**2/12))
            exact = -KAPPA*db*integral if slab == 1 else 0.
            actual = float(direct_face(t, u, slab, side, 2, plane, high_cell).sum())
            rows[f"{edge}_{side}"] = dict(expected=float(exact), measured=actual)
    return rows


def direct(t,u):
    m=t.model
    lower=int(np.argmin(abs(m.z[1]-t.shape.length+.12)))
    upper=int(np.argmin(abs(m.z[2]-t.shape.length-.12)))
    last=lambda arr:len(arr)-1
    shared={
        'former_lower_edge':((1,0,2,lower,False),(1,0,2,lower,True)),
        'former_upper_edge':((2,0,2,upper,False),(2,0,2,upper,True)),
        'terminal_disk':((1,0,2,last(m.z[1]),False),(2,0,2,0,True)),
        'entrance_open_disk':((0,0,2,last(m.z[0]),False),(1,0,2,0,True)),
        'intact_slab':((0,1,2,last(m.z[0]),False),(1,1,2,0,True)),
        'lateral_two_banks':((1,0,0,last(m.rr[0]),False),(1,1,0,0,True))}
    out={}
    for name,(left,right) in shared.items():
        a,b=direct_face(t,u,*left),direct_face(t,u,*right)
        ap,bp=direct_face(t,u,*left,degree=4),direct_face(t,u,*right,degree=4)
        out[name]=dict(positive_axis_integrals=[float(a.sum()),float(b.sum())],
            cellwise_unsigned_difference=float(abs(a-b).sum()),integrated_jump=float((b-a).sum()),
            quadrature_change=max(float(abs(a-ap).sum()),float(abs(b-bp).sum())))
    return out


def main(a):
    from mouth_continuous_order import ContinuousTrial
    from mouth_joined_front_audit import Front
    from mouth_last_rings_fluid import GRIDS
    from mouth_continuous_evidence import held_receipt
    receipt = held_receipt(a.held/'held')
    saved=a.held/'saved'
    rec=read(saved/'record.json')
    pin = a.held/'pin'
    for module, key in (('mouth_continuous_order.py', 'model'), ('mouth_held_continuous_order.py', 'code')):
        if digest(pin/module) != rec[key]:
            raise ValueError(f'Changed held model input {module}')
    sources = a.original/'sources.json'
    if digest(sources) != rec['saved_original_pin']:
        raise ValueError('Changed original source manifest')
    for module, sha in read(sources).items():
        if digest(a.original/'pin'/module) != sha:
            raise ValueError(f'Changed original input {module}')
    state=Path(read(a.original/'fine/invocation.json')['command'][5])
    if digest(state)!=rec['state']:
        raise ValueError('Changed preparation')
    base=Front(state,GRIDS['double'],1)
    out=dict(guard=digest(receipt),receipt_type=receipt.name,held=rec,generator=sys.argv,
             reader=digest(Path(__file__)),transmission={})
    for key in ('rest','wound'):
        p=saved/(key+'.npz')
        if digest(p)!=rec['solutions'][key]['sha256']:
            raise ValueError(f'Changed {key} field')
        with np.load(p) as f:
            if key=='rest':
                t=ContinuousTrial(base,state,f['x'],**rec['solid'],precondition=False)
                manufactured = manufactured_axial_traces(t)
                out['manufactured_axial_order_flux'] = manufactured
                write(a.output,out)
                if any(abs(v['measured']-v['expected']) > 1.e-9*max(1.,abs(v['expected']))
                       for v in manufactured.values()):
                    raise ArithmeticError('Manufactured one-sided axial order flux failed')
            u,q,load=f['field'],f['q'],f['load']
        out['transmission'][key]=dict(weak=weak(t,u,q,load),direct=direct(t,u))
        write(a.output,out)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for key in ('original','held','output'):
        p.add_argument(key,type=Path)
    a=p.parse_args()
    sys.path[:0]=[str(a.held/'pin'),str(a.original/'pin/software/stage1_solver/tools')]
    main(a)

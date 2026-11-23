"""Signed checks of the R4 finite strain, all free directions and front transfer."""
import argparse
import json
from pathlib import Path
from time import perf_counter

import numpy as np

from mouth_free_shape_solid import Solid
from mouth_joint_endpoint_audit import orientation
from mouth_last_ring_solid import source
from mouth_joined_front_audit import Front
from mouth_joint_polar_shape import Trial
from mouth_supercore_joined_front import file_sha


def run(args):
    sh, meridional, field, x = source(args.state,args.wall)
    with np.load(args.wall) as f:
        amplitude=float(f['internal_amplitude'])
    base=Front(args.state,(13,10,14),1)
    model=Trial(base,args.state,x,(meridional.shape[1]-1)//3,meridional.shape[0]-1,False).model
    print(json.dumps(dict(state=file_sha(args.state),wall=file_sha(args.wall),
          fields='wound order on the held native six-domain mesh',model_size=model.size)),flush=True)
    for degree in args.degree:
        b=Solid(args.state,sh,x[5],args.n,args.k,args.angles,model,degree=degree)
        q=b.seed(meridional,amplitude)
        grid_t,grid_x,grid_a=np.meshgrid(b.t,b.x,b.a,indexing='ij')
        perturb=np.zeros(b.chart_shape)
        for c in range(3):
            perturb[...,c]=.0018*(4*grid_x*(1-grid_x))**2*np.sin(np.pi*grid_t)**2*(
                np.cos(2*grid_a) if c!=2 else np.sin(2*grid_a))*(c+1)
        # Use a nonzero polar velocity and angularly varying axial/radial offsets.
        perturb[:,0,:,2]=perturb[:,1,:,2]
        q+=b.pack(perturb)
        dc=field[2].copy()
        if not args.common_only:
            tic=perf_counter()
            e,g,load,info=b.evaluate(q,dc,diagnose=True)
            elapsed=perf_counter()-tic
            off=b.evaluate(q,None,gradient=False,load=False)[0]
            transfer=load@dc
            orient=orientation(sh,x[5],b,b.unpack(q))
            print(json.dumps(dict(kind='state',degree=degree,grid=[args.n,args.k,args.angles],
                active=len(q),energy=e,forward_reverse_error=e-off-transfer,
                deposited=transfer,seconds=elapsed,**info,**b.residual(q,g),orientation=orient)),flush=True)
            if abs(e-off-transfer)>2e-8 or abs(transfer-info['transfer_energy'])>2e-10:
                raise ArithmeticError('Front transfer is not adjoint to material work')
            rng=np.random.default_rng(813)
            order_direction=rng.normal(size=model.size)*model.masks[2]
            h=2e-5
            v=[b.evaluate(q,dc+s*h*order_direction,gradient=False,load=False)[0] for s in (-1,1)]
            order_error=(v[1]-v[0])/(2*h)-load@order_direction
            print(json.dumps(dict(kind='conversion',degree=degree,analytic=load@order_direction,
                                  signed=(v[1]-v[0])/(2*h),error=order_error)),flush=True)
            if abs(order_error)>2e-5:raise ArithmeticError('Order finite work differs from deposition')
            for name in ('mixed','bond_radial','last_axial','angular','pole_radial'):
                d=np.zeros(b.chart_shape)
                if name=='mixed':
                    a=rng.normal(size=len(q))
                    a/=max(abs(a))
                else:
                    key={'bond_radial':(args.k//2,b.cut_node,1,0),'last_axial':(-2,-2,1,1),
                         'angular':(args.k//2,b.cut_node,1,2),'pole_radial':(args.k//2,b.cut_node,0,0)}[name]
                    d[key]=1.
                    a=b.pack(d)
                slope=g@a
                for h in (2e-6,1.e-6):
                    readings=[b.evaluate(q+s*h*a,dc,gradient=False,load=False)[0] for s in (-1,1)]
                    fd=(readings[1]-readings[0])/(2*h)
                    print(json.dumps(dict(kind='same_energy',degree=degree,mode=name,step=h,
                                  derivative=slope,signed=fd,difference=fd-slope)),flush=True)
                    if abs(fd-slope)>3e-4*max(1.,abs(slope)):raise ArithmeticError('Invalid material gradient')
        # Compare physically equivalent field basis against a refined mesh carrying
        # the exact piecewise-Q1 polynomial, with the same crystal material state.
        class Refined:
            pass
        from scipy.interpolate import RegularGridInterpolator
        fine=Refined();axes=[model.rr[0],model.a,model.z[1]]
        def split(a):return np.sort(np.r_[a,(a[1:]+a[:-1])/2])
        fx,fa,fz=[split(v) for v in axes]
        fine.a=fa;fine.rr=(fx,);fine.z=(None,fz)
        fine.size=len(fx)*len(fa)*len(fz);fine.ids={(1,0):np.arange(fine.size).reshape(len(fx),len(fa),len(fz))}
        exact=RegularGridInterpolator(axes,dc[model.ids[1,0]])(tuple(np.meshgrid(fx,fa,fz,indexing='ij'))).ravel()
        other=Solid(args.state,sh,x[5],args.n,args.k,args.angles,fine,degree=degree)
        ee,gg,ll,ii=other.evaluate(q,exact)
        if not args.common_only:
            print(json.dumps(dict(kind='mesh_only_remap',degree=degree,same_material_coordinates=True,
                        reference_volume_change=ii['reference_volume']-info['reference_volume'],
                        total_energy_change=ee-e,transfer_work_change=ll@exact-transfer,
                        peak_native_gradient_change=float(abs(gg-g).max()))),flush=True)
        common=Solid(args.state,sh,x[5],args.n,args.k,args.angles,model,degree=degree,
                     integration_axes=(fx,fa,fz))
        E,G,L,I=common.evaluate(q,dc)
        errors=dict(energy=ee-E, field_work=ll@exact-L@dc, volume=ii['reference_volume']-I['reference_volume'],
                    native_gradient_max=float(abs(gg-G).max()))
        print(json.dumps(dict(kind='common_crystal_partition',degree=degree,**errors)),flush=True)
        if max(abs(v) for v in errors.values())>2e-10:raise ArithmeticError('Crystal work changed on prolongation')


def benchmark(args):
    from types import SimpleNamespace
    from mouth_supercore_joined_front import mapped_nodes
    sh,d,field,x=source(args.state,args.wall)
    with np.load(args.wall) as f:
        amp=f['internal_amplitude'].item()
    # Trace topology and physical order grid of the doubled carrier. Timing the
    # isolated assembler needs only its duct trace and the held source chart.
    nr,na,nz=26,20,28
    rr,_=mapped_nodes(nr,sh.length,28.,22.)
    z=np.r_[np.linspace(0,sh.length-.12,nr+1)[:-1],np.linspace(sh.length-.12,sh.length,nz+1)]
    a=np.pi/2*np.linspace(0,1,na+1)**1.35
    size=len(rr[0])*len(z)*len(a)
    model=SimpleNamespace(rr=rr,z=(None,z),a=a,size=size,
                          ids={(1,0):np.arange(size).reshape(len(rr[0]),len(a),len(z))})
    from scipy.interpolate import RegularGridInterpolator
    base=Front(args.state,(13,10,14),1)
    original=Trial(base,args.state,x,(d.shape[1]-1)//3,d.shape[0]-1,False).model
    dc=RegularGridInterpolator((original.rr[0],original.a,original.z[1]),field[2,original.ids[1,0]])(
                              tuple(np.meshgrid(rr[0],a,z,indexing='ij'))).ravel()
    for n,k,nang in ((8,24,6),(11,32,10)):
        for degree,sub in ((2,1),(3,2)):
            b=Solid(args.state,sh,x[5],n,k,nang,model,degree=degree,subdivisions=sub)
            tic=perf_counter();q=b.seed(d,amp)
            en,g,load,info=b.evaluate(q,dc)
            print(json.dumps(dict(kind='timing',grid=[3*n,k,nang],degree=degree,subdivisions=sub,
                   seconds=perf_counter()-tic,nodes=len(q),energy=en,**info,**b.residual(q,g))),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--state',type=Path,default=Path('_scratch/s11c/supercore_joined_front01/adopt309/state.npz'))
    p.add_argument('--wall',type=Path,default=Path('_scratch/s11c/r4_fixed_wall01/stall/latest.npz'))
    p.add_argument('--n',type=int,default=4);p.add_argument('--k',type=int,default=10)
    p.add_argument('--angles',type=int,default=3)
    p.add_argument('--degree',type=int,nargs='+',default=[2,3])
    p.add_argument('--common-only', action='store_true')
    p.add_argument('--benchmark',action='store_true')
    args=p.parse_args()
    (benchmark if args.benchmark else run)(args)

"""Integrated finite-aperture readouts using the acoustic assembly quadrature."""
import numpy as np
from mouth_live_waist_fields import WIDTH


def disk_fluxes(op, omega, fields):
    """Lab and relative oscillatory current, including density at the moving cut.

    This is the Eulerian first perturbation on the reference disk plus the
    sampled base-current change from expanding or translating the section.
    The steady current changes sharply at the rim so quadrature and the value
    of its one-sided lip trace are deliberately exposed.
    """
    ids=op.capi
    length=op.shape.length
    cap_shape=np.array([length,0.,0.])
    integrated={key:np.zeros(3,complex) for key in
                ('through_minus','through_plus','relative_minus','relative_plus','layer_force')}
    quadrature=[]
    for r,area,N in op.cap_quad:
        w=np.full(r.shape,length)
        nminus,*_=op._fields(1,0,r,w)
        nplus,*_=op._fields(2,0,r,w)
        tr=np.stack([np.einsum('i,kij->kj',N,fields[inn]) for inn in
                     (ids[:,:2],ids[:,2:],op.density[ids[:,:2]],op.density[ids[:,2:]])],axis=1)
        disp=np.broadcast_to(cap_shape,(len(r),1,3))
        if op.sheet is None:
            flow=(tr[:,1]-tr[:,0])/(WIDTH/6)-1j*omega*cap_shape
            minus,plus=flow,flow
            normal=1j*omega*(tr[:,1]-tr[:,0])
        else:
            forces=op.sheet(omega,r)@np.concatenate([tr,disp],axis=1)
            minus,plus,normal=-forces[:,0],forces[:,1],forces[:,4]
        readings=dict(through_minus=minus,through_plus=plus,
                      relative_minus=minus+1j*omega*nminus[:,None]*cap_shape,
                      relative_plus=plus+1j*omega*nplus[:,None]*cap_shape, layer_force=normal)
        for key,value in readings.items():
            integrated[key]+=area@value
        quadrature.append(dict(r=r,area=area,nminus=nminus,nplus=nplus,**readings))
    # Eulerian sheet currents exclude the first variation of a moving domain of
    # integration. For constant background normal j both definitions differ by
    # 4 pi b^2 j(b) * delta b; the sharp R2 sheet supplies its own j(b−).
    rim_flux=op.base.disk(op.shape.exit).item() if op.base is not None else 0.
    area_change=np.array([0.,4*np.pi*op.shape.exit**2*rim_flux*op.exit_scale,0.])
    # Stationary div(j)=0 gives d/dw integral_{r<b} j_w dA = -4 pi b² j_r(b).
    # At the native zero-thickness junction these one-sided reconstructions can
    # disagree with the recovered cap trace, particularly at the transonic edge.
    shift={}
    for name,slab in [('minus',1),('plus',2)]:
        n,ur,*_=op._fields(slab,0,np.array([op.shape.exit]),np.array([op.shape.length]))
        shift[name]=np.array([-4*np.pi*op.shape.exit**2*(n*ur).item()*length,0.,0.])
    bins={key:sum(sample['area'][:,None]*sample[key] for sample in quadrature)/op.disk_weights[:,None]
          for key in integrated}
    return dict(**integrated,area_change=area_change,rim_base_flux=rim_flux,axial_sample_change=shift,
                material_section_minus=integrated['relative_minus']+area_change+shift['minus'],
                material_section_plus=integrated['relative_plus']+area_change+shift['plus'],
                bins=bins,quadrature=quadrature)


def profile(path):
    """Read a guarded tracked pole, reconstruct the driven field at its real frequency."""
    import hashlib
    import json
    from mouth_curved_sound import CurvedSound, mean_front
    from mouth_curved_steady import SteadyMouth, reference
    from mouth_poly_shallow_release import emit
    check=json.loads((path/'outcome.json').read_text())
    if check['exitCode'] or check['childOutcome']['exitCode'] or not check['limitsVerified']:
        raise ValueError('Unverified pole input')
    data=(path/'stdout').read_bytes()
    setup,measured,_=map(json.loads,data.splitlines())
    args=setup['args']
    shape,chi,native,scale,_,_,digest=reference(args['offset'])
    if digest!=setup['material_sha256']:
        raise ValueError('Mismatched saved preparation')
    st=SteadyMouth(shape,chi,native)
    j=st.info['terminal_current']/(4*np.pi*shape.exit**3/3)
    front=mean_front(j) if args['front']=='mean' else None
    op=CurvedSound(shape,scale,args['mesh'],args['box'],args['pml'],front,st,stretch=args['stretch'])
    w=measured['pole'][0]
    z=np.asarray(measured['mass_unit_vector']); v=z[:,0]+1j*z[:,1]
    _,X,_=op._solve(w)
    values=disk_fluxes(op,w,X)
    pair=lambda a: np.stack((np.real(a),np.imag(a)),axis=-1).tolist()
    exact={k:pair(values[k]@v) for k in ('through_minus','through_plus','relative_minus','relative_plus',
                                           'layer_force','material_section_minus','material_section_plus','area_change')}
    bins={k:pair(a@v) for k,a in values['bins'].items()}
    # Sum the core's phase-equation rows including both cap and feeding-end reactions.
    # Its row sum is i omega int delta_n across the fixed cavity. The side-wall
    # source is compiled separately (including advected original shape), so the
    # two projected opening currents, storage and moving-wall flow can be compared.
    core=op.grid[1,0]
    ids=np.unique(core)
    force=(op.K[2]+1j*w*op.G[2])@X
    edge=force-(op.B0-1j*w*op.Bv)
    upper=edge[core[-1]].sum(axis=0)@v
    inlet=edge[core[0]].sum(axis=0)@v
    side=(op.B0[ids]-1j*w*op.Bv[ids]).sum(axis=0)@v
    storage=force[ids].sum(axis=0)@v
    emit(kind='disk',material_sha256=digest,source_sha256=hashlib.sha256(data).hexdigest(),
         source=str(path),omega=w,exit_radius=shape.exit,mean_current=j,base_at_lip=values['rim_base_flux'],
         integrated=exact,bins=bins,radii=op.radii.tolist(),area=op.disk_weights.tolist(),
         moving_cut_axial_change={k:pair(a@v) for k,a in values['axial_sample_change'].items()},
         fixed_cavity=dict(storage_rate_iomega_dn=pair(storage),side_flux=pair(side),
                           feeding_end_outward_flux=pair(inlet),terminal_outward_flux=pair(upper),
                           closure=pair(storage-side-inlet-upper),
                           cap_current_crosscheck=pair(upper-values['through_minus']@v)))


def side_injection(path):
    """Receiver-side wall sweep, compared to the lab-frame flux at its disk."""
    import json
    from pathlib import Path
    from mouth_curved_steady import SteadyMouth,reference
    from mouth_poly_shallow_release import emit
    guard=json.loads((path/'outcome.json').read_text())
    if guard['exitCode'] or not guard['limitsVerified'] or guard['childOutcome']['exitCode']:
        raise ValueError('Unverified disk readout')
    f=json.loads((path/'stdout').read_text())
    src=Path(f['source']); inp,mode,_=map(json.loads,(src/'stdout').read_bytes().splitlines())
    args=inp['args']
    shape,chi,native,scale,_,_,digest=reference(args['offset'])
    if digest!=f['material_sha256']:
        raise ValueError('Mismatched reference')
    state=SteadyMouth(shape,chi,native)
    v=np.asarray(mode['mass_unit_vector'])@np.array([1.,1j])
    omega=mode['pole'][0]
    m=args['mesh']
    q,wt=np.polynomial.legendre.leggauss(3)
    t=(np.arange(m)[:,None]+(q+1)/2)/m
    r,s,_=shape.profile(t)
    b=state.base('bulk_beside',r,shape.length*t)[0]
    dx=scale*t*v[1]+4*shape.neck*t*(1-t)*v[2]-s*shape.length*t*v[0]
    geom=-1j*omega*4*np.pi*r**2*dx*shape.length/m*wt[None,:]/2
    side=(geom*b).sum(); cap=complex(*f['integrated']['through_plus'])
    pair=lambda z:[z.real,z.imag]
    emit(kind='injection',side_mass=pair(side),side_volume=pair(geom.sum()),
         disk_mass_lab=pair(cap),net_mean_surface_mass_source=pair(side+cap),
         fraction_of_wall_source_magnitude=abs(side+cap)/abs(side),
         net_over_sum_magnitudes=abs(side+cap)/(abs(side)+abs(cap)),
         omega=omega,mesh=m,material_sha256=digest,source=str(path))


if __name__=='__main__':
    import argparse
    from pathlib import Path
    p=argparse.ArgumentParser()
    p.add_argument('source_directory',type=Path)
    p.add_argument('--receiver-source',action='store_true')
    args=p.parse_args()
    (side_injection if args.receiver_source else profile)(args.source_directory)

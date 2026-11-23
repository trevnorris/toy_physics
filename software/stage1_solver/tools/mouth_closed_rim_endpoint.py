"""Audit the actual relaxed angular R4 states on matched fixed-lattice controls."""
from argparse import ArgumentParser
from pathlib import Path
from json import loads

import numpy as np

from mouth_closed_rim_work import Closed
from mouth_closed_rim_surface_review import Reader, quadrature
from mouth_fixed_wall_angular import AngularBonded
from mouth_joint_polar_shape import Trial
from mouth_joined_front_audit import Front
from mouth_polar_exit_readout import accepted
from mouth_supercore_joined_front import check_guard, emit, file_sha


def wall_reaction(reader, axial, depth0, n):
    """Oriented physical slit stress on the two held, coincident sealed faces."""
    m=reader.model; L=m.shape.length
    zs=m.z[1]; z,w=quadrature(np.r_[depth0,zs[(zs>depth0)&(zs<L)],L],n)
    if depth0<=.62*L:
        raise ValueError('Exposed wall overlaps polar distortion')
    a,v=quadrature(m.a,n)
    zz,aa=np.meshgrid(z,a,indexing='ij')
    b,s=m.shape.radius(zz)
    area=4*np.pi*b*b*np.sin(aa)*w[:,None]*v[None,:]*axial(zz)
    f=[]
    for bank,sign in ((0,1.),(1,-1.)):
        tr=reader.sample(1,bank,b,aa,zz,np.array([1.,0.,0.]))
        tz=reader.sample(1,bank,b,aa,zz,np.array([0.,0.,1.]))
        f.append(sign*(area.ravel()@(tr-s.ravel()[:,None]*tz)))
    return np.array(f)


def main(a):
    for p in (a.state,a.endpoint,a.controls): check_guard(p)
    prep,record,x,_,original=accepted(a.state,a.snapshot)
    read=lambda path:[loads(s) for s in (path.parent/'stdout').read_text().splitlines()]
    certified=lambda p:any(row.get('sha256')==file_sha(p) for row in read(p))
    if not (certified(a.endpoint) and certified(a.controls)):
        raise ValueError('Completed endpoint or matched controls not in physics readouts')
    with np.load(a.controls) as f:
        rest=f['rest'].copy()
        if f['checkpoint_sha'].item()!=record['sha256'] or not np.array_equal(f['wound'],original):
            raise ValueError('Original same-shape controls mismatched')
    with np.load(a.endpoint) as f:
        if (f['input_sha'].item()!=record['sha256'] or not np.array_equal(f['x'],x)
                or f['state_sha'].item()!=file_sha(a.state)):
            raise ValueError('Unmatched relaxed angular endpoint')
        d,u,s=f['displacement'],f['field'],f['internal_amplitude'].item()
    base=Front(a.state,tuple(prep['field_grid']),1)
    n,k=(d.shape[1]-1)//3,d.shape[0]-1
    tr=Trial(base,a.state,x,n,k,True)
    m=tr.model
    b=AngularBonded(a.state,n,k,shape=tr.shape,epsilon=x[5],full_depth=True)
    b.internal_amplitude=s
    m.material_load=b.get_material_load(d[b.mask]*b.scale,m)
    drive_error=np.max(abs(m.norm(m.vars(m.evaluate(u)[1]))))
    rest_start=np.max(abs(m.norm(m.vars(m.evaluate(rest)[1]))))
    if drive_error>2e-7:raise ArithmeticError('Drive changed since its endpoint')
    if rest_start>=2e-7:
        rest,metrics=m.solve(rest=True,initial=rest)
    rest_error=np.max(abs(m.norm(m.vars(m.evaluate(rest)[1]))))
    if rest_error>=2e-7 or max(abs(rest[1]))>1e-10:raise ArithmeticError('Unsettled same-lattice rest')
    np.savez_compressed(a.save,wound=u,rest=rest,checkpoint_sha=file_sha(a.endpoint),original_sha=record['sha256'])
    emit(kind='matched_control',state=file_sha(a.state),endpoint=file_sha(a.endpoint),
         input_control=file_sha(a.controls),controls_sha=file_sha(a.save),solid_cells=[3*n,k],
         ledger='wound drive minus wound rest on identical relaxed lattice',
         angular_amplitude=s,drive_error=drive_error,rest_initial=rest_start,rest_final=rest_error,
         field_change=float(np.max(abs(u-original))),fixed_rest_reused=bool(rest_start<2e-7))
    study=Closed(m,b,d,[u,rest])
    reader=Reader(m,base.c,u,rest)
    for ix in (2,3):
        ax,_,_,geo=study.lifts(ix,4)
        z0=m.shape.length-geo['lower_half_thickness'][-1]
        wall=wall_reaction(reader,ax,z0,6)
        other=wall_reaction(reader,ax,z0,10)
        if abs(wall-other).max()>2e-6:raise ArithmeticError('Slit stress needs quadrature refinement')
        for j in (4,5):
            row=study.audit(ix,j)
            emit(**row, held_slit_work_banks_mass_and_axial=other.tolist(),
                 sharp_slit_total_axial=other[:,1].sum(),
                 outer_plus_slit_axial=row['outward_momentum']+other[:,1].sum(),
                 slit_quadrature_difference=float(abs(wall-other).max()))


if __name__=='__main__':
    p=ArgumentParser(description=__doc__)
    for name in ('state','snapshot','controls','endpoint','save'):
        p.add_argument(name,type=Path)
    main(p.parse_args())

"""Conservative rim audit on the rejected radius-four wall at fixed material displacement."""
import argparse
import json
from pathlib import Path

import numpy as np

from mouth_closed_rim_endpoint import wall_reaction
from mouth_closed_rim_surface_review import Reader
from mouth_closed_rim_work import Closed
from mouth_disk_mesh_rematch import remap
from mouth_fixed_wall_angular import AngularBonded
from mouth_joined_front_audit import Front
from mouth_joint_polar_shape import Trial
from mouth_polar_joined import PolarCoupled
from mouth_supercore_joined_front import check_guard, emit, file_sha


GRIDS = {'native': (13, 10, 14), 'double': (26, 20, 28), 'local': (32, 20, 34)}


def bisect(nodes, toward_high):
    """Insert midpoints of six contiguous cells alongside one end of this axis."""
    mid = (nodes[1:] + nodes[:-1]) / 2
    add = mid[-6:] if toward_high else mid[:6]
    return np.sort(np.r_[nodes, add])


def local_axes(model):
    rr = (bisect(model.rr[0], True), bisect(model.rr[1], False))
    L = model.shape.length
    lower = np.searchsorted(model.z[1], L-model.h-1.e-13)
    upper = np.searchsorted(model.z[2], L+model.h-1.e-13)
    z = (bisect(model.z[0],True),
         np.r_[bisect(model.z[1][:lower+1],True)[:-1], bisect(model.z[1][lower:],True)],
         np.r_[bisect(model.z[2][:upper+1],False)[:-1], bisect(model.z[2][upper:],False)])
    return rr, z


class RimCarrier(PolarCoupled):
    """Start with the established structured topology, then assemble on the supplied axes.

    Joined numbers degrees of freedom by the integer cell indices and builds
    its angularly mobile order mask from the terminal-band cell counts. The
    common junctions and those counts are verified before assembling any FE
    weights. Only the physical coordinate of an otherwise identical node changes.
    """
    def __init__(self,*args,axes=None,**kwargs):
        self.nested_axes=axes
        super().__init__(*args,**kwargs)

    def matrix(self, order, elastic_density):
        if self.nested_axes is not None:
            rr,z=self.nested_axes
            for before,after in zip((*self.rr,*self.z),(*rr,*z)):
                if len(before)!=len(after) or np.any(np.diff(after)<=0) or not np.array_equal(before[[0,-1]],after[[0,-1]]):
                    raise ValueError('New axes alter structured topology or geometry')
            for i, plane in ((1,self.shape.length-self.h),(2,self.shape.length+self.h)):
                orig=np.flatnonzero(np.isclose(self.z[i],plane,atol=1.e-13,rtol=0))
                new=np.flatnonzero(np.isclose(z[i],plane,atol=1.e-13,rtol=0))
                if len(orig)!=1 or not np.array_equal(orig,new):
                    raise ValueError('Order band or native cut indices changed')
            self.rr,self.z=rr,z
        super().matrix(order,elastic_density)


class SameWindows(Closed):
    """Evaluate the original native piecewise-linear lift on a nested grid."""
    def __init__(self,*args,reference,**kwargs):
        super().__init__(*args,**kwargs)
        self.reference=reference

    def lifts(self,k,j):
        m=self.m
        ref=self.reference
        L=m.shape.length
        iz=np.flatnonzero(np.isclose(ref.z[1],L-m.h,rtol=0,atol=1e-12))[0]
        jo=np.flatnonzero(np.isclose(ref.z[2],L+m.h,rtol=0,atol=1e-12))[0]
        lo1,lo2=ref.z[1][iz-k:iz-k+2]
        hi1,hi2=ref.z[2][jo+k-1:jo+k+1]
        x1,x2=ref.rr[1][j:j+2]
        for knots, edges in ((m.z[1],[lo1,lo2]),(m.z[2],[hi1,hi2]),(m.rr[1],[x1,x2])):
            if any(np.min(abs(knots-edge))>2.e-13 for edge in edges):
                raise ArithmeticError('Enclosure left the measured nodal planes')
        def axial(w):
            return np.clip((w-lo1)/(lo2-lo1),0,1)*np.clip((hi2-w)/(hi2-hi1),0,1)
        a=axial(self.znodes)
        r=np.ones_like(a)
        for (slab,side),idx in m.ids.items():
            if side:
                r[idx]=np.broadcast_to(np.clip((x2-m.rr[1])/(x2-x1),0,1)[:,None,None],idx.shape)
        geo=dict(lower_half_thickness=[L-lo2,L-lo1],upper_half_thickness=[hi1-L,hi2-L],
                 lateral_radius_at_disk=(m.shape.exit+(28.-m.shape.exit)*np.array([x1,x2])).tolist())
        return axial,a,r,geo


def read_source(state, checkpoint):
    check_guard(state)
    p=checkpoint.parent
    guard=json.loads((p/'outcome.json').read_text())
    limits=json.loads((p/'limit-validation.json').read_text())
    if not guard['limitsVerified'] or not limits['verified'] or int(limits['actual']['memory.swap.max']) or guard['exitCode']!=1:
        raise ValueError('Requires the guarded incomplete fixed-wall trajectory')
    rows=[json.loads(s) for s in (p/'stdout').read_text().splitlines()]
    last=next(r for r in reversed(rows) if r.get('kind')=='fixed_wall_sweep')
    with np.load(checkpoint) as f:
        x,d,u,s=f['x'],f['displacement'],f['field'],float(f['internal_amplitude'])
        if (f['state_sha'].item()!=file_sha(state) or f['input_sha'].item()!=rows[0]['snapshot']
                or not np.array_equal(x,rows[0]['coordinates'])
                or json.loads(f['readout'].item())['energy']!=last['energy'] or last['sweep']!=2):
            raise ValueError('Rejected shape or latest fully saved sweep differs')
    return x,d,u,s


def run(a):
    x,d,old,s=read_source(a.state,a.checkpoint)
    base=Front(a.state,GRIDS['native'],1)
    tr=Trial(base,a.state,x,11,32,a.level=='native')
    ref=tr.model
    if a.level=='native':
        m=ref
    else:
        axes=None
        if a.level=='local':
            # Generate only inexpensive physical axes at the double level;
            # the uniform formulas and both terminal-band joins match Joined.
            from types import SimpleNamespace
            from mouth_exchange_carrier import mapped_nodes
            nr,_,nz=GRIDS['double']
            rr,zz=mapped_nodes(nr,tr.shape.length,28.,22.)
            L,h=tr.shape.length,ref.h
            below=np.r_[np.linspace(0,L-h,nr+1)[:-1],np.linspace(L-h,L,nz+1)]
            above=np.r_[np.linspace(L,L+h,nz+1)[:-1],L+h+(22.-L-h)*(np.sinh(2*np.linspace(0,1,nr+1))/np.sinh(2))**2]
            axes=local_axes(SimpleNamespace(rr=rr,z=(zz[0],below,above),shape=tr.shape,h=h))
        m=RimCarrier(tr.shape,base.c,1,GRIDS[a.level],polar_shift=x[5],axes=axes,
                     terminal_profile=ref.terminal_profile)
    for fine0,native in zip((*m.rr,m.a,*m.z),(*ref.rr,ref.a,*ref.z)):
        if any(np.min(abs(fine0-val))>2e-13 for val in native):
            raise ValueError('Missing native node on refined mesh')
    b=AngularBonded(a.state,11,32,shape=tr.shape,epsilon=x[5])
    b.internal_amplitude=s
    m.material_load=b.get_material_load(d[b.mask]*b.scale,m)
    common=dict(checkpoint_sha=file_sha(a.checkpoint),grid=GRIDS[a.level],level=a.level,
                inside_x=m.rr[0],outside_x=m.rr[1],feed_z=m.z[0],duct_z=m.z[1],receiver_z=m.z[2],angles=m.a)
    emit(kind='input',checkpoint=common['checkpoint_sha'],state=file_sha(a.state),level=a.level,
         grid=GRIDS[a.level],solid_cells=[33,32],wall_x=x.tolist(),target_depth=tr.shape.length,
         disk_radius=tr.shape.exit,node_count=m.size,held_internal_amplitude=s,
         local_last_inside_cell=tr.shape.exit*(1.-m.rr[0][-2]),
         local_first_outside_cell=(28.-tr.shape.exit)*m.rr[1][1],
         local_duct_depth_cell=m.z[1][-1]-m.z[1][-2],local_receiver_depth_cell=m.z[2][1]-m.z[2][0])
    pair=[]
    for rest in (False,True):
        init=old if a.level=='native' else remap(ref,m,old)
        updated,info=m.solve(rest=rest,initial=init)
        measure=m.measure(updated)
        emit(kind='matched',rest=rest,**measure,**info)
        if measure['port_error']>3e-5 or info['stationarity']>=2e-7:
            raise ArithmeticError('Unresolved pressure/phase/front or mass balance')
        pair.append(updated)
        if not rest:
            corr=updated[2,m.ids[1,0]]
            dst=a.folder/'solid_load.npz'
            np.savez_compressed(dst,duct_dc=corr,**common)
            emit(kind='solid_load',sha256=file_sha(dst))
    p=a.folder/'fields.npz'
    np.savez_compressed(p,wound=pair[0],rest=pair[1],duct_dc=corr,**common)
    emit(kind='fields',sha256=file_sha(p))
    # Sparse stationary solves no longer needed by weak/surface readers.
    m.solvers=[]
    study=SameWindows(m,b,d,pair,reference=ref)
    reader=Reader(m,base.c,*pair)
    for k in (2,3):
        taper,_,_,geo=study.lifts(k,4)
        low=tr.shape.length-max(geo['lower_half_thickness'])
        flux=wall_reaction(reader,taper,low,10)
        other=wall_reaction(reader,taper,low,6)
        if np.max(abs(flux-other))>2e-6:
            raise ArithmeticError('Wall force not quadrature resolved')
        for j in (4,5):
            row=study.audit(k,j)
            emit(**row,slit_banks_mass_axial=flux.tolist(),slit_force=flux[:,1].sum(),
                 closed_plus_slit=row['outward_momentum']+flux[:,1].sum())


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    for key in ('state','checkpoint','folder'):
        parser.add_argument(key,type=Path)
    parser.add_argument('level',choices=tuple(GRIDS))
    run(parser.parse_args())

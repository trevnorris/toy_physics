"""One-sided solid traction on the held terminal disk in the paired field audit."""
import argparse
import json
from pathlib import Path

import numpy as np

from mouth_bonded_front_adjoint import angular_weights
from mouth_closed_rim_surface_review import quadrature
from mouth_joined_front_audit import Front
from mouth_rim_shared_crystal import TerminalQuadrature
from mouth_supercore_joined_front import check_guard,emit,file_sha


def audit(state,wall,fields):
    check_guard(fields)
    records=[json.loads(l) for l in (fields.parent/'stdout').read_text().splitlines()]
    if not any(r.get('kind')=='fields' and r['sha256']==file_sha(fields) for r in records):
        raise ValueError('Unverified paired field')
    base=Front(state,(13,10,14),1)
    from mouth_exchange_carrier import mapped_nodes
    with np.load(fields) as f:
        if f['checkpoint_sha'].item()!=file_sha(wall):
            raise ValueError('Different material or terminal disk')
        ax,z=f['inside_x'],f['duct_z']
        material=TerminalQuadrature(state,wall,ax,z)
        # Native radial weight, retained under both field refinements.
        qx=mapped_nodes(13,material.L,28.,22.)[0][0]
        u=np.clip((qx-.90)/.10,0,1)
        q=1-u*u*(3-2*u)
        from mouth_last_ring_surface_review import FieldChart
        m=FieldChart(material.sh,base.c,1,f['grid'],precondition=False)
        disk=m.ids[1,0][:,:,-1]
        if not np.array_equal(disk,m.ids[2,0][:,:,0]):
            raise ValueError('Shared disk has non-identical field nodes')
        delta=(f['wound'][2,disk]-f['rest'][2,disk])@angular_weights(f['angles'])
        if np.max(abs(material.d[-1]))!=0:
            raise ValueError('Surface measure assumes the held, flat terminal disk')
        prep=(f['load_sha'].item() if 'load_sha' in f else 'archived native quadrature')
    # At t=1 both solid and registered triangle seams fall on their axial-row
    # vertices. Real material labels on the terminal face have zero displacement.
    ref=np.exp(-material.p*material.chart.a)
    knots=np.unique(np.r_[material.xt,ax,qx,ref[(ref>material.cut)&(ref<1.)]])
    emit(kind='material_trace_source',wall=file_sha(wall),fields=file_sha(fields),load=prep)
    for n in (4,8,16):
        x,wx=quadrature(knots,n)
        v=material._material(x,np.ones_like(x),np.ones_like(x))
        J,K=v['J'],v['K']
        if np.max(abs(v['x']-x))>1e-12 or np.max(abs(v['z']-material.L))>1e-12:
            raise ValueError('Invalid terminal trace')
        chi=np.interp(x,ax,delta)
        sigma=chi*(K[:,1,1]-1+.5*J*(J-1))/J
        area=4*np.pi*material.sh.exit**3*x*x*wx
        g=np.interp(x,qx,q)
        emit(kind='solid_face',quadrature=n,full_solid_traction=area@sigma,
             clamped_solid_traction=(area*g)@sigma,
             complementary_edge_solid_traction=(area*(1-g))@sigma,
             minimum_face_jacobian=J.min(),minimum_meridional_stretch=np.sqrt(np.linalg.eigvalsh(K)[:,0]).min())


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('state','wall','fields'):p.add_argument(name,type=Path)
    a=p.parse_args()
    audit(a.state,a.wall,a.fields)

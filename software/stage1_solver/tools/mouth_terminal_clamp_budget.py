"""Measure the complementary clamp work for comparing a bowed disk with whole-rim forces.

The companion variation moves the complementary disk profile and the surrounding
bank constraints. This is a virtual reaction measurement on the fixed shape.
"""
import argparse
import json
from itertools import product
from pathlib import Path

import numpy as np
from mouth_bonded_front_adjoint import angular_weights
from mouth_closed_rim_surface_review import Reader
from mouth_joined_front_audit import Front
from mouth_joint_polar_shape import Trial
from mouth_last_ring_surface_review import FieldChart
from mouth_last_rings_fluid import GRIDS,read_source
from mouth_rim_consistent_fields import checked_pair
from mouth_rim_shared_crystal import TerminalQuadrature
from mouth_terminal_disk_work import DiskMotion,crystal_work,carrier_work
from mouth_supercore_joined_front import file_sha


class RimExtension(DiskMotion):
    def __init__(self,ref,uniform=False):
        super().__init__(ref)
        self.uniform=uniform
        del self.ramps['short']
    def profile(self,r,z):
        if self.uniform:
            return np.ones_like(r),np.zeros_like(r),np.zeros_like(r)
        g,gr,gz=super().profile(r,z)
        return 1-g,-gr,-gz


def outside(reader,ref,motion):
    m=reader.model
    x1,x2=ref.rr[1][4:6]
    out={s:0. for s in motion.ramps}
    gauss=(1+np.array([-1,1])/np.sqrt(3))/2
    for slab in (1,2):
        ax,az=m.rr[1],m.z[slab]
        ii=np.flatnonzero(ax[:-1]<x2-1e-14)
        jj=np.flatnonzero((az[1:]>motion.ramps['far'][0])&(az[:-1]<motion.ramps['far'][-1]))
        I,A,J=np.meshgrid(ii,np.arange(len(m.a)-1),jj,indexing='ij')
        I,A,J=I.ravel(),A.ravel(),J.ravel()
        dx,da,dz=np.diff(ax)[I],np.diff(m.a)[A],np.diff(az)[J]
        for xx,aa,zz in product(gauss,repeat=3):
            x,a,z=ax[I]+xx*dx,m.a[A]+aa*da,az[J]+zz*dz
            R,slope=m.shape.radius(z);width=28.-R
            r=R+x*width
            wt=4*np.pi*r*r*np.sin(a)*width*dx*da*dz/8
            Tr=reader.sample(slab,1,r,a,z,np.array([1,0,0]))[:,1]
            Tz=reader.sample(slab,1,r,a,z,np.array([0,0,1]))[:,1]
            q=np.clip((x2-x)/(x2-x1),0,1)
            qr=-((x>x1)&(x<x2)).astype(float)/(width*(x2-x1))
            qz=-(1-x)*slope*qr
            for name,(lo,lmid,hmid,hi) in motion.ramps.items():
                rise=np.clip((z-lo)/(lmid-lo),0,1)*np.clip((hi-z)/(hi-hmid),0,1)
                grad=((z>lo)&(z<lmid))/(lmid-lo)-((z>hmid)&(z<hi))/(hi-hmid)
                out[name]-=wt@(Tr*rise*qr+Tz*(rise*qz+grad*q))
    return out


def run(state,wall,fields,disk):
    checked_pair(fields)
    x,_,_,_=read_source(state,wall)
    base=Front(state,GRIDS['native'],1)
    tr=Trial(base,state,x,11,32,False)
    digest=file_sha(fields)
    disk_data=json.loads(disk.read_text())
    if disk_data['fields']!=digest:raise ValueError('Unmatched bowed disk')
    with np.load(fields) as a:
        m=FieldChart(tr.shape,base.c,1,a['grid'],precondition=False,terminal_profile=tr.model.terminal_profile)
        m.rr=tuple(a[k] for k in ('inside_x','outside_x'))
        m.z=tuple(a[k] for k in ('feed_z','duct_z','receiver_z'))
        pair=[a['wound'],a['rest']]
        if not np.array_equal(m.a,a['angles']):raise ValueError('Wrong angles')
        dc=(pair[0][2]-pair[1][2])[m.ids[1,0]]
        delta=np.einsum('raj,a->jr',dc,angular_weights(m.a))
    reader=Reader(m,base.c,*pair)
    comp=RimExtension(tr.model)
    inside=carrier_work(reader,comp)
    exterior=outside(reader,tr.model,comp)
    quad=TerminalQuadrature(state,wall,m.rr[0],m.z[1],64)
    solid=crystal_work(quad,delta,comp)
    whole=crystal_work(quad,delta,RimExtension(tr.model,True))
    workers=[json.loads(l) for l in (fields.parent/'stdout').read_text().splitlines()]
    results={}
    for name,k in [('near',2),('far',3)]:
        F=next(r for r in workers if r['kind']=='closed_work' and r['axial_band']==k and r['lateral_band']==4)
        disk_row=disk_data['reactions'][name]
        complement_fluid=inside[name][5]+exterior[name]
        comp_material=solid[name][4]
        Wcomp=complement_fluid+comp_material
        Wdisk=disk_row['physical_actuator']
        uniform_solid=whole[name][4]
        uniform_fluid=F['signed_common_translation_work']['1e-05']
        closure=Wdisk+Wcomp-uniform_solid-uniform_fluid
        if abs(closure)>3e-5:raise ArithmeticError(f'Incompatible disk and rim virtual motions {closure}')
        results[name]=dict(disk=Wdisk,complement=Wcomp,complement_fluid=complement_fluid,
                           complement_solid=comp_material,full_terminal_virtual=uniform_fluid+uniform_solid,
                           full_solid=uniform_solid,external_virtual=uniform_fluid,sloping_banks=F['slit_force'],
                           measured_closed=F['closed_plus_slit'],split_error=closure,
                           closed_minus_full_terminal=F['closed_plus_slit']-(uniform_fluid+uniform_solid))
    return dict(fields=digest,measured=results,
                scope='Full virtual test moves all terminal disk sites and continues along the held banks. The bowed disk test fixes the banks. Free original interior crystal nodes stay fixed in both terminal tests; the common outer translation convects all of them.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for key in ('state','wall','fields','disk','output'):p.add_argument(key,type=Path)
    a=p.parse_args()
    a.output.write_text(json.dumps(run(a.state,a.wall,a.fields,a.disk),indent=2,allow_nan=False)+'\n')

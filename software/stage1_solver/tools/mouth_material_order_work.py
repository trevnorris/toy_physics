"""Signed work of independent material and carrier trial spaces on the unchanged R4 duct."""
import argparse
import json
from pathlib import Path
import numpy as np

from mouth_closed_rim_surface_review import Reader, quadrature
from mouth_disk_pillbox import dot
from mouth_rim_bank_work import Bank
from mouth_rim_disk_work_review import Disk
from mouth_rim_pinned_comparison import validate_guard, digest
from mouth_supercore_joined_front import ORDER_INERTIA, state_from


def measure(a):
    prior = json.loads(a.manifest.read_text())
    validate_guard(a.job)
    case = next(c for c in prior['cases'] if c['inputs']['grid'] == a.job.name)
    meta = case['inputs']
    if (digest(a.job/'zero.npz') != case['zero_state_sha256'] or digest(a.original) != meta['field_sha256'] or
        digest(a.state) != meta['state_sha256'] or digest(a.wall) != meta['wall_sha256']):
        raise ValueError('Geometry or field hash mismatch')
    disk = Disk(a.state,a.wall,a.original)
    with np.load(a.job/'zero.npz') as read:
        disk.reader = Reader(disk.model,state_from(a.state)[1],read['wound'],read['rest'])
    m = disk.model
    lo,mid,hi = meta['terminal_knots']
    b4,b8 = [Bank(disk,0,(lo,mid,hi),edge) for edge in meta['duct_bounds']]
    normalized_edge = b8.edge
    zaxis = m.z[1]
    readouts = []
    for nq in (2,4):
        xs,wx = quadrature(m.rr[0][m.rr[0]>=normalized_edge],nq)
        aa,wa = quadrature(m.a,nq)
        zz,wz = quadrature(np.unique(np.r_[lo,mid,hi,zaxis[(zaxis>lo)&(zaxis<hi)]]),nq)
        x,angle = np.meshgrid(xs,aa,indexing='ij')
        x,angle = x.ravel(),angle.ravel()
        weights = (wx[:,None]*wa).ravel()
        steps = np.array([2e-5,-2e-5,1e-5,-1e-5])
        signed = np.zeros((len(steps),3))
        for z,dz in zip(zz,wz):
            R=disk.sh.radius(z)[0]
            r=x*R
            eta,er,ez = np.subtract(b4.motion(r,z),b8.motion(r,z))
            basis = list(b4.physical(x,angle,z))
            volume = 4*np.pi*weights*R*r*r*np.sin(angle)*dz
            for j,h in enumerate(steps):
                jac=1+h*ez
                z1=z+h*eta
                x1=r/disk.sh.radius(z1)[0]
                shifted=list(b4.physical(x1,angle,z1))
                for sign, before, after in zip((1.,-1.),basis,shifted):
                    (_,p,c,n,_,_),amp,gamp,F = before
                    C = after[0][2]
                    def pull(g):
                        g1=g.copy()
                        g1[:,2] /= jac
                        g1[:,0] -= h*er*g1[:,2]
                        return g1
                    p1=pull(p)
                    v1=pull(gamp)
                    y=amp**2
                    interaction=F**10*(y-1)**2*(y**3+2*y**2+3*y+4)/20
                    def projection(gp,gc):
                        return -.5*ORDER_INERTIA*n*dot(gp,gc)**2/(n+ORDER_INERTIA*dot(gc,gc))
                    changes = [
                        .5*F**2*(jac*dot(v1,v1)-dot(gamp,gamp))+(jac-1)*interaction,
                        .5*n*(jac*dot(p1,p1)-dot(p,p)),
                        jac*projection(p1,C)-projection(p,c),
                    ]
                    signed[j] -= sign*(np.array(changes)@volume)/h
        for h, out in zip((steps[0],steps[2]),(signed[:2],signed[2:])):
            readouts.append(dict(gauss=nq,step=h,positive=out[0].tolist(),negative=out[1].tolist(),
                                 centered_sectors=out.mean(axis=0).tolist(),
                                 centered_total=out.mean(axis=0).sum()))
    return dict(level=a.job.name,original=meta,zero_sha256=case['zero_state_sha256'],
                sign='Minus driven-minus-resting energy difference per unit step.',
                sectors=['amplitude+native density','phase','phase/order projection'],
                order_well_and_gradient_work=0.,registered_crystal_work=0.,
                readouts=readouts)


if __name__ == '__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    for key in ('state','wall','original','job','manifest'):
        parser.add_argument(key,type=Path)
    print(json.dumps(measure(parser.parse_args()), allow_nan=False))

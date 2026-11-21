"""Find freely attached thin order fronts on saved four-dimensional carrier fields.

The stationary front satisfies nu*div(N)=u.N. Use t=r/a(w), shooting
outward from its center, and measure orthogonal contact at the wall.
"""
import argparse
from hashlib import file_digest
import json
from pathlib import Path

import numpy as np
from scipy.integrate import simpson

from mouth_curved_phase import quadratic, tensor
from mouth_flared_receiver import EXIT, EXIT_RADIUS, ROUNDING, FlaredDomain

NU = .4
ARCHIVE = Path("research/pde_ledger_v3/_measurements/mouth_order_and_flare.json")


class Carrier:
    def __init__(self, folder):
        field_file = folder / "field.npz"
        with field_file.open("rb") as stream:
            field_hash = file_digest(stream, "sha256").hexdigest()
        assert field_hash == json.loads(ARCHIVE.read_text())["output_field_sha256"][str(field_file)]
        row = json.loads((folder/"stdout").read_text().splitlines()[0])
        self.domain = FlaredDomain(row["n"], row["profile"], row["roof"])
        with np.load(field_file) as saved:
            np.testing.assert_array_equal(saved["points"], self.domain.points)
            self.field = saved["psi"]
        self.label = folder.name

    def radius(self, z):
        a, slope = self.domain.radius(z)
        low = z < ROUNDING
        if np.any(low):
            dz = z[low]-ROUNDING
            root = np.sqrt(ROUNDING**2-dz**2)
            a[low], slope[low] = self.domain.rc-root, dz/root
        return a, slope

    def eval(self, t, z):
        a, ap = self.radius(z)
        r = t*a
        f, fr, fw = np.empty((3, z.size), dtype=complex)
        low = z < ROUNDING
        if np.any(low):
            f[low], fr[low], fw[low] = self.domain.read(self.field, r[low], z[low])
        for block, mask in ((2,(z >= ROUNDING) & (z < EXIT)), (3,z >= EXIT)):
            if not np.any(mask):
                continue
            pg, vg, dofs = self.domain.blocks[block]
            p, v = z[mask], 1-np.broadcast_to(t,z.shape)[mask]
            i = np.clip(pg.searchsorted(p,side="right")-1,0,len(pg)-2)
            j = np.clip(vg.searchsorted(v,side="right")-1,0,len(vg)-2)
            hp, hv = (pg[i+1]-pg[i])/2, (vg[j+1]-vg[j])/2
            bp, dp = quadratic((p-pg[i])/hp-1)
            bv, dv = quadratic((v-vg[j])/hv-1)
            nodes = self.field[dofs[i*(len(vg)-1)+j]]
            f[mask] = np.sum(tensor(bp,bv)*nodes,axis=1)
            fp = np.sum(tensor(dp,bv)*nodes,axis=1)/hp
            fv = np.sum(tensor(bp,dv)*nodes,axis=1)/hv
            fr[mask] = -fv/a[mask]
            fw[mask] = fp-fr[mask]*(1-v)*ap[mask]
        return a, ap, f, (fr/f).imag, (fw/f).imag


def geometry(t, z, s, carrier):
    a, ap, f, ur, uw = carrier.eval(t,z)
    cw = np.sqrt(1-s*s)
    p = -s/cw
    dr = a/(1-t*ap*p)
    un = ur*s+uw*cw
    return a, ap, f, cw, p, dr, un


def rhs(t,z,s,carrier):
    a, _, _, _, p, dr, un = geometry(t,z,s,carrier)
    source = un/(3*NU) if t==0 else un/NU-2*s/(t*a)
    return p*dr, source*dr


def shoot(centers, carrier, steps, record=False):
    dt = 1/steps
    z, s = centers.copy(), np.zeros_like(centers)
    zs, ss = [z.copy()], [s.copy()]
    for k in range(steps):
        t=k*dt
        hz, hs = rhs(t,z,s,carrier)
        jz, js = rhs(t+dt/2,z+dt*hz/2,s+dt*hs/2,carrier)
        kz, ks = rhs(t+dt/2,z+dt*jz/2,s+dt*js/2,carrier)
        lz, ls = rhs(t+dt,z+dt*kz,s+dt*ks,carrier)
        z += dt*(hz+2*jz+2*kz+lz)/6
        s += dt*(hs+2*js+2*ks+ls)/6
        if record:
            zs.append(z.copy())
            ss.append(s.copy())
    _, ap = carrier.radius(z)
    boundary = ap*np.sqrt(1-s*s)-s
    if record:
        return z,boundary,np.array(zs),np.array(ss)
    return z,boundary


def measurements(c, steps):
    scan=np.linspace(.53,3.11,517)
    wall, residual=shoot(scan,c,steps)
    indices=np.where((residual[:-1]*residual[1:] < 0) & (wall[:-1]<EXIT) & (wall[1:]<EXIT))[0]
    lo, hi=scan[indices],scan[indices+1]
    left=residual[indices]
    for _ in range(21):
        mid=(lo+hi)/2
        _,res=shoot(mid,c,steps)
        hi=np.where(res*left<=0,mid,hi)
        lo=np.where(res*left>0,mid,lo)
    centers=(lo+hi)/2
    wall, mismatch, z, s=shoot(centers,c,steps,True)
    # The weighted area first variation along an interior-critical cap is
    # positive weight times boundary mismatch times displacement of wall contact.
    delta=.0005
    lower, rlow=shoot(centers-delta,c,steps)
    upper, rhigh=shoot(centers+delta,c,steps)
    results=[]
    ts=np.linspace(0,1,steps+1)
    for i,center in enumerate(centers):
        a, ap, f, cw, _, dr, un=geometry(ts,z[:,i],s[:,i],c)
        n=abs(f)**2
        area=4*np.pi*(ts*a)**2*dr/cw
        flux=simpson(area*n*un,x=ts)
        weight=np.exp(-np.angle(f)/NU)
        earea=simpson(area*weight,x=ts)
        coeff=un/n
        dpdz=(rhigh[i]-rlow[i])/(upper[i]-lower[i])
        row=dict(center=center,contact=wall[i],contact_radius=a[-1],
                 contact_orthogonality_error=mismatch[i]/np.hypot(1,ap[-1]),
                 min_normal_velocity=un.min(), integrated_current=flux,
                 work_coefficient_min=coeff.min(), work_coefficient_max=coeff.max(),
                 relative_span=(coeff.max()-coeff.min())/coeff.mean(),
                 effective_weighted_area=earea, boundary_slope=dpdz,
                 depth_class="minimum" if dpdz>0 else "maximum")
        results.append(row)
    print(json.dumps(dict(field=c.label, steps=steps, scan_center=[scan[0],scan[-1]],
                          scan_interval=scan[1]-scan[0], fronts=results)),flush=True)
    assert np.all(np.isfinite(residual)) and np.all(np.diff(wall)>0)
    assert all(abs(row["integrated_current"]-.25)<.0025 and
               abs(row["contact_orthogonality_error"])<.0001 for row in results)


if __name__ == "__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("folders",type=Path,nargs="+")
    parser.add_argument("--steps",type=int,default=96)
    args=parser.parse_args()
    for directory in args.folders:
        measurements(Carrier(directory),args.steps)

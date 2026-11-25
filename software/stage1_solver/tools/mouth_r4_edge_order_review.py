"""Independent direct 65-digit local energies for the three selected order probes."""
from pathlib import Path
import hashlib
import json
import sys
from tempfile import TemporaryDirectory

import mpmath as mp
import numpy as np

ROOT = Path(__file__).resolve().parents[3]
STUDY = ROOT / '_scratch/s11c/r4_newton_held01'
PIN = STUDY / 'pin/software/stage1_solver/tools'
sys.path.insert(0, str(PIN))
from mouth_continuous_order import I, KAPPA, WELL
from mouth_continuous_depth_mesh import axes_for_depth
from mouth_continuous_shape_relax import ShapeTrial
from mouth_r4_eulerian_solid import SpatialCrystal
from mouth_r4_smooth_fields import SmoothOrder
from mouth_r4_smooth_held import SmoothHeld

SITES = (36353, 45532, 62811)
STATE = {'100': ('rematched', '63a2de2ae45c6c2b197573ac70b2f035142f538d59908bb068ceda87b0570542'),
         '108': ('last_active', 'a2210361aa5b38712de7db618b5220985e12f39f3a089351dfd3767bef0a7178')}
mp.mp.dps = 65


def hfloat(x):
    return mp.mpf(float(x))


def held_row(matrix, row, field, registry, site, h):
    start, stop = matrix.indptr[row:row+2]
    return mp.fsum(hfloat(w)*(hfloat(field[j]) + (h if registry[j] == site else 0))
                   for w, j in zip(matrix.data[start:stop], matrix.indices[start:stop]))


def local_energy(t, fields, indices, site, h):
    m = t.model
    mask = t.number
    I0 = hfloat(I)
    # Use the minimum of a positive quadratic in the eliminated velocity.
    # This avoids differencing any factored perturbation formula.
    parts = []
    for row in indices:
        v = held_row(m.V, row, fields[2], mask, site, h)
        amp = held_row(m.V, row, fields[0], mask, -1, 0)
        n = (hfloat(m.D[row, 4])*amp)**2
        gc = [held_row(g, row, fields[2], mask, site, h) for g in m.G]
        gp = [held_row(g, row, fields[1], mask, -1, 0) for g in m.G]
        gg = mp.fdot(gc, gc)
        z = I0*mp.fdot(gc, gp)/(n+I0*gg)
        kinetic = (n*mp.fsum((p-z*c)**2 for p,c in zip(gp,gc)) + n*n*z*z/I0)/2
        density = kinetic + hfloat(KAPPA)*gg/2 + hfloat(WELL)*v*v*(1-v)**2 + hfloat(m.D[row,-1])*v
        parts.append(hfloat(m.hq[row])*density)
    return mp.fsum(parts)


def old_gauge_work(solver, chi, site, h):
    col = solver.Vg[:,site].tocoo()
    rows = col.row
    v = solver.Vg[rows]@chi
    smooth = lambda x: (1-x)**2*(1+2*x)
    return mp.fsum(hfloat(solver.work[i]) * (smooth(hfloat(x)+hfloat(w)*h)-smooth(hfloat(x)))
                   for i,w,x in zip(rows,col.data,v))


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    manifest=STUDY/'pinned.json'
    assert sha(manifest)=='22a06d53f4b64bd6d6be19fc2c8c2f8d038909ebd5f7357c9a327941d386820b'
    assert all(sha(PIN/k)==v for k,v in json.loads(manifest.read_bytes()).items())
    h=hfloat(1e-9)
    result={'command': sys.orig_argv, 'mpmath_digits': mp.mp.dps, 'increment': float(h),
            'selected_joined_nodes': list(SITES), 'step': {}, 'state_sha256': {s: v[1] for s,v in STATE.items()}}
    with TemporaryDirectory(dir=Path(__file__).parent) as work:
        driver = SmoothHeld(STUDY.parent, Path(work), 1e-8)
        t = ShapeTrial(driver.base, driver.prep, driver.x, 8,24,6,3,2, precondition=False,
                       axes=axes_for_depth(driver.shape(driver.x).length,driver.reference_depth))
        t.model.extended_energy_sum = True
        solid = SpatialCrystal(t.b)
        order_sites=t.number[t.b.ids.transpose(0,2,1)]
        previous_fields=None
        for step,(name,expected) in STATE.items():
            src=STUDY/'primary'/(name+'.npz')
            assert sha(src)==expected
            with np.load(src,allow_pickle=False) as saved:
                labels,fields=saved['labels'],saved['field']
                assert np.array_equal(saved['x'], driver.x) and int(saved['radius'])==4 and int(saved['winding'])==1
            if previous_fields is not None:
                assert np.array_equal(previous_fields,fields)
            previous_fields=fields
            chi=fields[2,t.first]
            gauge=SmoothOrder(t,labels,1e-8)
            probes={}
            for i in SITES:
                assert chi[i]==0 and t.free[i] and np.array_equal(chi[t.number],fields[2])
                rows=t.V[:,i].tocoo().row
                for G in t.G:
                    rows=np.union1d(rows,G[:,i].tocoo().row)
                delta_bath=local_energy(t,fields,rows,i,h)-local_energy(t,fields,rows,i,0)
                # At fixed labels the crystal energy is linear in order. Evaluate
                # its energy on the single joined order basis, omitting all
                # other coefficients, independently of the deposited source.
                unit_crystal=solid.evaluate(labels,(order_sites==i).astype(float),False)[0]
                ext=mp.fsum(hfloat(x) for x in t.ext[t.number==i])
                physical=delta_bath + h*(hfloat(unit_crystal)+ext)
                gauge_work=old_gauge_work(gauge,chi,i,h)
                probes[str(i)]=dict(touched_quadrature=len(rows), bath=float(delta_bath), exterior=float(h*ext),
                    interior_crystal=float(h*hfloat(unit_crystal)),physical=float(physical),
                    old_stencil_unit_gauge=float(gauge_work),
                    old_total={str(e):float(physical+hfloat(e)*gauge_work) for e in (1e-8,1e-10)})
            result['step'][step]=probes
    print(json.dumps(result,indent=2,allow_nan=False))


if __name__=='__main__':
    main()

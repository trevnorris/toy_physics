"""Check the phase boost and the base mass flux independently of the two drives."""
import argparse
import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np
from numpy.polynomial.legendre import leggauss

from mouth_carrier_free_center import motion_source
from mouth_carrier_lateral_response import Dipole, prepare, trace
from mouth_flared_receiver import EXIT
from mouth_occupied_push_baths import digest, load_bath
from mouth_steady_tow_partition import BASE, DATA, measure


def exact_mass_port(dom, f, block, along, fixed, axis):
    """Mass current integrated on the outer port directly from the base field."""
    pg, vg, _ = dom.blocks[block]
    edges = pg if along == 'p' else vg
    root, weights = leggauss(8)
    half = np.diff(edges)/2
    q = (edges[:-1, None]+half[:, None]*(1+root)).ravel()
    p, v = (q, np.full_like(q, fixed)) if along == 'p' else (np.full_like(q, fixed), q)
    r, _, rp, rv, wp, wv, val, fr, fw = trace(dom, block, p, v, f[:, None])
    tr, tw = (rp, wp) if along == 'p' else (rv, wv)
    nr, nw = tw, -tr
    sign = np.sign(nr if axis == 0 else nw)
    fn = sign*(nr*fr[:, 0]+nw*fw[:, 0])
    return float(4*np.pi*np.sum(r*r*(val[:, 0].conj()*fn).imag*
                                (half[:, None]*weights).ravel()))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('far', choices=(9,18), type=int)
    far = parser.parse_args().far
    if far == 9:
        base = Path('_scratch/s11c/mouth_flared_cubic_n34')
        dom, f = prepare(SimpleNamespace(source=base))
    else:
        base = BASE
        dom, f = load_bath(base)
    faces = [('feed',1,'v',dom.far,0), ('bulk_side',4,'p',1.,0),
             ('roof_inner',3,'v',dom.roof,1), ('roof_outer',4,'v',dom.roof,1)]
    flux = {key: exact_mass_port(dom,f,block,along,val,axis) for key,block,along,val,axis in faces}
    record = dict(far=far, field_sha256=digest(base/'field.npz'), base_outward_mass_traces=flux,
                  analytic_product_boost_medium=-sum(flux.values()))
    if far == 9:
        old = json.loads((DATA/'mouth_carrier_0pn.json').read_text())
        run = old['runs']['carrier_inertia_c34_01']
        assert run['exit_code'] == run['worker_exit'] == 0 and run['limits_verified']
        assert digest(Path(run['stdout_path'])) == run['stdout_sha256']
        rows = [json.loads(s) for s in run['stdout_lines']]
        assert rows[0]['sha256'] == record['field_sha256']
        record['old_stream_rows'] = [r for r in rows if r['test'] == 'steady_balance']
        record['uniform_nodal_boost'] = measure(dom,f,1j*dom.points[:,0]*f)
        packed = motion_source(dom,f)
        n = len(f)
        model = Dipole(dom,f)
        h, residual = model.solve(np.zeros_like(f), source=packed[:n]+1j*packed[n:])
        record['lab_chart_tow'] = measure(dom,f,h)
        record['lab_residual'] = residual
    print(json.dumps(record, allow_nan=False))


if __name__ == '__main__':
    main()

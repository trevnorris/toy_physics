"""One-sided single-site energy slopes as missing shear nodes reappear in a shallow pocket."""
import argparse
import hashlib
import json

import numpy as np

from mouth_exchange_carrier import BIAS
from mouth_live_waist_active_wave import ActiveWave
from mouth_live_waist_order import Order
from mouth_live_waist_reopen import reopen
from mouth_sealed_shape_relax import Meridian
from mouth_shallow_puncture_fields import ROOT, emit, read_case
from mouth_shallow_puncture_order import loads


def check(n):
    prep, _, _, _, _, _, source = read_case(1, n)
    folder = ROOT / f"_scratch/s11c/shallow_puncture_fields01/bounded_r1_m{n}"
    snapshot = folder/"center.npz"
    sha = hashlib.sha256(snapshot.read_bytes()).hexdigest()
    if sha not in (folder/"stdout").read_text():
        raise ValueError("Missing central state fingerprint")
    with np.load(snapshot) as data:
        if source != str(data["original_snapshot_sha256"].item()):
            raise ValueError("Order preparations differ")
        chi, vec, pole = data["order"], data["vector"], complex(data["pole"].item())
        shape = Meridian(1, data["coordinates"])
    order = Order(shape,n)
    mat = ActiveWave.install(loads(shape,n,prep["field_mesh"],(pole,vec),pole.real))
    wave = mat.wave
    fo, go = order.energy_gradient(chi)
    fc, gc, _ = mat.carrier.response(chi, BIAS)
    proj = wave.project.tocsc()
    weights = proj@chi.T.ravel()
    mass = wave.cut_m + wave.matrix(weights[:,None,None]*wave.mc)
    active = mass.diagonal()>0
    if np.any(vec[~active] != 0):
        raise ValueError("Nonzero displacement on absent wave degree of freedom")
    block = (wave.port.block(pole+2e-6)-wave.port.block(pole-2e-6))/4e-6
    edge = vec[wave.port.ids]
    den = 2*pole*(vec@(mass@vec))+edge@(block@edge)
    q = wave.kc-pole**2*wave.mc
    results = []
    for i,j in np.argwhere((chi==0)&order.free):
        col = j*(2*n+1)+i
        sl = slice(proj.indptr[col],proj.indptr[col+1])
        direct, correction, new, ext, _ = reopen(proj.indices[sl],proj.data[sl],wave,q,vec,active)
        results.append((go[i,j]+gc[i,j]+((direct+correction)/den).real,i,j,new,ext))
    lowest, i,j,new,ext = min(results,key=lambda r:r[0])
    v = vec.copy()
    v[new] = ext
    differences = []
    for h in (.002,.001):
        c = chi.copy()
        c[i,j] += h
        wave.seed = pole, v
        wave.target = pole.real
        ef,_,_ = mat.evaluate(c,1.,BIAS)
        differences.append((order.energy_gradient(c)[0]+ef-fo-fc-pole.real)/h)
    slopes, *_ = zip(*results)
    emit(kind="one_site_reopening",mesh=n,center_sha256=sha,zero_sites=len(slopes),
         newly_supported_cases=sum(len(r[3])>0 for r in results), minimum=lowest,
         maximum=max(slopes),below_minus_1e_6=sum(bool(y < -1e-6) for y in slopes),
         smallest_location=(int(i),int(j)),signed_positive_steps=differences)
    if abs(2*differences[1]-differences[0]-lowest)>2e-6:
        raise ArithmeticError("One-sided direct and extended-mode derivative disagree")
    if lowest < -1e-6:
        raise ArithmeticError("A single local reopening lowers energy")


if __name__ == "__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mesh",type=int,choices=(48,72),required=True)
    check(parser.parse_args().mesh)

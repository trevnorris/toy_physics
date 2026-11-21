"""Propagation checks for a torsional pulse on the prepared four-dimensional mouths.

Test the radial radiation channel of the feeding slab and the magnitude of
the optional material advection in the same order fields as the pulse.
"""
import argparse
from hashlib import file_digest
import json
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.optimize import brentq
from scipy.special import spherical_jn, spherical_yn

from mouth_free_contact import Carrier
from mouth_free_order_diffuse import read

MANIFEST = Path("research/pde_ledger_v3/_measurements/mouth_free_order_depth.json")


def emit(**data):
    print(json.dumps(data,allow_nan=False),flush=True)


def radiation():
    # X = U(r) (e_z x e_r). Natural radial shear traction from curl energy
    # is proportional to U' + U/r. The depth-uniform transverse mode in a
    # homogeneous feeding slab has U=h_1^(1)(kr) and exp(-i*omega*t).
    x = np.array([.05,.1,.25,.5,1.,2.,4.,8.])
    h = spherical_jn(1,x) + 1j*spherical_yn(1,x)
    dh = spherical_jn(1,x,True) + 1j*spherical_yn(1,x,True)
    impedance = (dh/h + 1/x) / 1j
    expected = (x*x + 1j/x)/(1+x*x)
    emit(test="feed_radiation",k_times_radius=x.tolist(),
         real_impedance_in_plane_wave_units=impedance.real.tolist(),
         imaginary_impedance=impedance.imag.tolist(),
         explicit_hankel_error=float(max(abs(impedance-expected))))
    np.testing.assert_allclose(impedance,expected,atol=2e-12,rtol=1e-12)


def evaluate_carrier(flow, z, eta):
    a,_=flow.radius(z)
    _,_,_,ur,uw=flow.eval(np.broadcast_to(eta,z.shape).copy(),z)
    return a,ur,uw


def local_trials(folder, source_hashes):
    rows=[json.loads(line) for line in (folder/"stdout").read_text().splitlines()]
    setup,rest=rows[0],rows[1:]
    flow=Carrier(Path(setup["source"]).parent)
    roots,weights=leggauss(32)
    rs=.24*(1+roots)
    zs=.25+.2*roots
    rr,ww=np.broadcast_arrays(rs[None,:],zs[:,None])
    static_shape=(1-((ww-.25)/.2)**2)**4*(1-(rr/.48)**2)**4
    weight=np.outer(weights,weights)*rr**4*static_shape**2
    z_axis=np.linspace(.04,2.75,601)
    _,_,uw=evaluate_carrier(flow,z_axis,0.)
    for row in rest:
        if row["seed"] not in (.44,1.55):
            continue
        path=folder/f"order_seed{row['seed']}.npz"
        with path.open("rb") as stream:
            assert file_digest(stream,"sha256").hexdigest()==source_hashes[str(path)]
        with np.load(path) as field:
            chi=field["chi"]
        # Entire initial velocity shape lies in the common rounded entrance.
        psi,pr,pw=flow.domain.read(flow.field,rr,ww)
        order=flow.domain.read(chi,rr,ww)[0]
        speed=np.hypot((pr/psi).imag,(pw/psi).imag)
        assert np.min(order)>.5
        mach=speed/np.sqrt(order)
        mean=np.sum(weight*mach)/np.sum(weight)
        on_axis=np.array([read(flow.domain,chi,0.,z)[0] for z in z_axis])
        # Constant transverse inertia: c_T²=chi in ordered speed units.
        # Near roots all chi are resolved positive. Avoid division by tail noise.
        diff=on_axis-uw**2
        indices=np.where((diff[:-1]>0)&(diff[1:]<=0)&(uw[:-1]>0))[0]
        crossings=[]
        for i in indices:
            def zero(z):
                _,_,v=evaluate_carrier(flow,np.array([z]),0.)
                return read(flow.domain,chi,0.,z)[0]-v[0]**2
            z=brentq(zero,z_axis[i],z_axis[i+1])
            _,_,v=evaluate_carrier(flow,np.array([z]),0.)
            crossings.append(dict(w=z,carrier_velocity=v[0],chi=read(flow.domain,chi,0.,z)[0]))
        emit(test="optional_advected_characteristics",folder=str(folder),seed=row["seed"],
             weak_kick_energy_weighted_mach=mean,weak_kick_max_mach=float(mach.max()),
             outward_axial_characteristic_stalls=crossings, peak_axis_velocity=float(uw.max()),
             order_proportional_inertia_maximum_axis_speed_ratio=float(uw[on_axis>1e-8].max()))
        assert 0<mean<1 and crossings


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("folders",nargs="+",type=Path)
    args=parser.parse_args()
    radiation()
    hashes=json.loads(MANIFEST.read_text())["diffuse_orders_sha256"]
    for directory in args.folders:
        local_trials(directory,hashes)

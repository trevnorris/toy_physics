"""Check what a fixed drain current determines about a sideways rim perturbation."""

from hashlib import sha256
from pathlib import Path

import numpy as np
import sympy as sp


def inspect_feed():
    seed = Path("_scratch/s11c/mouth_feed_depth01/fields/n32_restoring_slope_u0.npz")
    manifest = Path("docs/_measurements/mouth_material_readout.txt")
    pin = f"{sha256(seed.read_bytes()).hexdigest()}  {seed}"
    if pin not in manifest.read_text().splitlines():
        raise ValueError("Resting carrier packet differs from the audited preparation")
    with np.load(seed, allow_pickle=False) as packet:
        current, mouth, dx = packet["current"], packet["mouth"], float(packet["dx"])
        rate = dx**3 * current.sum()
        print(f"Fixed carrier packet: shape={packet['field'].shape} rate={rate:.16g} "
              f"reservoir_sideways_speed={float(packet['speed']):g} arrays={sorted(packet.files)}")
        if (packet["field"].shape != (65, 32, 32, 32) or float(packet["speed"]) != 0
                or not np.isclose(rate, .25, atol=1e-14, rtol=0)
                or not np.all(current[mouth] > 0) or np.any(current[~mouth] != 0)):
            raise ValueError("Inward carrier preparation or aperture support mismatch")


def check_mass_balance():
    # Work on the brane mass law in its own measure. No conversion of a carrier
    # number rate into a brane mass source is assumed.
    r, a, width, theta = sp.symbols("r a ell theta", positive=True)
    b = sp.Symbol("b", real=True)
    # A_phi=b*f(r)*sin(theta) in the band |r-a|<ell, zero outside.
    # ell is small enough for the band to exclude the origin and box faces.
    f = (r-a)*(1-((r-a)/width)**2)**4
    jr = 2*b*f*sp.cos(theta)/r
    jt = -b*sp.diff(r*f, r)*sp.sin(theta)/r
    divergence = sp.simplify(sp.diff(r**2*jr,r)/r**2
                             + sp.diff(sp.sin(theta)*jt,theta)/(r*sp.sin(theta)))
    edges = [sp.simplify(expr.subs(r,edge))
             for edge in (a-width, a+width) for expr in (f, sp.diff(f,r))]
    normal = sp.simplify(jr.subs(r,a))
    along_drive = sp.simplify(-jt.subs({r:a,theta:sp.pi/2}))
    print("Band divergence, rim normal mass-current change, equatorial tangential change:",
          divergence, normal, along_drive)
    print("Profile and first derivative on the two band edges:", edges)
    j, v, dj, dv, e = sp.symbols("j v delta_j delta_v epsilon")
    carried = sp.diff((j+e*dj)*(v+e*dv),e).subs(e,0)
    print("Linearized outward carried momentum; at fixed pointwise material loss:",
          carried, carried.subs(dj,0))
    if (divergence != 0 or normal != 0 or along_drive != b or any(edges)
            or sp.expand(carried - j*dv - dj*v) != 0):
        raise ArithmeticError("Continuity or exchange calculation failed")


if __name__ == "__main__":
    inspect_feed()
    check_mass_balance()

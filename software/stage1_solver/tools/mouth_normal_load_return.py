"""Read the feed's w-momentum at the aperture and resolve a radial return load."""

import json
from pathlib import Path

import numpy as np
import sympy as sp

from lateral_source_rate import verified_packet


def symbolic_checks():
    # Isotropic carrier stress in solver units, psi=u exp(i theta).
    w = sp.symbols("w")
    u, theta = (sp.Function(name)(w) for name in ("u", "theta"))
    n = u**2
    polar = n*theta.diff(w)**2 - n*sp.diff(sp.log(n), w, 2)/4
    amp_form = (u.diff(w)**2 - u*u.diff(w, 2))/2 + n*theta.diff(w)**2
    residual = sp.simplify(polar - amp_form)
    print("normal stress identity residual =", residual, flush=True)
    assert residual == 0

    # Between the two shells the homogeneous solutions are constant + 1/r.
    # The jump is the signed force applied to the slab in the +w direction.
    r, a, b, k = sp.symbols("r a b k", positive=True)
    f, reaction, c, d, e = sp.symbols("F R C D E", real=True)
    near, far = c + d/r, e/r
    def force(profile):
        return sp.factor(-4*sp.pi*k*r**2*sp.diff(profile, r))
    eqs = [force(near).subs(r, a)-f, (near-far).subs(r,b),
           (force(far)-force(near)).subs(r,b)+reaction]
    coeff = sp.solve(eqs, (c,d,e), dict=True)[0]
    near, far = (sp.factor(h.subs(coeff)) for h in (near,far))
    print("inside reaction shell =",near,"; outside =",far,flush=True)
    print("slope-inferred loads =",force(near),force(far),flush=True)
    assert sp.simplify(force(near)-f)==0 and sp.simplify(force(far)-f+reaction)==0
    assert far.subs(reaction,f)==0 and sp.simplify((near-far).subs(r,b))==0


def aperture(grid, hashes):
    path = Path(f"_scratch/s11c/lateral_bath_{grid}_01/fields/d0.npz")
    with verified_packet(path, hashes) as packet:
        if int(packet["d"]) != 0 or len(packet["currents"]) != 1:
            raise ValueError("Expected a stationary isolated mouth")
        f, j = packet["field"], packet["currents"][0]
        x, r, dz = packet["x"], packet["r"], float(packet["dz"])
    dx = x[1]-x[0]
    weight = np.maximum(1-x[:,None]**2-r**2, 0)**3
    mask = weight > 0
    face = f[-1]
    n = abs(face)**2
    measure = 2*np.pi*r*dx**2
    def integral(v):
        return np.sum((v*measure)[mask])

    if not np.all(n>0) or abs(integral(j)-0.25)>1e-12 or np.any(j[~mask]):
        raise ValueError("Invalid density, aperture or signed number current")
    pressure = np.expm1(5*np.log(n))/5
    carried = j**2/n
    # Independent face readouts. First use the one-sided curvature of log n.
    log_layers = np.log(abs(f[-1:-6:-1])**2)
    log_curvature = np.tensordot([35,-104,114,-56,11], log_layers, axes=1)/(12*dz**2)
    quantum_direct = -n*log_curvature/4

    # Then use the steady equation for psi_ww and the saved normal-current
    # law for psi_w. Only interior lateral cells enter the aperture integral.
    pad_x = np.pad(face,((1,1),(0,0)),mode="edge")
    lap_x = np.diff(pad_x,n=2,axis=0)/dx**2
    radial_flux = np.zeros((len(x),len(r)+1),complex)
    radial_flux[:,1:-1] = (r[:-1]+dx/2)*np.diff(face,axis=1)/dx
    lap_r = np.diff(radial_flux,axis=1)/(r*dx)
    if np.any(mask[[0,-1],:]) or np.any(mask[:,-1]):
        raise ValueError("Aperture touches lateral bath")
    fww = np.expm1(4*np.log(n))*face/2 - lap_x - lap_r
    fw = (-weight*(1-1/np.sqrt(n)) + 1j*j/n)*face
    total_pde = pressure + (abs(fw)**2 - (face.conj()*fww).real)/2
    readings = [integral(v) for v in (carried,pressure,quantum_direct,
                                       carried+pressure+quantum_direct,total_pde)]
    print(f"{grid} J={integral(j):.12g} carriage={readings[0]:+.12g} "
          f"pressure={readings[1]:+.12g} gradient={readings[2]:+.12g} "
          f"total_one_sided={readings[3]:+.12g} total_field_equation={readings[4]:+.12g}",
          flush=True)
    if not np.isfinite(readings).all():
        raise ValueError("Nonfinite normal momentum readout")


if __name__ == "__main__":
    symbolic_checks()
    runs = json.loads(Path("docs/_measurements/lateral_reservoir_pair.json").read_text())["runs"]
    hashes = {Path(p):h for run in runs for p,h in run.get("fields",{}).items()}
    for tag in ("e24_c4","e24_c6","e48_c4"):
        aperture(tag,hashes)

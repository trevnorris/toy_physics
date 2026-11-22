"""Small-strain registered 4D crystal with a shared static and transverse-wave shear modulus."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.integrate import simpson
from scipy.interpolate import PchipInterpolator
from scipy.linalg import eigh
from scipy.optimize import brentq

ROOT = Path(__file__).resolve().parents[3]
REFERENCE = ROOT / "research/pde_ledger_v3/_measurements/mouth_ordered_strip_bend.json"
MU = 1.       # Same pure-order inertia and speed normalization as the transverse benchmark.
RHO_T = 1.
K_REF = 1.    # Explicit comparison in carrier pressure units, not a measured lattice tangent.


def flat_sample():
    source = json.loads(REFERENCE.read_text())["runs"]["main"]["results"]
    native = source["native"]
    s, chi = [], []
    for i, region in enumerate(native["flat"]["stress_profiles"]):
        begin = int(i > 0)
        s.extend(region["depth"][begin:])
        chi.extend(region["fraction"][begin:])
    return PchipInterpolator(s, np.maximum(0., chi)), native["tangent"], source["comparison"]


def moduli(bulk):
    lam = bulk - MU / 2.  # Four spatial dimensions.
    normal = lam + 2 * MU
    # In a uniaxial experiment each of the three independent lateral tractions vanishes.
    poisson = lam / (2 * MU + 3 * lam)
    young = 2 * MU * (1 + poisson)
    one_free = normal - lam**2 / normal
    return dict(bulk=bulk, lame=lam, poisson_uniaxial=poisson, free_young=young,
                one_free_normal=one_free, all_transverse_held=normal)


def ramp_mode(omega, t, duration):
    """Response to unit half-cosine force ramp followed by a hold, scaled to static response."""
    omega = omega[:, None]
    drive = np.pi / duration
    during = np.minimum(t, duration)[None, :]
    den = omega**2 - drive**2
    resonant = abs(den) < 1.e-9 * drive**2
    ratio = omega**2 / np.where(resonant, 1., den)
    c, s = np.cos(omega * during), np.sin(omega * during)
    h = .5 * (1. - ratio * np.cos(drive * during) + (ratio - 1.) * c)
    v = .5 * (ratio * drive * np.sin(drive * during) - (ratio - 1.) * omega * s)
    h = np.where(resonant, .5 * (1. - c - .5 * omega * during * s), h)
    v = np.where(resonant, .25 * omega * (s - omega * during * c), v)
    after = np.maximum(0., t - duration)[None, :]
    return (1. + (h - 1.) * np.cos(omega * after) + v / omega * np.sin(omega * after),
            -(h - 1.) * omega * np.sin(omega * after) + v * np.cos(omega * after))


def setup(profile, refine, material):
    # Preserve every archived sample as a node on both grids.
    z = np.linspace(profile.x[0], profile.x[-1], (len(profile.x) - 1) * refine + 1)
    h = np.diff(z)
    q, wt = np.polynomial.legendre.leggauss(5)
    xq = z[:-1, None] + h[:, None] * (q + 1.) / 2.
    shell = h[:, None] * wt / 2. * profile(xq)
    int_chi = shell.sum(axis=1)
    I0, I2 = float(int_chi.sum()), float((shell * xq**2).sum())
    stiffness = np.zeros((len(z), len(z)))
    mass = np.zeros_like(stiffness)
    drive = np.zeros_like(z)
    basis = np.array([(1. - q) / 2., (1. + q) / 2.])
    for j, width in enumerate(h):
        band = np.s_[j:j + 2, j:j + 2]
        stiffness[band] += (material["all_transverse_held"] * int_chi[j] / width**2
                            * np.array([[1., -1.], [-1., 1.]]))
        mass[band] += RHO_T * ((basis * shell[j]) @ basis.T)
        drive[j:j + 2] += material["lame"] * int_chi[j] / width * np.array([-1., 1.])
    eigenvals, modes = eigh(stiffness, mass)
    if abs(eigenvals[0]) > 1.e-6 or eigenvals[1] <= 0. or abs(drive.sum()) > 1.e-10:
        raise ArithmeticError("Unbalanced normal force or unexpected translation modes")
    return z, I0, I2, np.sqrt(eigenvals[1:]), modes[:, 1:], drive


def hold(material, grid, amplitude, ramp, profile):
    z, i0, _, omega, modes, drive = grid
    t = np.linspace(0., ramp + 20., int((ramp + 20.) * 100) + 1)
    during = np.minimum(t, ramp)
    extension = amplitude * .5 * (1. - np.cos(np.pi * during / ramp))
    speed = np.where(t < ramp, amplitude * .5 * np.pi/ramp * np.sin(np.pi * during / ramp), 0.)
    h, v = ramp_mode(omega, t, ramp)
    g = modes.T @ drive
    equilibrium = -g / omega**2
    coordinates = amplitude * equilibrium[:, None] * h
    velocities = amplitude * equilibrium[:, None] * v
    bare = material["all_transverse_held"] * i0
    load = bare * extension + g @ coordinates
    total_energy = (.5 * bare * extension**2 + extension * (g @ coordinates)
                    + .5 * (velocities**2 + (omega[:, None] * coordinates)**2).sum(axis=0))
    work = simpson(load * speed, x=t)
    equilibrium_slope = bare + g @ equilibrium
    late = t >= ramp + 10.
    u = modes @ coordinates[:, -1]
    du = np.diff(u) / np.diff(z)
    centers = (z[1:] + z[:-1]) / 2.
    local_axial = profile(centers) * (material["all_transverse_held"] * amplitude
                                      + material["lame"] * du)
    analytic = material["one_free_normal"] * i0
    sample = np.arange(0, len(t), 100)
    result = dict(amplitude=amplitude, ramp=ramp, hold=20., times=t[sample].tolist(),
                  end_load=load[sample].tolist(), end_energy=total_energy[-1],
                  actuator_work=float(work), late_mean_slope=float(np.mean(load[late])/amplitude),
                  late_mean_ratio=float(np.mean(load[late])/(analytic * amplitude)),
                  late_max_relative_ripple=float(np.max(abs(load[late]/(analytic * amplitude)-1.))),
                  relative_energy_work_error=float(abs(work/total_energy[-1]-1.)),
                  hold_energy_range_relative=float(np.ptp(total_energy[t>=ramp])/total_energy[-1]),
                  static_slope=equilibrium_slope, stress_at_end={
                      str(s): float(np.interp(s, centers, local_axial)) for s in (-.6,-.5,-.25,0.,.25,.5,.6)})
    return result


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--refine", type=int, choices=(1, 2), required=True)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    profile, front, prior = flat_sample()
    baseline = moduli(K_REF)
    grid = setup(profile, args.refine, baseline)
    _, i0, i2, *_ = grid
    support90 = 2.*brentq(lambda h: profile.integrate(-h, h) - .90*i0, 0., profile.x[-1])
    alternatives = []
    for bulk in (.001, .25, K_REF, 2.):
        row = moduli(bulk)
        for fixture in ("free_young", "one_free_normal", "all_transverse_held"):
            row["bending_tangent_" + fixture] = front + row[fixture] * i2
        alternatives.append(row)
    trials = [hold(baseline, grid, a, 6., profile) for a in (-.002, -.001, .001, .002)]
    trials.append(hold(baseline, grid, .002, 1., profile))
    readings = dict(checkpoint="d0ecb87d", profile_source=str(REFERENCE.relative_to(ROOT)),
                    profile_source_sha256=hashlib.sha256(REFERENCE.read_bytes()).hexdigest(),
                    geometry="flat 4D strip; 2 tangential transverse strains held, incremental normal crystal traction free",
                    profile_law="mu=chi and rho_T=chi, K=K0*chi; fixed flat order in reference coordinates at linear order",
                    shear_modulus_ordered=MU, transverse_inertia_ordered=RHO_T,
                    transverse_speed=np.sqrt(MU/RHO_T), reference_material=baseline,
                    incremental_lattice_bulk_is_assumed=True, refinements=args.refine,
                    force_thickness=i0, centered_width_90pct_static_load=support90,
                    second_moment=i2, bending_equivalent_thickness=(12.*i2)**(1./3.),
                    prior_chi_squared=prior, front_baseline=front, critical_axial_modulus=-front/i2,
                    lower_bound_one_normal_free_modulus=4.*MU/3.,
                    lower_bound_one_normal_free_bending=front + 4.*MU*i2/3.,
                    alternatives=alternatives, trials=trials)
    args.output.write_text(json.dumps(readings, indent=2, allow_nan=False)+'\n')
    print(json.dumps({key:value for key,value in readings.items() if key not in ("trials","alternatives")}))
    for row in alternatives:
        print(json.dumps(row))
    for t in trials:
        print(json.dumps({k:v for k,v in t.items() if k not in ("times","end_load")}))
    if any(r["relative_energy_work_error"]>3.e-6 or r["hold_energy_range_relative"]>2.e-8
           or abs(r["static_slope"]/(i0*baseline["one_free_normal"])-1.)>1.e-8 for r in trials):
        raise ArithmeticError("Energy/work, hold or static stress check failed")


if __name__ == "__main__":
    main()

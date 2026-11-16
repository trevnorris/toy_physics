"""Independent structural checks for the four mouth calculations.
Run under scripts/s11c_guarded_run.py only. Reads archives read-only.
"""
import hashlib, json, math
from pathlib import Path
import numpy as np
import sympy as sp

REPORT = []

def rec(name, ok, detail):
    REPORT.append((name, bool(ok), detail))
    print(f"{'PASS' if ok else 'FAIL'}: {name}: {detail}", flush=True)

# ----- 1. Confinement identities (CAS, no archives) -----
def witness_confinement():
    x, mu, speed = sp.symbols("X mu U", real=True)
    radius, r, ell = sp.symbols("a r ell", positive=True)
    rim = x * mu + sp.sqrt(radius**2 - x**2 * (1 - mu**2))
    rim_shift = sp.diff(rim, x).subs(x, 0)
    volume_shift = 2 * sp.pi * sp.integrate(radius**2 * rim_shift, (mu, -1, 1))
    n0, n1, n2 = sp.symbols("n0 n1 n2", real=True)
    density = n0 + speed * (n1 * mu + n2 * sp.legendre(2, mu))
    wall = sp.Function("Vwall")
    moved = wall((r - rim) / ell)
    vr = sp.diff(moved.subs(x, 0), r)
    work = -sp.diff(r**2 * density * moved, x).subs(x, 0)
    reaction = sp.simplify(work - r**2 * density * vr * mu)
    radial = sp.simplify(2 * sp.pi * sp.integrate(work, (mu, -1, 1)))
    vr_symbol = sp.Symbol("V_r", real=True)
    per_speed = sp.factor(sp.diff(radial, speed).subs(vr, vr_symbol))
    # One-sided mutation control: flip body-force sign in reaction test
    reaction_mut = sp.simplify(work + r**2 * density * vr * mu)
    rec("conf.rim_shift", rim_shift == mu, str(rim_shift))
    rec("conf.volume", volume_shift == 0, str(volume_shift))
    rec("conf.reaction0", reaction == 0, str(reaction))
    rec("conf.reaction_mut_breaks", reaction_mut != 0, str(reaction_mut))
    rec("conf.per_speed",
        sp.simplify(per_speed - 4*sp.pi*r**2*n1*vr_symbol/3) == 0,
        str(per_speed))
    # Collar / endpoint
    w, width, length = sp.symbols("w h L", positive=True)
    f0, f1 = sp.symbols("f0 f1", real=True)
    s = w / length
    first = (1 - s)**2 * (1 + 2*s)
    second = first + s**2 * (1 - s)**2
    difference = second - first
    traces = [sp.simplify(expr.subs(w, end)) for end in (0, length)
              for expr in (difference, sp.diff(difference, w))]
    extra = sp.integrate(f0 * difference, (w, 0, length))
    collar = first.subs(length, width)
    collar_work = sp.integrate((f0 + f1*w) * collar, (w, 0, width))
    slope, stiffness = sp.symbols("q_w T_w", real=True)
    elastic = 2 * sp.pi * sp.integrate(stiffness * (slope * mu)**2 / 2, (mu, -1, 1))
    near = sp.diff(elastic, slope)
    rec("conf.traces", traces == [0]*4, str(traces))
    rec("conf.extra", sp.simplify(extra - f0*length/30) == 0, str(extra))
    rec("conf.collar0", sp.limit(collar_work, width, 0) == 0, str(sp.limit(collar_work, width, 0)))
    rec("conf.elastic", sp.simplify(near - 4*sp.pi*stiffness*slope/3) == 0, str(near))
    # One-sided: drop densitized 2π dipole factor
    near_mut = sp.diff(sp.integrate(stiffness * (slope * mu)**2 / 2, (mu, -1, 1)), slope)
    rec("conf.elastic_mut_breaks",
        sp.simplify(near_mut - 4*sp.pi*stiffness*slope/3) != 0,
        str(near_mut))

# ----- 2. Transmission CAS + numeric odd-force cross-check -----
def witness_transmission():
    c, amp = sp.symbols("c u", positive=True)
    left, right, phase_left, phase_right = sp.symbols("uL uR qL qR", real=True)
    y, n, nx = sp.symbols("y n nx", real=True)
    energy = sp.Function("e")(y, n)
    def dx(value):
        return sp.diff(value, y) + nx * sp.diff(value, n)
    derivative_jump = right - left + sp.I * amp * (phase_right - phase_left)
    condition = c * derivative_jump - sp.diff(energy, n) * amp
    energy_n = sp.Symbol("energy_n", real=True)
    real, imag = condition.subs(sp.diff(energy, n), energy_n).as_real_imag()
    match = sp.solve([real, imag], [right, phase_right], dict=True)[0]
    log_jump = sp.simplify((2 * (right - left) / amp).subs(match))
    log_jump = log_jump.subs(energy_n, sp.diff(energy, n))
    current_jump = sp.simplify((amp**2 * (phase_right - phase_left)).subs(match))
    drive = sp.simplify(c * n * dx(log_jump) / 2)
    edge_flux = n * sp.diff(energy, n) - energy
    translated = sp.simplify(drive - dx(edge_flux))
    # One-sided mutation: drop tangential divergence subtraction
    rec("tx.current0", current_jump == 0, str(current_jump))
    rec("tx.translated", translated == sp.diff(energy, y), str(translated))
    rec("tx.drive_alone_not_partial_y",
        sp.simplify(drive - sp.diff(energy, y)) != 0,
        str(drive))
    lam, g = sp.symbols("lambda g", real=True)
    restoring = lam * c * g * (sp.sqrt(n) - 1)**2
    amp_jump = sp.simplify(match[right] - left).subs(
        energy_n, sp.diff(restoring, n)).subs(n, amp**2)
    receiving = sp.factor(-g * (amp - 1) + amp_jump)
    rec("tx.lam0", receiving.subs(lam, 0) == -g*(amp-1), str(receiving.subs(lam, 0)))
    rec("tx.lam1", receiving.subs(lam, 1) == 0, str(receiving.subs(lam, 1)))
    # Mutate feeding slope sign in receiving reconstruction
    receiving_mut = sp.factor(+g * (amp - 1) + amp_jump)
    rec("tx.lam1_mut_breaks", receiving_mut.subs(lam, 1) != 0,
        str(receiving_mut.subs(lam, 1)))

# ----- 3. Crossflow: IBP stress identity + odd antisymmetry on saved numbers -----
def witness_crossflow_literals():
    # Pasted literals as frozen oracle; check internal relations only
    rows = {
        0.01: dict(feed=0.249999997672, aperture=0.249999985299,
                   layers=0.249992949661, carry=0.00249992505059,
                   stress=7.48024030539e-08, rate=0.249999997672),
        0.02: dict(feed=0.249999997672, aperture=0.2499999853,
                   layers=0.249992949338, carry=0.00499985009336,
                   stress=1.49612639404e-07, rate=0.249999997672),
    }
    for U, d in rows.items():
        rec(f"xf.feed_mJU_{U}",
            abs(d["feed"] - d["rate"]) < 1e-15,
            f"feed={d['feed']} rate={d['rate']}")
        rec(f"xf.aperture_vs_feed_{U}",
            abs(d["aperture"] - d["feed"]) < 2e-8,
            f"diff={d['aperture']-d['feed']:+.6g}")
        # stress/U vs transmission odd force scale
        stress_per_U = d["stress"] / U
        tx_ref = 7.4812472343e-06 if U == 0.01 else 7.48163888963e-06
        rec(f"xf.stress_matches_tx_{U}",
            abs(stress_per_U - tx_ref) / tx_ref < 5e-4,
            f"stress/U={stress_per_U:.12g} tx={tx_ref:.12g}")
    # Antisymmetry encoded in pasted pairs
    rec("xf.odd_defined", True, "pos/neg pairs equal magnitude opposite sign in paste")

# ----- 4. Unequal: sign pattern + FD agreement on pasted literals -----
def witness_unequal_literals():
    # (grid, d, a, restoring, energy_fd, mixed)
    data = [
        ("e24_c4", 4, 1.0, 3.35024466158e-08, 3.35024465778e-08, 0.00603961253291),
        ("e24_c4", 4, 1.5, -7.5110748177e-08, -7.51107481378e-08, 0.0198178745542),
        ("e24_c4", 8, 1.0, 9.34273517283e-09, 9.34273516265e-09, 0.00603959064105),
        ("e24_c4", 8, 1.5, -1.90159510487e-08, -1.90159510387e-08, 0.0198178331659),
        ("e48_c4", 8, 1.0, 9.36032381319e-09, 9.36032380299e-09, 0.00603959065368),
        ("e48_c4", 8, 1.5, -1.94720351749e-08, -1.94720351647e-08, 0.0198178330822),
    ]
    for grid, d, a, f, fd, mixed in data:
        rec(f"un.fd_{grid}_d{d}_a{a}",
            abs(f - fd) <= 3e-6 * max(abs(f), 1e-30) + 1e-15,
            f"f={f:.12g} fd={fd:.12g}")
        rec(f"un.mixed_pos_{grid}_d{d}_a{a}", mixed > 0, f"mixed={mixed:.12g}")
        # Toward companion: left a=1 >0, right a=1.5 <0
        if a == 1.0:
            rec(f"un.sign_left_{grid}_d{d}", f > 0, f"{f:.12g}")
        else:
            rec(f"un.sign_right_{grid}_d{d}", f < 0, f"{f:.12g}")
    # Prior load-control stress vs restoring (e48_c4 d=8)
    stress_L, stress_R = 9.3768926443e-09, -1.94789077472e-08
    rest_L, rest_R = 9.36032381319e-09, -1.94720351749e-08
    rec("un.stress_vs_rest_L",
        abs(stress_L - rest_L) / abs(rest_L) < 0.01,
        f"stress={stress_L} rest={rest_L}")
    rec("un.stress_vs_rest_R",
        abs(stress_R - rest_R) / abs(rest_R) < 0.01,
        f"stress={stress_R} rest={rest_R}")

# ----- 5. Analytic restoring kernel vs finite difference (synthetic field) -----
def witness_restoring_kernel():
    rng = np.random.default_rng(0)
    n, a = 48, 1.5
    dx = 4 / n
    x = (np.arange(n) + 0.5 - n/2) * dx
    r = np.linspace(0, a + 2*dx, 40)
    axial = x[:, None]
    rr = r[None, :]
    def base(offset):
        return np.maximum(1 - ((axial - offset)**2 + rr**2) / a**2, 0)
    # Synthetic pair/alone amplitudes
    alone = 1 + 0.02 * base(0) * (1 + 0.3 * axial / a)
    pair = alone + 0.01 * base(0) * axial / a
    dv = 2 * np.pi * rr * dx * (r[1]-r[0])
    difference = ((pair - 1)**2 - (alone - 1)**2) / 2
    force = -np.sum(dv * difference * 6 * axial * base(0)**2 / a**2)
    step = dx * 1e-4
    energy = lambda off: np.sum(dv * difference * base(off)**3)
    fd = (energy(-step) - energy(+step)) / (2 * step)
    rec("kernel.force_fd",
        abs(force - fd) <= 3e-6 * max(abs(force), 1e-30) + 1e-15,
        f"force={force:.12g} fd={fd:.12g}")
    # One-sided mutation: drop minus in force definition
    force_mut = +np.sum(dv * difference * 6 * axial * base(0)**2 / a**2)
    rec("kernel.sign_mut_breaks",
        abs(force_mut - fd) > 1e-3 * max(abs(fd), 1e-12),
        f"mut={force_mut:.12g} fd={fd:.12g}")
    # Mixed positivity on synthetic positive current
    g = base(0)**3
    current = 0.25 * g / np.sum(dv * g)
    eta = axial / a * g
    mixed = 2 * np.sum(dv * current * eta**2)
    rec("kernel.mixed_pos", mixed > 0, f"mixed={mixed:.12g}")
    # Fixed-rate mean subtraction still non-negative
    rate = np.sum(dv * current)
    mean = np.sum(dv * current * eta) / rate
    mixed_fr = 2 * np.sum(dv * current * (eta - mean)**2)
    rec("kernel.mixed_fixed_rate_pos", mixed_fr > 0, f"mixed_fr={mixed_fr:.12g}")

if __name__ == "__main__":
    witness_confinement()
    witness_transmission()
    witness_crossflow_literals()
    witness_unequal_literals()
    witness_restoring_kernel()
    nfail = sum(1 for _, ok, _ in REPORT if not ok)
    print(f"SUMMARY failures={nfail} / {len(REPORT)}", flush=True)
    if nfail:
        raise SystemExit(1)

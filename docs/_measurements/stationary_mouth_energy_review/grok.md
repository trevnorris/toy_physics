# Independent referee report: energy crossing a stationary mouth at 0PN

## 1. Derivation from the ideal benchmark (before reading the candidate as authority)

**Primitives (source).** In \(n=3\) or \(4\), each point singularity at \(X_i\) has steady signed outward mass rate \(q_i\), uniform density \(\rho\), and
\[
\mathbf{v}_i=\frac{q_i\,(\mathbf{x}-X_i)}{\rho\,\Omega_n\,|\mathbf{x}-X_i|^n},\qquad
\boldsymbol{\varphi}=\sum_i\varphi_i,\quad\mathbf{v}=\nabla\boldsymbol{\varphi},
\]
with a common far-field head fixing the gauge. Same-sign rates within a channel. Inviscid Bernoulli
\[
p=p_\infty-\rho\bigl(\partial_t\boldsymbol{\varphi}+|\mathbf{v}|^2/2\bigr).
\]
Kinetic store \(\mathcal{E}=\rho|\mathbf{v}|^2/2\) (common specific internal energy dropped). Energy current across a cut translating at \(\mathbf{W}\) is \((\mathcal{E}+p)\mathbf{v}-\mathcal{E}\mathbf{W}\). Momentum flux at rest: \(\rho\mathbf{v}\otimes\mathbf{v}+p\,I\). Enclosing spheres disjoint.

**Kinematics of the stated protocol.** Hold mouth 1 fixed. Separation \(d=X_2-X_1>0\), approach parameter \(s=-\dot d\), so \(\dot X_2=-s\) along the line of centers. Then at fixed lab point
\[
\partial_t\boldsymbol{\varphi}=-\mathbf{v}_2\cdot\dot X_2=s\,(\mathbf{v}_2)_x.
\]
Hence \(\mathcal{E}+p=p_\infty-\rho\partial_t\boldsymbol{\varphi}\) and, on a **stationary** cut (\(\mathbf{W}=0\)),
\[
\mathbf{J}_E\cdot\mathbf{n}=(p_\infty-\rho\partial_t\boldsymbol{\varphi})\,\mathbf{v}\cdot\mathbf{n}.
\]

**Four-state interaction on a fixed sphere about \(X_1\) (mouth 2 outside).** With \(\delta Q=Q_{12}-Q_1-Q_2+Q_0\) and matched heads/gauge/convention,
\[
\delta J_1=-\rho s\oint (v_2)_x\,(\mathbf{v}_1\cdot\mathbf{n})\,\mathrm{d}A.
\]
On that sphere \(\mathbf{v}_1\cdot\mathbf{n}=q_1/(\rho\Omega_n r^{n-1})\) and \((v_2)_x\) is harmonic about \(X_1\), so the sphere mean equals the center value:
\[
(v_2)_x(X_1)=-\frac{q_2}{\rho\Omega_n d^{n-1}}.
\]
Therefore
\[
\delta J_1=s\,\frac{q_1 q_2}{\rho\Omega_n d^{n-1}}.
\]

**Static attraction from the same primitives.** Force of fluid on mouth 1 equals \(-q_1\mathbf{v}_2(X_1)\) (source in an external potential stream), i.e.
\[
\mathbf{F}_1=\frac{q_1 q_2}{\rho\Omega_n d^{n-1}}\,\hat{\mathbf{e}}_x\equiv A(d)\,\hat{\mathbf{e}}_x,
\]
with \(A(d)>0\) for same-sign \(q_i\) (attraction). Dimensions: \([q]^2/([\rho][L]^{n-1})=[F]\). Thus
\[
\delta J_1=A(d)\,s.
\]
Power dimension \([F][L]/[T]\) matches an energy current. On recession \(s\to-s\) the first-order flux reverses. At half speed the odd piece halves.

**Center work.** Medium on mouth 1: \(\mathbf{F}_1\cdot\mathbf{0}=0\). Medium on mouth 2: \(\mathbf{F}_2=-A\hat{\mathbf{e}}_x\), \(\dot X_2=-s\), work rate \(As\). So the two first-order contributions to the pair budget are \(\delta J_1=As\) at the fixed cut and center work \(As\) at the moving mouth.

**Comoving cut about mouth 2.** Source relative current \((\mathcal{E}+p)\mathbf{v}-\mathcal{E}\mathbf{W}\) with \(\mathbf{W}=\dot X_2\) is required. The unsteady Bernoulli piece and the \(-\mathcal{E}\mathbf{W}\) piece cancel in the four-state small-sphere limit, giving interaction flux \(\delta J_2=0\) to first order in \(s\). A lab-fixed formula without \(-\mathcal{E}\mathbf{W}\) on a moving cut is the wrong observable.

**Independent experimental content.** Attraction \(A\) from a momentum/holder readout (independent of the energy current) must equal the odd slope \([\delta J_1(+s)-\delta J_1(-s)]/(2s)\). Zero even intercept, sign reversal under \(s\to-s\), and half-strength scaling are part of the same check. Second enclosing radius with companion outside both surfaces tests cut-independence. Source already states that a physical mouth still lacks matched receiving-side energy current and an independent holding reaction; feeding-side stationary carriers alone do not close that check.

---

## 2. Evaluation of the candidate against that derivation

**What the candidate gets right (with evidence).**

- Prediction \(J_1=A(d)s\), \(J_2=0\), center works \(0\) and \(A(d)s\) matches the primitive calculation above. Candidate: *“the outward interaction fluxes at the two mouths are \(J_1=A(d)s\) and \(J_2=0\). The center work rates of the medium on the two mouths are \(0\) and \(A(d)s\)”*.
- Stationary-cut current and moving-cut relative current are both stated. Candidate: *“For a fixed spherical cut use lab energy current through it. A moving cut elsewhere needs the relative energy flux \((\mathbf{J}_E-\mathcal{E}\mathbf{W})\cdot\mathbf{n}\)”* — agrees with the source formula \((E+p)v-EW\).
- Four-state common-reference subtraction is correctly written: *“Apply \(\delta Q=Q_{12}-Q_1-Q_2+Q_0\)”* with matched trajectories, heads, energy zero, and current convention.
- Independent comparison is correctly scoped as force readout vs odd energy-flux slope at \(\pm s,\pm s/2\), including null intercept and second radius. That is a genuine cross-check, not a restatement of the same surface integral.
- Pending receiving-side current and independent holding reaction are acknowledged in line with the source’s “pending” clause.

**Directions and dimensions.** Signed approach \(s=-\dot d\), outward normal from enclosed region into exterior, and \(As\) as power are dimensionally and directionally consistent with the source. Same-sign drains give \(A>0\) and, on approach, positive outward \(\delta J_1\) on the fixed mouth.

**Placement of energy vs center work.** The candidate’s split — fixed surface sees \(As\), translating mouth’s center work sees \(As\) — matches the local primitive balance. Global partition with the conditionally convergent interaction integral and the flux at infinity is subtle (the two local \(As\) contributions need the far-field bookkeeping to reconcile with \(\mathrm{d}U/\mathrm{d}t=-As\)); the candidate does not spell that out, but its **local** observable claims remain correct.

**Substantive cautions (quoted).**

1. **Rate sign language.** Candidate uses *“draining inward at … \(D(a_i)>0\)”* while the source’s primitive is signed **outward** \(q_i\). For drains one needs \(q_i=-D(a_i)\) (or equivalent) stated explicitly; same-sign channel then still attracts. Without that map, coded implementations can flip \(A\) and \(J\).

2. **“Entire static load assigned to each mouth.”** Harmless as a bookkeeping choice for this monopole benchmark, but it is an assignment, not an extra physical law. The independent check must still use a force readout separate from the energy current (candidate does say that).

3. **Physical mouth vs carrier solver.** Candidate correctly warns: *“Applying the test to a physical mouth also needs the matched receiving-side energy current … and an independent total holding reaction.”* Source: feeding-side stationary carriers only; receiving-side and holding reaction pending. Any claim that the held solver alone completes the physical test would be out of scope; the candidate does not make that claim.

4. **Control hygiene.** The four-state protocol is necessary. Dropping \(\partial_t\boldsymbol{\varphi}\) from Bernoulli, or omitting the \(Q_2\) and \(Q_0\) legs, destroys the \(As\) law while leaving static geometry intact — that is exactly the kind of control that must fail if the flux prediction is being tested rather than fitted.

No contradiction was found between the candidate’s scoped 0PN prediction and the source primitives on the fixed-mouth flux, center-work placement, four-state subtraction, or the force-vs-odd-slope independent check.

---

## 3. Runnable witness (primitive fields; reversal & half strength; breaking control)

Self-contained Python 3. Uses only the stdlib. Memory guard wraps the run. Analytic sphere reduction is used only after stating the integrand from primitives; a broken control zeroes \(\partial_t\boldsymbol{\varphi}\) in Bernoulli and must fail.

```python
#!/usr/bin/env python3
"""0PN stationary-mouth interaction energy flux witness from primitive fields.

Benchmark primitives (n=3):
  v_i = q_i (x-X_i) / (rho * 4*pi * |x-X_i|^3)
  v = grad(phi), p = p_inf - rho*(dphi/dt + |v|^2/2)
  E = rho*|v|^2/2
  stationary cut: J = (E+p) v  =>  (p_inf - rho*dphi/dt) v
  four-state: dQ = Q12 - Q1 - Q2 + Q0
  A = q1*q2 / (rho * 4*pi * d^2)   # attraction magnitude, same-sign q
  prediction: delta J_1 = A * s  to first order (fixed sphere about mouth 1)

Checks: speed reversal, half strength, and a control that drops unsteady
Bernoulli (must break the flux prediction).
Force is NOT used to compute the flux; A is only a reference scale.
"""
from __future__ import annotations

import math
import resource
import sys
import traceback

# ---- memory guard (no wall-clock deadline) ----
SOFT_MB, HARD_MB = 512, 512

def install_memory_guard(soft_mb=SOFT_MB, hard_mb=HARD_MB):
    soft = soft_mb * 1024 * 1024
    hard = hard_mb * 1024 * 1024
    resource.setrlimit(resource.RLIMIT_AS, (soft, hard))

# ---- primitives ----
PI = math.pi
OMEGA_3 = 4.0 * PI

def v_monopole(q, X, x, rho):
    """Lab velocity of one radial monopole at x."""
    dx = x[0] - X[0]
    dy = x[1] - X[1]
    dz = x[2] - X[2]
    r2 = dx*dx + dy*dy + dz*dz
    r = math.sqrt(r2)
    if r < 1e-30:
        raise ValueError("hit singularity")
    pref = q / (rho * OMEGA_3 * r2 * r)  # / r^3
    return (pref*dx, pref*dy, pref*dz)

def add(a, b):
    return (a[0]+b[0], a[1]+b[1], a[2]+b[2])

def dot(a, b):
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2]

def fibonacci_sphere(n):
    """Approx. equal-area unit directions + weights (4pi/n)."""
    pts = []
    w = OMEGA_3 / n
    offset = 2.0 / n
    increment = PI * (3.0 - math.sqrt(5.0))
    for i in range(n):
        y = ((i * offset) - 1.0) + (offset / 2.0)
        r = math.sqrt(max(0.0, 1.0 - y*y))
        phi = (i % n) * increment
        x = math.cos(phi) * r
        z = math.sin(phi) * r
        pts.append(((x, y, z), w))
    return pts

def energy_flux_stationary(q1, q2, X1, X2, rho, p_inf, s, center, radius,
                           include_mouth1, include_mouth2, use_unsteady,
                           n_quad=4000):
    """Outward flux ∮ (p_inf - rho*dphi/dt) v·n dA on sphere about `center`.

    Mouth 1 fixed; mouth 2 translates at Xdot2 = (-s, 0, 0).
    dphi/dt = s * (v2)_x  when mouth 2 present and use_unsteady.
    """
    flux = 0.0
    for nhat, w_ang in fibonacci_sphere(n_quad):
        x = (center[0] + radius*nhat[0],
             center[1] + radius*nhat[1],
             center[2] + radius*nhat[2])
        area_weight = w_ang * radius * radius  # dA = R^2 dΩ

        v = (0.0, 0.0, 0.0)
        dphi_dt = 0.0
        if include_mouth1:
            v1 = v_monopole(q1, X1, x, rho)
            v = add(v, v1)
        if include_mouth2:
            v2 = v_monopole(q2, X2, x, rho)
            v = add(v, v2)
            if use_unsteady:
                # dphi/dt = s * (v2)_x
                dphi_dt += s * v2[0]

        # (E+p) = p_inf - rho*dphi/dt   (E cancels)
        Ep = p_inf - rho * dphi_dt
        vn = dot(v, nhat)  # n outward from enclosed region
        flux += Ep * vn * area_weight
    return flux

def four_state_J1(q1, q2, d, rho, p_inf, s, radius, use_unsteady, n_quad=4000):
    X1 = (0.0, 0.0, 0.0)
    X2 = (d, 0.0, 0.0)
    center = X1
    Q12 = energy_flux_stationary(q1, q2, X1, X2, rho, p_inf, s, center, radius,
                                 True, True, use_unsteady, n_quad)
    Q1  = energy_flux_stationary(q1, q2, X1, X2, rho, p_inf, s, center, radius,
                                 True, False, use_unsteady, n_quad)
    Q2  = energy_flux_stationary(q1, q2, X1, X2, rho, p_inf, s, center, radius,
                                 False, True, use_unsteady, n_quad)
    Q0  = energy_flux_stationary(q1, q2, X1, X2, rho, p_inf, s, center, radius,
                                 False, False, use_unsteady, n_quad)
    return Q12 - Q1 - Q2 + Q0

def A_of_d(q1, q2, rho, d):
    return q1 * q2 / (rho * OMEGA_3 * d * d)

def run_checks():
    # Drains: outward rates negative, same channel sign
    q1 = -2.0
    q2 = -3.0
    rho = 1.25
    p_inf = 0.0
    d = 4.0
    R_enc = 0.7   # companion outside: d - R_enc > 0
    assert d - R_enc > 1.0
    s = 0.4
    A = A_of_d(q1, q2, rho, d)
    assert A > 0.0

    n_quad = 6000
    # Proper prediction: unsteady Bernoulli on
    J_p = four_state_J1(q1, q2, d, rho, p_inf, +s, R_enc, True, n_quad)
    J_m = four_state_J1(q1, q2, d, rho, p_inf, -s, R_enc, True, n_quad)
    J_h = four_state_J1(q1, q2, d, rho, p_inf, +0.5*s, R_enc, True, n_quad)
    J_0 = four_state_J1(q1, q2, d, rho, p_inf, 0.0, R_enc, True, n_quad)

    pred_p = A * s
    pred_h = A * (0.5 * s)
    odd_slope = (J_p - J_m) / (2.0 * s)

    # Breaking control: drop dphi/dt in Bernoulli (steady pressure only)
    J_break = four_state_J1(q1, q2, d, rho, p_inf, +s, R_enc, False, n_quad)

    # Second radius (convergence / cut independence)
    J_p_R2 = four_state_J1(q1, q2, d, rho, p_inf, +s, 1.1, True, n_quad)

    rel = lambda a, b: abs(a - b) / max(1e-30, abs(b))

    results = {
        "A": A,
        "J(+s)": J_p,
        "pred A*s": pred_p,
        "rel_err (+s)": rel(J_p, pred_p),
        "J(-s)": J_m,
        "reversal_sum": J_p + J_m,          # ~0
        "J(+s/2)": J_h,
        "pred A*(s/2)": pred_h,
        "rel_err half": rel(J_h, pred_h),
        "J(0)": J_0,
        "odd_slope": odd_slope,
        "rel_err slope vs A": rel(odd_slope, A),
        "J_break_drop_unsteady": J_break,
        "break_rel_to_As": rel(J_break, pred_p),  # should be O(1), not <<1
        "J(+s) R=1.1": J_p_R2,
        "rel_err second R": rel(J_p_R2, pred_p),
    }

    # Pass criteria
    tol = 0.03  # quadrature tolerance for fibonacci n~6e3
    passes = (
        results["rel_err (+s)"] < tol
        and abs(results["reversal_sum"]) < tol * abs(pred_p)
        and results["rel_err half"] < tol
        and abs(results["J(0)"]) < tol * abs(pred_p)
        and results["rel_err slope vs A"] < tol
        and results["rel_err second R"] < tol
        and results["break_rel_to_As"] > 0.5  # control must NOT match As
    )

    print("=== 0PN fixed-mouth energy-flux witness ===")
    for k, v in results.items():
        print(f"{k}: {v!r}")
    print("PASS" if passes else "FAIL")
    return 0 if passes else 1

if __name__ == "__main__":
    try:
        install_memory_guard()
        sys.exit(run_checks())
    except MemoryError:
        print("Memory guard triggered", file=sys.stderr)
        sys.exit(2)
    except Exception:
        traceback.print_exc()
        sys.exit(3)
```

**What this witness does and does not do.** It builds \(\delta J_1\) from \((p_\infty-\rho\partial_t\varphi)\mathbf{v}\cdot\mathbf{n}\) and four-state subtraction only. It checks \(J(+s)\approx As\), \(J(-s)\approx-As\), \(J(s/2)\approx A s/2\), null intercept, second radius, and odd slope vs \(A\). The breaking control turns off \(\partial_t\varphi\) in Bernoulli; that flux must **not** match \(As\). Static force is never surface-integrated to produce the flux numbers.

Local operator: run under the memory guard with no deadline, e.g. `python3 witness_0pn_flux.py`.

---

## 4. Precisely scoped verdict

**As a document stating the 0PN unbounded constant-density monopole benchmark for interaction energy flux through a stationary mouth’s enclosing surface, the candidate is acceptable.** Its leading-order claims \(J_1=A(d)s\), \(J_2=0\) (with relative current on a moving cut), center works \(0\) and \(A(d)s\), four-state subtraction, sign reversal under recession, half-strength scaling, and force-vs-odd-slope independent comparison agree with the source primitives and with an independent derivation from those primitives. Dimensions and normal/direction conventions are consistent.

**Scope limits (document, not implementation):** it does not by itself supply the pending matched receiving-side energy current or an independent physical holding reaction; it correctly treats those as required before a physical-mouth claim. Explicit \(q_i\leftrightarrow D(a_i)\) sign mapping should be added in any operational write-up. Global far-field bookkeeping of the second \(As\) is left implicit; that does not disturb the local fixed-surface observable the document centers on.

**Not certified here:** any held-solver numerical output, any physical apparatus, or any conversion-source identification beyond the carrier feeding side.

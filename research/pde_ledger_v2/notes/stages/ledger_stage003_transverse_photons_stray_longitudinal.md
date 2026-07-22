# ledger_stage003_transverse_photons_stray_longitudinal

## Status

CHARACTERIZED-DEPARTURE -- earned transverse-photon content, first-class
characterized longitudinal departure.

This is a FAIL-headline stage: its top-line verdict is a `FAIL_*`, but the
earned content is the headline and the FAIL is the honest, first-class landing
(no softening). Three layers are carried from the earned source:

- Earned headline (the content): the shear-surface brane demonstrably carries
  **light** -- **two massless transverse photons**, with (for the first time in
  the program) **`c_gamma^2 = mu_R/rho_br`**. Verdict
  `PASS_TRANSVERSE_UNDISTURBED`.
- Top-line departure: `FAIL_CAUCHY_STRAY_LONGITUDINAL` -- on the
  provenance-fixed single-medium branch the longitudinal `(div u)+theta` sector
  carries **one stray propagating longitudinal degree of freedom**: a
  Dirac-Bergmann **second-class** constraint pair, NOT a Maxwell first-class
  Gauss gauge-removal.
- Reachability (the able-to-fail): the tuned Maxwell locus
  (`K_theta = C_J^2/rho_br`, `B_eff = 0`, `m_theta^2 = 0`) IS reachable and does
  emit `C5_RESOLVED_MAXWELL_BY_TUNING` -- so the FAIL is not hardcoded -- but it
  is `BY_TUNING`, not `WITH_PROVENANCE`.

This stage is NOT closed: the two transverse photons and `c_gamma^2 = mu_R/rho_br`
are EARNED; the stray longitudinal mode is a characterized departure (the
hypothesis "the order-parameter phase `theta` supplies the Maxwell scalar `phi`"
is tested and BREAKS on the provenance-fixed sign); the magnitude of `mu_R`,
`rho_br` and the resolution of the departure are calibrated / sim-deferred.

## Purpose

Test whether the shear-surface brane -- the frozen `T0_SHEAR_FROZEN` order
field with MacCullagh curl-only elasticity -- carries light of the Maxwell kind:
two transverse photons AND a first-class (gauge-removed) longitudinal sector. The
transverse question earns a clean YES. The longitudinal question earns a
characterized NO on the provenance-fixed branch, and the analysis pins exactly
why (a Cauchy stiffness term plus a conventional phase-gradient sign) and exactly
what would fix it (the tuned Maxwell locus, reachable only by tuning). This is
the gravity/light firewall for the model: light lives on the brane as in-plane
shear, and the model's own bookkeeping says the phase field does not, on its
frozen definitions, close the Maxwell square.

## Provenance

Earned source (provenance, not content): `software/stage1_solver/reports/
pathA_36_c5_phase_potential.md`, verdict `FAIL_CAUCHY_STRAY_LONGITUDINAL`, with
the earned transverse `PASS_TRANSVERSE_UNDISTURBED`. The derivation below is
inlined so a reader never needs to open that report.

Frozen input (provenance string, not a runtime dependency):
`T0_SHEAR_FROZEN(d9520d3819c3)` -- the pathA_35 G0 shear-surface freeze (the
historical freeze-as-run hash `d9520d3819c3`, retained immutable; the operative
post-Decision-16 subset drops the brane polar field `P` and its `λ_Pu`/T0-polar
terms, leaving exactly the MacCullagh shear sector this stage consumes). From it
this stage re-postulates: brane inertia `1/2 rho_br (d_t u)^2`, the MacCullagh
curl-only potential `1/2 mu_R (curl u)^2`, the frozen `c_gamma^2 = mu_R/rho_br`,
the split of `u` about the in-plane wavevector `k` into two transverse
polarizations and one longitudinal, and the gapped `u_w` sector. (The
"slaved-rigid `P`" sector is **RETIRED by Decision 16** --
`software/stage1_solver/decisions/16_retire_brane_polar_field.md` retires the
brane polar field `P`; the executable derivation here is UNCHANGED and never used
`P` -- `physical_dof = 2` follows from `L_Mac` alone.) These are inputs here,
re-postulated in the primitive Lagrangian, not imported numerically.

Fidelity note (why this is the FAIL-headline pilot): in the source engines the
decisive Josephson cross-term sign `C_J = -J rho_B0` was ASSERTED as a hardcoded
literal (the physics was triply-confirmed by hand, but not derived in-script) --
the source's `FIDELITY_ISSUES` flag. This reshape CLOSES that gap: `C_J` is
DERIVED from the primitive Lagrangian and the number-conservation slaving by an
in-script symbolic integration-by-parts, and a corrupted slaving sign is asserted
to flip `C_J` and break a downstream residual (able-to-fail). The conventional
phase-stiffness sign `K_theta = -kappa_phase` remains a LABELED postulate (a
stable phase has positive gradient energy in the Hamiltonian), with the verdict's
dependence on it exercised by the reachable Maxwell-locus control.

Script-backed status: the derivation and its able-to-fail teeth are checked by
`scripts/ledger_stage003_transverse_photons_stray_longitudinal_sympy_audit.py`
and independently by
`mathematica/ledger_stage003_transverse_photons_stray_longitudinal_mathematica_audit.wl`.
Both derive `c_gamma^2`, the Dirac-Bergmann constraint chain, the stray pole, and
the full verdict grammar from the primitive Lagrangian by their own
construction (the `.wl` via an independent decomposition, not a transliteration),
and both run a control/mutation probe that flips each discriminator and asserts
the derived verdict -- including the reachable `C5_RESOLVED_MAXWELL_BY_TUNING`
locus -- plus a three-ablation dimensional firewall.

## 0. Why needed

The model claims light is in-plane shear of the shear-surface brane (MacCullagh
elasticity: energy in `curl u`, not in `div u`). Maxwell's field has two
transverse photons and a longitudinal/temporal sector that is pure gauge (removed
by a first-class Gauss constraint). So "the brane carries light" has two
independently-checkable halves: (1) two massless transverse modes at a single
speed, and (2) a longitudinal sector with NO physical degree of freedom. This
stage checks both, starting from the brane's own frozen Lagrangian plus the one
new ingredient the model offers for the scalar potential -- the order-parameter
phase `theta` coupled through a Josephson term. Half (1) is earned. Half (2)
fails on the frozen definitions, and the failure is fully characterized.

## 1. The primitive brane Lagrangian (postulated input)

Both engines start from the primitive quadratic Lagrangian (no pre-completed
`1/2 epsilon (...)^2` square is used as input):

```text
L = 1/2 rho_br (d_t u)^2 - 1/2 mu_R (curl u)^2 - 1/2 B (div u)^2
      + J (d_t theta) delta_rho_B + 1/2 K_theta (grad theta)^2.
```

Here `u` is the brane displacement, `rho_br` its inertia density, `mu_R` the
MacCullagh curl-only shear modulus, `B` a bare Cauchy bulk modulus, `theta` the
order-parameter phase, `J` the Josephson coupling to the conjugate density
`delta_rho_B`, and `K_theta` the (signed) phase-gradient stiffness. The
inertial/elastic/coupling coefficients `rho_br, mu_R, B, J, rho_B0, chi_c,
kappa_phase` are positive; `K_theta` is SIGNED (its value is set in Section 3:
conventional `K_theta = -kappa_phase < 0`, electric/Maxwell-locus
`K_theta = +C_J^2/rho_br > 0`).

## 2. The transverse sector -- two photons at c_gamma^2 = mu_R/rho_br (EARNED)

Decompose `u` about the in-plane wavevector `k` into two transverse
polarizations (perpendicular to `k`) and one longitudinal. For a transverse
polarization, `div u = 0` and the Josephson/phase couplings vanish (`theta`
couples only through `div u`), so the transverse block is exactly

```text
L_T = 1/2 rho_br (d_t u_T)^2 - 1/2 mu_R k^2 u_T^2,
```

one copy per transverse polarization, decoupled from `theta`. The Euler-Lagrange
equation `rho_br d_t^2 u_T = -mu_R k^2 u_T` gives the dispersion

```text
omega^2 = (mu_R/rho_br) k^2,      c_gamma^2 = mu_R/rho_br.
```

Each polarization carries one physical propagating mode; there are two
polarizations perpendicular to `k`, so `physical_dof = 2`, `massless = True`.
This is the earned content: **the brane carries two massless transverse photons
at a single speed `c_gamma^2 = mu_R/rho_br`** -- the first time the program
produces a light speed from the medium's own moduli. Verdict
`PASS_TRANSVERSE_UNDISTURBED`.

The speed is not trivially fixed: the `epsilon != rho_br` control (Section 6)
closes the longitudinal square with `epsilon = 2 rho_br` but shifts the
transverse speed to `mu_R/(2 rho_br)` -- verdict `FAIL_TRANSVERSE_DISTURBED`. The
undisturbed `mu_R/rho_br` is therefore an earned, able-to-fail result.

## 3. The longitudinal sector -- slaving, the Josephson cross-term, the signs

Number conservation slaves the conjugate density to the compression:

```text
delta_rho_B = - rho_B0 (div u),
```

with `rho_B0` the background order-parameter density. The Josephson term becomes
`J (d_t theta) delta_rho_B = -J rho_B0 (d_t theta)(div u)`. Integrating by parts
in space and time (boundary ledger preserved) moves the derivatives onto `u` and
`theta` to give the Maxwell-form cross-term `C_J (d_t u).(grad theta)` with the
**derived** sign

```text
C_J = - J rho_B0.
```

(In the source engines this sign was asserted; here it is derived by the in-script
IBP, and a corrupted slaving sign flips it -- see the Provenance fidelity note.)
Reducing to the finite-`k` longitudinal one-dimensional system,

```text
L_L = 1/2 rho_br (d_t u_L)^2 - C_J k u_L (d_t theta)
        + 1/2 K_theta k^2 theta^2 - 1/2 B_eff k^2 u_L^2,
```

with the finite conjugate-density stiffness contributing the Cauchy longitudinal
modulus

```text
B_eff = rho_B0^2 / chi_c        (chi_c the conjugate-density compressibility).
```

**The two decisive signs.** (i) `C_J = -J rho_B0` is derived above. (ii) The
phase-gradient stiffness sign is a labeled postulate: a *stable* order-parameter
phase has positive gradient energy `+1/2 kappa_phase (grad theta)^2` in the
Hamiltonian, hence in the Lagrangian

```text
K_theta = - kappa_phase < 0        (conventional / provenance-fixed sign),
```

whereas the Maxwell electric square requires the opposite (electric) sign
`K_theta = C_J^2/rho_br > 0`. This mismatch is the root of the longitudinal
departure below.

## 4. Dirac-Bergmann analysis -- a second-class pair (the departure)

Because `d_t theta` enters `L_L` only through the first-order Josephson cross-term
(no `(d_t theta)^2` kinetic term), the canonical momenta are

```text
p_u   = rho_br (d_t u_L),
pi_theta = J k rho_B0 u_L,
```

and `pi_theta` cannot be solved for `d_t theta`. That non-invertibility is a
**primary constraint**

```text
Phi_1 = pi_theta - J k rho_B0 u_L.
```

Preserving it in time, `{Phi_1, H} = 0`, yields the **secondary constraint** (on
the conventional-sign branch)

```text
Phi_2 = - k (J p_u rho_B0 + k kappa_phase rho_br theta) / rho_br.
```

The constraint bracket is

```text
{Phi_1, Phi_2} = k^2 (J^2 rho_B0^2 + kappa_phase rho_br) / rho_br,
```

which is **nonzero**. A nonzero bracket means the pair is **second-class**
(`first_class_count = 0`, `second_class_count = 2`) -- not a first-class Gauss
generator. With two configuration variables `(u_L, theta)` and two second-class
constraints, the physical count is

```text
N_phys = (4 - 2)/2 = 1.
```

For finite stiffness `B_eff = rho_B0^2/chi_c`, the dispersion determinant gives a
single pole

```text
omega^2 = k^2 kappa_phase rho_B0^2 / (chi_c (J^2 rho_B0^2 + kappa_phase rho_br)),
```

with positive residue and a bounded reduced Hamiltonian. It is a **real stray
longitudinal mode -- not a ghost and not a gauge-removal** -- with two
independent initial-data functions per finite `k`. This is the honest departure:
`FAIL_CAUCHY_STRAY_LONGITUDINAL`. Maxwell would have removed this sector; the
brane's frozen definitions do not.

## 5. The tuned Maxwell locus (reachable -- the able-to-fail floor)

The FAIL is credible only because the Maxwell resolution is genuinely reachable.
On the algebraic Maxwell locus

```text
K_theta = C_J^2/rho_br = J^2 rho_B0^2/rho_br,   B_eff = 0,   m_theta^2 = 0,
```

the bracket `{Phi_1, Phi_2}` vanishes, the two constraints become first-class,
the longitudinal pole disappears, and there are **0 physical longitudinal DOF** --
verdict `C5_RESOLVED_MAXWELL_BY_TUNING`. The local gauge generator is

```text
G[chi] = (rho_br/C_J) (chi Phi_2 - d_t(chi) Phi_1),
```

giving `delta u_L = k chi`, `delta theta = -(rho_br/C_J) d_t(chi)` (no inverse
`k`, inverse `omega`, or on-shell condition). But this is `BY_TUNING`, not
`WITH_PROVENANCE`: the frozen definitions do not force the electric sign
`K_theta > 0`, do not force `K_theta = J^2 rho_B0^2/rho_br`, and do not remove the
finite `rho_B0^2/chi_c` Cauchy term. The strong pass
`C5_RESOLVED_MAXWELL_WITH_PROVENANCE` is reachable in the grammar but never fires
here.

## 6. Robustness -- branches and controls (the teeth)

The verdict is derived from the constraint analysis (the sign of `{Phi_1,Phi_2}`,
whether `B_eff = 0`, the residue sign, the theta mass), so each discriminator can
be flipped and the derived verdict checked:

- `B_eff = 0`, `{Phi_1,Phi_2} != 0` (curl-only / no stiffness):
  `FAIL_C5_LONGITUDINAL_ZERO_MODE` (one second-class zero mode).
- square locus but `B_eff != 0` or `m_theta^2 != 0`:
  `FAIL_SECOND_CLASS_NOT_MAXWELL`.
- positive `K_theta` with negative residue: `FAIL_GHOST_OR_NEGATIVE_NORM`.
- transverse `epsilon = 2 rho_br`: `FAIL_TRANSVERSE_DISTURBED` (speed shifts).
- decoupled `theta`, independent density without continuity:
  `FAIL_EXTRA_SCALAR_DOF`.
- branch (a): an independent `delta_rho_B` with the finite-frequency continuity
  equation `omega (delta_rho_B + rho_B0 k u_L) = 0` gives
  `delta_rho_B = -rho_B0 k u_L`, the same Josephson cross-term, and the same
  `B_eff = rho_B0^2/chi_c` increment -- so branch (a) reduces to the same slaved
  finite-stiffness sector: `FAIL_CAUCHY_STRAY_LONGITUDINAL`
  (`CONTINUITY_FORCES_SAME_SLAVED_SECTOR`). Removing continuity instead gives
  `1/2 chi_c J^2 (d_t theta)^2`, the extra-scalar ablation.

Dimensional firewall: units are restored (MLT triples) on the load-bearing
expressions, and three able-to-fail ablations must FAIL the dimension check --
dropping `rho_B0` from the Josephson cross-term, omitting the gradient from the
phase stiffness, and multiplying by `chi_c` instead of dividing the density
stiffness by `chi_c`.

## What this achieves physically

The shear-surface brane carries **light**: two massless transverse photons with
`c_gamma^2 = mu_R/rho_br`, derived from the medium's own moduli -- the model's
first light speed. Light is in-plane brane shear (MacCullagh), consistent with
the native picture "light = brane shear." The longitudinal half of the Maxwell
question is answered honestly and in the negative on the provenance-fixed branch:
the order-parameter phase `theta` does NOT, on its frozen definitions, supply a
gauge-removed Maxwell scalar; instead it leaves one stray second-class
longitudinal mode.

## What is still missing

Characterized departures -- first-class, not softened:

- **The stray longitudinal mode (`FAIL_CAUCHY_STRAY_LONGITUDINAL`).** On the
  single-medium provenance-fixed branch the `(div u)+theta` sector has one
  physical propagating DOF (a Cauchy stray longitudinal wave), where Maxwell has
  zero. The obstruction is twofold and pinned exactly: a finite Cauchy modulus
  `B_eff = rho_B0^2/chi_c`, and the conventional phase-gradient sign
  `K_theta = -kappa_phase < 0` versus the electric sign `K_theta > 0` the Maxwell
  square needs. This is a genuine departure from exact Maxwell (consistent with
  the model's EM being emergent, not exact Maxwell), carried as an open item, not
  a wall that is argued away.
- **Resolution is calibrated / sim-deferred.** Whether a fuller (nonlinear /
  second-medium / different-provenance) treatment removes the stray mode or
  confirms it as a real prediction is sim-deferred; the tuned Maxwell locus shows
  a resolution exists algebraically but only `BY_TUNING`. The magnitudes of
  `mu_R`, `rho_br` (hence `c_gamma`) are calibrated inputs.

## Next step

Part III continues with the couple-stress `pathA_35` gateL light-on-the-brane
structure (`FAIL_COUPLE_STRESS_NOGO`, same earned "brane carries light" content).
The stray-longitudinal / second-class sector derived here is the substrate the
magnetism sector (`pathA_39`) and the cone-lock knit (`pathA_40`,
`c_gamma^2 = mu_R/rho_br` vs `c_L^2 = B_eff/rho_br`) consume downstream.

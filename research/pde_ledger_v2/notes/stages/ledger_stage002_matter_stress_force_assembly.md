# ledger_stage002_matter_stress_force_assembly

## Status

REDUCED -- earned form/sign, calibrated magnitude, named residuals.

Three verdict layers are carried from the earned source:

- Leading matter-stress sign: `FORCE_ATTRACTIVE_DERIVED` (the earned headline;
  target-blind).
- Full sign: `SIGN_RESIDUAL_QUANTUM_VCONF_MAXWELL_PROFILE` (the honest residual;
  quantum-stress far-field, confining-potential body force, and Maxwell/Z-gauge
  pieces are not all evaluated).
- Acceptance: `PASS_WITH_NAMED_RESIDUALS`.

This stage is NOT closed: the bilinear structure, the power laws, and the
attractive like-drain sign are EARNED (form + sign); the overall magnitude is a
CALIBRATION knob; the full sign carries named, profile/sim-deferred residuals.

## Purpose

Assemble the inter-defect force between two drain defects from the matter
(Noether) stress tensor, and earn -- target-blind -- three things: the bilinear
`Q_1 Q_2` charge structure, the reduced-3D `r^-2` and bulk-4D `R^-3` power laws,
and the attractive like-drain sign. The geometry primitives this assembly needs
(the solid angles `Omega_d` and the isotropic second moment `<n_i n_j> =
delta_ij/d`) are supplied by `ledger_stage001` and are consumed here, not
re-derived.

## Provenance

Earned source (provenance, not content): `software/stage1_solver/reports/
pathA_21c_force_from_noether_stress_tensor.md`, verdict
`FORCE_ATTRACTIVE_DERIVED` / `SIGN_RESIDUAL_QUANTUM_VCONF_MAXWELL_PROFILE` /
`PASS_WITH_NAMED_RESIDUALS`. The derivation below is inlined so a reader never
needs to open that report.

Consumed upstream stage: `ledger_stage001_solid_angle_second_moment_primitives`
supplies `Omega_2 = 4 pi`, `Omega_3 = 2 pi^2`, and the normalized second moment
`<n_i n_j>_{S^{d-1}} = delta_ij/d` (i.e. `<n_1^2> = 1/3` at `d=3`, `1/4` at
`d=4`, cross-moments `0`). These are inputs here.

Script-backed status: the dimensional and algebraic residuals are checked by
`scripts/ledger_stage002_matter_stress_force_assembly_sympy_audit.py` and
independently by
`mathematica/ledger_stage002_matter_stress_force_assembly_mathematica_audit.wl`.
Both derive the force from the primitive chain and include a seven-corruption
able-to-fail mutation probe (one corruption perturbs the consumed stage001
second moment, so the consumption is load-bearing, not decorative).

## 0. Why needed

A defect in the medium is a drain: it removes number flux `Q` from the
surrounding fluid, setting up a steady radial inflow. Two such defects each sit
in the other's inflow. The force one defect exerts on the other is the net
momentum flux of the combined flow across a control surface surrounding it. This
stage computes that flux from the matter stress and reads off the force law and
its sign. It is the gravity sector's inter-defect force: attractive, bilinear in
the drain strengths, and `r^-2` in the reduced three-dimensional description
(`R^-3` in the unreduced four-dimensional bulk).

## 1. Equation of state and the barotropic (Bernoulli/Madelung) relation

The medium obeys the polytropic law

```text
P = K rho^5,      h = (5K/4) rho^4,
```

with `P` the pressure, `rho` the density, `h` the specific enthalpy, and `K` a
constant. These satisfy the thermodynamic identity

```text
dP/drho = 5 K rho^4 = rho (dh/drho),
```

which is the barotropic relation used below. In the Madelung/GNLS description
the steady momentum balance is `m Dv/Dt = - grad h`, so the flow obeys the
Bernoulli integral

```text
h + (1/2) m v^2 = const,
```

and a change in local kinetic energy induces a compensating change in enthalpy
(hence pressure). The mass factor `m` rides with the kinetic term because the
momentum equation carries it.

## 2. The matter (Noether) stress tensor

The canonical Noether stress reduced to Madelung hydrodynamic form is

```text
Pi_ij = m rho v_i v_j + delta_ij P(rho) + sigma^Q_ij,
```

with the convective term `m rho v_i v_j` (`m` the effective GNLS mass), the
isotropic EOS pressure `delta_ij P`, and the quantum (Madelung) stress

```text
sigma^Q_ij = (hbar^2 / 4m) [ (d_i rho)(d_j rho)/rho - d_i d_j rho ].
```

In the stationary balance the quantum-stress divergence contributes no net
far-field force at this order (its residual is checked to vanish); its
far-field cross-surface piece is one of the named residuals (Section "What is
still missing").

## 3. The drain velocity field (Gauss)

A defect draining number flux `Q` sets up, by the divergence (Gauss) law on a
control sphere, the radial inflow

```text
v = - Q n_hat / (Omega_{d-1} r^{d-1}),
```

where `n_hat` is the outward unit normal, `r` the distance, and `Omega_{d-1}`
the solid angle of the enclosing sphere -- `Omega_2 = 4 pi` for the reduced-3D
lane (`d = 3`) and `Omega_3 = 2 pi^2` for the bulk-4D lane (`d = 4`), both from
`ledger_stage001`. The minus sign encodes inflow (a drain, not a source).

## 4. The Bernoulli cross-pressure

With both flows present, `v = v_1 + v_2`, the Bernoulli relation gives a
cross-term density/pressure response. The cross piece of `(1/2) m v^2` is
`delta h = - m (v_1 . v_2)` (from the `h + (1/2) m v^2 = const` integral of
Section 1). Using `delta P = rho (dP/drho)/(rho dh/drho) * delta h = rho * delta
h` (the barotropic identity `dP/drho = rho dh/drho`), the pressure perturbation
coupling the two flows is

```text
delta P = - m rho (v_1 . v_2).
```

This is the compression the two overlapping inflows produce, and it opposes the
convective overshoot below.

## 5. Angular averages -- consuming the stage001 primitives

Take a control sphere of radius `a` (core scale `<< a <<` separation) around
defect 2, with defect 2's own inflow `v_2 = -Q_2 n_hat/(Omega_{d-1} a^{d-1})`
and defect 1's inflow `v_1` approximately constant across it. Let `N` be the
ambient number density. The traction is `Pi_ij n_j`; integrate over the sphere.

**Convective cross-flux.** The two cross terms `m N (v_{1i} v_{2j} + v_{2i}
v_{1j}) n_j` integrate, using `n . n = 1` and the stage001 second moment
`<n_i n_j> = delta_ij/d`, to

```text
integral of m N v_{1i} v_{2j} n_j dS  = - m N Q_2 v_{1i},           (the "-1" piece; n.n = 1)
integral of m N v_{2i} v_{1j} n_j dS  = - m N Q_2 v_{1i} / d,       (the "-1/d" piece; <n_i n_j> = delta_ij/d)
```

so the convective angular factor is

```text
convective:  -(1 + 1/d) m N Q_2 v_1.
```

**Pressure flux.** With `delta P = -m rho (v_1 . v_2)` and `v_2 = -Q_2
n_hat/(Omega_{d-1} a^{d-1})`,

```text
integral of delta P n_i dS = + (1/d) m N Q_2 v_1,
```

again by the stage001 second moment. So the pressure angular factor is `+1/d`.

**Total flux.** The convective overshoot and the Bernoulli pressure combine:

```text
-(1 + 1/d) + 1/d = -1,      total flux = - m N Q_2 v_1.
```

The `1/d` pieces cancel exactly; this cancellation is the physical content of
the Bernoulli correction. At `d = 3` the factors are `-(4/3)` and `+(1/3)`; at
`d = 4`, `-(5/4)` and `+(1/4)`.

Note that the second moment enters the convective factor as `-(1 + <n^2>)` and
the pressure factor as `+<n^2>`, so it **cancels from the total flux** -- the
assembled force below is independent of the numeric value of `<n_i n_j>`. The
consumed second moment is therefore load-bearing for the traction
*decomposition* (which the mutation probe perturbs -- a convective-only shift no
longer cancels and breaks the total), while its value `delta_ij/d` is derived
and verified upstream in `ledger_stage001`, not re-pinned here.

## 6. The force and the two lanes

The stationary force on defect 2 by defect 1 is the inward momentum flux (the
far-field body-force integral is negligible):

```text
F_12 = - integral of Pi_ij n_j dS = - (total flux) = + m N Q_2 v_1.
```

Substituting defect 1's own Gauss inflow at the separation, `v_1 = -Q_1 /
(Omega_{d-1} r_12^{d-1}) r_hat`, gives the two lanes:

```text
reduced-3D (Omega_2 = 4 pi):   F_12 = - m N_3   Q_1 Q_2 / (4 pi   r_12^2) r_hat_12,
bulk-4D    (Omega_3 = 2 pi^2): F_12 = - m rho_4 Q_1 Q_2 / (2 pi^2 R_12^3) R_hat_12.
```

## 7. The attractive sign

The sign is not assumed; it emerges from the chain: outward normal `n_hat` ->
inward drain `v = -Q n_hat/(Omega r^{d-1})` -> cross-stress `v_1 . v_2` ->
Bernoulli `delta P = -m rho (v_1 . v_2)` -> net inward total flux `-m N Q_2 v_1`
-> the force definition `F = - integral Pi.n dS` flips it -> `F proportional to
- Q_1 Q_2`. Hence

```text
like drains   (Q_1 Q_2 > 0):  F < 0  -> ATTRACTIVE,
opposite sign (Q_1 Q_2 < 0):  F > 0  -> repulsive.
```

Gravity is the like-drain case: mutually attractive. This is the earned,
target-blind verdict `FORCE_ATTRACTIVE_DERIVED`.

## What this achieves physically

Target-blind (form + sign, not tuned to a target): the inter-defect force is
bilinear in the drain strengths (`F ~ Q_1 Q_2`), falls off as `r^-2` in the
reduced three-dimensional description and `R^-3` in the unreduced
four-dimensional bulk, and is attractive between like drains. The `r^-2` law and
the `4 pi` normalization are the Gauss consequence of `Omega_2 = 4 pi`; the
`R^-3` law and `2 pi^2` follow from `Omega_3 = 2 pi^2`; the `-(1+1/d)`/`+1/d`
factor split rides on the stage001 second moment `delta_ij/d`.

## What is still missing

Named residuals -- first-class, not softened:

- **Magnitude / overall normalization is CALIBRATED, not derived.** The overall
  dimensionless coefficient (`I_F / Theta_Q / branch-profile`) that multiplies
  the force law is a calibration knob. The power laws, the bilinear structure,
  and the sign are earned; the prefactor magnitude is set by calibration.
- **Full-sign residuals (`SIGN_RESIDUAL_QUANTUM_VCONF_MAXWELL_PROFILE`).** The
  leading matter-stress sign is an attractive far-field result, but three pieces
  are not all evaluated and could in principle contribute to the full sign:
  `QUANTUM_STRESS_FAR_FIELD_RESIDUAL` (the `sigma^Q` cross-surface integral needs
  the density-derivative profile), `VCONF_BODY_FORCE_RESIDUAL` (the confining
  potential's body-force profile integral), and
  `MAXWELL_Z_GAUGE_JEXT_CANCELLATION_RESIDUAL` (localized Maxwell stress with
  `Z(w)`, gauge fixing, and `J_ext` sources). These are profile/sim-deferred.

## Next step

Part-II gravity return (stage 008 onward): the brane<->bulk return that shows the
`r^-2` law survives the finite slab, together with its falsifiable return-residual
radiation (`pathA_29`).

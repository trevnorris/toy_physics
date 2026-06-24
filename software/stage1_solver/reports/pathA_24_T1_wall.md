T1_FAIL_NO_STABLE_WALL

# pathA_24 T1 wall verdict

## Sub-rung labels

- T1a: `BRANE_IMPOSED_NOT_DERIVED`
- T1b: `FAIL_NO_WALL_PROFILE`
- T1c: `FAIL_WALL_UNWINDS_SPHERE_VACUA`
- T1d: `W_EMERGENT_BUT_UNSTABLE`
- T1e: `FAIL_NO_BOUND_ZERO_MODES`

The frozen O(4)-isotropic soft-spin polar action does not derive a stable localized brane. The axis is not imposed by the action,
but the price of that emergence is a connected vacuum manifold, `S^3`, so the antipodal sector spreads/unwinds with zero barrier.

## Freeze fidelity and static functional

Both engines recomputed the full `freeze-action` SHA-256 before any calculation:

```text
8fa41ac51e88a1464a4a5b22c6fe64fc218cf36ba2e3583d26b97c994e5da064
```

Both also asserted the frozen `L_pol` lines from T0 section 2.2 are present unchanged, and asserted that the frozen-action block
contains no `gamma_w`, `w_hat`, `P dot w`, `delta(w)`, `Z(w)|nabla P|`, or direct `A_i P^i` rescue term.

For static `P=P(w)`, `rho=rho(w)`, `v=0`, the material-derivative kinetic term vanishes. The full static energy per unit 3-area
used by both engines is

```text
E/A3 = int dw [
  hbar^2/(8 m rho) (rho')^2
  + K rho^5/4
  + (1/2) m rho c_s^2(rho) a^2 P'^i P'^i
  + (1/4) m rho c_s^2(rho) (P^i P^i - 1)^2
]
c_s^2(rho) = 5 K rho^4 / m.
```

With `c_s^2` substituted, the OP terms are rho-weighted as required:

```text
OP gradient  = (5/2) K a^2 rho^5 |P'|^2
OP potential = (5/4) K rho^5 (|P|^2 - 1)^2
```

Dimension check: every static energy-density term has `(M,L,T)=(1,-2,-2)`. The surface tension dimension is
`(1,-1,-2) = M L^-1 T^-2`, i.e. energy per 3-area.

## T1a

`U(rho)=K rho^5/4` has

```text
dU/drho = 5 K rho^4/4
d2U/drho2 = 5 K rho^3
```

For nonnegative density the only stationary scalar endpoint is `rho=0`; there is no pair of scalar vacua and no scalar wall.
The existing stack still imposes the brane ingredients through `V_conf`, `Z(w)`, `W(w)`, `B_ell(w)`, and `k_w` rather than deriving
them from the scalar GNLS sector.

## T1b

The decisive admissible sequence keeps `rho=rho0`, `|P|=1`, and rotates on the vacuum sphere over the finite interval
`w in [-Lhalf,Lhalf]`:

```text
P_L(w) = (sin(pi (w+Lhalf)/(2 Lhalf)), 0, 0, -cos(pi (w+Lhalf)/(2 Lhalf))).
```

This finite-box representative has an in-plane midpoint, `P_L(0)=(1,0,0,0)`, but it is not a localized wall. Both engines compute

```text
K_P0 = m rho0 c_s0^2 a^2 = 5 K a^2 rho0^5
sigma_L = K_P0 pi^2/(4 Lhalf)
        = 5 pi^2 K a^2 rho0^5/(4 Lhalf)
d sigma_L/dLhalf < 0
lim_{Lhalf -> infinity} sigma_L = 0.
```

Thus the minimizing sequence spreads to infinite width, and no natural finite wall thickness or positive finite `sigma_brane` is
derived. A finite-tension radial stationary saddle exists,

```text
P = (0,0,0,tanh(w/(sqrt(2) a)))
sigma_saddle = (10 sqrt(2)/3) K a rho0^5,
```

but it is not the unconstrained minimum, and its core is amplitude collapse `P=0`, not a flat orientational core.

## T1c

Prerecorded long-lived threshold before the scan:

```text
tau_unwind must exceed 1000 * tau_Hubble.
```

The frozen vacuum manifold is

```text
{P in R^4 | P.P = 1} = S^3
pi0 = 0
```

so the `+w` and `-w` vacua are connected. The delocalization channel has

```text
DeltaE_unwind = 0
tau_unwind = not computed: no metastable local minimum against delocalization
```

The clamped spectrum check is not allowed to mask this. The finite-box connected-vacuum control has a clean clamped spectrum
(lowest eigenvalue about `-1.26e-7`, treated as zero at tolerance `1e-3`) while `pi0=0` still flags unwinding. The other controls
also passed: the phi4 kink returned no negative mode plus `pi0=Z2`, and the known saddle returned a negative mode near `-0.987`.
The soft-spin radial saddle's transverse operator also has a negative mode, with analytic
`omega^2 = -c_s0^2/(2 a^2) = -5 K rho0^4/(2 a^2 m)`.

## T1d and T1e

The frozen action is O(4)-isotropic and contains no hidden fixed background vector. The axis is selected only by the boundary
sector, so `w` is emergent in the weak sense, but it is unstable: `W_EMERGENT_BUT_UNSTABLE`.

There is no localized transverse potential well. The spread sequence gives

```text
lim_{Lhalf -> infinity} c_s0^2 pi^2/(4 Lhalf^2) = 0,
```

so the confinement gap vanishes. Existing `V_conf`/`Z`/`W`/`B_ell`/`k_w` remain imposed effective machinery, not derived from this
baseline polar wall.

## Honest provenance and diagnosis

Derived from the frozen postulate: the static functional, coefficient closure, spreading sequence, `sigma_L`, the connected
vacuum manifold, the zero unwinding barrier, the spectrum controls, and the absence of a confinement gap. Assumptions used only
as sensitivity: the fixed-background `rho=rho0` sequence is an admissible upper-bound/minimizing sequence; allowing `rho(w)` to
respond cannot restore a positive localized minimum because this admissible sequence already drives the wall tension to zero.

Diagnosis, not a rescue: a stable antipodal wall would need disconnected vacua, for example an easy-axis or component double-well
that makes `pi0=Z2`. That structure is explicitly pruned by T0 for this baseline and was not run.

Three-way status: emergent `w` survives only as symmetry selection; stable wall fails; a light-capable flat core does not survive
as a localized derived wall.

## Artifacts and commands

JSON:

- `software/stage1_solver/_scratch/pathA_24_T1_wall_mathematica.json`
- `software/stage1_solver/_scratch/pathA_24_T1_wall_sympy.json`

Run commands from `/var/projects/toy_physics`:

```text
timeout 600 python3 software/stage1_solver/tools/pathA_24_T1_wall_sympy.py
timeout 600 math -script software/stage1_solver/tools/pathA_24_T1_wall.wl
```

Both commands exited 0.

Engine agreement: PASS. Mathematica and SymPy agree on all T1 labels, the freeze hash, `sigma_L -> 0`, `sigma_saddle`, `pi0=0`,
`DeltaE_unwind=0`, zero confinement-gap limit, and the spectrum fixture outcomes to the stated tolerance.

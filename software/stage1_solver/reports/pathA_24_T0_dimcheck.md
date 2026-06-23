# pathA_24 T0 Dimcheck

## Verdict

`T0_DIMCHECK: PASS`

Both frozen T0 calculated claims are machine-verified:

- Claim 1, dimensional homogeneity with units restored: PASS by SymPy.
- Claim 2, bulk OP mode speed: PASS by SymPy and independently by Mathematica.
- Engine agreement: PASS. Both engines compute `c_OP = c_s0`, one longitudinal gapped mode, and three transverse modes.

This report does not edit the frozen `freeze-action` block or its hash.

## Artifacts

- SymPy primary checker: `software/stage1_solver/tools/pathA_24_T0_dimcheck_sympy.py`
- Mathematica cross-check: `software/stage1_solver/tools/pathA_24_T0_op_modes.wl`
- SymPy JSON: `software/stage1_solver/_scratch/pathA_24_T0_dimcheck_sympy.json`
- Mathematica JSON: `software/stage1_solver/_scratch/pathA_24_T0_op_modes_mathematica.json`

Both scripts include a fidelity guard requiring the frozen `L_pol` lines from
`reports/pathA_24_T0_freeze.md` section 2.2 to be present unchanged.

## Run Commands

From `/var/projects/toy_physics`:

```text
timeout 600 python3 software/stage1_solver/tools/pathA_24_T0_dimcheck_sympy.py
timeout 600 math -script software/stage1_solver/tools/pathA_24_T0_op_modes.wl
```

Both commands exited 0.

The combined agreement check printed:

```text
T0_DIMCHECK: PASS
term_dimensions_MLT: gradient=(1, -2, -2) kinetic=(1, -2, -2) potential=(1, -2, -2)
c_OP_sympy: sqrt(5*K*rho0**4/m)
c_OP_mathematica: Sqrt[(5*K*rho0^4)/m]
longitudinal_gap: sqrt(10*K*rho0**4/(a**2*m)) / Sqrt[(10*K*rho0^4)/(a^2*m)]
transverse_count: 3
```

## Claim 1: Dimensions

The SymPy checker represents dimensions as exponent triples `(M,L,T)` and
substitutes

```text
c_s^2(rho) = 5 K rho^4 / m
```

from the base dimensions in section 2.3:

```text
[m]   = (1, 0, 0)
[rho] = (0, -4, 0)
[K]   = (1, 18, -2)
[a]   = (0, 1, 0)
[P]   = (0, 0, 0)
```

Computed:

```text
[c_s^2] = (0, 2, -2) = L^2 T^-2
```

Term dimensions computed from the frozen `L_pol` factors:

```text
kinetic   (1/2) m rho a^2 (D_t^v P^i)(D_t^v P^i)                  -> (1, -2, -2)
gradient -(1/2) m rho c_s^2(rho) a^2 (partial_j P^i)(partial_j P^i) -> (1, -2, -2)
potential-(1/4) m rho c_s^2(rho) (P^i P^i - 1)^2                  -> (1, -2, -2)
```

Thus every term has action-density dimension

```text
[L_pol] = (1, -2, -2) = M L^-2 T^-2.
```

The measure gives:

```text
[dt d^4X L_pol] = (0, 4, 1) + (1, -2, -2) = (1, 2, -1) = M L^2 T^-1.
```

## Claim 2: Bulk OP Modes

Both engines linearize the frozen density about

```text
rho = rho0 const, v = 0, P0 = (1,0,0,0), |P0| = 1.
```

The fluctuation split is

```text
P = (1 + sigma, pi1, pi2, pi3),
```

where `sigma` is longitudinal/amplitude and `pi1..pi3` are transverse.

The computed quadratic density is equivalently

```text
L2 =
  1/2 I_P (sigma_t^2 + pi1_t^2 + pi2_t^2 + pi3_t^2)
  - 1/2 K_P (|grad sigma|^2 + |grad pi1|^2 + |grad pi2|^2 + |grad pi3|^2)
  - m rho0 c_s0^2 sigma^2,
```

with

```text
I_P = m rho0 a^2
K_P = m rho0 c_s0^2 a^2
c_s0^2 = 5 K rho0^4 / m.
```

The potential Hessian/mass matrix is

```text
diag(2 m rho0 c_s0^2, 0, 0, 0).
```

Therefore the transverse wave operator is

```text
-I_P omega^2 + K_P |k|^2 = 0,
```

so

```text
omega_T^2 = (K_P/I_P) |k|^2
c_OP^2 = K_P/I_P = c_s0^2 = 5 K rho0^4 / m
c_OP = c_s0.
```

The longitudinal/amplitude mode is gapped:

```text
omega_L^2 = c_s0^2 |k|^2 + 2 c_s0^2/a^2
gap^2 = 2 c_s0^2/a^2 = 10 K rho0^4/(a^2 m)
gap = sqrt(10 K rho0^4/(a^2 m)).
```

Mode count:

```text
longitudinal amplitude modes: 1
transverse Goldstone/spin-wave modes: 3
```

# Actual Branch Protocol V1

This protocol is the next simulation handoff implied by the unincorporated
5PN, barrier, atom, lepton, and `P22` notes.  It is not a candidate generator
yet.  It freezes what a future physical moving-throat exporter must declare
before any target residuals are evaluated.

## Status

The current executable simulation bundle shows that the reduced target-blind
families and the manufactured nonlinear readiness family do not hit the target.
The notes-intake guard in `diagnose_notes_intake.py` records the reason this is
not a license to tune the current packets after the fact:

- the 5PN handoff says the reduced theorem side is essentially closed;
- the remaining task is an actual PDE-selected orbit-lock and outgoing branch;
- support/source enhancement is not the active bottleneck in the reduced stack;
- scalar support or scalar `P0` flux cannot be treated as direct `P22` or
  `N2/N4` moment-shape control;
- atom/lepton notes provide useful finite-throat source physics, but not a
  calibrated actual-branch map into `D0`, `C`, `P0`, `N2`, and `N4`.

## Frozen Inputs

A future actual-branch exporter must freeze these items before target
evaluation:

| Item | Requirement |
|---|---|
| parent status | Declare strict `S_Sigma` promotion or explicitly label the wall equation as an effective closure |
| geometry | Open finite-throat branch with `R_exit > 0`; no hard cap |
| gauge and boundary class | Fixed before solve, with outgoing-port convention declared |
| wall/support basis | Wall coordinates, BdG support basis, Maxwell/mixed ports, and extraction formulas |
| stability gates | Wall positivity, stable BdG/Krein certificate, mixed-sector positivity, `D0 > 0`, and non-dark `N0` |
| source physics | Any leakage/work, `U/V` dressing, microscopic export, finite-throat `P0/P2`, scalar flux, and intrinsic `P22` terms included before residuals |
| freeze hash | Complete protocol hash and source hashes attached to every exported packet |

## Required Packet

The exporter must emit one or more target-blind packets with these fields.

### Packet A

For each grouped lane `A in {20, 21, 22}`:

```text
D_A0, D_A2, D_A4
N_A0, N_A2, N_A4
```

The packet must also include:

```text
mhat_0
chi_Q or N_Q
parent_action_status
boundary_protocol
stability_certificate
source_hashes
freeze_hash
```

### Packet B

The same realized branch must export one equivalent orbit-lock representation:

```text
m_T, m_K, m_mu
```

or:

```text
R_tr, R_nt, R_eta
```

or:

```text
q_tr, q_nt, q_eta
```

## Target Quantities

The simulation aliases are:

```text
D0 = K - B0 - Z0
A  = -D2 = M + B2 + Z2
C  = -D4 = B4 + Z4
```

The one-pole surface is:

```text
D0*C/(3*A^2) = 1
```

The outgoing prefactor moments are:

```text
P0 = N0/D0
P2 = (D0*N2 + 2*A*N0)/D0^2
P4 = (D0^2*N4 + 2*D0*(A*N2 + C*N0) + 3*A^2*N0)/D0^3
```

The constant-prefactor branch requires `P2 = 0` and `P4 = 0`, equivalently:

```text
N2 = -2*A*N0/D0
N4 = N0*(A^2 - 2*D0*C)/D0^2
```

On the one-pole surface:

```text
N4 = -5*A^2*N0/D0^2
```

## Work Packages

| Work package | Output | Guard |
|---|---|---|
| WP0 freeze protocol | Protocol JSON and source hashes | No target modules imported |
| WP1 stationary branch | Open stationary branch state | Refinement and stability certificate |
| WP2 Packet A export | D/N moments, `mhat_0`, `N_Q` or `chi_Q` | No canonical projection after solve |
| WP3 weak-axisymmetric tangent | Tangent data on the same branch | Same branch and same freeze hash |
| WP4 orbit-lock packet | `dln_R_tr`, `dln_R_target`, `dln_epsilon_eta` or equivalent | Support enhancement not mixed into orbit residuals |
| WP5 frozen verdict | Existing V2-22B -> V2-22A -> V2-21 residual report | Post-hoc evaluation only |

## Non-Rescue Rules

- Do not project a realized `chi_Q != 1` branch back to the canonical outgoing
  branch.
- Do not use `zeta` or support enhancement to explain an orbit-lock miss.
- Do not use scalar `P0` mouth hammering as direct `P22` moment-shape control.
- Do not mutate source support, boundary class, gauge convention, port
  normalization, or packet extraction formulas after seeing residuals.
- Do not report an algebraically projected zero-residual packet as a physical
  target-blind hit.

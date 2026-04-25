# Stage V2-25 - Actual Branch Protocol and Notes Intake

## Purpose

This stage records the mathematical handoff created after the first
target-blind simulation miss and the intake of the unincorporated 5PN, barrier,
atom/lepton, and `P22` notes.

The point is to keep the next calculation from being lost in the notes:

- the reduced target algebra is already sharp;
- the current simulated branch families miss the target;
- the unincorporated notes add useful source and mouth physics;
- but those notes do not yet provide a calibrated physical exporter.

So the next load-bearing object is an actual moving-throat branch packet, not a
post-hoc retuning of the present simulation coefficients.

The executable companion for this note is
`simulation/diagnose_notes_intake.py`, which verifies the source anchors and
writes `simulation/output/notes_intake_report.json`.

---

## 1. Status after the target-blind simulation miss

The current simulation bundle freezes candidate packets before target
evaluation, then runs the existing V2-22B -> V2-22A -> V2-21 chain.

The current referee run gives:

- `0/192` target-passing reduced frozen candidates;
- `0/3` target-passing manufactured nonlinear candidates;
- reduced open-stable one-pole ratio
  `D0*C/(3*A^2)` between `0.0033775383274364888` and
  `0.1353664855760648`;
- reduced median required `C` or `D0` multiplier
  `16.30132163440465`;
- reduced median required `P0` multiplier `171.65261223353198`;
- projection stress shows that one-pole support alone and uniform outgoing
  amplitude scaling are both insufficient.

This is evidence against the current reduced and manufactured simulation
families.  It is not a theorem against every possible nonlinear moving-throat
branch.

---

## 2. Coefficient map

On the isotropic grouped real `P2` branch, write the conservative denominator
and outgoing numerator as

```text
D(omega) = D0 + D2 omega^2 + D4 omega^4 + O(omega^6),
N(omega) = N0 + N2 omega^2 + N4 omega^4 + O(omega^6).
```

The full-bundle aliases are

```text
D0 = K - B0 - Z0,
D2 = -(M + B2 + Z2),
D4 = -(B4 + Z4).
```

Equivalently,

```text
A = M + B2 + Z2 = -D2,
C = B4 + Z4 = -D4.
```

Here:

- `K,M` are wall/worldtube stiffness and inertia data;
- `B0,B2,B4` are stable BdG support moments;
- `Z0,Z2,Z4` are conservative Maxwell/mixed moments;
- `N0,N2,N4` are outgoing-transfer moments.

The normalized conservative response is

```text
Y(omega) = D0 / D(omega)
         = 1 + u2 omega^2 + u4 omega^4 + O(omega^6),
```

with

```text
u2 = A/D0,
u4 = (A^2 + D0*C)/D0^2.
```

The one-pole condition is `u4 = 4 u2^2`, so

```text
D0*C = 3*A^2,
```

or

```text
D0*C/(3*A^2) = 1.
```

That is the one-pole ratio reported by the simulation diagnostics.

---

## 3. Outgoing prefactor and moment-shape conditions

The outgoing prefactor is

```text
P(omega) = D0*N(omega) / D(omega)^2.
```

Expanding through `O(omega^4)` gives

```text
P0 = N0/D0,
P2 = (D0*N2 + 2*A*N0)/D0^2,
P4 = (D0^2*N4 + 2*D0*(A*N2 + C*N0) + 3*A^2*N0)/D0^3.
```

The constant-prefactor branch requires

```text
P2 = 0,
P4 = 0.
```

Equivalently,

```text
N2 = -2*A*N0/D0,
N4 = N0*(A^2 - 2*D0*C)/D0^2.
```

On the one-pole surface `D0*C = 3*A^2`, this reduces to

```text
N4 = -5*A^2*N0/D0^2.
```

This is why the projection-stress diagnostic matters.  Raising `C` or `D0` to
the one-pole surface is not enough.  Scaling the outgoing amplitude to fix `P0`
is also not enough.  The branch must realize the outgoing moment shape through
`N2` and `N4` on the same frozen branch.

---

## 4. Actual branch packets

The 5PN handoff reduces the remaining computation to the actual PDE-selected
branch.  The symbolic target surface is not the missing piece; the missing
piece is the realized branch data.

### Packet A

For each grouped lane `A in {20, 21, 22}`, export

```text
D_A0, D_A2, D_A4,
N_A0, N_A2, N_A4.
```

The same packet must also include

```text
mhat_0,
N_Q or chi_Q,
parent_action_status,
boundary_protocol,
stability_certificate,
source_hashes,
freeze_hash.
```

### Packet B

The same realized branch must export one equivalent orbit-lock representation:

```text
m_T, m_K, m_mu
```

or

```text
R_tr, R_nt, R_eta
```

or

```text
q_tr, q_nt, q_eta.
```

The four finish-line conditions are

```text
dln_R_tr = 0,
dln_R_target = 0,
dln_epsilon_eta = 0,
N_Q = 1.
```

The last condition may also be written as `chi_Q = 1` on the natural source-map
branch.

---

## 5. What the unincorporated notes add

The note intake found useful physical ingredients, but not a completed
exporter.

Promotable reduced ingredients:

- leakage/work lane: useful support-side source/work accounting;
- non-rigid `U/V` dressing: useful branch component for mouth/dressing and
  orbit-side response;
- microscopic export kernel: useful odd passive/export term in the active `V`
  equation;
- finite-throat atomic `P0/P2` source: useful replacement for point-source
  forcing;
- open/radiative scalar `P0` flux: possible outgoing-normalization hook;
- intrinsic `P22` bracing: useful mouth-shape source if the half-flux/mixed
  closure is realized.

But the missing calibrated maps remain:

```text
source physics -> D0, C,
source physics -> P0,
source physics -> N2, N4.
```

In particular:

- support/source enhancement is not the active bottleneck in the reduced 5PN
  stack;
- scalar `P0` hammering does not linearly drive the area-preserving `P22`
  mouth mode;
- outgoing transfer moments are support-blind in the explicit BdG pair, so
  support tuning alone cannot supply `N2/N4` control.

---

## 6. No-refit rules

The next exporter must preserve the V2-16 freeze discipline.

Before target residuals are evaluated, freeze:

- parent action status or effective closure declaration;
- gauge convention;
- open-exit boundary class;
- wall/interface action;
- support basis;
- mixed/outgoing port list;
- stability gates;
- extraction formulas;
- source terms admitted from the atom/lepton/barrier notes;
- source hashes and protocol hash.

After residuals are known, do not change:

- support cardinality;
- boundary class;
- gauge convention;
- source normalization;
- port normalization;
- outgoing branch;
- extraction formulas.

Do not project a realized `chi_Q != 1` branch back to `chi_Q = 1`.  Do not use
support enhancement to explain an orbit-lock miss.  Do not report an
algebraically projected zero-residual packet as a physical target-blind hit.

---

## 7. Next calculation

The next honest calculation is:

1. solve or continue a stationary open moving-throat branch;
2. linearize wall, BdG, Maxwell/mixed, and outgoing-port sectors;
3. extract `K,M`, `B_n`, `Z_n`, and `N_n`;
4. export Packet A and Packet B with freeze/source hashes;
5. run the unchanged target-blind guard and post-hoc residual chain;
6. report `R_pole`, `R_norm`, `R_P2`, and `R_P4` without refitting.

That is the point at which the unincorporated notes become testable branch
physics instead of promising but uncalibrated reduced ingredients.

---

## Current claim status

```text
Reduced target algebra: passed.
Target-blind reduced simulation: target miss.
Manufactured nonlinear readiness: target miss.
Notes intake: passed as source/protocol evidence.
Physical actual-branch exporter: not yet implemented.
Post-hoc retuning allowed: no.
Next required artifact: actual moving-throat branch packet.
```

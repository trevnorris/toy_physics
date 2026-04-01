# Local 4PN master/referee assembly — result ledger

## What this step accomplished

This step closes the **local** 4PN chain end to end.

The new referee audit rebuilds the full fixed-chart local 4PN ordinary target from the imported local ADM Hamiltonian target using the exact quartic Hamiltonian-to-ordinary compiler, then verifies that one exact generic-frame ordinary representative of the **entire local sector** exists after combining:

1. the natural generic-frame local seed,
2. the generic-frame lift of the aligned-seed correction,
3. the already solved canonical comparable-mass local residual.

So the local-first 4PN program is now no longer blocked by local existence.

---

## Exact ingredients used

The master audit reuses the exact stage results already frozen in the current 4PN chain:

- `4pn_local_hamiltonian_to_ordinary_audit.py`
  - exact 21-slot quartic COM map,
  - exact local Hamiltonian target,
  - exact ordinary target,
  - strict one-body gate,
  - natural local self/static seed.

- `4pn_generic_frame_ordinary_translation_audit.py`
  - exact sign-flip theorem for the canonical comparable-mass residual,
  - exact aligned ordinary seed,
  - exact ordinary residual slots.

- `4pn_hamiltonian_chart_generic_frame_lift_audit.py`
  - COM reduction map,
  - block-slot extractor,
  - exact image-matrix machinery for the generic-frame lift.

---

## The three decisive checks

### 1. The free block is already fixed by the aligned seed

The aligned ordinary seed reproduces the entire free block of the local ordinary target exactly:

- slot 1 matches the full translated free coefficient,
- slots 2–6 remain zero.

So the free local 4PN ordinary sector does **not** need a separate generic-frame rescue.

### 2. The canonical local residual translates exactly by sign flip

The Hamiltonian-derived generic-frame comparable-mass local residual satisfies

\[
\Delta L_{4,\mathrm{loc}}^{\mathrm{can}} = -\Delta H_{4,\mathrm{loc}}^{\mathrm{can}}.
\]

This removes the earlier ordinary-chart worry about the true comparable-mass residual: the only remaining ordinary obstruction comes from **seed alignment**, not from the comparable-mass local lift itself.

### 3. The aligned-seed correction is exactly generic-frame liftable

The seed-alignment correction

\[
\delta L_{4,\mathrm{seed}} = L_{4,\mathrm{seed}}^{(\mathrm{aligned})} - L_{4,\mathrm{seed}}^{(\mathrm{natural})}
\]

is liftable in the minimal structured seed spaces

\[
Q,
\qquad
T \oplus (pq)T,
\qquad
S \oplus (pq)S,
\qquad
U \oplus (p^2 q^2)U,
\qquad
W=0.
\]

That is the exact generic-frame repair of the old ordinary-chart obstruction.

---

## Final theorem statement for the local sector

Define

\[
L_{4,\mathrm{loc}}^{(\mathrm{gen})}
=
L_{4,\mathrm{seed}}^{(\mathrm{natural,gen})}
+
\delta L_{4,\mathrm{seed}}^{(\mathrm{gen})}
+
\Delta L_{4,\mathrm{loc}}^{(\mathrm{can,gen})}.
\]

Then the master audit verifies that the COM reduction of this generic-frame candidate reproduces the exact fixed-chart local ordinary 4PN target slot by slot in all interaction blocks:

- `Q` block: 5 slots,
- `T` block: 4 slots,
- `S` block: 3 slots,
- `U` block: 2 slots,
- `W` block: 1 slot.

Equivalently, the **entire local 4PN ordinary interaction sector is now assembled exactly**.

---

## Canonical structured-lift sparsity

The final master audit also records one canonical sparse realization of the aligned-seed lift with the following structured basis sizes / nonzero-direction counts:

- `Q`: basis size 52, nonzero directions 13,
- `T`: basis size 92, nonzero directions 14,
- `S`: basis size 58, nonzero directions 12,
- `U`: basis size 20, nonzero directions 8,
- `W`: no lift required.

These are not yet uniqueness theorems for the aligned-seed sector, but they show that the local 4PN lift is concretely realizable rather than merely formal.

---

## What remains open after this step

This step does **not** close the full conservative 4PN theorem.

What remains open is:

1. the sharper uniqueness / interpretation question for the aligned-seed sector,
2. the **tail / hereditary bridge**, which stays separate from the local theorem by construction.

So the next mathematically clean step is to turn back to the tail side and build the quadrupole bridge on top of the now-closed local scaffold.

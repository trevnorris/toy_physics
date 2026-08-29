# Independent build review — S11c-b ADJUDICATION layer (the cross-engine reconcile that decides the number)

## Artifact
`research/pde_ledger_v3/scripts/S11c_b_adjudicated_comparison.py` (Codex-built, uncommitted; imports the
committed comparator `S11c_b_cross_engine_comparator as C` and v1 reconcile `S11c_b_handcoded_comparison as H`).
It applies ONE source-verified algebraic identity (Bridge A: `bRho ↦ B_rho_3/W_0`), routes each joined case by
operand sort, and re-checks zero. Its run reported: **ALGEBRAIC 38 MATCH / 36 FLAG; CONTAINER 0 MATCH / 4 FLAG
(KINETIC); STRUCTURE_INCOMPLETE 69; COVERAGE 84; 12 NAMESPACE_INCOMPLETE; 231 JET_CONSERVED / 0 JET_LOST**.
Directive: `directives/S11c_b_adjudication_build_directive.md`. This layer decides whether the S11c-b coupling
kernel + brane operator AGREE cross-engine, so the ONE failure that matters is a **FALSE MATCH** — a case
reported as reducing to zero that only does so because the layer over-collapsed (a blanket fold, a subset
container match, a Boolean routed as arithmetic, a quotient representative identified, or Bridge A over-reaching).
The pre-adjudication reconcile had only 2 MATCH; this reports 38. Your job is to prove each of those 38 is
GENUINE (renames + Bridge A + integral linearity only) or find the ones that are not.

## Method + constraints (read first)
- Derive/ablate independently; SAVE every ablation script + its literal stdout to named absolute paths under
  your scratch dir and report those paths. Prose "I checked" is discarded.
- ⛔ Copy anything you ablate to /tmp and ablate the COPY; NEVER modify the working tree; NEVER commit or `git
  stash`.
- ⭐ Run every ablation in the FOREGROUND; do NOT background anything and do NOT spawn a monitor/poll loop.
- The layer run is ~11 min (668 s); each ablation re-runs it. Budget for it; run with a generous
  `--residual-leaf-budget` (the residual serialization is timing-dependent — a pre-existing comparator quirk).
- Pure Python (reads committed `.out` transcripts). ⛔ Do NOT spawn Mathematica/wolframscript.

## Required checks — a finding only if it catches a way the comparison/physics could be wrong

1. **NO FALSE MATCH — the critical check, via ablation (MANDATORY, this is a FORM ablation).**
   - `--drop-bridge-a`: re-run and confirm the MATCHes that depend on the bRho normalization move to FLAG (i.e.
     Bridge A is LOAD-BEARING, not decorative, and the base residuals genuinely carry the `W_0` scale). If
     dropping Bridge A leaves the same 38 MATCH, Bridge A is either a no-op or the MATCHes were already zero —
     report which.
   - `--drop-rename <WLname>` on a rename that occurs in a MATCHing family: confirm it is surgical (turns that
     family's MATCH into FLAG) — proving the rename is the same variable, not two different quantities forced
     together.
   - `--collapse-jet <token>=<base>` on a real background-gradient/inertia jet (e.g. `w1_profile_d1`, `u_1_tt`):
     confirm it computes `JET_LOST` and changes the residual — proving the jets are load-bearing and the layer
     is NOT already collapsing them.
   - Then: pick 2-3 of the 38 ALGEBRAIC MATCH cases and INDEPENDENTLY reduce operand_A − operand_B yourself
     (apply only the cited renames + `bRho→B_rho_3/W_0` + integral linearity; nothing else) and confirm zero.
     If any MATCH needs a fold NOT in {renames, Bridge A, integral linearity} to reach zero, it is a FALSE MATCH
     — BLOCKING.

2. **Routing soundness.** (a) ALGEBRAIC entered only for arithmetic `sp.Expr` — confirm no `Equivalent`/`Not`/
   relational/`TextAtom` case is in the ALGEBRAIC route (check the admissibility Boolean operands). (b) The
   CONTAINER route yields a verdict ONLY via a total bijection — verify the KINETIC adapter pairs ALL four
   leaves by cited labels (PY tuple pos0=`U_MOMENTUM_ROWS`/pos1=`THICKNESS_ROW` PY:1573 ↔ WL:851), with no
   dropped/duplicate/unlabeled leaf and no differently-labelled pairing. (c) STRUCTURE_INCOMPLETE cases are
   genuine sort mismatches (not algebraically comparable families wrongly deferred to hide a residual).

3. **Protected one-engine physics.** Confirm `07/10`, `gamma-DivGrad` (08/11), and every ENERGY_BASIS quotient
   representative are NEVER folded or adapted into a MATCH — no substitution or container adapter identifies
   them. A MATCH that depended on equating a protected representative is BLOCKING.

4. **Accounting integrity.** Confirm the classified case-ID multiset EQUALS the comparator-emitted
   `join+py_only+wl_only` multiset (no silent drop, duplicate, or invented case). Confirm the per-route counts
   sum to the total. Try to make the accounting miss a case (e.g. a family with a duplicate key) and confirm it
   raises rather than silently dropping.

5. **Bridge A correctness + reach.** Re-verify from source that `B_rho_3 = bRho·W_0` (θ² coefficient alignment
   WL:472 / PY:1130 / spec:102), and that the SYMBOL substitution rewrites every `bRho` occurrence (bare and
   `bRho·W_bg`) without over-reaching onto a non-compression-modulus use of `bRho`.

6. **Script hygiene.** No assert on a measured payload; the layer PRINTS operands/residuals/verdicts; no
   `PASS/FAIL/VERDICT/target`. Fixtures assert only on synthetic inputs.

## Report
Numbered findings (BLOCKING vs non-blocking), each with file+line, the ablation command + its literal stdout
path, and a concrete correction. Then an overall verdict (SOUND / NOT SOUND) naming what you verified with
which command, and — critically — your independent confirmation of whether the 38 ALGEBRAIC MATCH are genuine
(no false MATCH) or which are not.

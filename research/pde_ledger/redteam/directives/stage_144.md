---
unit_id: 144
batch: IV.5
created_at: 2026-05-29T00:00:00Z
supersedes: 2026-05-27T00:00:00Z
findings_count: 1
stop_cold: null
applied: true
applied_at: 2026-05-30T03:30:23Z
findings_applied: 1
findings_blocked: 0
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 144 (REWRITE, paper-grounded, post-review)

## What this rewrite does

This directive REPLACES the stale 2026-05-27 directive (which predated the Codex
read-only review and split the issue across F1/F2/F3/F4 with two unresolved
"Resolve before fix_loop" blocks). The 2026-05-27 F1 (banner relabel) and F2
(SymPy numerical-target asserts) were ALREADY APPLIED and the Codex review
(`redteam/codex_reviews/stage_144.md`) confirms them PASS and non-tautological —
those are carried as **tainted-applied / reconcile** below (verify literals are
sourced, no re-application).

The review collapses the live blocking issue to ONE finding:

- **R1 (mathematica_transliteration, `.wl:64-71`):** the Mathematica branch/target
  block — in particular the `Pi_*` drift check against the literal
  `1.50882951349316` — is a line-by-line MIRROR of the SymPy block: same `gPi`,
  `gMinus`, same `FindRoot`/`nsolve` seeds, same target constants. A
  transcendental root claim is the load-bearing quantity here, and per the
  MIRROR POLICY a load-bearing transcendental root is **not** a sanctioned mirror.
  The second engine currently gives NO independent evidence for the `Pi_*`
  numerical branch claim: a shared wrong target copied into both engines would
  pass in both.

The R1 fix follows the batch-4 §131-F3 / batch-5 §142 precedent EXACTLY: give the
Mathematica engine an INDEPENDENT route to `Pi_*` via a **cleared-denominator
bracketed `FindRoot`** on a polynomial-in-`(piM, Exp[piM])` residual, then anchor
the independently-found root to the OWNING-stage (131) target. That converts
"engines agree" from a degenerate copy into two genuinely independent root-finds
(SymPy `nsolve` of the rational `gPi - gminus` vs Mathematica `FindRoot` of a
CLEARED-DENOMINATOR residual).

Apply R1 only. After applying, append an `## Applied: R1` block with
`files_changed`, `summary`, and `deviation` (or "none").

Do NOT introduce new features, refactors, or stylistic changes beyond the named
edit. Do NOT touch paper.tex, notes/, or any prose documents. Edit only the
Mathematica `.wl` (the SymPy `.py` is correct as-is — see Reconcile below).

After editing, RUN the Mathematica script (`math -script <path>`, under
`timeout 600`) and iterate until it exits 0 with all in-file checks passing.
Getting it to run cleanly is your job; the orchestrator independently re-runs
afterward and refreshes the committed `output/*.txt`.

---

## R1 — Mathematica transliteration: independent cleared-denominator `Pi_*` route

**File:**
`/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage144_unique_regular_canonical_branch_mathematica_audit.wl`

**Current anchors (confirmed against the live file 2026-05-29):**
- `piStar` is currently found at **line 46** by a NON-bracketed transliteration of
  the SymPy seed:
  ```
  piStar = piM /. FindRoot[gPi == gMinus, {piM, 1.5}, WorkingPrecision -> 80, AccuracyGoal -> 30, PrecisionGoal -> 30, MaxIterations -> 100];
  ```
- The transliterated target block is **lines 64-71**; the load-bearing offending
  line is **line 67**:
  ```
  If[!(Abs[N[piStar - 1.50882951349316`30, 30]] < tol), fail["Pi_* drift", N[piStar, 30]], pass["Pi_* matches notes target"]];
  ```

### Issue
`piStar` is obtained from the SAME root equation (`gPi == gMinus`), SAME seed
(`1.5`), and then compared to the SAME literal target as the SymPy script. Both
engines share `gPi`'s exact symbolic form. There is no independent corroboration
of the `Pi_*` numerical branch claim.

### Required change — two parts

**(R1a) Replace the line-46 `FindRoot` with an INDEPENDENT cleared-denominator
bracketed `FindRoot`.**

Derivation of the residual (from 144's ACTUAL `Pi_*` defining equation
`gPi(piStar) == gMinus`, with `gPi` as defined at line 35):

```
gPi = 2*piM*(2*piM*Exp[piM] + Pi) / ((4*piM^2 + Pi^2)*(Exp[piM] - 1))
```

The denominator `(4*piM^2 + Pi^2)*(Exp[piM] - 1)` is strictly positive for
`piM > 0` (both factors positive: `4*piM^2 + Pi^2 > 0`, and `Exp[piM] - 1 > 0`).
Multiplying `gPi == gMinus` through by it gives the polynomial-in-`(piM,Exp[piM])`
residual (SIGN: keep `(Exp[piM] - 1)`, the positive factor — do NOT flip to
`(1 - Exp[piM])`; that is the §131-F3 sign-error trap that shipped a residual
~6366 instead of ~0):

```
gThresholdResidual[p_] := 2*p*(2*p*Exp[p] + Pi) - gMinus*(4*p^2 + Pi^2)*(Exp[p] - 1)
```

(Here `Pi` is Mathematica's symbol pi = 3.14159..., `p`/`piM` is the variable Pi
of the paper; this matches the verified §131 `gThresholdResidual` up to the
overall nonzero `20*Pi` scale that §131 needed only because its `gMinus` carried
a `/(20 Pi)` factor — 144's `gMinus` is the simplified `r - Sqrt[1+r^2]/2`, so the
direct form above is correct.)

Replace line 46 with:

```
(* INDEPENDENT Pi_* route (NOT a transliteration of SymPy's nsolve on the *)
(* rational gPi == gMinus): clear the strictly-positive denominator *)
(* (4 piM^2 + Pi^2)(Exp[piM]-1) so the root equation is polynomial-in- *)
(* (piM, Exp[piM]), and isolate the unique positive root with a bracketing *)
(* seed pair. gPi is monotone on (0, Infinity), so the bracket robustly fixes *)
(* the root. Precedent: batch-4 stage 131-F3, batch-5 stage 142. *)
gThresholdResidual[p_] := 2*p*(2*p*Exp[p] + Pi) - gMinus*(4*p^2 + Pi^2)*(Exp[p] - 1);
piStar = piM /. FindRoot[gThresholdResidual[piM] == 0, {piM, 1.4, 1.6},
  WorkingPrecision -> 80, AccuracyGoal -> 30, PrecisionGoal -> 30, MaxIterations -> 100];
```

Note `gMinus` is defined at line 33 as `FullSimplify[r - Sqrt[1 + r^2]/2, ...]`, so
it is already a concrete real before this line — fine to use inside the residual.
The bracket `{piM, 1.4, 1.6}` straddles the known root `~1.5088` and matches the
§131 bracket.

**(R1b) Re-anchor the `Pi_*` check (line 67) to the OWNING-stage target with an
independence label, and add an explicit cleared-denominator residual-near-zero
assert so the bracketed root is provably a root (not just self-consistent).**

The canonical `Pi_*` literal MUST be the **owning stage 131** value, NOT 144's own
solver output. Source: `scripts/output/moving_throat_pde_stage131_parent_mouth_threshold_sympy_audit.txt:2`
(stage 131 nsolve, prec 80):
```
Pi_* = 1.50882951349315558300555075595
```
Do NOT anchor to 144's transliterated `1.50882951349316` (that is the truncated
notes display value) nor to 144's own nsolve `...5274704351177` (which diverges at
~digit 16 — the §142 trap). Stage 144's own Mathematica FindRoot already lands on
`...5830055...` (see the current `.txt:12`), so anchoring to 131's full value at
tol `10^-12` is satisfiable.

Replace line 67:
```
If[!(Abs[N[piStar - 1.50882951349316`30, 30]] < tol), fail["Pi_* drift", N[piStar, 30]], pass["Pi_* matches notes target"]];
```
with:
```
(* (R1, INDEPENDENT root anchor) piStar above came from the cleared-denominator *)
(* FindRoot, NOT from gPi == gMinus; assert the residual is genuinely ~0 there, *)
(* then anchor to STAGE 131's owning Pi_* (its independent nsolve, prec 80): *)
(* scripts/output/moving_throat_pde_stage131_parent_mouth_threshold_sympy_audit.txt:2 *)
piStarResidual = Chop[N[gThresholdResidual[piStar], 30], 10^-30];
If[!(piStarResidual === 0 || Abs[piStarResidual] < 10^-25),
  fail["Pi_* cleared-denominator residual not ~0", piStarResidual],
  pass["Pi_* solves cleared-denominator residual (independent root)"]];
piStarOwner = N[Rationalize[1.50882951349315558300555075595, 0], 30];
If[!(Abs[N[piStar - piStarOwner, 30]] < tol), fail["Pi_* drift vs stage 131 owner", N[piStar, 30]], pass["Pi_* matches stage-131 owning value (independent route)"]];
```

Leave lines 64-65 (`g_+^F1 > 1`, `2/Pi < g_-^F1 < 1` brackets), line 66
(`tol = 10^(-12)`), and lines 68-71 (`That_*`, `Sigma0_*`, `Pi_match`,
`That_match` checks) as-is — those are not the load-bearing transcendental-root
claim and the SymPy F2 block already corroborates them (see Reconcile).

### Why this is independent (anti-tautology)
- The SymPy root comes from `nsolve(Eq(gPi, gminus), 1.5)` on the **rational**
  `gPi - gminus` (denominator left in).
- The Mathematica root now comes from `FindRoot` on the **cleared-denominator
  polynomial-in-`(piM,Exp[piM])`** residual with a **bracketing** seed pair.
- The compared target is the **owning stage 131** value, located by stage 131's
  OWN independent nsolve — not copied from 144's own solver.
- The new `piStarResidual === 0` (Chop'd) assert can FAIL if the bracket converges
  to a non-root or the residual sign is wrong (the §131-F3 trap), so it exercises
  the claim rather than restating it.

### Expected new PASS lines (Mathematica)
Two `Pi_*` lines now appear in place of the single old one:
```
PASS: Pi_* solves cleared-denominator residual (independent root)
PASS: Pi_* matches stage-131 owning value (independent route)
```
plus the unchanged `upper branch g_+^F1 > 1`, `lower branch bracket ...`,
`That(Pi_*) ...`, `Sigma0(Pi_*) ...`, `Pi_match ...`, `That(Pi_match) ...`.

**Verification command:**
After Codex applies, the verifier runs `redteam exec-mathematica 144` and confirms
it exits 0; that `grep -n "1.50882951349316\`30" <wl>` returns no match for the
OLD truncated literal in the `Pi_*` assert (the owning value `...5830055075595`
is now used); and that the transcript shows BOTH new `Pi_*` PASS lines above. If
the residual assert FAILS (sign/bracket problem), Codex must NOT relax tolerances
— `## Blocked: R1` with the residual value and the question (this is the §131-F3
sign-check gate).

---

## RESOLVED (Claude+Codex consult Q2, `redteam/codex_reviews/_consult_batch6.md` — CONCUR, non-conceptual)

This was a MATH-COVERAGE call (which engine route corroborates the load-bearing
root), not a conceptual change. The read-only consult CONCURRED with no conceptual
concern, so `needs_user_resolution` is now **cleared (false)**. Resolution:

- **Adopt the independent-`Pi_*` route (R1 above).** Reject "accept transliteration
  for this stage" (stale-directive F4 option (a)): a load-bearing transcendental
  root is explicitly NOT a sanctioned mirror under the MIRROR POLICY, and §131/§142
  already set the precedent that THIS root gets a cleared-denominator independent
  Mathematica `FindRoot`. The cleared-denominator residual gives the Mathematica
  engine a genuinely non-transliterated derivation of `Pi_*`; anchoring to stage
  131's owning value (not 144's own solver) closes the §142 "self-anchor" trap.
- Codex confirmed (Q2): `gThresholdResidual` is the correct cleared-denominator form
  with the SIGN correct (`(Exp[p]-1)`, the positive factor — NOT the §131-F3
  `(1-Exp[p])` trap), and anchoring to the stage-131 owning value is the right
  non-self-anchor.
- This subsumes the stale F4 entirely (F4's options (b)/(c) — alternative `Sigma_0`
  fixed-point or `Limit[]` endpoint checks — are not needed once the root itself is
  independently corroborated).

---

## Paper-card downgrade status (stale F3) — RESOLVED-no-fix (already applied IV.5)

The stale directive's F3 asked to downgrade the paper-card `Checks` items
"(i) gain pair `(M_s,M_q)` vs outlet consistency" and "(ii) self-matched
susceptibility closure" to forward-carry citations. **This was ALREADY DONE in an
earlier batch (IV.5).** Confirmed in `paper/stages/stage_144.tex:21-25`:
```
\item Outlet consistency of the gain pair \((M_s,M_q)\) is verified upstream at \ref{stage:135}; this stage carries that result forward.
\item Self-matched susceptibility closure is verified upstream at \ref{stage:140}; this stage carries that result forward.
\item Check numerical fixed points are recorded as numerically located, not closed-form constants.
```
Both items (i) and (ii) already cite `\ref{stage:135}` / `\ref{stage:140}`. **No
paper edit is prescribed here**, and in any case paper edits are Codex-applied +
Claude-reviewed in the paper loop, never in this scripts-only fix loop. Record as
RESOLVED-no-fix.

---

## Reconcile (tainted-applied — verify, do NOT re-apply)

The following were applied 2026-05-27 (orchestrator-direct era) and the Codex
review confirms them PASS and non-tautological. Codex should VERIFY their target
literals are sourced as cited below (cite file:line) and that they still pass
after the R1 edit; do NOT modify them.

- **SymPy banner relabel** (`.py:16`, `.py:70`): `STAGE 144 ...` / `STAGE 144
  LEDGER` — present, correct.
- **Mathematica banner relabel** (`.wl:27`, `.wl:73`): `STAGE 144 ...` / `STAGE
  144 LEDGER` — present, correct.
- **SymPy `Sigma0` reporting** (`.py:40-43`): prints `Sigma0(Pi_*)` /
  `Sigma0(Pi_match)` — print-only, fine.
- **SymPy branch inequalities** (`.py:53-56`): `g_+^F1 > 1` and
  `2/pi < g_-^F1 < 1` — can fail, exercises branch exclusion/bracketing. Targets
  are structural (no literal). Sourced: notes §1 (`g_+^F1 ~ 2.79795...`,
  `g_-^F1 ~ 0.758035...`) and §2 (`2/pi ~ 0.636619...`).
- **SymPy numerical-target asserts** (`.py:57-68`): `Pi_*`, `That_*`, `Sigma0_*`,
  `Pi_match`, `That_match`. VERIFY each literal is sourced from the notes file
  `notes/stages/moving_throat_pde_stage144_unique_regular_canonical_branch.md`,
  NOT fabricated:
  - `Pi_* = 1.50882951349316` — notes line 63.
  - `That_* = 0.901484054174205` — notes line 78.
  - `Sigma0_* = 1.80594111095636` — notes line 77.
  - `Pi_match = 1.90848600654854` — notes line 84.
  - `That_match = 1.01132972803599` — notes line 86.
  These are the truncated notes display values at tol `10^-12`; they remain valid
  for the SymPy side (the SymPy script is the engine the notes literals were read
  into). Only the Mathematica `Pi_*` ANCHOR is upgraded to the owning stage 131
  full value (R1b) because that is the cross-engine corroboration line.

## Anti-X-X guard

The new R1 checks must not reuse the compared quantity's own derivation primitive.
- The `piStarResidual === 0` assert evaluates the CLEARED-DENOMINATOR residual at
  `piStar`; `piStar` was found by `FindRoot` on that SAME residual, so this assert
  is a convergence/sign witness (it CAN fail on the §131-F3 sign error), NOT a
  paper-claim test — that is its only job and it is labelled as such.
- The load-bearing paper-claim test is the SECOND assert: `piStar` (cleared-denom
  bracketed FindRoot, Mathematica) vs stage-131's OWNING `Pi_*` (independent SymPy
  nsolve at prec 80). The two roots come from DIFFERENT equations (cleared-denom
  polynomial vs rational) in DIFFERENT engines, so equality is non-vacuous.
- Never anchor 144's Mathematica `Pi_*` to 144's own nsolve/FindRoot output.

---

## For orchestrator/Codex consult (feasibility uncertainty)

1. **Bracketed `FindRoot` convergence.** The §131 verified script uses the IDENTICAL
   bracket `{piM, 1.4, 1.6}` and the same residual shape on the same `gPi`, and the
   current 144 FindRoot already lands on `1.5088295134931555830055...` — so
   convergence at `WorkingPrecision -> 80` is expected. If Mathematica's bracketed
   `FindRoot` is fussy about the `Exp[piM]` in the polynomial residual at WP 80,
   fall back to the verified §131 form `40*Pi*p*(2*p*Exp[p]+Pi) - 20*Pi*gMinus*
   (4*p^2+Pi^2)*(Exp[p]-1)` (same root, overall `20*Pi` scale) — both are
   cleared-denominator and independent of the SymPy rational route. Do NOT widen
   the cap beyond `timeout 600`.
2. **Residual-zero tolerance.** `Chop[..., 10^-30]` then `Abs[...] < 10^-25` is the
   near-zero gate at WP 80; if `Chop` to exact `0` is flaky on the transcendental
   residual, the `< 10^-25` fallback covers it. If the residual is NOT ~0 (e.g.
   ~6366), that is the §131-F3 SIGN ERROR — `## Blocked: R1`, do not relax.
3. **`needs_user_resolution` stays true** pending the consult CONCUR; orchestrator
   finalizes.

## Applied: R1

- files_changed:
  - `mathematica/moving_throat_pde_stage144_unique_regular_canonical_branch_mathematica_audit.wl`
- summary: Replaced the Mathematica `Pi_*` solve with a cleared-denominator bracketed `FindRoot` and anchored the resulting root to the stage-131 owner value with an explicit residual check.
- deviation: none

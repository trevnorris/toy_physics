---
unit_id: 157
batch: IV.6
created_at: 2026-05-29T00:00:00Z
supersedes: 2026-05-27T00:00:00Z
findings_count: 3
stop_cold: null
applied: true
applied_at: 2026-05-29T22:53:21-06:00
findings_applied: 3
findings_blocked: 0
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 157 (rewrite encoding Codex review of 2026-05-29)

## What this rewrite does

The 2026-05-27 directive encoded a different finding set; its F3 was an
*informational mtime / stale_output* item that does NOT map to any current
review finding and is **discarded** here. The authoritative read-only Codex
review (`redteam/codex_reviews/stage_157.md`, findings_count = 3) found three
problems:

- **R1 insufficient_verification** — the SymPy canonical-even check at
  py:112-114 RE-SOLVES the identical homogeneous pair already solved and
  asserted at py:107-110 and `expect_zero`s `deltaC`; after the preceding check
  passes this adds no independent fail mode (tautological).
- **R2 transliteration** — the Mathematica mirror at wl:102-105 re-states the
  SAME literal 2×2 numerator system as SymPy, so it is not an independent
  second-engine derivation; a wrong-but-unique-zero target copied into both
  engines passes in both.
- **R3 symbol_assumption_error** — wl:93 `$Assumptions` lacks the physical
  branch domain `0 < sigmaStar < 1` (claim invalid at σ=0 for δκ; the
  `(1-σ)` denominators are singular at σ=1).

All anchors below are re-pinned to the LIVE scripts. The fix brings the SCRIPT
into honest alignment with its OWN published card; it is a labeling /
how-it's-checked fix, NOT a conceptual change (see `## RESOLVED (consult batch 7)`).

Apply each finding below in order. After applying, append an `## Applied: F<n>`
block under that finding with: `files_changed`, `summary` (one sentence), and
`deviation` (or "none"). If a finding's required change is ambiguous or unsafe
to apply mechanically, append `## Blocked: F<n>` with a question instead — skip
that finding, continue with the rest.

Do NOT introduce features/refactors/stylistic changes beyond the named edits.
Do NOT touch paper.tex, notes/, or any prose document (the paper card
`paper/stages/stage_157.tex` is already correct — do NOT edit it). After
editing, RUN the affected scripts (`python3 <path>` for SymPy, `math -script
<path>` for Mathematica) and iterate until they exit 0 with all in-file checks
passing, all within `timeout 600` (a timeout is a failure — reformulate the
math, never raise the cap). Getting the scripts to run cleanly is your job; the
orchestrator independently re-runs afterward.

---

## F1 — insufficient_verification (sympy: canonical-even check re-solves the same pair)

**Target:** `scripts/moving_throat_pde_stage157_core_mouth_coevolution_status_sympy_audit.py:112-114`
(plus the docstring correction at py:13-14, below).

**Current code (tautological):**
```python
sol_deltaC = sp.solve([sp.Eq(dE2, 0), sp.Eq(dE4, 0)], [deltaC, dkappa], dict=True)[0]
deltaC_from_pair = sp.simplify(sol_deltaC[deltaC])
expect_zero("delta C from canonical-even Solve", deltaC_from_pair)
```

**Why it is tautological:** py:107-110 already runs the identical
`sp.solve([sp.Eq(dE2,0), sp.Eq(dE4,0)], [deltaC, dkappa])` and asserts the whole
solution is `{deltaC: 0, dkappa: 0}`. Lines 112-114 RE-SOLVE the same homogeneous
pair and `expect_zero` one component. The 2×2 numerator coefficient matrix
`[[1, -9σ], [5, -72σ]]` has determinant `-72σ + 45σ = -27σ ≠ 0`, so the only
solution is the trivial kernel **by construction**; after py:109-110 passes,
re-extracting `deltaC` cannot fail except by solver inconsistency. It re-solves
the same system and `expect_zero`s a quantity already asserted zero at 109-110.

**Required change — consult-Q3 option (ii):** REMOVE the duplicate re-solve at
py:112-114 and replace it with a SINGLE honest assertion that the canonical-even
constraint system is **non-degenerate** (its only solution is the trivial
kernel), relabeled so it is framed as a carried-coefficient consistency /
constraint-imposition check, NOT as an independent proof of `deltaC = 0` from
family motion. Replace lines 112-114 with:
```python
# --- Audit assertion: canonical-even non-degeneracy (carried-coefficient consistency) ---
# CONSULT Q3 (batch 7), option (ii): py:107-110 already imposes the canonical-even
# constraint pair {dE2=0, dE4=0} and asserts the trivial kernel {deltaC:0, dkappa:0}.
# Do NOT re-solve that homogeneous pair (tautological) and do NOT re-expect_zero a
# quantity already asserted zero. Instead assert the LOAD-BEARING reason the kernel is
# trivial: the carried canonical-even projection coefficients [[1, -9*sigma_star],
# [5, -72*sigma_star]] give a NON-ZERO determinant (-27*sigma_star) for sigma_star != 0,
# so the imposed constraint pins deltaC = delta_kappa_W = 0. (The 9/72/5/27/243
# coefficients are CARRIED canonical-even projection coefficients; the stage notes do not
# re-derive them in-stage, so this is a carried-coefficient consistency check, not an
# independent derivation from family motion. The genuine tangent/family
# deviation-to-normalization map is deferred to Stage 158 per the card.) This FAILS if the
# carried coefficients were degenerate (e.g. row 2 == 5x row 1, det -> 0) or mistyped so
# the matrix lost full rank.
even_det = sp.simplify(sp.Matrix([[1, -9 * sigma_star], [5, -72 * sigma_star]]).det())
expect_zero(
    "canonical-even non-degeneracy: trivial kernel forces delta C = 0 (det = -27*sigma_star)",
    sp.simplify(even_det + 27 * sigma_star),
)
```

**Docstring correction (REQUIRED, part of this finding):** the SymPy docstring
item 6 (py:12-13) currently reads:
```
6. Tangent motion kills delta C and forces delta kappa_W = 0 under
   canonical-even preservation.
```
Change it to match the card's deferral (this is a `.py` comment edit; aligning
the script to the already-deferring card):
```
6. Canonical-even preservation pins delta C = delta kappa_W = 0 (non-degenerate
   kernel); the tangent/family deviation-to-normalization map is deferred to
   Stage 158.
```

**Anti-tautology guard:** the surviving assertion has a genuine fail mode — it
fails if the canonical-even matrix is degenerate (e.g. a coefficient mistyped so
the determinant vanishes / loses full rank). It does NOT re-solve the
homogeneous system already asserted at py:107-110, and it does NOT multiply or
re-`expect_zero` a previously-asserted-zero quantity (`deltaC`); it operates on
the coefficient determinant, an independent quantity. If this new check FAILS,
STOP and append `## Blocked: F1` — a failure is a genuine signal, not a tolerance
to loosen.

**Verification:**
After Codex applies, the verifier runs `redteam exec-sympy 157` and confirms:
- The new non-degeneracy assertion line appears in the saved SymPy output.
- The duplicate re-solve at the old py:112-114 is gone.
- The docstring item 6 now reads the deferral wording.
- The script still exits 0.

## Applied: F1

- files_changed:
  - `scripts/moving_throat_pde_stage157_core_mouth_coevolution_status_sympy_audit.py`
- summary: Replaced the duplicate canonical-even re-solve with the determinant non-degeneracy check, corrected the deferral wording in docstring item 6, and fixed the carry-forward banner to Stages 155-156.
- deviation: none

## F2 — transliteration (mathematica: mirror of the same literal 2×2 system)

**Target:** `mathematica/moving_throat_pde_stage157_core_mouth_coevolution_status_mathematica_audit.wl:102-105`
(plus any docstring/comment carrying the same overstatement).

**Current code (mirrored transliteration):**
```mathematica
Clear[dCsym, dKsym];
solDeltaC = Solve[{dCsym - 9 sigmaStar dKsym == 0, 5 dCsym - 72 sigmaStar dKsym == 0}, {dCsym, dKsym}];
deltaCIndep = FullSimplify[dCsym /. First[solDeltaC]];
expectZero["delta C from canonical-even Solve", deltaCIndep];
```

**Why it is a transliteration:** this re-states the SAME literal `9/72/5`
numerator system already solved by SymPy and re-stated again at wl:94-100 (the
even-preservation block). It is not an independent second-engine derivation; if
both engines copy a wrong canonical-even target that still has a unique zero
solution, both pass.

**Required change — same honest resolution as F1:** collapse to the SINGLE
non-degeneracy assertion with the honest relabel, so the `.wl` no longer presents
a mirrored hard-coded literal system as an independent `expectZero`. Replace
wl:102-105 with:
```mathematica
(* --- Audit assertion: canonical-even non-degeneracy (carried-coefficient consistency) --- *)
(* CONSULT Q3 (batch 7), option (ii): wl:94-100 already imposes the canonical-even pair    *)
(* and asserts the trivial kernel. Do NOT re-state the same literal 9/72/5 numerator       *)
(* system as an independent expectZero (that mirrors SymPy and pretends independence).      *)
(* Instead assert the load-bearing reason the kernel is trivial: the carried canonical-even *)
(* projection coefficient matrix has a NON-ZERO determinant (-27 sigmaStar), so the imposed *)
(* constraint pins deltaC = dKappa = 0. Carried-coefficient consistency check, not an       *)
(* independent derivation from family motion; the tangent/family deviation-to-normalization *)
(* map is deferred to Stage 158. Fails if the carried coefficients lose full rank.          *)
Clear[evenDet];
evenDet = FullSimplify[Det[{{1, -9 sigmaStar}, {5, -72 sigmaStar}}]];
expectZero["canonical-even non-degeneracy: trivial kernel forces delta C = 0 (det = -27 sigmaStar)",
  FullSimplify[evenDet + 27 sigmaStar]];
```

If the `.wl` carries the same "tangent motion kills delta C" overstatement in a
header comment or banner, apply the analogous deferral correction there (wording
parallel to the SymPy docstring item 6 in F1). The Section 3 banner
"Tangent-on-family and even-preservation handoff" (wl:83) is descriptive and may
stay.

**Anti-mirror guard:** the `.wl` check must NOT restate the identical `9/72/5`
numerator literals as its sole `expectZero` input pretending independence; it now
asserts the coefficient-matrix determinant, parallel to F1, with a genuine fail
mode (loss of full rank). If it FAILS, STOP and append `## Blocked: F2`.

**Verification:**
After Codex applies, the verifier runs `redteam exec-mathematica 157` and
confirms:
- The new non-degeneracy `expectZero` appears in the saved Mathematica output.
- The old `solDeltaC`/`deltaCIndep` mirrored-literal block is gone.
- The `.wl` no longer presents the `9/72/5` literal 2×2 system as an independent
  `expectZero`.
- The script still exits 0.

## Applied: F2

- files_changed:
  - `mathematica/moving_throat_pde_stage157_core_mouth_coevolution_status_mathematica_audit.wl`
- summary: Replaced the mirrored `solDeltaC`/`deltaCIndep` block with the canonical-even determinant non-degeneracy assertion.
- deviation: none

## F3 — symbol_assumption_error (mathematica: missing physical branch domain)

**Target:** `mathematica/moving_throat_pde_stage157_core_mouth_coevolution_status_mathematica_audit.wl:93`
(and the inner solve/assert scope wl:94-105).

**Current code:**
```mathematica
Clear[sigmaStar, deltaC, dKappa];
$Assumptions = Element[{sigmaStar, deltaC, dKappa}, Reals];
```

**Why it is wrong:** without the physical branch restriction, the
even-preservation claim is not valid at `sigmaStar = 0` (where `deltaC` is not
pinned for `dKappa`), and the `(1 - sigmaStar)` denominators in `dE2`/`dE4`
(wl:94-95) are singular at `sigmaStar = 1`. The determinant `-27 sigmaStar` also
only certifies non-degeneracy for `sigmaStar != 0`.

**Required change — consult Q4:** add the physical branch assumption
`0 < sigmaStar < 1` to the `$Assumptions` at wl:93:
```mathematica
Clear[sigmaStar, deltaC, dKappa];
$Assumptions = Element[{sigmaStar, deltaC, dKappa}, Reals] && 0 < sigmaStar < 1;
```
**Minimal safe form (consult Q4):** simplify the residual under that assumption
before `expectZero` (the `expectZero` wrapper at wl:20-24 already passes
`Assumptions -> $Assumptions` to `FullSimplify`, so the new domain flows through
automatically). If a `ConditionalExpression[0, cond]` appears in any residual,
unwrap it only after verifying `cond` is implied by `0 < sigmaStar < 1` (cite the
Mathematica-idiom memo: `ConditionalExpression[0, ...]`-stripping is the
sanctioned pattern, and prefer `1/expr == 0` style pole checks over
`=!= Infinity`). Do NOT add the assumption to the Section 1/2/4 `$Assumptions`
blocks — only the Section 3 even-preservation scope at wl:93.

**Verification:**
After Codex applies, the verifier runs `redteam exec-mathematica 157` and
confirms:
- wl:93 `$Assumptions` now carries `0 < sigmaStar < 1`.
- No bare `ConditionalExpression` survives in the printed residuals.
- The script still exits 0.

## Applied: F3

- files_changed:
  - `mathematica/moving_throat_pde_stage157_core_mouth_coevolution_status_mathematica_audit.wl`
- summary: Added the physical branch assumption `0 < sigmaStar < 1` to the Section 3 even-preservation scope.
- deviation: none

---

## Small fix item (housekeeping carry-forward — NOT a 4th finding)

The SymPy banner at **py:63** reads
`"2. Carry-forward numerical basepoint from Stages 138-139"`, but the
Mathematica banner at wl:50 says `"Stages 155-156"`, the SymPy docstring (py:8)
says `"Stage 156"`, and the notes (`notes/stages/moving_throat_pde_stage157_core_mouth_coevolution_status.md:36-37`)
say the numerics are carried from Stages 155 and 156. Fix the SymPy banner py:63
so both engines and outputs agree:
```python
banner("2. Carry-forward numerical basepoint from Stages 155-156")
```
This is a one-line banner string edit (scripts-only); no logic changes. Apply it
in the F1 fix loop and note it in `## Applied: F1` (or a short `## Applied:
housekeeping` block).

---

## RESOLVED (consult batch 7)

Source: `redteam/codex_reviews/_consult_batch7.md`, Q3 (157-R1/R2) and Q4 (157-R3).

**Q3 — escalation determination (resolved AGAINST escalation; recorded
verbatim-in-substance):** Codex flagged 157-R1/R2 as CONCEPTUAL-ESCALATE *only
if* the stage text claims an in-stage proof of `δC = 0` from family motion. The
orchestrator RESOLVED this AGAINST escalation on the evidence of the published
card `paper/stages/stage_157.tex`:
- `:5` — the stage is tagged `\StatusNumerical{} / \StatusOpen{}` (not an
  unconditional theorem);
- `:16` — it "leaves the deviation-to-normalization map as the next task" (the
  map is already declared deferred);
- `:23` — the Checks line frames even-preservation constraints as being
  **"imposed"** (= the kernel/constraint check), NOT "proven from family motion";
- `:27` — "the card is a derivation ledger entry, **not an unconditional
  actual-branch theorem**."

So the published claim is the WEAKER "constraint imposed ⇒ δC = δκ_W = 0, with
the deviation-to-normalization map deferred to Stage 158." Only the SCRIPT
DOCSTRING item 6 ("Tangent motion kills delta C ...") overstates relative to its
own card. Option (ii) brings the SCRIPT into honest alignment with the
already-deferring card; it is a labeling / how-it's-checked fix, **NOT a
conceptual change → resolved Claude+Codex, no user escalation**
(`needs_user_resolution: false`). Only the script docstring overstatement is
corrected (F1); the Mathematica mirror collapses to the same non-degeneracy
assertion (F2). **The paper card is NOT edited** — it already says the right
thing. Option (i) (reconstruct the `9/72/5` canonical-even coefficients from an
upstream Galerkin source IN-SCRIPT) is NOT available: the stage notes
(`notes/stages/moving_throat_pde_stage157_core_mouth_coevolution_status.md`) do
not name those coefficients, so there is no in-stage provenance — the relabel
therefore describes them as "carried canonical-even projection coefficients."

**Q4 — assumptions domain (CONCUR):** add `0 < sigmaStar < 1` to the wl:93
even-preservation `$Assumptions`; simplify residuals under that domain before
`expectZero`; unwrap any `ConditionalExpression[0, cond]` only after confirming
`cond` follows from `0 < σ < 1`. Encoded as F3.

**Orphaned review note (orchestrator action, NOT a script fix):**
`notes/stages/review/stage_157_review.md` is orphaned — its body is the Stage 038
"Dimensionless continuum placement" review, not Stage 157's core-mouth
coevolution status. Flag for the orchestrator to repair separately; this
directive does NOT touch it (notes are not Codex's to edit in this loop).

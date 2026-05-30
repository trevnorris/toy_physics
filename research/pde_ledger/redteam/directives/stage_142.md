---
stage: 142
findings_count: 2
applied: true
needs_user_resolution: false
applied_at: 2026-05-29T16:52:19-06:00
findings_applied: 2
findings_blocked: 0
---

# Codex directive — stage 142 (REWRITTEN from authoritative codex_review)

This directive supersedes the prior 2026-05-27 directive, which PRE-DATED the
`redteam/codex_reviews/stage_142.md` review and prescribed two "independence"
checks that the review then flagged as still tautological / transliterated. Only
TWO findings remain open (R1, R2). The F2 external-decimal-target anchors that
the prior directive added are CORROBORATED by the saved outputs and are KEPT
(see `## RESOLVED: F2-kept`). F5 (banner mismatch) is already fixed in the live
scripts (see `## RESOLVED: F5`).

Apply F1 then F2 in order. After applying each, append an `## Applied: F<n>`
block under that finding with `files_changed`, `summary` (one sentence), and
`deviation` (or "none"). If a finding's required change is ambiguous or unsafe to
apply mechanically, append `## Blocked: F<n>` with a question instead — skip that
finding, continue with the rest. Do NOT introduce new features, refactors, or
stylistic changes beyond what each finding names. Do NOT touch paper.tex, notes/,
or any prose documents.

After editing, RUN the affected scripts (`python3 <path>` for SymPy,
`math -script <path>` for Mathematica) and iterate until they exit 0 with all
in-file checks passing. Getting the scripts to run cleanly is your job; the
orchestrator independently re-runs afterward.

## Anti-fabrication rule (load-bearing for BOTH findings)

Do NOT introduce any numeric literal or series coefficient invented out of thin
air, and NEVER produce an "expected" value by evaluating the script's own `gPi`
at runtime — that is exactly the transliteration R2 flagged. Every literal you
add must be one of:
- a value already present and corroborated in the saved outputs / notes (the
  external decimal targets), OR
- the result of an INDEPENDENT symbolic derivation the verifier can re-check by
  hand (the closed-form integral of the upstream mouth-source law, see F2).

---

## F1 — tautological_check (R1)

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage142_selfconsistent_mouth_branch_sympy_audit.py:75-78`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage142_selfconsistent_mouth_branch_mathematica_audit.wl:109`

**Issue:**
`Rq_star_residual = abs(float(Rq_star - sp.Rational(1,4)))` (py) and
`expectApprox["R_q(Pi_*) numeric = 1/4", rQStar, 1/4, 10^-20]` (wl) are
determined by the SAME construction they claim to verify. `Pi_*` is solved from
`gPi(Pi_*) = g_-`, and `R_q = (gPi - r)^2/(1+r^2)` with
`g_- = r - sqrt(1+r^2)/2`, so `R_q(Pi_*) = ((g_- - r)^2)/(1+r^2) = 1/4`
identically AT WHATEVER point the solve lands. For any self-consistent-but-wrong
`r`, or a copied-wrong `gPi` that still self-solves, this still yields exactly
`1/4`. It exercises only the solver/numeric residual, not the paper claim. (The
Mathematica side also uses tol `10^-20`, which is below SymPy's actual `nsolve`
residual `1.945e-18`; the conceptual tautology is independent of the tolerance.)

**Required change:**
Do NOT delete the existing `R_q(Pi_*) = 1/4` check — KEEP it, but RELABEL it as
what it actually is: a redundant solver-consistency check. Then ADD the
genuinely non-tautological anchor the review prescribed: evaluate `R_q` at
**Stage 131's independently-derived** `Pi_*` (`1.50882951349315558300555075595`,
from `scripts/output/moving_throat_pde_stage131_parent_mouth_threshold_sympy_audit.txt:2`,
corroborated to 30 digits by 131's Mathematica engine; Stage 131 found this parent
mouth-threshold bias by a structurally different route — the cleared-denominator
`FindRoot`, batch-4-verified), NOT 142's OWN nsolve output and NOT the point solved
from the same equations. (Claude+Codex consult Q6: CONCUR; 142's canonical point IS
the same lower Family-1 point that 131 owns.) This is non-tautological because
`R_q(Pi_ext)` lands on `1/4` ONLY IF the hardcoded `gPi` closed form genuinely passes
through `g_-` at that externally-fixed, independently-derived bias; perturb `gPi` or
`r` and `R_q(Pi_ext)` drifts off `1/4`.

In **sympy** — replace the current block at lines 75-78 with:

```python
# (R1, solver-consistency) R_q(Pi_*) = 1/4 holds IDENTICALLY at the solved
# point because Pi_* was found from gPi(Pi_*) = g_-, and R_q(g_-) = 1/4 for any
# r. So this is ONLY a redundant solver-residual check, NOT a paper-claim test.
Rq_star_residual = abs(float(Rq_star - sp.Rational(1,4)))
print(f"R_q(Pi_*) - 1/4 (solver-consistency) = {Rq_star - sp.Rational(1,4)}")
# 1e-15 tracks nsolve's actual precision here (residual ~1.945e-18); do NOT
# tighten to 1e-20 (that was too tight and was already loosened in a prior batch).
if Rq_star_residual > 1e-15:
    raise AssertionError(f"R_q(Pi_*) does not equal 1/4 at nsolve'd Pi_* (residual {Rq_star_residual}).")

# (R1, NON-tautological anchor) Evaluate R_q at STAGE 131's INDEPENDENTLY-DERIVED Pi_*
# (NOT 142's own nsolve output, NOT the self-solved point). Stage 131 found this parent
# mouth-threshold bias by a structurally different route (cleared-denominator FindRoot,
# batch-4-verified); 142's gPi crosses g_- there ONLY IF the hardcoded gPi closed form is
# right, so R_q lands on 1/4. A typo in gPi or a wrong r breaks this.
Pi_ext = sp.Float("1.50882951349315558300555075595", 30)  # Stage 131 Pi_* (independent)
Rq_at_ext = sp.N(Rq.subs(Pi, Pi_ext), 30)
Rq_ext_residual = abs(float(Rq_at_ext - sp.Rational(1,4)))
print(f"R_q(Pi_ext) - 1/4 (independent anchor) = {Rq_at_ext - sp.Rational(1,4)}")
if Rq_ext_residual > 1e-12:
    raise AssertionError(f"R_q at external Pi_*={Pi_ext} is not 1/4 (residual {Rq_ext_residual}); gPi or r is off.")
```

In **mathematica** — replace line 109 with both a relabeled consistency check
AND the independent-anchor check:

```mathematica
(* (R1, solver-consistency) Identical at the solved point; keep as a redundant
   solver residual, NOT a paper-claim test. Tol tracks the SymPy nsolve gap. *)
expectApprox["R_q(Pi_*) numeric = 1/4 (solver-consistency)", rQStar, 1/4, 10^-15];
(* (R1, NON-tautological anchor) R_q at STAGE 131's independently-derived Pi_* *)
(* (cleared-denominator FindRoot, batch-4-verified) — NOT 142's own nsolve output. *)
piExt = N[Rationalize[1.50882951349315558300555075595, 0], 30];
rQAtExt = N[rQ /. piM -> piExt, 30];
expectApprox["R_q(Pi_ext) = 1/4 (independent anchor)", rQAtExt, 1/4, 10^-12];
```

**Anti-tautology guard:**
The new anchor is independent ONLY because `Pi_ext` is **Stage 131's** independently-derived
value (a different stage, a structurally different root-finder), NOT 142's own nsolve output
and NOT the output of the same `gPi(Pi)=g_-` solve. Do NOT replace `Pi_ext` with `Pi_star` /
`piStar` (that reintroduces the tautology), and do NOT substitute 142's own line-81 literal
(`...5274704351177`, which diverges from 131's value at digit ~16 because it is 142's own
nsolve output). Do NOT compute `Pi_ext` at runtime from `gPi`; it must be the literal Stage-131
target. The `1e-12` tolerance is comfortably loose: 131's `Pi_*` is accurate to ~30 digits, so
`R_q(Pi_ext)` matches `1/4` to ~1e-20+; the falsifiability comes from the structure (a `gPi`
sign typo shifts `R_q(Pi_ext)` by O(1), not O(1e-12)), not from the tolerance.

**Verification command:**
Verifier runs `redteam exec-sympy 142` and `redteam exec-mathematica 142`;
confirm BOTH new "independent anchor" lines appear and PASS, the relabeled
solver-consistency lines remain, and both scripts exit 0. Verifier should also
mentally confirm the anchor would FAIL under a perturbed `gPi` (e.g. a sign typo
in the numerator) — that is the whole point.

## Applied: F1

- files_changed:
  - `scripts/moving_throat_pde_stage142_selfconsistent_mouth_branch_sympy_audit.py`
  - `mathematica/moving_throat_pde_stage142_selfconsistent_mouth_branch_mathematica_audit.wl`
- summary: Relabeled the self-solved R_q(Pi_*) checks as solver-consistency checks and added the Stage 131 external Pi_* anchor in both engines.
- deviation: none

---

## F2 — transliteration (R2)

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage142_selfconsistent_mouth_branch_mathematica_audit.wl:64-73`
  (the `gPiSeries = Normal[Series[gPi,...]]` block and its sample-residual loop)

**Issue:**
`gPiSeries = Normal[Series[gPi, {piM, 0, 4}]]` then comparing `gPi` against
`gPiSeries` at `piM = {0.1, 0.2, 0.3}` is NOT an independent derivation. A typo
copied into `gPi` is copied verbatim into its own Taylor series, so the residuals
only show that an expression is close to its own truncated series near 0. This
catches nothing in the closed-form encoding.

**Required change:**
REMOVE the self-series check (lines 64-73, i.e. the `gPiSeries` assignment, the
`gPiSampleVals` table, the residual `Print`, the `Do[...fail...]` loop, and the
trailing `pass["g_Pi closed-form/series consistency at small piM"]`) and REPLACE
it with a genuinely INDEPENDENT re-derivation of `gPi` from the UPSTREAM
mouth-source law, then compare that to the hardcoded closed form.

The upstream owner of `gPi` is the Stage 130 / Stage 129 mouth-source law
(`notes/stages/moving_throat_pde_stage130_mouth_bias_map.md §1`, building on
`...stage129_mouth_boundary_layer.md §2`). There, `gPi` is DEFINED — not as the
closed form — but as the projection integral of the normalized exponential
mouth-source profile against the first D/N derivative shape `cos(pi z / (2L))`,
with `L = 1` (the closed form already assumes the normalized interval):

```
sigma_Pi(z) = Pi * Exp[-Pi z] / (1 - Exp[-Pi])          (Stage 129 §2, L=1)
gPi := Integrate[ sigma_Pi(z) * Cos[Pi_geom z / 2], {z, 0, 1} ]   (Stage 130 §1)
```

where `Pi_geom` is the geometric constant pi (3.14159...), distinct from the
bias variable `piM`. Symbolically integrating this and simplifying MUST reproduce
the hardcoded closed form `2 piM (2 piM Exp[piM] + pi)/((4 piM^2 + pi^2)(Exp[piM]-1))`.
This is a true independent route: the integral is built ONLY from `sigma_Pi`
(the source law) and the `Cos` projection shape — it does NOT reuse the hardcoded
`gPi` primitive, so a typo in the hardcoded `gPi` is NOT shared by the integral.

Insert, in place of the removed block (after the definitions block, before
`subbanner["Core-to-mouth reduction"]`):

```mathematica
subbanner["Independent re-derivation of g_Pi from the mouth-source law (Stage 130 §1)"];

(* g_Pi is DEFINED upstream (Stage 129 §2 source law + Stage 130 §1 projection)
   as the projection of the normalized exponential mouth source against the first
   D/N derivative shape Cos[pi z / 2] on the normalized interval z in [0,1] (L=1).
   Re-derive it here by direct symbolic integration and confirm it equals the
   hardcoded closed form. This route is built only from sigma_Pi and Cos, so it
   does NOT share a typo with the hardcoded gPi closed form. *)
Clear[zVar];
sigmaPi = piM*Exp[-piM*zVar]/(1 - Exp[-piM]);              (* Stage 129 §2, L=1 *)
gPiFromSource = FullSimplify[
    Integrate[sigmaPi*Cos[Pi*zVar/2], {zVar, 0, 1},
        Assumptions -> piM > 0],
    Assumptions -> $Assumptions];
gPiDerivResidual = FullSimplify[gPiFromSource - gPi, Assumptions -> $Assumptions];
Print["g_Pi (source integral) = ", fmt[gPiFromSource]];
expectZero["g_Pi closed form = integral of mouth-source law (Stage 130 §1)", gPiDerivResidual];
```

Notes for Codex while iterating:
- `Pi` in this Mathematica file is the GEOMETRIC pi (the bias variable is `piM`),
  matching the existing `Sqrt[4107 - 100*Pi^2]` usage. So `Cos[Pi*zVar/2]` is
  `cos(pi z / 2)` — the correct projection shape. Do NOT use `piM` in the cosine.
- If `Integrate` returns a `ConditionalExpression`, strip it: the integral is
  valid for `piM > 0` (already in `$Assumptions`); take the unconditional branch
  via `FullSimplify` under the assumption, or wrap with the project's
  `ConditionalExpression[expr, _] -> expr` strip idiom before `expectZero`.
- If `FullSimplify` of the difference is slow but the symbolic forms are clearly
  equal, you MAY instead assert equality of `Together`-combined numerators, but
  the FIRST choice is the direct symbolic `expectZero` on the difference. Keep
  total runtime under the 10-minute cap; if `Integrate`+`FullSimplify` hangs,
  see the Blocked path below rather than raising any cap.
- Do NOT add this to the SymPy script — R2 is Mathematica-only (the SymPy script
  has no `gPiSeries` check; its independence comes from F1's external anchor and
  the F2-kept decimal targets). Leave the SymPy file untouched for this finding.

**Anti-tautology guard:**
The independence is REAL only because `gPiFromSource` is rebuilt from `sigmaPi`
(the exponential source law) and the `Cos` projection — NOT from `gPi`. Do NOT
let Codex shortcut this by defining `gPiFromSource := gPi` or by seeding the
integrand from `gPi`. Do NOT replace the integral with `Series[gPi]` again (that
is the exact transliteration being removed). The residual must be a symbolic
zero, not a numeric small-piM coincidence.

**If genuinely blocked:**
If symbolic `Integrate` of `sigmaPi*Cos[Pi*zVar/2]` cannot be made to close to
the exact hardcoded form within the time cap (e.g. it returns an irreducible
hypergeometric/conditional form Mathematica won't simplify), append
`## Blocked: F2` describing exactly what `Integrate` returned, and STOP — do not
fall back to the self-series check. (Orchestrator note: this is the
"needs Claude+Codex consult" branch flagged in the audit summary.)

**Verification command:**
Verifier runs `redteam exec-mathematica 142`; confirm the `gPiSeries` self-check
is GONE, the new `g_Pi closed form = integral of mouth-source law` `expectZero`
PASSES, the SymPy script is unchanged for this finding, and Mathematica exits 0.
Verifier should confirm the integrand is built from `sigmaPi` + `Cos`, NOT from
`gPi`.

## Applied: F2

- files_changed:
  - `mathematica/moving_throat_pde_stage142_selfconsistent_mouth_branch_mathematica_audit.wl`
- summary: Replaced the self-series Mathematica check with a symbolic source-law projection integral that independently re-derives g_Pi.
- deviation: none

---

## RESOLVED: F2-kept (the F2 external-decimal anchors from the prior directive — KEEP)

The prior directive's F2 added five external-decimal-target `expectApprox` /
`expect_close` checks in both engines:
- `g_-^{F1} value` = 0.7580350789446628269196808904 (tol 1e-25)
- `Pi_* value`      = 1.5088295134931555274704351177 (tol 1e-12)
- `S_q(Pi_*) value` = 0.6580759376054292719303153134 (tol 1e-12)
- `Sigma_0(Pi_*) value` = 1.8059411109563538072179672471 (tol 1e-12)
- `That(Pi_*) value` = 0.9014840541742040227024016887 (tol 1e-12)

These are sympy lines 80-84 and Mathematica lines 110-114. The codex_review
table marks all of these PASS ("Yes, against external decimal target"), and the
saved sympy output (`...sympy_audit.txt` lines 28-32) confirms residuals
~1e-29. These are the genuinely non-tautological anchors of the stage. KEEP them
unchanged — do NOT remove, re-tolerance, or rewrite them. They are recorded as
`tainted-applied` (kept) from the prior batch.

## RESOLVED: F5 — paper_misalignment, direction (a), ALREADY APPLIED

The prior directive held F5 (script banners read `STAGE 125` while the canonical
stage is `142`) for user resolution. Direction (a) is selected (pre-resolved by
the orchestrator): the banners were stale labels from when this stage was
numbered 125 (the `142 - 17 = 125` stale-offset class); replace `STAGE 125` with
`STAGE 142` in the scripts; no paper-side change.

HOWEVER, this is already done in the live scripts. A grep of both files finds
ZERO `STAGE 125` occurrences; all banners now read `STAGE 142`:
- sympy line 28: `banner("STAGE 142 — SELF-CONSISTENT MOUTH-BRANCH LAW")`
- sympy line 86: `banner("STAGE 142 LEDGER")`
- wl line 39: `banner["STAGE 142 — SELF-CONSISTENT MOUTH-BRANCH LAW"];`
- wl line 116: `banner["STAGE 142 LEDGER"];`
- wl line 128 already: `Print["Stage 142 Mathematica audit passed."];`

No edit is required for F5. Codex must NOT re-touch the banners. Recording this
as RESOLVED + applied. (No Mathematica comment-terminator pitfall to worry about
since no `STAGE 125` text remains.)

## Tolerance note (settled — do NOT revert)

The `1e-15` tolerance on the F1 solver-consistency check is CORRECT and must be
kept. SymPy's `nsolve` precision here is ~`1.945e-18` (see sympy output line 27),
so the prior directive's `1e-20` was too tight and was already loosened to
`1e-15` in an earlier batch. The directive DOCUMENTS `1e-15` as the
nsolve-precision gap via the inline comment in F1; do NOT revert to `1e-20`. On
the Mathematica side, F1 likewise sets the solver-consistency tol to `10^-15`
(replacing the stale `10^-20`) for the same reason.

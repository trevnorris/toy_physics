---
unit_id: 142
batch: IV.5
created_at: 2026-05-27T00:00:00Z
findings_count: 5
stop_cold: null
applied: false
verification_status: pending
needs_user_resolution: true
---

# Codex directive — unit 142

Apply each non-paper_misalignment finding below in order (F1, F2, F3, F4). After applying each, append an `## Applied: F<n>` block under that finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

F5 is a `paper_misalignment` finding — DO NOT touch it; the orchestrator is holding for user resolution. Do not edit paper.tex, notes/, or scripts to "fix" F5 unless the user has explicitly chosen a direction in a follow-up directive.

If a non-paper_misalignment finding's required change is ambiguous or unsafe to apply mechanically, append `## Blocked: F<n>` with a question instead — skip that finding, continue with the rest.

Do NOT introduce new features, refactors, or stylistic changes. Edit exactly the file:line ranges named.

Do NOT run python or mathematica. Only edit files.

Do NOT touch paper.tex, notes/, or any prose documents.

## F1 — tautological_check

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage142_selfconsistent_mouth_branch_sympy_audit.py:41`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage142_selfconsistent_mouth_branch_mathematica_audit.wl:61`

**Issue:**
The symbolic assertion `R_q(g_minus) - 1/4 == 0` is algebraically guaranteed for any `r` and so does not exercise the specific value of `r_F1`. The check passes whether `r` is `sqrt(4107-100*pi^2)/(10*pi)` or anything else, because `R_q(g_minus) = ((1/2)sqrt(1+r^2))^2/(1+r^2) = 1/4` identically.

**Required change:**
Keep the existing assertion (it confirms the algebraic identity `R_q(g_-) = 1/4` from `g_- = r - (1/2)sqrt(1+r^2)`), but ADD a numeric assertion that pins `R_q(Pi_*) = 1/4` where `Pi_*` is the nsolve'd canonical point. This is a non-tautological test because `Pi_*` is determined by `gPi(Pi_*) = g_-` with the specific `r_F1`; if `r_F1` were perturbed, `Pi_*` would shift and `R_q(Pi_*)` would no longer be exactly `1/4` (or rather, since the structural identity still holds at the nsolve'd point, swap the check to: `gPi(Pi_*) - g_- = 0` as a symbolic-into-numeric residual which DOES depend on `r_F1`).

Concretely, in **sympy** (`scripts/moving_throat_pde_stage142_selfconsistent_mouth_branch_sympy_audit.py`), after the existing line 60 (the AssertionError on solver convergence), insert:

```python
Rq_star_residual = abs(float(Rq_star - sp.Rational(1,4)))
print(f"R_q(Pi_*) - 1/4 = {Rq_star - sp.Rational(1,4)}")
if Rq_star_residual > 1e-20:
    raise AssertionError(f"R_q(Pi_*) does not equal 1/4 at nsolve'd Pi_* (residual {Rq_star_residual}).")
```

In **mathematica** (`mathematica/moving_throat_pde_stage142_selfconsistent_mouth_branch_mathematica_audit.wl`), after the existing line 78 (`expectApprox["Pi_* compensation solve", ...]`), insert:

```mathematica
expectApprox["R_q(Pi_*) numeric = 1/4", rQStar, 1/4, 10^-20];
```

**Verification:**
After Codex applies, the verifier will run `redteam exec-sympy 142` and `redteam exec-mathematica 142` and confirm both new assertions appear AND both scripts exit 0.

## F2 — insufficient_verification

**Target:** `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage142_selfconsistent_mouth_branch_sympy_audit.py:57-60` and `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage142_selfconsistent_mouth_branch_mathematica_audit.wl:77-78`

**Issue:**
The notes pin five canonical-point numerical targets that the script computes but does not anchor: `g_-^{F1} ≈ 0.758035078944663`, `Pi_* ≈ 1.50882951349316`, `S_q(Pi_*) ≈ 0.658075937605429`, `Sigma_0(Pi_*) ≈ 1.80594111095636`, `That(Pi_*) ≈ 0.901484054174205`. None is asserted.

**Required change:**
In **sympy**, after the existing assertion block (after line 60), insert (using the same `abs(float(... - target)) > tol` style as the existing solver-convergence check, since `expect_close` is not in this file):

```python
def expect_close(name, value, target, tol):
    res = abs(float(sp.N(value, 30) - sp.N(target, 30)))
    print(f"{name} residual = {res}")
    if res > tol:
        raise AssertionError(f"{name} off by {res} > tol {tol}")

expect_close("g_-^{F1} value", gminus, sp.Float("0.7580350789446628269196808904", 30), 1e-25)
expect_close("Pi_* value",      Pi_star, sp.Float("1.5088295134931555274704351177", 30), 1e-12)
expect_close("S_q(Pi_*) value", Sq_star, sp.Float("0.6580759376054292719303153134", 30), 1e-12)
expect_close("Sigma_0(Pi_*) value", Sigma_star, sp.Float("1.8059411109563538072179672471", 30), 1e-12)
expect_close("That(Pi_*) value", That_star, sp.Float("0.9014840541742040227024016887", 30), 1e-12)
```

(The `expect_close` def should be inserted near the top with the other helpers, around line 13 next to `expect_zero`, NOT inline in the canonical-point block. If a more idiomatic helper already exists in this file family, reuse it; otherwise define this one local.)

In **mathematica**, after the existing `expectApprox["Pi_* compensation solve", ...]` at line 78, insert:

```mathematica
expectApprox["g_-^{F1} value",      N[gMinus, 30],   N[Rationalize[0.7580350789446628269196808904, 0], 30], 10^-25];
expectApprox["Pi_* value",          piStar,          N[Rationalize[1.5088295134931555274704351177, 0], 30], 10^-12];
expectApprox["S_q(Pi_*) value",     sQStar,          N[Rationalize[0.6580759376054292719303153134, 0], 30], 10^-12];
expectApprox["Sigma_0(Pi_*) value", sigmaStar,       N[Rationalize[1.8059411109563538072179672471, 0], 30], 10^-12];
expectApprox["That(Pi_*) value",    tHatStar,        N[Rationalize[0.9014840541742040227024016887, 0], 30], 10^-12];
```

Tolerances rationale: `g_-` is purely algebraic in `r_F1` so engines should match to ~25 digits; the other four depend on numerically solving for `Pi_*`, where SymPy's nsolve is configured for 30 digits and Mathematica is configured for `WorkingPrecision -> 80` — both engines agreed to ~14 digits in the prior transcript, so `1e-12` is a safe shared tolerance.

**Verification:**
After Codex applies, verifier confirms five new assertion-pass lines per engine and both scripts exit 0.

## F3 — mathematica_transliteration

**Target:** `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage142_selfconsistent_mouth_branch_mathematica_audit.wl:44-49` (definitions block)

**Issue:**
The Mathematica script transcribes SymPy's symbolic forms for `gPi`, `Sq`, and `r` verbatim. Two engines transcribing the same expressions cannot independently catch a typo in either expression.

**Required change:**
Add an independent numerical cross-check block to the Mathematica script ONLY (do not edit the SymPy script for this finding). After the definitions block (current line 49) and before the `subbanner["Core-to-mouth reduction"]` (current line 51), insert:

```mathematica
subbanner["Independent numerical cross-checks (Mathematica)"];

(* Independent series expansion of gPi near Pi = 0. The closed form
   gPi = 2 Pi (2 Pi exp(Pi) + pi) / ((4 Pi^2 + pi^2)(exp(Pi)-1))
   expanded to O(Pi^3) should agree with a direct Series[] of the same
   expression, cross-checking the algebraic encoding. *)
gPiSeries = Normal[Series[gPi, {piM, 0, 4}]];
gPiSampleVals = Table[
    {pVal, N[gPi /. piM -> pVal, 30] - N[gPiSeries /. piM -> pVal, 30]},
    {pVal, {1/10, 2/10, 3/10}}];
Print["g_Pi closed-form vs series residuals at piM={0.1,0.2,0.3}: ", fmt[gPiSampleVals[[All, 2]]]];
Do[
    If[Abs[gPiSampleVals[[i, 2]]] > 10^-3,
        fail["g_Pi closed-form vs series small-piM disagreement", gPiSampleVals[[i, 2]]]],
    {i, 1, Length[gPiSampleVals]}];
pass["g_Pi closed-form/series consistency at small piM"];

(* Independent r_F1 cross-check: r_F1 = sqrt(4107 - 100 pi^2)/(10 pi)
   should satisfy 100 pi^2 (1 + r^2) = 4107. *)
rSquared = 1 + r^2;
rIdentity = FullSimplify[100*Pi^2*rSquared - 4107];
expectZero["r_F1 satisfies 100 pi^2 (1+r^2) = 4107", rIdentity];
```

This gives two independent checks: a small-`piM` series-vs-closed-form sanity check on `gPi`, and an algebraic identity check on `r_F1` that does not appear in the SymPy script.

**Verification:**
After Codex applies, the Mathematica script contains the two new checks in a new `subbanner` block; the SymPy script is unchanged; Mathematica exits 0.

## F4 — hardcoded_result

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage142_selfconsistent_mouth_branch_sympy_audit.py:24,26`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage142_selfconsistent_mouth_branch_mathematica_audit.wl:44,46`

**Issue:**
`r_F1` and the closed form of `S_q(Pi)` are inserted without a comment naming the upstream stage that derives them.

**Required change:**
In **sympy** (`scripts/moving_throat_pde_stage142_selfconsistent_mouth_branch_sympy_audit.py`), insert a comment above line 24 and above line 26:

Above line 24 (the `r = sp.sqrt(...)` line) insert:
```python
# r_F1: Family-1 reduced mixed-core ratio. Carried forward from upstream
# (see notes/stages/moving_throat_pde_stage142_selfconsistent_mouth_branch.md
# section 1; original derivation is in the upstream "shell/mixed core" block
# referenced by paper/stages/stage_142.tex Inputs field).
```

Above line 26 (the `Sq = ...` line) insert:
```python
# S_q(Pi) closed form: carried forward from the self-matched mouth-susceptibility
# closure (Stage 242 / Sigma_0 = (20/9) That_m^2). The closed form here is
# S(Pi, pi/2), evaluated at the fixed second argument pi/2.
```

In **mathematica**, mirror those two comments above lines 44 and 46 using Mathematica comment syntax `(* ... *)`.

**Verification:**
After Codex applies, four new comment blocks appear (two per script) explaining provenance. No assertion changes; no output changes.

## F5 — paper_misalignment

**Subtype:** target_mismatch

**Paper side:**
- `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_142.tex:1` quote: `\section[Stage 142]{Stage 142: Self-Consistent Mouth-Branch Law}`

**Script side:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage142_selfconsistent_mouth_branch_sympy_audit.py:22` quote: `banner("STAGE 125 — SELF-CONSISTENT MOUTH-BRANCH LAW")`
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage142_selfconsistent_mouth_branch_sympy_audit.py:62` quote: `banner("STAGE 125 LEDGER")`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage142_selfconsistent_mouth_branch_mathematica_audit.wl:39` quote: `banner["STAGE 125 — SELF-CONSISTENT MOUTH-BRANCH LAW"];`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage142_selfconsistent_mouth_branch_mathematica_audit.wl:80` quote: `banner["STAGE 125 LEDGER"];`

## Resolve before fix_loop

The scripts banner themselves as **"STAGE 125"** in four places (sympy lines 22, 62; mathematica lines 39, 80) while the paper card, the filenames, and the notes file all unambiguously identify this as **Stage 142**. (The Mathematica script's very last line, line 92, correctly says "Stage 142 Mathematica audit passed." — so even the Mathematica script is internally inconsistent.)

Possible directions (the user picks one):
- (a) Banners are stale labels from when this stage was numbered 125 in an earlier draft — replace all four occurrences of `STAGE 125` with `STAGE 142` in the script files; no paper-side change.
- (b) There is a separate, legitimate "Stage 125" that this script was originally intended to audit, and stage_142.tex was mis-pointed at the wrong script files — verify which script corresponds to which stage and possibly rename or repath; this requires user judgment.
- (c) Stage 125 and Stage 142 are intentionally the same stage in different parts of the paper draft (unlikely; the orchestrator should flag the original numbering decision).

Direction (a) is the strongly likely answer (filenames are `stage142_*`, paper card is `stage_142.tex`, math itself matches what the notes describe for Stage 142), but the audit policy does not allow Codex to auto-resolve a `paper_misalignment` without user confirmation.

The orchestrator will not invoke Codex on F5 until the user has chosen a direction. F1-F4 can be applied independently.

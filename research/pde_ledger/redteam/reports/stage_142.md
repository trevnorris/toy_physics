---
unit_id: 142
batch: IV.5
auditor_model: claude-opus-4-7-1m
audit_date: 2026-05-27T00:00:00Z
verdict: findings
stop_cold: null
findings_count: 5
paper_alignment: partial
scripts_checked:
  sympy: present
  mathematica: present
  engines_agree: true
  outputs_fresh: true
docs_read:
  paper_stage_tex: present
  notes_stage_files:
    - /var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage142_selfconsistent_mouth_branch.md
  paper_appendix: present
---

# Audit unit 142 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_142.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage142_selfconsistent_mouth_branch.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part04.tex` (only contains `\input{stages/stage_142}`; no narrative row)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage142_selfconsistent_mouth_branch_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage142_selfconsistent_mouth_branch_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage142_selfconsistent_mouth_branch_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage142_selfconsistent_mouth_branch_mathematica_audit.txt`

## What the paper claims

Stage 142 is the "Self-Consistent Mouth-Branch Law" stage. The paper card's body quotes (Output): `\[\Pi=\Sigma_0[1-R_q(\Pi)\mathcal S_q(\Pi)]\]` and `\stagefield{Checks}` enumerates three: (i) gain pair `(M_s, M_q)` against outlet consistency, (ii) self-matched susceptibility closure before using one-scalar branch law, (iii) numerical fixed points recorded as numerically located. The notes are richer and authoritative: they boxed-formulate `R_q(\Pi)=(g_\Pi-r_{F1})^2/(1+r_{F1}^2)`, `\Sigma_0(\Pi)=\Pi/(1-R_q(\Pi)S_q(\Pi))`, `That(\Pi)=\sqrt{9\Pi/(20[1-R_q(\Pi)S_q(\Pi)])}`, identify the canonical Family-1 point `g_-^{F1}=r_{F1}-(1/2)\sqrt{1+r_{F1}^2}\approx 0.758035078944663` with corresponding `Pi_* \approx 1.50882951349316`, and quote `R_q(Pi_*)=1/4`, `S_q(Pi_*)\approx 0.658075937605429`, `\Sigma_0(Pi_*)\approx 1.80594111095636`, `That(Pi_*)\approx 0.901484054174205`. Six numerical canonical-point deliverables and the underlying boxed identities are paper-side claims that the script must touch.

## What the script claims to verify

The SymPy script defines symbolic `r`, `gPi`, `Sq`, `Rq`, `Sigma0`, `That`, then performs exactly one symbolic assertion: `expect_zero("R_q(g_minus)-1/4", Rq_minus - 1/4)` at line 41, where `gminus = r - (1/2)sqrt(1+r^2)` (line 39). It also numerically locates `Pi_*` via `nsolve(Eq(gPi, gminus), 1.5)` (line 44) and asserts `abs(g(Pi_*) - g_minus) < 1e-12` (line 59) as a solver-convergence check. The numerical values of `R_q(Pi_*)`, `S_q(Pi_*)`, `Sigma_0(Pi_*)`, `That(Pi_*)`, and `Pi_*` itself are printed but never compared to the notes-stated targets. The Mathematica `.wl` mirrors this structure line-by-line with one additional `expectApprox` on the same `g(Pi_*) - g_minus` quantity.

## Paper ↔ script cross-check

| Paper-side deliverable | Script-side check | Status |
|---|---|---|
| `Pi = Sigma_0[1 - R_q(Pi) S_q(Pi)]` (Output identity) | `Sigma0 = Pi/(1 - Rq*Sq)` (line 28) — defines `Sigma_0` from this identity | match (definitional; not a falsifiable check, but consistent) |
| `R_q(Pi) = (g_Pi - r_F1)^2/(1+r_F1^2)` (boxed in notes) | Line 27 `Rq = simplify((gPi-r)**2/(1+r**2))` | match (definitional) |
| `R_q(Pi_*) = 1/4` | Line 41 `expect_zero("R_q(g_minus)-1/4", Rq_minus - 1/4)` | **mismatch** — checks `R_q(g_minus)=1/4` algebraically (tautological in `r`), not `R_q(Pi_*)=1/4` numerically |
| `Pi_* ≈ 1.50882951349316` | Line 44 `nsolve(Eq(gPi, gminus), 1.5)`; line 59 `abs(g(Pi_*)-g_minus) < 1e-12` | partial — numerically located but not anchored against the notes' value |
| `S_q(Pi_*) ≈ 0.658075937605429` | Only printed, no assertion | missing |
| `Sigma_0(Pi_*) ≈ 1.80594111095636` | Only printed, no assertion | missing |
| `That(Pi_*) ≈ 0.901484054174205` | Only printed, no assertion | missing |
| `g_-^{F1} ≈ 0.758035078944663` | Line 39 `gminus = simplify(r - (1/2)sqrt(1+r^2))`; printed but not asserted against value | partial (formula match; numeric not anchored) |
| Check (i) gain pair `(M_s, M_q)` against outlet consistency | Not checked | missing |
| Check (ii) self-matched susceptibility closure `Sigma_0 = (20/9) That_m^2` | Used as definition `That = sqrt((9/20)Sigma_0)` (line 29), not verified | missing |
| Check (iii) numerical fixed points recorded as numerically located | Pi_* recorded but not anchored to a notes-stated value | partial |

Paper alignment: `partial` — main identity present, but multiple numerical anchors and two of three card-listed checks are missing.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 41 | `expect_zero("R_q(g_minus)-1/4", Rq_minus - 1/4)` (i.e., `simplify(...) == 0`) | Intended to test `R_q(Pi_*)=1/4` | no — algebraically guaranteed for any `r`, see F1 |
| A2 | sympy | 59 | `abs(float(g_star - N(gminus,30))) > 1e-12` (raises) | Confirms `nsolve` converged | partial — convergence only, not the canonical value |
| A3 | mathematica | 61 | `expectZero["R_q(g_minus)-1/4", rQMinus - 1/4]` | Same as A1 | no — algebraically guaranteed, see F1 |
| A4 | mathematica | 78 | `expectApprox["Pi_* compensation solve", gStar, N[gMinus,30], 10^-12]` | Same as A2 | partial — convergence only |

There are no assertions on `Pi_*`, `S_q(Pi_*)`, `Sigma_0(Pi_*)`, `That(Pi_*)` themselves. There is no assertion that verifies the susceptibility-closure constant 20/9 quoted in the notes.

## Findings

### F1 — tautological_check

**Severity:** high
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage142_selfconsistent_mouth_branch_sympy_audit.py:39-41`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage142_selfconsistent_mouth_branch_mathematica_audit.wl:59-61`

**What's wrong:**
The single load-bearing symbolic assertion `R_q(g_minus) - 1/4 == 0` is algebraically guaranteed for ANY value of `r`. Substituting `g_minus = r - (1/2)sqrt(1+r^2)` into `R_q = (g_Pi - r)^2/(1+r^2)`:

```
R_q(g_minus) = (g_minus - r)^2/(1+r^2)
             = (r - (1/2)sqrt(1+r^2) - r)^2/(1+r^2)
             = ((1/2)sqrt(1+r^2))^2/(1+r^2)
             = (1/4)(1+r^2)/(1+r^2)
             = 1/4   identically.
```

This identity holds whether `r = sqrt(4107-100*pi^2)/(10*pi)` or `r = 42`. The assertion never exercises the specific value of `r_F1` carried in from upstream. The script could redefine `r` to anything and this check would still pass.

**Why this matters:**
The stage card and notes attach physical meaning to `R_q(Pi_*) = 1/4` *as a property of the canonical Family-1 compensation point* (i.e., the property that fixes the gain ratio at the specific `Pi_*` for the specific `r_F1`). The script's check verifies a category-theoretic algebraic identity about a manifestly identically-1/4 expression, not the physical claim. If `r_F1` were stale or wrong (e.g., a typo in the radical), nothing in this audit would notice.

**Required change:**
Either (preferred) add an additional numeric assertion that `R_q(Pi_*)` evaluated at the nsolve'd `Pi_*` equals `1/4` to high precision — `expect_close(N(Rq.subs(Pi, Pi_star), 30), Rational(1,4), 1e-25)` — which is a non-trivial property of the specific `Pi_*` solving `g_Pi(Pi_*) = g_minus`, OR replace the symbolic check with the equivalent symbolic identity that does depend on `r`, e.g., `expect_zero("(1+r^2) R_q(g_minus) - (1+r^2)/4", (1+r**2)*Rq_minus - (1+r**2)/4)` — same warning. The numeric route is the right one because it ties the check to the actual `r_F1` value.

**Verification:**
After fix: a new line (sympy ~ line 50, mathematica ~ line 70) appears that evaluates `R_q` at the *numerical* `Pi_*` and compares to `1/4` with tight tolerance. Removing or perturbing `r` would now break the test.

### F2 — insufficient_verification

**Severity:** high
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage142_selfconsistent_mouth_branch_sympy_audit.py:44-60`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage142_selfconsistent_mouth_branch_mathematica_audit.wl:64-78`

**What's wrong:**
The notes quote five explicit canonical-point numerical values (`g_-^{F1} ≈ 0.758035078944663`, `Pi_* ≈ 1.50882951349316`, `S_q(Pi_*) ≈ 0.658075937605429`, `Sigma_0(Pi_*) ≈ 1.80594111095636`, `That(Pi_*) ≈ 0.901484054174205`). None of these is asserted against the notes' value; all five are merely printed. The only numerical "assertion" is `abs(g_star - gminus) < 1e-12`, which only confirms that `nsolve` converged to *some* solution of the equation it was handed — it is not a test that the canonical point matches the paper-stated value.

**Why this matters:**
If `r_F1` (the upstream-carried input) drifted, or if `gPi`/`Sq` were slightly mis-encoded, the script would still PASS while reporting a wrong `Pi_*`, wrong `Sigma_0`, wrong `That`. The script silently rebrands whatever it finds as "the canonical point". The paper card's `Checks` line (iii) — "numerical fixed points are recorded as numerically located, not closed-form constants" — is honored in form (they are numerical) but the script does not actually pin them.

**Required change:**
Add four `expect_close`/`expectApprox` style numeric assertions after the nsolve block:
- `g_-^{F1}` value against `0.7580350789446628269196808904` to ≥ 25-digit tolerance.
- `Pi_*` against `1.5088295134931555` to ≥ 12-digit tolerance.
- `S_q(Pi_*)` against `0.6580759376054293` to ≥ 12-digit tolerance.
- `Sigma_0(Pi_*)` against `1.8059411109563538` to ≥ 12-digit tolerance.
- `That(Pi_*)` against `0.9014840541742040` to ≥ 12-digit tolerance.

Mirror in Mathematica with `expectApprox`. Tolerances should be loose enough to absorb engine numerical jitter (the engines agree to ~25 digits already).

**Verification:**
After fix: five new numeric assertions in each script, each tied to a notes-quoted value. Output transcript should show each PASS line.

### F3 — mathematica_transliteration

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage142_selfconsistent_mouth_branch_sympy_audit.py:19-49`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage142_selfconsistent_mouth_branch_mathematica_audit.wl:41-69`

**What's wrong:**
The Mathematica script is a line-by-line port of the SymPy script's algebra, not an independent derivation. Corresponding sections:

SymPy 24-29:
```
r = sp.sqrt(sp.Integer(4107) - 100*pi**2)/(10*pi)
gPi = 2*Pi*(2*Pi*sp.exp(Pi)+pi)/((4*Pi**2+pi**2)*(sp.exp(Pi)-1))
Sq = Pi*(((pi/2)*sp.tanh(pi/2)) + Pi*(sp.exp(-Pi)*sp.sech(pi/2)-1))/((1-sp.exp(-Pi))*(((pi/2)**2)-Pi**2))
Rq = sp.simplify((gPi-r)**2/(1+r**2))
Sigma0 = Pi/(1-Rq*Sq)
That = sp.sqrt(sp.Rational(9,20)*Sigma0)
```

Mathematica 44-49:
```
r = Sqrt[4107 - 100*Pi^2]/(10*Pi);
gPi = 2*piM*(2*piM*Exp[piM] + Pi)/((4*piM^2 + Pi^2)*(Exp[piM] - 1));
sQ = piM*(((Pi/2)*Tanh[Pi/2]) + piM*(Exp[-piM]*Sech[Pi/2] - 1))/((1 - Exp[-piM])*((Pi/2)^2 - piM^2));
rQ = (gPi - r)^2/(1 + r^2);
sigma0 = piM/(1 - rQ*sQ);
tHat = Sqrt[(9/20)*sigma0];
```

Same constant `4107 - 100*pi^2`, same expression tree for `gPi`, same expression tree for `Sq`, same definitions, same single assertion `R_q(g_minus) - 1/4 == 0`. The Mathematica script does not derive `r_F1`, `gPi`, or `Sq` from any independent physical premise — it transcribes the symbolic forms.

**Why this matters:**
Second-engine policy requires both engines to derive the result independently so that a bug in symbolic encoding (e.g., a sign flip in `Sq`'s numerator) cannot pass both engines. Here a typo carried between the two `.py` and `.wl` files would propagate identically and both would PASS.

**Required change:**
The Mathematica script should derive at least one of `gPi`, `Sq`, or `r_F1` from a different starting form (e.g., series expansion + closed-form recognition for `gPi`, or alternate integral form for `Sq`) and then `FullSimplify` to the form used in the algebra. Alternatively, introduce an independent numerical cross-check: e.g., construct `gPi` from a finite-Pi integral / series and verify it agrees with the closed form to 25 digits at several `Pi` values.

A lighter-weight option that still satisfies independence: in Mathematica, define `gPi` and `Sq` symbolically as in SymPy, but additionally introduce *one* independent numerical verification (e.g., compare `Sq` at `Pi = 1.5` against a direct numerical integral of whatever Stage 242 defined `S` to be — if Stage 242's integral form is not available here, the closed-form re-derivation route above is the only option).

**Verification:**
After fix: Mathematica script contains at least one block that re-derives one of `gPi`, `Sq`, or `r_F1` from a different starting form OR runs an independent numerical sanity check against the closed form. The line count of the Mathematica file should grow; the SymPy file does not need to change.

### F4 — hardcoded_result

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage142_selfconsistent_mouth_branch_sympy_audit.py:24,26`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage142_selfconsistent_mouth_branch_mathematica_audit.wl:44,46`

**What's wrong:**
Two closed-form expressions are inserted without in-stage derivation OR a comment naming the upstream stage they came from:
- `r = sqrt(4107 - 100*pi^2)/(10*pi)` — this is `r_F1`, marked in notes as a Family-1 reduced-mixed-core input "imported" from upstream stages, but the script gives no provenance comment.
- `Sq = Pi*(((pi/2)*tanh(pi/2)) + Pi*(exp(-Pi)*sech(pi/2)-1))/((1-exp(-Pi))*(((pi/2)**2)-Pi**2))` — the notes describe `S_q(Pi) = S(Pi, pi/2)` but never spell out this closed form; the closed form is presumed to come from Stage 242. No comment.

**Why this matters:**
A future reader (or auditor) has no way to confirm these closed forms without leaving the stage. If the upstream stage's `S_q` or `r_F1` is ever re-derived to a corrected form, this stage will silently disagree.

**Required change:**
Add inline comments above lines 24 and 26 (sympy) and lines 44 and 46 (mathematica) naming the upstream stage and the exact symbol carried forward, e.g.:
```
# r_F1: imported reduced Family-1 ratio. Carried forward from <upstream-stage>; see notes/stages/<file>.md
# S_q(Pi) closed form: carried forward from Stage 242 self-matched susceptibility closure.
```
This finding is `informational + paper card alignment` — the paper card's `\stagefield{Inputs}` mentions "shell/mixed core, the mouth source law, outlet consistency, core-to-mouth gain maps, and self-matched susceptibility closure" so the upstream provenance exists; what is missing is the script-side label tying these literals to that text.

**Verification:**
After fix: two new comment lines in each script. No assertion change; no output change.

### F5 — paper_misalignment (subtype: target_mismatch)

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage142_selfconsistent_mouth_branch_sympy_audit.py:22,62`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage142_selfconsistent_mouth_branch_mathematica_audit.wl:39,80`
- output transcripts repeat the same banner text

**What's wrong:**
Both scripts banner themselves as **"STAGE 125 — SELF-CONSISTENT MOUTH-BRANCH LAW"** (sympy line 22, mathematica line 39) and `"STAGE 125 LEDGER"` (sympy line 62, mathematica line 80). The paper card is `paper/stages/stage_142.tex` (`Stage 142: Self-Consistent Mouth-Branch Law`). Only the very last printout in Mathematica (`"Stage 142 Mathematica audit passed."`, line 92) uses the correct number. The output `.txt` transcripts inherit the wrong banner.

Paper quote (stage_142.tex:1):
> `\section[Stage 142]{Stage 142: Self-Consistent Mouth-Branch Law}`

Script quote (sympy line 22):
> `banner("STAGE 125 — SELF-CONSISTENT MOUTH-BRANCH LAW")`

**Why this matters:**
The math is unaffected, but the verification transcript carries the wrong stage label, which makes audit-by-grep brittle and could cause downstream pipelines that key off the transcript banner to attribute the result to a different stage. This is a `target_mismatch` in *labeling*, not in identity verified.

**Required change:**
Replace `"STAGE 125"` with `"STAGE 142"` in:
- sympy script line 22 banner
- sympy script line 62 banner
- mathematica script line 39 banner
- mathematica script line 80 banner

After fix, regenerate the output transcripts (which the verifier will do as a matter of course when rerunning).

**Verification:**
After fix: `grep "STAGE 125" scripts/moving_throat_pde_stage142_*` and `grep "STAGE 125" mathematica/moving_throat_pde_stage142_*` both return no matches; output transcripts say "STAGE 142" in banners.

## Independent-derivation check (Mathematica)

The Mathematica script is structurally a transliteration of the SymPy script, as documented in F3. Both scripts:
1. Hardcode `r = sqrt(4107 - 100*pi^2)/(10*pi)` in identical form.
2. Encode `gPi`, `Sq` as identical syntactic expressions.
3. Use the same `(gPi-r)^2/(1+r^2)` form for `Rq`.
4. Perform the same single symbolic assertion `R_q(g_minus) - 1/4 == 0`.
5. Solve `gPi == gminus` numerically with the same starting guess `1.5`.

The Mathematica script has no independent path. Its presence verifies that two engines can agree on transcription, not that two engines can independently confirm the physics.

## Engine cross-check

Both engines agree where they overlap. Selected residuals:

| Quantity | SymPy (30 digits) | Mathematica (30 digits) |
|---|---|---|
| `g_-^{F1}` | 0.758035078944662826919680890414 | 0.758035078944662826919680890414... |
| `Pi_*` | 1.50882951349315552747043511772 | 1.50882951349315558300555075595... |
| `R_q(Pi_*)` | 0.250000000000000001945042992303 | 0.250000000000000000000000000000... |
| `S_q(Pi_*)` | 0.658075937605429271930315313436 | 0.658075937605429274616601849160... |
| `Sigma_0(Pi_*)` | 1.80594111095635380721796724713 | 1.80594111095635387237367290920... |
| `That(Pi_*)` | 0.901484054174204022702401688674 | 0.901484054174204038964512711141... |

Both agree to ~14 decimal places on `Pi_*` (limited by SymPy's nsolve precision configured at 30 digits, which the script does not push tighter); Mathematica with `WorkingPrecision -> 80, AccuracyGoal -> 30, PrecisionGoal -> 30` gets ~25 digits. No `engine_disagreement` finding — the two engines agree to within stated precision.

## Verdict justification

The stage card's main identity is *defined* into existence by the script (`Sigma0 = Pi/(1 - Rq*Sq)`, then `That = sqrt((9/20)Sigma0)`) rather than verified. The single symbolic assertion (`R_q(g_minus)=1/4`) is algebraically identically true for any `r`, so it does not exercise the specific `r_F1`. The numerical canonical-point values quoted in the notes are never anchored. The Mathematica script is a structural mirror of the SymPy script. The stage banners use the wrong stage number (125). Engines agree where they compute, but they're computing the same hand-encoded expressions, not independently verifying them.

Verdict: **findings**. Five findings: one tautological_check (high), one insufficient_verification (high), one mathematica_transliteration (medium), one hardcoded_result (medium), one paper_misalignment / target_mismatch (low). `stop_cold: null` — F1 and F2 fix into the same script with new assertions; F3 fix adds Mathematica-side independence; F4 is informational; F5 is a label fix. None of these would propagate destructively to downstream units before being resolved, because the *math* used by Stages 146-153 is the same closed-form `Sigma_0(\Pi)`/`That(\Pi)` formulas the notes already box — those formulas are physically correct as written; the script just doesn't verify them, it asserts them.

## Self-test notes

- For F1, mentally substituted `g_minus = r - (1/2)sqrt(1+r^2)` into `R_q = (g-r)^2/(1+r^2)` and confirmed the residual reduces to `1/4` identically, regardless of `r`'s value — so the existing check is tautological.
- For F2, confirmed the script's `Pi_*` is found by `nsolve(gPi == gminus, 1.5)`; the only assertion against it is `g(Pi_*) == g_minus` which is solver-convergence, not value-anchor. New numeric assertions compare against notes-quoted decimal strings; tolerances are looser than engine residuals so the checks won't be too tight.
- For F3, the Mathematica file uses `piM` where SymPy uses `Pi`, but every other identifier, expression, and assertion is a direct one-to-one transcription — confirmed by reading both files end-to-end.
- For F5, confirmed the file path uses `stage142` and the paper section name is "Stage 142" — the banner literal `"STAGE 125"` is unambiguously wrong.

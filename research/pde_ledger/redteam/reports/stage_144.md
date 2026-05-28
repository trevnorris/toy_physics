---
unit_id: 144
batch: IV.5
auditor_model: claude-opus-4-7-1m
audit_date: 2026-05-27T00:00:00Z
verdict: findings
stop_cold: null
findings_count: 4
paper_alignment: partial
scripts_checked:
  sympy: insufficient
  mathematica: insufficient
  engines_agree: true
  outputs_fresh: true
docs_read:
  paper_stage_tex: present
  notes_stage_files:
    - /var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage144_unique_regular_canonical_branch.md
  paper_appendix: present
---

# Audit unit 144 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_144.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage144_unique_regular_canonical_branch.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part04.tex` (only the `\input{stages/stage_144}` line at L1322; no narrative row in this appendix)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage144_unique_regular_canonical_branch_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage144_unique_regular_canonical_branch_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage144_unique_regular_canonical_branch_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage144_unique_regular_canonical_branch_mathematica_audit.txt`

## What the paper claims

The paper card (stage_144.tex) is a "coupled mouth fixed point and gain selection ledger step" and quotes the bottom-line claim: "Lower compensated branch is the unique finite-bias, finite-traction branch in the explicit positive exponential closure." The paper's `\stagefield{Checks}` block requires three explicit verifications: (i) the gain pair `(M_s, M_q)` against outlet consistency, (ii) the self-matched susceptibility closure before using the one-scalar branch law, and (iii) that the numerical fixed points are recorded as numerically located (not closed-form). The notes elaborate four concrete deliverables: (a) the upper branch `g_+^{F1} ≈ 2.79795 > 1` is impossible since `0 ≤ g_c ≤ 1`; (b) bracketing `2/π ≈ 0.6366 < g_-^{F1} ≈ 0.7580 < 1` so a unique `Π_* ≈ 1.50882951349316` exists with `g_{Π_*} = g_-^{F1}`; (c) at that point `Σ_0(Π_*) ≈ 1.80594111095636` and `T̂_m(Π_*) ≈ 0.901484054174205` are finite/moderate; (d) the comparison point `Π_match ≈ 1.90848600654854` (where `g = π/4`) and `T̂_m(Π_match) ≈ 1.01132972803599` lies past `Π_*`, establishing the canonical branch is reached before the self-matched derivative point.

## What the script claims to verify

Both engines compute the symbolic constants `r = √(4107 − 100π²)/(10π)`, the two compensated branches `g_∓^{F1} = r ∓ (1/2)√(1+r²)`, the bias function `g_Π = 2Π(2Π·e^Π + π)/((4Π² + π²)(e^Π − 1))`, an auxiliary susceptibility kernel `s_Q`, the residual `r_Q = (g_Π − r)²/(1+r²)`, the susceptibility `Σ_0 = Π/(1 − r_Q·s_Q)`, and the traction `T̂ = √((9/20)·Σ_0)`. They then numerically locate `Π_*` from `g_Π = g_-^{F1}` (seed 1.5) and `Π_match` from `g_Π = π/4` (seed 1.9), and print numerical values. The only enforced assertion is the ordering check: `Π_* > 0` AND `Π_match > Π_*`. No assertion checks the printed numerical values against the paper/notes targets; no assertion checks `g_+^{F1} > 1`; no assertion checks the bracketing `2/π < g_-^{F1} < 1`; the gain pair `(M_s, M_q)` and the self-matched susceptibility closure (paper-card checks i and ii) are not exercised at all.

## Paper ↔ script cross-check

| Paper-side deliverable | Script-side check | Status |
|---|---|---|
| Upper branch impossible: `g_+^{F1} > 1` | computes `g_+^F1` and prints it; no assertion that `> 1` | partial |
| Bracket: `2/π < g_-^{F1} < 1` | computes and prints all three values; no assertion | partial |
| Unique `Π_*` ≈ 1.50882951349316 such that `g_{Π_*} = g_-^{F1}` | `nsolve`/`FindRoot` from seed 1.5; no value-target assertion | partial |
| `Σ_0(Π_*) ≈ 1.80594111095636` (notes §3) | `Σ_0` symbol defined but never printed or asserted | missing |
| `T̂_m(Π_*) ≈ 0.901484054174205` | prints `That(Pi_*)`; no value-target assertion | partial |
| `Π_match ≈ 1.90848600654854` and `T̂_m(Π_match) ≈ 1.01132972803599` | prints values; ordering check `Π_match > Π_*` is the only assert | partial |
| Paper-check (i): gain pair `(M_s, M_q)` vs outlet consistency | not present | missing |
| Paper-check (ii): self-matched susceptibility closure | not present (script computes `Σ_0` but never tests closure) | missing |
| Paper-check (iii): fixed points recorded numerically | yes — values printed at full 30-digit precision | match |

Dominant pattern is `partial`/`missing`: `paper_alignment: partial`.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 43-44 | `if not (float(Pi_star) > 0 and float(Pi_match) > float(Pi_star)): raise AssertionError(...)` | partial — only the ordering aspect of notes §3 | partial |
| A2 | mathematica | 57 | `If[!(N[piStar,30] > 0 && N[piMatch,30] > N[piStar,30]), fail[...]]` | same as A1 | partial |

Both engines have exactly one substantive guard. Every other "check" is a `Print[...]` of a derived number with no assertion attached. Per the v2 rubric, prints are not checks: a script whose only assertion is "`Π_match > Π_* > 0`" cannot fail unless the bracketing is grossly wrong, and in particular it cannot detect a wrong numerical value of `Π_*` or `T̂_m` (e.g., a unit-conversion error in `g_Π` or a sign flip in `s_Q`).

## Findings

### F1 — paper_misalignment

**Subtype:** target_mismatch (banner / docstring labelling)
**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage144_unique_regular_canonical_branch_sympy_audit.py:16,46`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage144_unique_regular_canonical_branch_mathematica_audit.wl:27,59`

**What's wrong:**
Both scripts (and both saved output transcripts) print `STAGE 127 — UNIQUE REGULAR CANONICAL MOUTH BRANCH` and `STAGE 127 LEDGER` for a stage that the filename, paper card, and notes call stage 144. The paper card label is `\section[Stage 144]{Stage 144: Unique Regular Canonical Mouth Branch}` and the notes title is `## Moving-Throat PDE — Stage 144: Unique Regular Canonical Mouth Branch`. The mismatched "127" inside both engines is identical, indicating a copy-paste source from a different unit's scaffolding.

**Why this matters:**
Transcripts go into the audit ledger; a script whose own banner claims it is verifying a different stage is a provenance defect. It will also surface as a discrepancy in any downstream tooling that greps transcripts by stage number.

**Required change:**
Replace the two banner strings in each script to read `STAGE 144 — UNIQUE REGULAR CANONICAL MOUTH BRANCH` and `STAGE 144 LEDGER`, then refresh the saved output by re-running the engines (verifier action, not a code edit).

**Verification:**
Both transcripts should print `STAGE 144` headers; no `127` substring remains in either script or saved output.

### F2 — insufficient_verification

**Severity:** high
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage144_unique_regular_canonical_branch_sympy_audit.py:43-44`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage144_unique_regular_canonical_branch_mathematica_audit.wl:57`

**What's wrong:**
The only enforced assertion in either engine is the ordering `Π_* > 0 ∧ Π_match > Π_*`. The paper card claims the stage proves "Lower compensated branch is the unique finite-bias, finite-traction branch in the explicit positive exponential closure." The notes spell out four concrete numerical deliverables: `g_+^{F1} ≈ 2.798 > 1` (upper branch impossible), `2/π < g_-^{F1} < 1` (lower branch bracketed by the `g_Π` range), `Π_* ≈ 1.50882951349316`, `Σ_0(Π_*) ≈ 1.80594111095636`, `T̂_m(Π_*) ≈ 0.901484054174205`, `Π_match ≈ 1.90848600654854`, and `T̂_m(Π_match) ≈ 1.01132972803599`. None of these targets is asserted. In particular:
- `g_+^{F1} > 1` is never asserted (the entire "upper branch impossible" leg of the argument).
- `2/π < g_-^{F1} < 1` is never asserted (the bracketing that gives uniqueness).
- `Π_* − 1.50882951349316`, `T̂_m(Π_*) − 0.901484054174205`, etc., are never bounded.
- `Σ_0(Π_*) − 1.80594111095636` is never asserted; in fact `Σ_0` is computed but never even printed.

A typo in `g_Π` that left the qualitative shape intact (e.g., a numeric factor error somewhere in `s_Q` or `r_Q`) could shift `Π_*` and `T̂_m` away from the notes values, but the ordering `Π_match > Π_* > 0` would still hold and the assertion would still pass. This is the textbook "prints look right, assert is vacuous" pattern v2 is designed to catch.

**Why this matters:**
The paper card downstream-uses this stage in Stages 146-153 and the finite mouth-profile correction. If the numerical fixed point is silently off, those downstream stages will inherit a wrong constant. The stage card's own `\stagefield{Checks}` block lists three concrete checks; the script implements at most one (numerical recording).

**Required change:**
Add closed-tolerance numerical assertions tying the computed quantities to the notes-anchored targets. Specifically add the following non-tautological checks to both engines:

1. `assert N(g_+^{F1}) > 1` (upper branch impossible) — and `assert N(g_-^{F1}) < 1 and N(g_-^{F1}) > 2/π` (bracketing).
2. `assert abs(N(Π_*) - 1.50882951349316) < 1e-12` (canonical-point value).
3. `assert abs(N(Σ_0(Π_*)) - 1.80594111095636) < 1e-12` (susceptibility at canonical point — also requires adding a print of `Σ_0(Π_*)` so the transcript records it).
4. `assert abs(N(T̂_m(Π_*)) - 0.901484054174205) < 1e-12` (canonical traction).
5. `assert abs(N(Π_match) - 1.90848600654854) < 1e-12` and `assert abs(N(T̂_m(Π_match)) - 1.01132972803599) < 1e-12` (matched-derivative point).

All targets come directly from the notes file under `## 2. Lower branch is uniquely reachable at finite bias` and `## 3. Finite regular traction at the canonical point`. The tolerances `1e-12` are conservative against the 30-digit precision both engines already use; the assertions will fail if any algebraic constant in `r`, `g_Π`, `s_Q`, `r_Q`, `Σ_0`, or `T̂` is mistyped.

**Verification:**
After the fix, both transcripts should print explicit `assert`/`pass` lines for items 1-5 above, and the saved outputs should include a printed `Sigma_0(Pi_*)` line. Exit codes remain 0 only if all five numerical guards hold.

### F3 — insufficient_verification

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage144_unique_regular_canonical_branch_sympy_audit.py:1-56` (whole script)
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage144_unique_regular_canonical_branch_mathematica_audit.wl:1-73` (whole script)

**What's wrong:**
The paper card's `\stagefield{Checks}` block enumerates three items. Item (i) is "Check the gain pair `(M_s, M_q)` against outlet consistency." Item (ii) is "Check the self-matched susceptibility closure before using the one-scalar branch law." Neither engine references `M_s`, `M_q`, an outlet consistency relation, or a closure identity for the self-matched susceptibility. The script computes `Σ_0 = Π/(1 − r_Q s_Q)` but never tests an independent closure (e.g., that `Σ_0` satisfies some self-consistency equation that the construction is meant to enforce). The script also never names `M_s` or `M_q`, so paper-check (i) has zero engine coverage.

**Why this matters:**
Two of the three explicit paper-card checks have no script counterpart. The "Inputs" line of the card states the stage "imports the shell/mixed core, the mouth source law, outlet consistency, core-to-mouth gain maps, and self-matched susceptibility closure" — these are the inputs that ought to be exercised. If the upstream gain maps or susceptibility closure are off, this stage's script cannot detect it.

**Why this is `insufficient_verification` rather than `paper_misalignment`:**
The script's own docstring/comments don't make the claim that `(M_s, M_q)` and the self-matched susceptibility closure are tested, so this isn't "script claims X, paper claims Y" disagreement. It's that the script under-covers the paper's claim list. The relevant subtype could equally be argued as `paper_misalignment / script_missing_paper_claim`; we record it here because Codex can mechanically add the checks once the upstream symbol names are inlined into this script. However, if the user prefers to treat it as paper-misalignment (because `M_s, M_q` are defined in upstream stages whose results should be carried forward by reference rather than re-tested here), see the directive's `## Resolve before fix_loop` note for F3.

**Required change:**
The minimal "carry-forward" fix is to add into both engines (a) symbols/values for `M_s, M_q` (carried forward from the upstream gain-map stage — the script's introductory comment should cite which stage) and a numeric assertion that the outlet-consistency identity that defines them is satisfied at `Π_*`, and (b) the self-matched susceptibility closure relation (typically of the form `Σ_0(Π) − Π/(1 − r_Q s_Q) == 0`, but the actual closure identity must come from the upstream notes; the present script *defines* `Σ_0` by that formula, so the equation is tautological — the closure check needs a second, independent expression for `Σ_0` and asserts equality). Because the second expression is not derivable from material inside this stage alone, this finding's directive routes the resolution decision to the user (see `## Resolve before fix_loop`).

**Verification:**
After resolution, both transcripts should print explicit `M_s = …, M_q = …, outlet-consistency residual = 0` and `Σ_0 (independent form) = …, susceptibility-closure residual = 0` lines, each guarded by an assertion.

### F4 — mathematica_transliteration

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage144_unique_regular_canonical_branch_mathematica_audit.wl:32-50`
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage144_unique_regular_canonical_branch_sympy_audit.py:18-36`

**What's wrong:**
The `.wl` is a line-by-line transliteration of the `.py`. Compare:

SymPy (lines 18-25):
```
r = sp.sqrt(sp.Integer(4107) - 100*pi**2)/(10*pi)
gminus = sp.simplify(r - sp.Rational(1,2)*sp.sqrt(1+r**2))
gplus  = sp.simplify(r + sp.Rational(1,2)*sp.sqrt(1+r**2))
gPi = 2*Pi*(2*Pi*sp.exp(Pi)+pi)/((4*Pi**2+pi**2)*(sp.exp(Pi)-1))
Sq = Pi*(((pi/2)*sp.tanh(pi/2)) + Pi*(sp.exp(-Pi)*sp.sech(pi/2)-1))/((1-sp.exp(-Pi))*(((pi/2)**2)-Pi**2))
Rq = (gPi-r)**2/(1+r**2)
Sigma0 = Pi/(1-Rq*Sq)
That = sp.sqrt(sp.Rational(9,20)*Sigma0)
```

Mathematica (lines 32-39):
```
r = Sqrt[4107 - 100*Pi^2]/(10*Pi);
gMinus = FullSimplify[r - Sqrt[1 + r^2]/2, Assumptions -> $Assumptions];
gPlus = FullSimplify[r + Sqrt[1 + r^2]/2, Assumptions -> $Assumptions];
gPi = 2*piM*(2*piM*Exp[piM] + Pi)/((4*piM^2 + Pi^2)*(Exp[piM] - 1));
sQ = piM*(((Pi/2)*Tanh[Pi/2]) + piM*(Exp[-piM]*Sech[Pi/2] - 1))/((1 - Exp[-piM])*((Pi/2)^2 - piM^2));
rQ = (gPi - r)^2/(1 + r^2);
sigma0 = piM/(1 - rQ*sQ);
tHat = Sqrt[(9/20)*sigma0];
```

The seeds for root-finding (1.5 and 1.9), the symbolic forms, the ordering of derivations, the printed quantities, and even the ledger prose are identical. There is no independent re-derivation in Mathematica from physical premises (e.g., neither engine derives `g_Π` from a mouth-layer ansatz; both just paste the closed form). Per the project's second-engine policy, this defeats the purpose of running two engines: a transliteration error in `s_Q` or `r_Q` would be present in both, and both engines would "agree" on the wrong answer.

**Why this matters:**
The "engines agree" signal that the verifier reports becomes meaningless when both engines execute the same algebraic recipe. The independent-derivation policy exists so that engine-disagreement can flag mistakes; transliteration neutralises that signal.

**Required change:**
This is the same systemic pattern called out elsewhere in the project. The minimal mitigation that does not require re-deriving anything is to have the Mathematica script obtain `g_Π`, `s_Q`, and `Σ_0` by a route that is *not* a transcription of the SymPy definitions — for example, by writing `g_Π` via an integral form that Mathematica evaluates symbolically rather than as the closed-form ratio, and computing `Σ_0` via a self-consistent fixed point of the closure equation rather than as the explicit ratio `Π/(1 − r_Q s_Q)`. The auditor cannot prescribe the precise alternative form because the upstream derivation of `g_Π` and `s_Q` is not in scope for this stage; the directive therefore routes this to the user under `## Resolve before fix_loop` rather than to Codex.

**Verification:**
After resolution, the Mathematica script should not contain a string match for the SymPy definitions of `g_Π`, `s_Q`, `r_Q`, `Σ_0`, `T̂` (allowing for syntactic differences); the two engines should still agree numerically at `Π_*` and `Π_match` to 12+ digits.

## Independent-derivation check (Mathematica)

The Mathematica script is a transliteration of the SymPy script. Same variable choreography (`r → gMinus/gPlus → gPi → sQ → rQ → sigma0 → tHat`), identical algebraic forms for every defined quantity, identical numerical seeds for `FindRoot` (1.5 and 1.9) matching the SymPy `nsolve` seeds, identical printed lines, identical ledger prose. See F4 for the side-by-side quotation. This is `mathematica_transliteration`.

## Engine cross-check

Both engines produce 15+ digits of agreement on every printed numerical quantity:

| Quantity | SymPy | Mathematica |
|---|---|---|
| `g_-^{F1}` | 0.758035078944662826919680890414 | 0.758035078944662826919680890414... |
| `g_+^{F1}` | 2.79795199200529341011158893417 | 2.79795199200529341011158893417... |
| `2/π` | 0.636619772367581343075535053490 | 0.636619772367581343075535053490... |
| `Π_*` | 1.50882951349315552747043511772 | 1.5088295134931555830055507559... (diverges at the 17th digit) |
| `T̂(Π_*)` | 0.901484054174204022702401688674 | 0.9014840541742040389645127111... (diverges at the 17th digit) |
| `Π_match` | 1.90848600654854538838378630317 | 1.9084860065485453947610306052... (diverges at the 17th digit) |
| `T̂(Π_match)` | 1.01132972803599475860454058210 | 1.0113297280359947602555983007... (diverges at the 17th digit) |

The ~16-digit agreement is the expected behaviour of SymPy's `nsolve` (defaults to ~15-digit precision) against Mathematica's 30-digit `FindRoot`. No engine_disagreement. (The agreement is, however, partly a side-effect of F4 — both engines execute the same recipe.)

## Verdict justification

`findings`. The numerical results that the script prints match the notes targets to displayed precision and the two engines agree, so nothing is mathematically wrong on the printed surface. But the script-side verification is thin: a single ordering assertion, no numerical-target assertions, no exercise of paper-card checks (i) and (ii), wrong stage banner, and the Mathematica script is a transliteration of the SymPy. The stage's downstream consumers (146-153) inherit `Π_*, T̂_m(Π_*)` as constants, so silent numeric drift would propagate — `insufficient_verification` (F2) is the high-severity finding. F1 (banner mislabelling) is cosmetic but trivial to fix. F3 (missing `(M_s, M_q)` and self-matched susceptibility closure checks) and F4 (transliteration) need user direction on scope; the directive lists `## Resolve before fix_loop` blocks for those. No `stop_cold`: the math holds up; the verification scaffolding is the weak point. Paper card and notes were read before the scripts; the algebraic constants in the scripts (e.g., `4107 − 100π²`, `(9/20)`, the bracket `(2/π, 1)`) match what the notes use.

## Self-test notes

I confirmed: (1) the trivial-case behaviour of `g_Π` — at `Π → 0⁺` it limits to `2/π` and at `Π → ∞` it limits to `1`, matching the notes' range claim; (2) the bracketing inequality `2/π < g_-^{F1} < 1` is numerically true at the printed values (0.6366 < 0.7580 < 1); (3) the proposed F2 assertions are non-tautological — they would each catch a sign flip, a missing factor of 2, or a wrong rational coefficient in `r`, `s_Q`, `r_Q`, or `T̂`; (4) the proposed `1e-12` tolerances are loose against the 30-digit precision both engines already use, so legitimate calculation will pass cleanly. No derivative-of-constant traps in this stage's proposed checks because the new assertions are scalar value-comparisons, not derivative-based.

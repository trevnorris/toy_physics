---
unit_id: 072
batch: III.3
auditor_model: claude-opus-4-7-1m
audit_date: 2026-05-26T00:00:00Z
verdict: findings
stop_cold: null
findings_count: 1
paper_alignment: aligned
scripts_checked:
  sympy: present
  mathematica: present
  engines_agree: true
  outputs_fresh: true
docs_read:
  paper_stage_tex: present
  notes_stage_files:
    - /var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage072_explicit_branch_thresholds.md
  paper_appendix: present
---

# Audit unit 072 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_072.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage072_explicit_branch_thresholds.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part03.tex` (row at line 122)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage072_explicit_branch_thresholds_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage072_explicit_branch_thresholds_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage072_explicit_branch_thresholds_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage072_explicit_branch_thresholds_mathematica_audit.txt`

## What the paper claims

The paper card `\stagefield{Output}{A threshold-surface map in the canonical branch variables.}` quotes the bottom-line identity
`Upsilon_w * Lambda_ell^2  <=>  Pe_req / Delta_{inf,0}(kappa(chi_s, Lambda_ell), eta=Lambda_ell)`
with `kappa = 4*chi_s^2 + (4/5)*Lambda_ell^2`, `eta = Lambda_ell`. The notes flesh this out: (i) explicit threshold surfaces `Upsilon_fail = Pe_req/(Lambda_ell^2 Delta_inf)` and `Upsilon_suff = Pe_req/(Lambda_ell^2 Delta_0)`; (ii) physical wall-amplitude thresholds `V0_fail^2 = hbar^2 c_sw^2 Upsilon_fail / (4 rho_w^2)` and `V0_suff^2 = hbar^2 c_sw^2 Upsilon_suff / (4 rho_w^2)`; (iii) shell-gradient regime `(4/5)*Lambda_ell^2 >> 4*chi_s^2` giving `Upsilon_fail ~ 2 Pe_req/(sqrt(5) Lambda_ell)` and `Upsilon_suff ~ (4/5)(1 + 2/sqrt(5)) Pe_req`; (iv) compression regime `4*chi_s^2 >> (4/5)*Lambda_ell^2` giving `Upsilon_fail ~ 2 Pe_req chi_s/Lambda_ell^2` and `Upsilon_suff ~ 4 Pe_req chi_s^2 (Lambda_ell+2 chi_s)/Lambda_ell^3`. The appendix row at part03.tex:122 summarizes this as "Fail/succeed surfaces in (chi_s, Lambda_ell, Upsilon_w)" with status ExactClosure.

## What the script claims to verify

Both engines define `kappa = 4*chi_s^2 + (4/5)*Lambda_ell^2`, `eta = Lambda_ell`, `alpha = sqrt(kappa)`, the closed forms `Delta_0 = eta*(cosh(alpha)-1)/(alpha^2*(alpha*sinh(alpha)+eta*cosh(alpha)))` and `Delta_inf = (cosh(alpha)+(eta/alpha)*sinh(alpha)-1)/(alpha*sinh(alpha)+eta*cosh(alpha))`. They compute `Upsilon_fail = Pe_req/(Lambda_ell^2 * Delta_inf)`, `Upsilon_suff = Pe_req/(Lambda_ell^2 * Delta_0)`, and the physical `V0_fail^2`, `V0_suff^2`. They then construct hand-built shell forms (`Delta0_shell`, `DeltaInf_shell`) by substituting `alpha = c*Lambda_ell` with `c = 2/sqrt(5)`, and hand-built compression forms (`Delta0_comp`, `DeltaInf_comp`) by substituting `alpha = 2*chi_s`. The substantive assertions are: (a) `Delta_0.subs(chi_s, 0)/Delta0_shell -> 1` as `Lambda_ell -> oo` and similarly for `DeltaInf`; (b) `Pe_req/(Lambda_ell^2 * DeltaInf_shell) = 2*Pe_req/(sqrt(5)*Lambda_ell)` and the suff analog; (c) `Delta_0/Delta0_comp -> 1` as `chi_s -> oo` and similarly for `DeltaInf`; (d) `Pe_req/(Lambda_ell^2 * DeltaInf_comp) = 2*Pe_req*chi_s/Lambda_ell^2` and the suff analog.

## Paper <-> script cross-check

| Paper deliverable | Script-side coverage | Result |
|---|---|---|
| `Upsilon_fail = Pe_req/(Lambda^2 Delta_inf)` | sympy line 42; mathematica line 46 | match |
| `Upsilon_suff = Pe_req/(Lambda^2 Delta_0)` | sympy line 43; mathematica line 47 | match |
| `V0_fail^2 = hbar^2 c_sw^2 Upsilon_fail/(4 rho_w^2)` | sympy line 52; mathematica line 56 | match |
| `V0_suff^2 = hbar^2 c_sw^2 Upsilon_suff/(4 rho_w^2)` | sympy line 53; mathematica line 57 | match |
| Shell `Upsilon_fail ~ 2 Pe_req/(sqrt(5) Lambda)` | sympy 62, 79-80; mathematica 64, 85 | match |
| Shell `Upsilon_suff ~ (4/5)(1 + 2/sqrt(5)) Pe_req` | sympy 63, 81-82; mathematica 65, 86 | match |
| Shell-regime leading-order anchor for Delta_0, Delta_inf | sympy 69-78; mathematica 72-81 | match |
| Compression `Upsilon_fail ~ 2 Pe_req chi/Lambda^2` | sympy 106; mathematica 110 | match |
| Compression `Upsilon_suff ~ 4 Pe_req chi^2 (Lambda+2 chi)/Lambda^3` | sympy 107-108; mathematica 111-113 | match |
| Compression-regime leading-order anchor for Delta_0, Delta_inf | sympy 93-100; mathematica 97-106 | match |

`paper_alignment: aligned` — every paper-side deliverable maps to a non-tautological script-side assertion, and the script tests nothing extraneous to the paper.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 75-76 | `expect_zero(Delta0_shell_ratio - 1)` (limit Lambda_ell -> oo of full Delta_0 / Delta0_shell at chi_s=0) | shell-regime leading order anchor | yes |
| A2 | sympy | 77-78 | `expect_zero(DeltaInf_shell_ratio - 1)` | shell-regime leading order anchor | yes |
| A3 | sympy | 79-80 | `expect_zero(Pe_req/(L^2 * DeltaInf_shell) - Upsilon_fail_shell)` | shell `Upsilon_fail` value | yes |
| A4 | sympy | 81-82 | `expect_zero(Pe_req/(L^2 * Delta0_shell) - Upsilon_suff_shell)` | shell `Upsilon_suff` value | yes |
| A5 | sympy | 97-98 | `expect_zero(Delta0_comp_ratio - 1)` (limit chi_s -> oo of full Delta_0 / Delta0_comp) | compression leading order anchor | yes |
| A6 | sympy | 99-100 | `expect_zero(DeltaInf_comp_ratio - 1)` | compression leading order anchor | yes |
| A7 | sympy | 106 | `expect_zero(Upsilon_fail_comp - 2 Pe_req chi_s/L^2)` | compression `Upsilon_fail` value | yes |
| A8 | sympy | 107-108 | `expect_zero(Upsilon_suff_comp - 4 Pe_req chi_s^2 (L+2 chi_s)/L^3)` | compression `Upsilon_suff` value | yes |
| B1 | mathematica | 78-79 | `expectZero[delta0ShellRatio - 1]` | shell Delta_0 anchor | yes |
| B2 | mathematica | 80-81 | `expectZero[deltaInfShellRatio - 1]` | shell Delta_inf anchor | yes |
| B3 | mathematica | 85 | `expectZero[peReq/(L^2 * deltaInfShell) - upsilonFailShell]` | shell `Upsilon_fail` value | yes |
| B4 | mathematica | 86 | `expectZero[peReq/(L^2 * delta0Shell) - upsilonSuffShell]` | shell `Upsilon_suff` value | yes |
| B5 | mathematica | 103-104 | `expectZero[delta0CompRatio - 1]` | compression Delta_0 anchor | yes |
| B6 | mathematica | 105-106 | `expectZero[deltaInfCompRatio - 1]` | compression Delta_inf anchor | yes |
| B7 | mathematica | 110 | `expectZero[upsilonFailComp - 2 peReq chiS/L^2]` | compression `Upsilon_fail` value | yes |
| B8 | mathematica | 111-113 | `expectZero[upsilonSuffComp - 4 peReq chiS^2 (L+2 chiS)/L^3]` | compression `Upsilon_suff` value | yes |

Notes:
- A1/A2/A5/A6 (and the Mathematica counterparts) are the substantive anchors: they take limits of the full closed-form `Delta_0, Delta_inf` and check the ratio against the hand-built leading-order forms is 1. Without these anchors, A3/A4/A7/A8 would be algebraic identities between hand-built forms with no physics content; with these anchors, A3/A4/A7/A8 establish that the paper's quoted asymptotic values follow from the full closed forms in the asserted regimes.
- All eight assertions exit 0 in the saved transcripts (sympy `expect_zero` produces `= 0` lines; mathematica produces `PASS:` lines).
- Manually verified the algebra: `c = 2/sqrt(5)` gives `c^2 = 4/5`, `c+1 = 1 + 2/sqrt(5)`, so `DeltaInf_shell = (1 + 1/c)/((c+1)*Lambda) = (1 + sqrt(5)/2)/((1+2/sqrt(5))*Lambda) = sqrt(5)/(2*Lambda)` after rationalization; thus `Pe_req/(Lambda^2 * DeltaInf_shell) = 2 Pe_req/(sqrt(5) Lambda)`. Compression `DeltaInf_comp = (1 + Lambda/(2 chi))/(2 chi + Lambda) = 1/(2 chi)` identically (this one is unconditional algebra, but the connection to the full `Delta_inf` is established by B6/A6). Both match paper.

## Findings

### F1 — paper_misalignment

**Subtype:** notes_contradicts_script

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage072_explicit_branch_thresholds_sympy_audit.py:3, 5, 24, 110`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage072_explicit_branch_thresholds_mathematica_audit.wl:26`

**What's wrong:**
The script docstring/banner strings still reference the old "Stage 55" / "Stage 055" identifier although the filename, paper card, notes file, and appendix row all label this unit Stage 072.

SymPy script:
- Line 3 (docstring): `moving_throat_pde_stage55_explicit_branch_thresholds_sympy_audit.py`
- Line 5 (docstring): `SymPy audit for Stage 55:`
- Line 24 (banner): `banner("STAGE 55 — EXPLICIT BRANCH THRESHOLD SURFACES")`
- Line 110 (banner): `banner("STAGE 55 THEOREM LEDGER")`

Mathematica script:
- Line 26 (banner): `banner["STAGE 055 — EXPLICIT BRANCH THRESHOLD SURFACES"]`
- Line 117 already says `Stage 072 Mathematica audit passed.` — inconsistent with the line-26 banner.

Paper card:
- `/paper/stages/stage_072.tex:1` quote: `\section[Stage 072]{Stage 072: Explicit Branch Placement Map and Threshold Surfaces}`

Notes file: `notes/stages/moving_throat_pde_stage072_explicit_branch_thresholds.md` (filename "stage072").

These artifacts of an earlier stage renumbering bleed into the saved output transcripts: the saved `.txt` outputs contain `STAGE 55 — EXPLICIT BRANCH THRESHOLD SURFACES` and `STAGE 55 THEOREM LEDGER` headers (sympy output lines 3, 31) and `STAGE 055 — EXPLICIT BRANCH THRESHOLD SURFACES` (mathematica output line 3), which contradict the surrounding paper labelling.

**Why this matters:**
Cosmetic but real. A verifier diffing the saved transcripts against the paper card sees "Stage 55"/"Stage 055" printouts from a file named `stage072_…`. The physics is unaffected; only the labelling lags behind the renumbering. This is `notes_contradicts_script` per the subtype taxonomy — strictly speaking the docstring contradicts the .tex paper card, not the notes; the direction of correction is unambiguous (the file path, the paper card, the notes file, and the appendix row all agree on 072; the script's internal strings are the outdated artifact).

**Required change:**
- In SymPy script, line 3 replace `moving_throat_pde_stage55_explicit_branch_thresholds_sympy_audit.py` with `moving_throat_pde_stage072_explicit_branch_thresholds_sympy_audit.py`.
- Line 5 replace `SymPy audit for Stage 55:` with `SymPy audit for Stage 072:`.
- Line 24 replace `STAGE 55 — EXPLICIT BRANCH THRESHOLD SURFACES` with `STAGE 072 — EXPLICIT BRANCH THRESHOLD SURFACES`.
- Line 110 replace `STAGE 55 THEOREM LEDGER` with `STAGE 072 THEOREM LEDGER`.
- In Mathematica script, line 26 replace `STAGE 055 — EXPLICIT BRANCH THRESHOLD SURFACES` with `STAGE 072 — EXPLICIT BRANCH THRESHOLD SURFACES`.

This is a pure label change; no math/assertion edits are involved. The direction is unambiguous and does not require user resolution (the paper-side label has been Stage 072 for some time; the script labels are stale artifacts of a prior renumbering).

**Verification:**
After the patch, `grep -i "Stage 55\|Stage 055" scripts/moving_throat_pde_stage072_*.py mathematica/moving_throat_pde_stage072_*.wl` returns no matches. Re-running the scripts refreshes the saved `.txt` outputs so their headers no longer mention Stage 55/055. All eight `expect_zero`/`expectZero` checks still pass (they are unaffected by string changes).

## Independent-derivation check (Mathematica)

The Mathematica script's central blocks (lines 33-47 defining `kappa`, `eta`, `alpha`, `delta0`, `deltaInf`, `upsilonFail`, `upsilonSuff`) are a transliteration of the SymPy script's central block (lines 29-43). Identical variable sequence, identical right-hand sides, identical `cShell = 2/Sqrt[5]` choice. However, both engines independently execute their own native asymptotic-limit machinery to anchor the leading-order forms — SymPy uses `sp.limit(..., Lambda_ell, sp.oo)` and `sp.limit(..., chi_s, sp.oo)`, while Mathematica uses `Limit[..., lambdaEll -> Infinity]` and `Limit[..., chiS -> Infinity]`. The CAS-level execution of those limits proceeds via different simplification routines (SymPy's hyperbolic-rewrite-and-Gruntz vs Mathematica's Limit applied to `Cosh`/`Sinh`). A sign or factor error in the hand-typed `Delta_0`/`Delta_inf` closed forms would propagate identically through both engines, but the asymptotic-extraction step gives partial independence. Per the orchestrator's prior resolution on this unit, the transliteration concern is accepted as mitigated by the F1-style ratio checks (which were applied in a previous audit cycle and are now present in both scripts at sympy 69-78, mathematica 72-81 for shell and sympy 93-100, mathematica 97-106 for compression).

No new `mathematica_transliteration` finding raised in this cycle: the prior cycle's F1 was applied, the prior cycle's F2 was orchestrator-resolved as won't-fix-here (the closed forms are inherited from upstream and re-deriving them inside this unit's audit script is out of scope).

## Engine cross-check

Both outputs agree at the level claimed:

| Quantity | SymPy (.txt) | Mathematica (.txt) |
|---|---|---|
| `kappa` | `4*Lambda_ell**2/5 + 4*chi_s**2` | `(4*(5*chiS^2 + lambdaEll^2))/5` |
| `Delta_0` | hyperbolic form via `cosh(2*sqrt(5)*sqrt(L^2+5*chi^2)/5)` | algebraically equivalent form via `Sinh[...]^2` and `Cosh[...]` |
| `Upsilon_fail_shell` | `2*sqrt(5)*Pe_req/(5*Lambda_ell)` | `(2*peReq)/(Sqrt[5]*lambdaEll)` |
| `Upsilon_suff_shell` | `4*Pe_req*(2*sqrt(5)+5)/25` | `(4*(5+2*Sqrt[5])*peReq)/25` |
| `Upsilon_fail_comp` | `2*Pe_req*chi_s/Lambda_ell**2` | `(2*chiS*peReq)/lambdaEll^2` |
| `Upsilon_suff_comp` | `4*Pe_req*chi_s**2*(Lambda_ell+2*chi_s)/Lambda_ell**3` | `(4*chiS^2*(2*chiS+lambdaEll)*peReq)/lambdaEll^3` |
| Shell `Delta_0` ratio | `1` (after `sp.simplify`) | `1` (mathematica output line 21) |
| Shell `Delta_inf` ratio | `1` (after `sp.simplify`) | surface form `2/Sqrt[5] + (5+2*Sqrt[5])^{-1}` which rationalizes to 1 (mathematica output line 22; `FullSimplify[ratio-1]` returns `0` so `expectZero` passes) |
| Compression ratios | `1`, `1` | `1`, `1` |

Note on the mathematica shell `DeltaInf` ratio surface form: `2/Sqrt[5] + 1/(5+2*Sqrt[5])`. Rationalizing the second term: `1/(5+2*sqrt(5)) = (5-2*sqrt(5))/((25)-(20)) = (5-2*sqrt(5))/5 = 1 - 2/sqrt(5)`. Sum = `2/sqrt(5) + 1 - 2/sqrt(5) = 1`. So this is `1` in an unreduced surface form, and `FullSimplify[ratio - 1]` correctly returns `0`. The Mathematica `expectZero` PASS is not a false pass.

Output transcripts mtime:
- sympy script: 2026-05-22 20:12:27; sympy output: 2026-05-22 21:38:42 (output newer, fresh)
- mathematica script: 2026-05-22 20:12:27; mathematica output: 2026-05-22 21:39:04 (output newer, fresh)

`outputs_fresh: true`. No `stale_output` finding.

## Verdict justification

`verdict: findings` with one low-severity finding.

What holds up against the paper:
- All ten paper deliverables (two threshold surfaces, two physical wall-amplitude thresholds, four asymptotic regime values, two leading-order anchors for `Delta_0`/`Delta_inf` per regime) have corresponding non-tautological script assertions in both engines.
- The ratio-limit anchors (A1/A2/A5/A6 and Mathematica B1/B2/B5/B6) genuinely exercise the full closed-form `Delta_0`/`Delta_inf` against the hand-built leading-order forms via each engine's native limit machinery. These were added in a prior audit cycle and are now present in the current scripts.
- Hand-built shell/comp forms substituted into the threshold formula reproduce the paper-quoted asymptotic values exactly.
- The compression-`DeltaInf_comp = 1/(2 chi_s)` simplification (which is an unconditional algebraic identity) is now anchored to the full closed form via A6/B6, so its meaning is no longer dependent on regime context.
- Mathematica surface forms agree with SymPy after simplification; engines agree at the level they claim.
- Outputs are fresh.

What doesn't hold up:
- F1 (paper_misalignment / notes_contradicts_script, low): the script docstring and banner strings still say "Stage 55"/"Stage 055" while every paper-side reference uses "Stage 072". Direction of fix is unambiguous (update scripts to say "Stage 072"). Does not require user resolution.

Attacks I tried that failed:
- Numeric-constant disagreement: paper says `2/sqrt(5)`, `4/5`, `(1 + 2/sqrt(5))`, `2 chi_s/Lambda^2`, `4 chi_s^2 (Lambda+2 chi_s)/Lambda^3`; script uses `sp.Rational(2)/sp.sqrt(5)`, `sp.Rational(4,5)`, `(1 + sp.Rational(2)/sp.sqrt(5))`, `2*chi_s/Lambda_ell**2`, `4*chi_s**2*(Lambda_ell+2*chi_s)/Lambda_ell**3`. All match exactly.
- Tautology: the four ratio-limit anchors (A1/A2/A5/A6, B1/B2/B5/B6) exercise the full closed-form `Delta_0`/`Delta_inf` against independently constructed shell/comp forms; they are not `x == x`.
- Missing branch: paper enumerates two asymptotic regimes; both have leading-order anchors and value-recovery checks.
- Symbol-domain error: all symbols declared `positive=True, real=True` in SymPy and `Element[..., Reals]` plus positivity in Mathematica. Consistent.
- Engine disagreement: residuals match algebraically.
- Stale output: both `.txt` files newer than their scripts.

Stop-cold: neither `UNFIXABLE` nor `CRITICAL_DOWNSTREAM`. F1 is a string-labelling fix; correcting it does not propagate.

`stop_cold: null`; orchestrator may invoke Codex on F1.

## Self-test notes

- Walked through `(1 + 1/c)/((c+1)*Lambda)` with `c = 2/sqrt(5)` and confirmed it rationalizes to `sqrt(5)/(2*Lambda)`. So `Pe_req/(Lambda^2 * DeltaInf_shell) = Pe_req * 2*Lambda/(sqrt(5)*Lambda^2) * Lambda = 2*Pe_req/(sqrt(5)*Lambda)`. Matches `Upsilon_fail_shell` literal. (No hidden factor.)
- Walked through `DeltaInf_comp = (1 + Lambda/(2 chi))/(2 chi + Lambda)` -> common factor `(2 chi + Lambda)/(2 chi)` cancels, leaving `1/(2 chi)`. So `Pe_req/(Lambda^2 * DeltaInf_comp) = 2 chi * Pe_req / Lambda^2`. Matches.
- F1 is purely a string replacement and cannot introduce a math regression. No assertion lines edited.
- Confirmed F1 fix would not introduce a new `paper_misalignment`: the new strings exactly match the paper card identifier "Stage 072" and the filename.
- Variable-independence trap (n/a for F1 — no diffs added).
- No new branch claim added.

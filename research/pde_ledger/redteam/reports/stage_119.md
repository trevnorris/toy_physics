---
unit_id: 119
batch: IV.3
auditor_model: claude-opus-4-7-1m
audit_date: 2026-05-27T00:00:00Z
verdict: findings
stop_cold: null
findings_count: 2
paper_alignment: aligned
scripts_checked:
  sympy: present
  mathematica: present
  engines_agree: true
  outputs_fresh: true
docs_read:
  paper_stage_tex: present
  notes_stage_files:
    - /var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage119_parent_balance_family.md
  paper_appendix: present
---

# Audit unit 119 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_119.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage119_parent_balance_family.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part04.tex` (stage 119 is included via `\input{stages/stage_119}` at line 1272; no further per-stage prose beyond the card itself)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage119_parent_balance_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage119_parent_balance_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage119_parent_balance_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage119_parent_balance_mathematica_audit.txt`

## What the paper claims

The stage card is a derivation-ledger entry whose load-bearing quoted claim is the boxed identity:

> Normalized ratios `(r, g)` obey `1 + r^2 = 4(g - r)^2`.

The card frames this as the reduction of the compensated outlet branch to a relation between two normalized parent ratios. The accompanying notes elaborate: define `r := lambda/sqrt(K_s K_q)` and `g := g_q sqrt(K_s)/(g_s sqrt(K_q))`; show that the Stage-98 core-balance theorem `g_s^2 (K_s K_q + lambda^2) = 4 (K_s g_q - lambda g_s)^2` collapses to the boxed quadratic; solve for the two-branch law `g = r +/- (1/2) sqrt(1+r^2)`; derive the D/N tube length `L_W = (pi a / 2) sqrt((1 + r^2)/3)` from `kappa_c = 1/3, r_c = r^2`; supply explicit parent-action expressions for `r`, `g_s`, `g_q`, `K_q`; and finally invert to a traction law `T_m = sqrt(2 Z_q K_s)/(J_s c_s sqrt(mu_0 L_W)) * 1/(r +/- (1/2) sqrt(1+r^2))`. The paper card itself quotes only the boxed reduction; the notes carry the auxiliary deliverables.

## What the script claims to verify

Both scripts perform four matched sections. Section I substitutes `lam -> rhat sqrt(K_s K_q)` and `g_q -> ghat g_s sqrt(K_q)/sqrt(K_s)` into the Stage-98 core-balance form and asserts the reduced expression equals `1 + rhat^2 - 4(ghat - rhat)^2`. Section II solves the dimensionless quadratic for `ghat` and asserts both branches satisfy the original equation. Section III sets `kappa0 := 4 L_W^2/(pi^2 a^2)`, solves `kappa0 == (1 + rc)/3` for `L_W`, and asserts the solution equals `pi a sqrt((1+rc)/3)/2`. Section IV substitutes the explicit `K_q` and `g_q` formulas into `ghat = g_q sqrt(K_s)/(T_m J_s sqrt(K_q))` and asserts the simplification equals `sqrt(2 Z_q K_s)/(T_m J_s c_s sqrt(mu_0 L_W))`; it then solves the two-branch matching condition for `T_m` and prints the result (no assertion).

## Paper ↔ script cross-check

| Paper deliverable | Script-side coverage |
|---|---|
| Boxed `1+r^2 = 4(g-r)^2` from Stage-98 form | match — Section I substitutes and asserts the reduction |
| Two-branch `g = r +/- (1/2) sqrt(1+r^2)` | match — Section II solves and double-checks both branches |
| `L_W = (pi a/2) sqrt((1+r^2)/3)` from `kappa_c=1/3, r_c=r^2` | partial — Section III only confirms the algebraic rearrangement of `4 L_W^2/(pi^2 a^2) = (1+rc)/3`; it does not show that `kappa_c=1/3` and the operator-eigenvalue identification produce `(1+rc)/3` |
| `ghat = sqrt(2 Z_q K_s)/(T_m J_s c_s sqrt(mu_0 L_W))` | match — Section IV substitutes the parent-action formulas for `K_q`, `g_q` and asserts the simplification |
| Boxed `T_m = sqrt(2 Z_q K_s)/(J_s c_s sqrt(mu_0 L_W)) * 1/(r +/- (1/2) sqrt(1+r^2))` (notes §5; not in paper card) | partial — Section IV computes `T_m` by solving but only prints; no assertion against the notes' boxed form |

The paper card's quoted boxed claim is fully exercised. Notes-level deliverables III and V are covered only partially (no false statements, just no anchor assertion).

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 36 | `expect_zero("dimensionless law", law_red - (1 + rhat**2 - 4*(ghat-rhat)**2))` | claim 1 (boxed reduction) | yes |
| A2 | sympy | 42 | `expect_zero("positive branch check", 1 + rhat**2 - 4*(ghat_sol[0]-rhat)**2)` | claim 2 (two-branch law) | yes |
| A3 | sympy | 43 | `expect_zero("negative branch check", ...)` | claim 2 | yes |
| A4 | sympy | 52 | `expect_zero("tube-length law", L_sel - sp.pi*a*sp.sqrt((1+rc)/3)/2)` | claim 3 (L_W) | partial (see F1) |
| A5 | sympy | 61-64 | `expect_zero("ghat explicit simplification", ghat_expr - sqrt(2*Zq*K_s)/(Tm*J_s*c_s*sqrt(mu0*L_W)))` | claim 4 (ghat explicit) | yes |
| A6 | sympy | 66-69 | print only (no assertion) | claim 5 (T_m traction law) | no (see F2) |
| B1 | mathematica | 40 | `expectZero["dimensionless law", lawRed - (1 + rHat^2 - 4*(gHat - rHat)^2)]` | claim 1 | yes |
| B2 | mathematica | 47 | `expectZero["first branch check", ...]` | claim 2 | yes |
| B3 | mathematica | 48 | `expectZero["second branch check", ...]` | claim 2 | yes |
| B4 | mathematica | 59 | `expectZero["tube-length law", lSel - (Pi*a*Sqrt[(1 + rC)/3])/2]` | claim 3 | partial (see F1) |
| B5 | mathematica | 73-76 | `expectZero["gHat explicit simplification", ...]` | claim 4 | yes |
| B6 | mathematica | 78-82 | print only (no assertion) | claim 5 | no (see F2) |

## Findings

### F1 — tautological_check

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage119_parent_balance_sympy_audit.py:48-52`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage119_parent_balance_mathematica_audit.wl:52-59`

**What's wrong:**
Section III of both scripts defines `kappa0 = 4*L_W**2/(pi**2*a**2)` and then `solve(kappa0 == (1+rc)/3, L_W)`, asserting the solution equals `pi*a*sqrt((1+rc)/3)/2`. This is purely algebraic: given the equation `4 L_W^2/(pi^2 a^2) = (1+rc)/3`, solving for positive `L_W` necessarily gives that closed form — the assertion cannot fail no matter what physics says. The notes' substantive content — that `kappa_c = 1/3` (D/N eigenvalue selection) plus `r_c = r^2` (parent-ratio identification) imply the equation `4 L_W^2/(pi^2 a^2) = (1+r^2)/3` in the first place — is never exercised in the script. The script effectively inputs the boxed formula in implicit form and confirms its explicit form. The variable name `kappa0` is misleading: it is not the κ_c eigenvalue defined elsewhere; it is just the dimensionless combination `4 L_W^2/(pi^2 a^2)` that the script asserts equals `(1+rc)/3`.

**Why this matters:**
If the upstream identification `kappa_c = 1/3` or `r_c = r^2` were wrong, this section would still pass. The check provides no defense against errors in the D/N-tube derivation chain.

**Required change:**
Either (a) strengthen the check to derive the implicit equation `4 L_W^2/(pi^2 a^2) = (1+r^2)/3` from a more upstream source (e.g., substitute `rc -> rhat**2` and confirm the result reads `L_W = pi a sqrt((1+rhat^2)/3)/2`, demonstrating consistency between sections I/II and the L_W formula), OR (b) annotate the assertion as algebraic-rearrangement-only so future readers do not over-credit it. Option (a) is preferable. Concretely, add (after the existing assertion in both scripts) a substitution-form check: substitute `rc -> rhat**2` into both `L_sel` and the boxed formula and assert they agree, anchoring this section's reduction step to the same `rhat` symbol used by sections I-II. This converts the section from "algebra check" to "consistency with the rest of the family parameter".

**Verification:**
After the fix, the .py and .wl scripts each contain an additional `expect_zero`/`expectZero` line at the bottom of Section III that substitutes `rc -> rhat**2` and confirms `L_sel == pi*a*sqrt((1+rhat**2)/3)/2`. The output files show the new check name and value `0`.

### F2 — insufficient_verification

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage119_parent_balance_sympy_audit.py:66-69`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage119_parent_balance_mathematica_audit.wl:78-82`

**What's wrong:**
The Section IV banner reads "Explicit traction law" / "IV. EXPLICIT TRACTION LAW", and the script solves the two-branch matching condition `ghat_expr == rhat +/- sqrt(1+rhat^2)/2` for `T_m`, but the result is only `print`-ed; no `expect_zero` / `expectZero` assertion follows. The notes (§5) supply an explicit boxed form `T_m = sqrt(2 Z_q K_s)/(J_s c_s sqrt(mu_0 L_W)) * 1/(rhat +/- (1/2) sqrt(1+rhat^2))` that the script never tests against. The script's own stated topic ("Explicit traction law") is therefore not actually verified.

**Why this matters:**
Either of the two `T_m` outputs could disagree with the notes' boxed formula and the script would still exit 0. The script's own banner overstates what it has confirmed.

**Required change:**
Append two assertions immediately after the two `print` lines in each script:
- assert `Tm_sol_plus - 2*sqrt(2)*sqrt(K_s)*sqrt(Zq)/(J_s*sqrt(L_W)*c_s*sqrt(mu0)*(2*rhat + sqrt(1+rhat**2))) == 0`
- assert `Tm_sol_minus - 2*sqrt(2)*sqrt(K_s)*sqrt(Zq)/(J_s*sqrt(L_W)*c_s*sqrt(mu0)*(2*rhat - sqrt(1+rhat**2))) == 0`

(The factor-of-2 rearrangement `1/(r +/- (1/2)*sqrt(1+r^2)) = 2/(2 r +/- sqrt(1+r^2))` matches the current printed output exactly.) For the Mathematica script, the analogous assertions go in the Section IV closing block; the Mathematica solver returns `ConditionalExpression[...]`, so the expected-form expression must be wrapped in the same form (or the `ConditionalExpression` head stripped by `expectZero` — see the established Mathematica idiom for this pipeline).

**Verification:**
After the fix, both engines emit two new check lines with value `0`/`PASS` for the `T_m` branches, and the script names contain a recognizable token (e.g., `T_m (+ branch) match` / `T_m (- branch) match`).

## Independent-derivation check (Mathematica)

The `.wl` script mirrors the `.py` script section-by-section with identical algebraic choreography. For a purely algebraic identity audit of this nature (substitution and `simplify`/`FullSimplify`), the canonical path is essentially fixed and an "independent derivation" cannot diverge in form without re-stating the same algebra. The Mathematica script uses native idioms (`FullSimplify[..., Assumptions -> ...]`, `Solve[..., Reals]`, `ConditionalExpression`), and the assumptions list is independently re-declared per section. I do not flag this as `mathematica_transliteration`; the parallel structure reflects the parallel algebra, not blind transliteration.

## Engine cross-check

Sections I, II, III, IV: both engines report residual `0` for all assertions. Sections I/II reduced laws match: SymPy `−4 ghat² + 8 ghat rhat − 3 rhat² + 1`, Mathematica `1 − 4 gHat² + 8 gHat rHat − 3 rHat²` (identical). Section III L_W: SymPy `pi a sqrt(3 rc + 3)/6`, Mathematica `(a pi sqrt(1 + rC))/(2 sqrt(3))` (identical after rationalization). Section IV `gHat explicit`: SymPy `sqrt(2) sqrt(K_s) sqrt(Zq)/(J_s sqrt(L_W) Tm c_s sqrt(mu0))`, Mathematica `(sqrt(2) kS zQ)/(cSound jS tM sqrt(kS lW mu0 zQ))` (identical after pulling `sqrt(kS zQ)` out of the radical). `T_m` branch outputs algebraically agree (Mathematica's `−` branch surface-form has a leading `−` over `(−2 rHat + sqrt(1 + rHat^2))`, which is the same as SymPy's `+` over `(2 rhat − sqrt(1 + rhat^2))`). No `engine_disagreement`.

## Verdict justification

The paper card's load-bearing quoted equation `1 + r^2 = 4(g − r)^2` is verified non-tautologically by both engines (Section I substitutes the dimensionless ratios into the Stage-98 form and the residual is `0`). Section II's two-branch resolution and Section IV's `ghat` explicit-form simplification are also non-tautological substantive checks. Section III is mildly tautological (F1) and Section IV's `T_m` derivation lacks a final assertion against the notes' boxed form (F2). Both findings are low severity: the paper card's strict quoted deliverable is fully covered, and the notes-level deliverables that are not asserted produce printed expressions that happen to match the notes when compared by hand. The audit is `findings` with both items localized to clearly-named blocks; nothing propagates downstream.

Attacks tried that did not succeed: (i) checked whether section IV's substitution of `g_s -> Tm*J_s` matches the notes — it does (notes line `g_s = T_m J_s`); (ii) checked sign consistency between SymPy and Mathematica T_m outputs — they agree algebraically despite differing surface forms; (iii) checked whether the script's `kappa0` matches a κ_c definition — it does not match κ_c (it is just `4 L_W^2/(pi^2 a^2)`), which feeds F1 but is not itself a contradiction; (iv) checked the `gq` and `Kq` symbolic forms against the notes' parent-action formulas — they match exactly including the `mu_0` denominators; (v) checked all positivity assumptions for compatibility with the physical setup — all consistent.

## Self-test notes

I mentally executed the F1 fix: substituting `rc -> rhat**2` into `pi*a*sqrt((1+rc)/3)/2` yields `pi*a*sqrt((1+rhat**2)/3)/2`, which matches the boxed L_W formula in the notes; the new assertion's residual is identically zero, and the check is non-tautological because it links section III's `rc` to sections I-II's `rhat`. For F2, I mentally substituted `rhat -> 0` into the proposed `T_m` assertion: `Tm_sol_plus(rhat=0) = 2 sqrt(2) sqrt(K_s) sqrt(Zq)/(J_s sqrt(L_W) c_s sqrt(mu0) * 1) = 2 sqrt(2 K_s Z_q)/(J_s c_s sqrt(mu_0 L_W))`, matching the notes' `sqrt(2 Z_q K_s)/(J_s c_s sqrt(mu_0 L_W)) * 1/(0 + 1/2 * 1) = 2 sqrt(2 Z_q K_s)/(J_s c_s sqrt(mu_0 L_W))` ✓. No variable-independence or parity traps apply — the section IV expressions only contain symbolic constants, no derivatives or integrals. Paper round-trip: both proposed fixes add assertions consistent with the notes; neither introduces a constant the paper or notes does not state.

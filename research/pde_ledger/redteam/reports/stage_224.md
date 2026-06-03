---
unit_id: 224
batch: VII.1
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-06-02T00:00:00-06:00
verdict: findings
stop_cold: null
findings_count: 2
paper_alignment: aligned
scripts_checked:
  sympy: present
  mathematica: missing
  engines_agree: n/a
  outputs_fresh: true
docs_read:
  paper_stage_tex: present
  notes_stage_files: [moving_throat_pde_stage224_pde_branch_packet_compiler_weak_axisymmetric_ceiling_transport_and_first_actual_branch_kill_test_sympy_audit.md]
  paper_appendix: present
---

# Audit unit 224 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_224.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage224_pde_branch_packet_compiler_weak_axisymmetric_ceiling_transport_and_first_actual_branch_kill_test_sympy_audit.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part07.tex` (rows at lines 60, 518, 599-679, 1437)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage224_pde_branch_packet_compiler_weak_axisymmetric_ceiling_transport_and_first_actual_branch_kill_test_sympy_audit.py`
- mathematica: (missing)
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage224_pde_branch_packet_compiler_weak_axisymmetric_ceiling_transport_and_first_actual_branch_kill_test_sympy_audit.txt`
- mathematica output: (missing)

## What the paper claims

Stage 224 (`\stagefield{Output}`): "Actual-branch kill test: the branch survives only if the transported quantity $(\Delta_{\rm norm}+T_{\rm quad})(1+|\varepsilon\Xi_1|)/\hat m_0^2$ stays below the critical prefactor budget." The stage compiles the final 5PN branch packet into the three lane prefactors via the grouped inverse map $P_{20}=\bar P_0+4a_{P_0}$, $P_{21}=\bar P_0-a_{P_0}+b_{P_0}$, $P_{22}=\bar P_0-a_{P_0}-b_{P_0}$ with $\bar P_0=(\Delta_{\rm norm}+T_{\rm quad})/\hat m_0^2$; states the isotropic ceiling $\Delta_{\rm norm}\le\hat m_0^2 P_{\rm crit}-T_{\rm quad}$; the weak-axisymmetric signature $(\lambda_{20},\lambda_{21},\lambda_{22})=(1,\tfrac12,-1)$ giving $P_A=\bar P_0(1+\varepsilon\lambda_A\Xi_1)$, hence $a_{P_0}=\varepsilon\bar P_0\Xi_1/4$, $b_{P_0}=3\varepsilon\bar P_0\Xi_1/4$, so $b_{P_0}=3a_{P_0}$; the robust all-lane collapse $\bar P_0(1+|\varepsilon\Xi_1|)\le P_{\rm crit}$ (= $\bar P_0+4|a_{P_0}|$); the calibrated lower bounds $\hat m_0^2\ge T_{\rm quad}/P_{\rm crit}$ and $\hat m_0^2\ge T_{\rm quad}(1+|\varepsilon\Xi_1|)/P_{\rm crit}$; and four explicit headroom budgets at the Stage-240 (script comments say Stage-223) compatibility point $\bar P_0\approx0.002069792318062885$. Notes section 9 enumerates these as the eight things the script verifies. The appendix Theorem (`thm:app-part07-actual-branch-ceiling`) states the transported necessary condition. Claim status: "exact for ... packet compiler"; numerical headroom is `\StatusNumerical{}` and the primitive-family ceilings are a transported reduced test, not an unconditional full-PDE theorem.

## What the script claims to verify

The SymPy script verifies, symbolically: (1) the grouped inverse map round-trips $(\bar P_0,a_{P_0},b_{P_0})\leftrightarrow(P_{20},P_{21},P_{22})$ (lines 30-36); (2) the isotropic ceiling rearrangement $(\bar P_0-P_{\rm crit})\hat m_0^2 = \Delta_{\rm norm}-(\hat m_0^2 P_{\rm crit}-T_{\rm quad})$ (line 52); (3) the weak-axisymmetric lane law produces $a_{P_0}=\varepsilon\bar P_0\Xi_1/4$, $b_{P_0}=3\varepsilon\bar P_0\Xi_1/4$, $b_{P_0}=3a_{P_0}$, and that those defects re-expand to the WA lanes (lines 76-86); (4) the robust max-lane collapse for both signs of $\varepsilon\Xi_1$ equals $\bar P_0(1+|\varepsilon\Xi_1|)=\bar P_0+4|a_{P_0}|$ (lines 106-120). It then numerically (lines 133-159) hardcodes a compatibility point and four ceilings and asserts the derived budgets $P_{\rm crit}/\bar P_0-1$ and $(P_{\rm crit}-\bar P_0)/4$ equal four further hardcoded budget literals via `assert_close(tol=1e-12)`.

## Paper ↔ script cross-check

| Paper deliverable | Script check | Status |
|---|---|---|
| Grouped inverse map $\bar P_0=(\Delta_{\rm norm}+T_{\rm quad})/\hat m_0^2$, lanes (eq app-part07-prefactor-lanes) | lines 25-36 round-trip | match |
| Isotropic ceiling $\Delta_{\rm norm}\le\hat m_0^2 P_{\rm crit}-T_{\rm quad}$ | line 52 | match |
| WA signature $(1,\tfrac12,-1)$, $P_A=\bar P_0(1+\varepsilon\lambda_A\Xi_1)$ | lines 65-71 | match |
| $a_{P_0}=\varepsilon\bar P_0\Xi_1/4$, $b_{P_0}=3\varepsilon\bar P_0\Xi_1/4$, $b_{P_0}=3a_{P_0}$ | lines 73-78 | match |
| Robust collapse $\bar P_0(1+|\varepsilon\Xi_1|)=\bar P_0+4|a_{P_0}|$ | lines 99-120 | match |
| Calibrated bounds $\hat m_0^2\ge T_{\rm quad}/P_{\rm crit}$; $\ge T_{\rm quad}(1+|\varepsilon\Xi_1|)/P_{\rm crit}$ | lines 54, 122 (computed, not asserted) | partial (printed, not assert-checked) |
| Four headroom budgets at compat point | lines 133-159 | match in form, but self-referential (see F2) |

`paper_alignment: aligned` — every symbolic deliverable has a faithful, non-tautological script-side check; the WA lane law, signature, and robust collapse match the appendix equations exactly. The calibrated-bound expressions (line 54, 122) are computed and printed but never `assert`-ed; they are trivially `T_quad/Pcrit` by construction so a check would be tautological anyway, so this is not flagged as `insufficient`.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 34 | `simplify(inv_bar - Pbar) == 0` | inverse map $\to\bar P_0$ | yes |
| A2 | sympy | 35 | `simplify(inv_a - aP) == 0` | inverse map $\to a_{P_0}$ | yes |
| A3 | sympy | 36 | `simplify(inv_b - bP) == 0` | inverse map $\to b_{P_0}$ | yes |
| A4 | sympy | 52 | `expand((Pbar-Pcrit)*mhat0**2)==expand(Delta_norm-rhs)` | isotropic ceiling | yes |
| A5 | sympy | 76 | `simplify(aP_wa - eps*Pbar*Xi1/4)==0` | $a_{P_0}$ WA compiler | yes |
| A6 | sympy | 77 | `simplify(bP_wa - 3*eps*Pbar*Xi1/4)==0` | $b_{P_0}$ WA compiler | yes |
| A7 | sympy | 78 | `simplify(bP_wa - 3*aP_wa)==0` | $b_{P_0}=3a_{P_0}$ | yes |
| A8 | sympy | 84-86 | `simplify(P2x_from_ab - P2x_wa)==0` | WA lanes re-expand | yes |
| A9 | sympy | 106-107 | pos-branch lane gaps vanish | robust collapse (eps Xi1>0) | yes |
| A10 | sympy | 112-113 | neg-branch lane gaps vanish | robust collapse (eps Xi1<0) | yes |
| A11 | sympy | 119-120 | `robust_*_form.subs(...) - Pbar*(1+zabs)==0` | robust form identity | yes |
| A12 | sympy | 152-159 | `assert_close(budget, <literal>)` x8 | headroom budgets | partial — self-referential (F2) |

## Findings

### F1 — missing_verification_script

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/` (no `.wl` for unit 224)
- `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_224.tex:11` ("Mathematica audit: none yet.")

**Subtype:** `missing_mathematica`

**What's wrong:**
No Mathematica audit script exists for this unit. Under the dual-engine rule ("is it POSSIBLE for Mathematica to independently verify," not "is it necessary"), Mathematica can trivially re-derive every claim here: all of the symbolic content is elementary linear algebra over a small symbol set (grouped inverse map, weak-axisymmetric lane law, robust max-lane collapse) and the numerical content is exact rational/high-precision decimal arithmetic. There is no transcendental closure, no special function, no numerical root-finding done *in this stage* (the bisection lives upstream in stage 223; stage 224 only divides the carried-forward literals). It is fully within Mathematica's reach. The card itself records "none yet," confirming this is a gap, not an impossibility.

**Why this matters:**
This is a checkpoint-class non-status unit (manifest `is_status_only_candidate: False`, `is_checkpoint: False` per prompt). The dual-engine policy requires a second, independently-derived engine wherever possible. A single-engine symbolic identity can hide a SymPy-specific `simplify`/`Abs` artifact (the script leans on the fragile `subs({Abs(eps*Xi1): zabs, ...})` trick at lines 117-120); an independent Mathematica derivation that never forms those `Abs` substitutions guards against that.

**Required change:**
Create the Mathematica audit script (claim manifest M1-M6 in the directive). It must derive the claims natively (Mathematica primitives, a different decomposition than the SymPy round-trip), NOT transliterate the `.py`.

**Verification:**
`math -script <target>` exits 0 with all `expectZero`/`If[...,Exit[1]]` checks passing; the new `.wl` and its committed `output/*.txt` appear; the file is an independent derivation, not a line-by-line port of the SymPy script.

### F2 — hardcoded_result

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/...stage224...sympy_audit.py:133-159`

**What's wrong:**
The numerical headroom block hardcodes BOTH the inputs and the expected outputs, then checks the outputs against themselves. `barP0_compat = Decimal("0.002069792318062885")` (line 133) and the four `ceilings` literals (lines 135-140) are hardcoded; the budgets are computed as `Pcrit_val/barP0_compat - 1` and `(Pcrit_val - barP0_compat)/4` (lines 144-145); then lines 152-159 assert those computed budgets equal eight further hardcoded literals (`0.367930328492646`, etc.). Because the expected literals are just the pre-evaluated arithmetic of the same two hardcoded inputs, the assertions cannot fail unless decimal arithmetic itself is broken — they verify nothing about the upstream physics. Separately, the carried literals do not exactly match the upstream stage-223 source values: stage 223's saved output gives the compat point `0.002069792318062883` and the 10% walls `0.0028313316855593336` / `0.003596510589684656`, whereas stage 224 hardcodes `...62885`, `0.0028313316855593175`, `0.0035965105896846573`. The differences are at the ~13th-16th significant figure, inside stage 223's own `tol=5e-13` and stage 224's `tol=1e-12`, i.e. floating-point/bisection tail noise — NOT a genuine value conflict (so not a `paper_misalignment`), but the literals are not the exact upstream strings.

**Why this matters:**
The block contributes no independent verification: it is a self-consistent arithmetic tautology over hardcoded data. It also silently re-types upstream numbers with altered trailing digits, which weakens traceability to stage 223 and the notes (which print the same headroom values, themselves derived from the same hardcoded inputs).

**Required change:**
Replace the self-referential `assert_close(budget, <literal>)` pair-set with a check that ties the budgets to the upstream-carried ceiling/compat data rather than to a second copy of the answer. Concretely: keep the four ceiling literals and the compat point as the carried inputs (cite stage 223 in a comment), compute the budgets, and assert the *defining relations* hold — e.g. `assert_close((eps_xi_budget + 1) * barP0_compat, Pcrit_val)` and `assert_close(barP0_compat + 4*a_budget, Pcrit_val)` — so the assertion exercises the budget formula (eq app-part07-xi-prefactor-ceiling) against the ceiling, not against a pre-baked copy of itself. Do NOT add or change any physical constant; do not alter the carried literal values (they are within tolerance of upstream). This is a low-severity hardening; if Codex judges the round-trip rewrite ambiguous it may instead simply add the inverse-check assertions alongside the existing ones.

**Verification:**
The new assertions reference `Pcrit_val` and `barP0_compat` on both sides of the budget relation (not a third hardcoded budget literal); script exits 0; printed headroom numbers unchanged.

## Independent-derivation check (Mathematica)

N/A — no `.wl` exists. See F1 (`missing_mathematica`).

## Engine cross-check

N/A — single engine present.

## Verdict justification

The SymPy script's symbolic core is sound and faithfully aligned with the paper: the grouped inverse round-trip (A1-A3), isotropic ceiling rearrangement (A4), weak-axisymmetric compiler and $b_{P_0}=3a_{P_0}$ (A5-A8), and the two-sign robust max-lane collapse (A9-A11) each non-tautologically reproduce the exact appendix equations (`app-part07-prefactor-lanes`, `app-part07-ap0-bp0-xi`, `app-part07-xi-prefactor-ceiling`) and notes section 9. Attacks tried that failed: (a) checked the inverse-map coefficients $(1,2,2)/5$, $(2,-1,-1)/10$ actually conspire with the lane definitions — they do, and the round-trip is not guaranteed by construction; (b) checked the WA signature $(1,\tfrac12,-1)$ propagates correctly through the trace/anomaly projector to give exactly $a_{P_0}=\varepsilon\bar P_0\Xi_1/4$ — it does; (c) checked the `Abs`-substitution trick (lines 117-120) is valid given $\bar P_0>0$ (positive symbols) so $\mathrm{Abs}(\varepsilon\bar P_0\Xi_1/4)=\bar P_0\,\mathrm{Abs}(\varepsilon\Xi_1)/4$ — valid, and the saved output confirms the asserts passed. Two findings remain: a medium `missing_mathematica` (the dual-engine rule requires a second engine and this is fully Mathematica-tractable), and a low `hardcoded_result` (the numerical headroom block checks pre-baked answers against themselves and re-types upstream literals with altered trailing digits). Neither is stop-cold; the symbolic physics holds. Verdict: `findings`.

## Self-test notes

I checked (1) variable independence — no `diff` calls in this stage, so the zero-derivative trap does not apply; (2) parity/symmetry — no integrals; (3) trivial-case substitution for the inverse map and WA compiler (worked the algebra by hand above; the round-trips are genuine, not construction-guaranteed); (4) the `Abs` positivity assumption (`Delta_norm,T_quad,mhat0` positive => `Pbar` positive, so the `subs` of `Abs(eps*Pbar*Xi1/4)` resolves cleanly — the saved output confirms it passed). For F2's proposed inverse-relation assertions I verified `(eps_xi_budget+1)*barP0_compat == Pcrit_val` and `barP0_compat + 4*a_budget == Pcrit_val` reduce to identities of the budget definitions, so they are non-tautological w.r.t. the budget formula yet true for the carried data. Paper round-trip: F2 introduces no new constant and leaves the carried literals/printed values unchanged, so it cannot create a new misalignment.

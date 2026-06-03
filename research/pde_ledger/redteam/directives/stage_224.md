---
unit_id: 224
batch: VII.1
created_at: 2026-06-02T00:00:00-06:00
findings_count: 2
stop_cold: null
applied: true
applied_at: 2026-06-02T16:34:54-06:00
findings_applied: 2
findings_blocked: 0
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 224

Apply each finding below in order. After applying, append an `## Applied: F<n>` block under that finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

If a finding's required change is ambiguous or unsafe to apply mechanically, append `## Blocked: F<n>` with a question instead — skip that finding, continue with the rest.

Do NOT introduce new features, refactors, or stylistic changes. Edit exactly the file:line ranges named.

After editing, RUN the affected scripts (`python3 <path>` for SymPy, `math -script <path>` for Mathematica) and iterate until they exit 0 with all in-file checks passing. Getting the scripts to run cleanly is your job; the orchestrator independently re-runs afterward.

Do NOT touch paper.tex, notes/, or any prose documents.

## F1 — missing_verification_script (subtype: missing_mathematica)

**Target:** `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage224_pde_branch_packet_compiler_weak_axisymmetric_ceiling_transport_and_first_actual_branch_kill_test_mathematica_audit.wl`

**Issue:** No Mathematica audit script exists for unit 224. Under the dual-engine rule, Mathematica CAN independently verify every claim in this stage (all elementary symbolic algebra over a small symbol set, plus exact arithmetic), so a `.wl` is required. The new script must be an INDEPENDENT re-derivation, NOT a transliteration of the SymPy file. Anti-transliteration guard: use native Mathematica primitives (`Solve`/`LinearSolve` for the grouped inversion rather than echoing the SymPy `(P20+2*P21+2*P22)/5` projector; `Simplify`/`FullSimplify` with `Assumptions`; `Abs` resolved via `Simplify[..., Assumptions -> {Pbar > 0, ...}]` or by case-splitting on the sign of `eps*Xi1` with `Reduce`, NOT by an `Abs`-substitution dictionary). Choose a different decomposition than the SymPy round-trip wherever practical (e.g. derive the inverse map by solving the 3x3 forward system symbolically rather than asserting pre-chosen projector coefficients).

**Required change:**
Write the `.wl` at the exact Target path above. It must verify, with hard checks (`If[Simplify[lhs - rhs] =!= 0, Print["FAIL ..."]; Exit[1]]` or an `expectZero` helper that strips `ConditionalExpression[0, ...]`):

**Claim manifest:**

- **M1 (grouped inverse map):** With `Pbar = (Delta_norm + T_quad)/mhat0^2`, `P20 = Pbar + 4 aP`, `P21 = Pbar - aP + bP`, `P22 = Pbar - aP - bP`, solving the forward 3-vector system recovers exactly `Pbar`, `aP`, `bP`. Derive the inverse by `Solve[{P20 == Pbar + 4 aP, P21 == Pbar - aP + bP, P22 == Pbar - aP - bP}, {Pbar, aP, bP}]` (independent route), then confirm the solution equals the input symbols. (LaTeX: $\bar P_0=(P_{20}+2P_{21}+2P_{22})/5,\ a_{P_0}=(2P_{20}-P_{21}-P_{22})/10,\ b_{P_0}=(P_{21}-P_{22})/2$.)
- **M2 (isotropic ceiling rearrangement):** `Simplify[(Pbar - Pcrit)*mhat0^2 - (Delta_norm - (mhat0^2*Pcrit - T_quad))] == 0`. (LaTeX: $\Delta_{\rm norm}\le\hat m_0^2 P_{\rm crit}-T_{\rm quad} \iff \bar P_0\le P_{\rm crit}$.)
- **M3 (weak-axisymmetric compiler):** With `lam20=1, lam21=1/2, lam22=-1` and `PAwa = Pbar (1 + eps lamA Xi1)`, the trace/anomaly projector gives `aPwa = eps Pbar Xi1 / 4`, `bPwa = 3 eps Pbar Xi1 / 4`, hence `bPwa == 3 aPwa`. Verify each with `Simplify[... ] == 0`. (LaTeX: $a_{P_0}=\varepsilon\bar P_0\Xi_1/4$, $b_{P_0}=3\varepsilon\bar P_0\Xi_1/4$, $b_{P_0}=3a_{P_0}$.)
- **M4 (WA lanes re-expand):** Substituting `aPwa, bPwa` back into `Pbar + 4 aP`, `Pbar - aP + bP`, `Pbar - aP - bP` reproduces `P20wa, P21wa, P22wa`. (LaTeX: $P_A=\bar P_0(1+\varepsilon\lambda_A\Xi_1)$.)
- **M5 (robust max-lane collapse):** For `eps Xi1 > 0` the max lane is `P20 = Pbar(1+|eps Xi1|)`; for `eps Xi1 < 0` the max lane is `P22 = Pbar(1+|eps Xi1|)`; and `Pbar(1+|eps Xi1|) == Pbar + 4 |aPwa|` given `Pbar > 0`. Establish this by sign-case `Reduce`/`Refine` (e.g. `Refine[Abs[eps Xi1], eps Xi1 > 0]`), NOT by an `Abs`-substitution dictionary. (LaTeX: $\max\{P_{20},P_{21},P_{22}\}=\bar P_0(1+|\varepsilon\Xi_1|)=\bar P_0+4|a_{P_0}|$.)
- **M6 (headroom budget relations at the carried compatibility point):** Using the carried `barP0compat = 0.002069792318062885` and the four carried ceilings (`both_10 = 0.0028313316855593175`, `one_10 = 0.0035965105896846573`, `both_30 = 0.00817339430971383`, `one_30 = 0.0116633929790174`; cite stage 223 in a comment), compute `epsXiBudget = Pcrit/barP0compat - 1` and `aBudget = (Pcrit - barP0compat)/4`, and verify the DEFINING relations `(epsXiBudget + 1)*barP0compat == Pcrit` and `barP0compat + 4*aBudget == Pcrit` to high precision for each ceiling (use `SetPrecision`/exact `Rationalize` so the check is not floating-point-fragile). Do NOT assert the budgets against a second hardcoded budget literal. (LaTeX: budgets satisfy $\bar P_0(1+|\varepsilon\Xi_1|)=P_{\rm crit}$ and $\bar P_0+4|a_{P_0}|=P_{\rm crit}$ at the ceiling.)

**Verification command:**
After Codex applies, the verifier will run `redteam exec-mathematica 224` and confirm the `.wl` exists at the exact Target path, all checks pass, and the script exits 0. The committed `mathematica/output/...` transcript must be generated.

## Applied: F1

- files_changed:
  - `mathematica/moving_throat_pde_stage224_pde_branch_packet_compiler_weak_axisymmetric_ceiling_transport_and_first_actual_branch_kill_test_mathematica_audit.wl`
  - `mathematica/output/moving_throat_pde_stage224_pde_branch_packet_compiler_weak_axisymmetric_ceiling_transport_and_first_actual_branch_kill_test_mathematica_audit.txt`
- summary: Added an independent Mathematica audit for M1-M6 and generated its saved transcript.
- deviation: none

## F2 — hardcoded_result

**Target:** `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage224_pde_branch_packet_compiler_weak_axisymmetric_ceiling_transport_and_first_actual_branch_kill_test_sympy_audit.py:152-159`

**Issue:** The numerical headroom block checks pre-baked answers against themselves. Lines 133-140 hardcode the compat point and four ceilings; lines 144-145 compute the budgets `Pcrit_val/barP0_compat - 1` and `(Pcrit_val - barP0_compat)/4`; lines 152-159 then assert those budgets equal eight further hardcoded budget literals. The assertions cannot fail except on broken arithmetic — they exercise no physics. They also re-type upstream stage-223 numbers with altered trailing digits (within tolerance, so not a value conflict, but the literals are not the exact upstream strings).

**Required change:**
Replace the eight self-referential `assert_close(float(budgets[...][i]), <literal>)` lines (152-159) with checks that tie each budget to its defining relation against the carried ceiling and compat point, so the assertion exercises the budget formula (paper eq `app-part07-xi-prefactor-ceiling`) rather than a second copy of the answer. For each `key` in `ceilings`, after computing `(eps_xi_budget, a_budget)`, assert:

```python
assert_close((eps_xi_budget + Decimal(1)) * barP0_compat, Pcrit_val) ...
```
but note `assert_close` takes floats — either add a `Decimal`-aware comparison (e.g. `abs((eps_xi_budget + 1) * barP0_compat - Pcrit_val) <= Decimal("1e-30")`) or cast to float with the existing `assert_close`. Concretely, replace lines 152-159 with a loop over `budgets`/`ceilings`:

```python
    for key, (eps_xi_budget, a_budget) in budgets.items():
        Pcrit_val = ceilings[key]
        assert_close(float((eps_xi_budget + Decimal(1)) * barP0_compat), float(Pcrit_val))
        assert_close(float(barP0_compat + Decimal(4) * a_budget), float(Pcrit_val))
```

Add a one-line comment citing the upstream source of the ceilings/compat point: `# ceilings and compat point carried from Stage 223 dynamic-survival-window audit`. Do NOT change any of the hardcoded ceiling or compat-point values; do NOT add or alter any physical constant. The printed headroom numbers (lines 148-150) must remain unchanged.

If Codex judges replacing the existing checks risky, the acceptable fallback is to ADD the two inverse-relation assertions above alongside the existing eight literal checks (keeping both), which still removes the "no independent check" gap.

**Verification command:**
The verifier will run `redteam exec-sympy 224` and confirm the new inverse-relation assertions appear (referencing `Pcrit_val`/`barP0_compat` on both sides, not a third budget literal) and the script exits 0 with printed headroom values unchanged.

## Applied: F2

- files_changed:
  - `scripts/moving_throat_pde_stage224_pde_branch_packet_compiler_weak_axisymmetric_ceiling_transport_and_first_actual_branch_kill_test_sympy_audit.py`
- summary: Replaced literal headroom budget checks with inverse-relation assertions against the carried ceiling and compatibility point.
- deviation: none

## Authorized notes renumber (USER-AUTHORIZED 2026-06-02)

The user authorized renumbering stale VII.1 notes-prose stage labels to canonical. The audit logged a notes-side renumbering drift here (notes reference pre-renumber stage numbers; the card/appendix/scripts use canonical). Notes-only cleanup in THIS fix loop (Codex applies notes/*.md; Claude reviews):
- In `notes/stages/moving_throat_pde_stage224_..._sympy_audit.md`, renumber every stale stage-number reference to match the canonical numbering used in this stage's SymPy script comments + the paper card (self-reference → Stage 224; cited upstream stages → the numbers the .py comments use). Math/content unchanged.
- Do NOT touch scripts, paper.tex, or appendix. Acceptance: notes stage labels match the .py comments + card. Append `## Applied: notes-renumber`.

## Applied: notes-renumber

- files_changed:
  - `notes/stages/moving_throat_pde_stage224_pde_branch_packet_compiler_weak_axisymmetric_ceiling_transport_and_first_actual_branch_kill_test_sympy_audit.md`
- summary: Renumbered stale Stage 240/241 notes labels to canonical Stage 223/224 references without changing math or values.
- deviation: none

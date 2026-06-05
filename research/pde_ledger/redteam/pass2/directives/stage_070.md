---
unit_id: 070
batch: III.3
created_at: 2026-06-05T00:00:00Z
findings_count: 2
stop_cold: null
applied: true
applied_at: 2026-06-05T19:49:26Z
findings_applied: 2
findings_blocked: 0
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 070

Apply F1 (math de-tautologization, both engines) and F2 (label-only self-labels), then run both scripts. After applying each, append an `## Applied: F<n>` block with `files_changed`, `summary`, `deviation` (or "none").

After editing, RUN `timeout 600 python3 scripts/moving_throat_pde_stage070_gnls_wall_shell_sympy_audit.py` and `timeout 600 math -script mathematica/moving_throat_pde_stage070_gnls_wall_shell_mathematica_audit.wl`; iterate until both exit 0 with all checks passing. Do NOT touch paper.tex or notes/.

## F1 — tautological_check (de-tautologize the sech-moment anchor)

**Target:**
- `scripts/moving_throat_pde_stage070_gnls_wall_shell_sympy_audit.py` (around line 56)
- `mathematica/moving_throat_pde_stage070_gnls_wall_shell_mathematica_audit.wl` (around lines 85–86)

**Issue:** The `I_1/J_1 = 4πa²ℓ` check is a self-cancelling tautology — `(4πa²ℓ·If/Hw)/(If/Hw) = 4πa²ℓ` identically (the `If/Hw` factor cancels), independent of the sech moment value, so the sech integration is inert for it. The genuinely informative quantity (the sech profile yields `I_f = 2/3` and `I_g = 14/15`) is only printed, never asserted. Make the sech moments load-bearing.

⚠️ **Correct constants (orchestrator-verified — the prior draft's `8/15` was WRONG):**
- `I_f = ∫_{-∞}^{∞} (f')² dξ = ∫ sech²ξ·tanh²ξ dξ = 2/3` (substitute `u=tanh ξ`: `∫_{-1}^{1} u² du = 2/3`).
- `I_g = ∫_{-∞}^{∞} (f'')² dξ = 2 − 16/3 + 64/15 = **14/15**` (NOT 8/15; `f'' = sech ξ(1 − 2 sech²ξ)`, and `∫sech²=2, ∫sech⁴=4/3, ∫sech⁶=16/15`; numerically ≈ 0.9333).

**Required change:**

SymPy — leave the existing `I_1/J_1` ratio line and the "expected 2/3" print in place (documentation), and ADD, right after `If_sym = sp.integrate(...)` (line 56), a load-bearing assertion:
```python
expect_zero("sech-profile moment I_f = 2/3", sp.simplify(If_sym - sp.Rational(2, 3)))
```

Mathematica — inside the `Module` (the `INDEPENDENT NUMERIC PROFILE CROSS-CHECK` block):
1. Correct the wrong print annotation at line 86: change `(analytic 8/15 = ", N[8/15, 30]` to `(analytic 14/15 = ", N[14/15, 30]`.
2. Right after that I_g print (line 86), ADD numeric assertions using the already-computed `IfNum`/`IgNum` and the in-scope `tol` (`10^-10`):
```
If[Abs[IfNum - 2/3] < tol, pass["sech-profile moment I_f = 2/3"],
   fail["sech-profile moment I_f = 2/3", IfNum - 2/3]];
If[Abs[IgNum - 14/15] < tol, pass["sech-profile moment I_g = 14/15"],
   fail["sech-profile moment I_g = 14/15", IgNum - 14/15]];
```
Do NOT remove the existing `I_1/J_1` `expectZero` (it documents the structural normalization). Do NOT touch the symbolic `kappa`/`W_wall`/`Xi` assertions. These deliverables are unchanged → `material_change: false`.

**Self-test:** Both added asserts depend on `If_sym`/`IfNum`/`IgNum`, which genuinely depend on `xi` (no zero-derivative trap). The constants are independently verified above (symbolic + numeric). No new paper-side constant (the sech moments are the script's own anchor; Stage 071 separately uses the tanh profile with `I_f=1/3`). If SymPy's `Ig` were ever added it would be `sp.Rational(14,15)`.

**Verification:** SymPy output contains `sech-profile moment I_f = 2/3 = 0`; Mathematica output contains `PASS: sech-profile moment I_f = 2/3` and `PASS: sech-profile moment I_g = 14/15`, and its I_g print annotation reads `(analytic 14/15 = ...)`. The kappa/W_wall/Xi checks still pass; both scripts exit 0.

## F2 — stale self-labels (numbering)

**Target:** `scripts/moving_throat_pde_stage070_gnls_wall_shell_sympy_audit.py:3` and `:5`

**Issue:** The docstring carries two stale +17 pre-renumber SELF-labels: filename-style (line 3, `stage53`) and prose (line 5, `Stage 53`). Canonical is 070.

⚠️ SCOPE GUARD: Do NOT touch the cross-references at py:57 (`Stage-47 integral`), py:59/wl:87/wl:88 (`Stage-48`/`Stage-47` normalization), or the variable names `J1_stage48` (py) / `J1Stage48` (wl) — all DEFERRED to the dedicated numbering plan. Do NOT pad the `STAGE 70` ledger banner.

**Required change:**
- Line 3: `moving_throat_pde_stage53_gnls_wall_shell_sympy_audit.py` → `moving_throat_pde_stage070_gnls_wall_shell_sympy_audit.py` (`stage53`→`stage070`, 3-digit)
- Line 5: `SymPy audit for Stage 53:` → `SymPy audit for Stage 70:` (`53`→`70`, 2-digit)
- Change ONLY the number tokens; preserve everything else.

**Verification:** `git diff` shows the two docstring lines changed (cross-refs and variable names untouched), identical except the number tokens.

## F3 — stale_output (orchestrator-refreshed)

Both committed transcripts are stale (predate the sech/numeric-profile blocks, carry `STAGE 53`/`STAGE 053` banners). No source edit beyond F1/F2 — the re-run regenerates them with canonical `STAGE 70`/`STAGE 070` banners and the new assertion lines. Handled by the orchestrator's independent re-run.

## Applied: F1

- files_changed:
  - `scripts/moving_throat_pde_stage070_gnls_wall_shell_sympy_audit.py`
  - `mathematica/moving_throat_pde_stage070_gnls_wall_shell_mathematica_audit.wl`
- summary: Added load-bearing sech-profile moment assertions for `I_f = 2/3` and `I_g = 14/15`, and corrected the Mathematica `I_g` analytic annotation.
- deviation: none

## Applied: F2

- files_changed:
  - `scripts/moving_throat_pde_stage070_gnls_wall_shell_sympy_audit.py`
- summary: Updated the SymPy audit docstring self-labels from stage 53 to canonical stage 070/70.
- deviation: none

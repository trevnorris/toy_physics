---
unit_id: 203
batch: VI.1
created_at: 2026-06-01T00:00:00-06:00
findings_count: 1
stop_cold: null
applied: true
applied_at: 2026-06-02T12:20:31-06:00
findings_applied: 1
findings_blocked: 0
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 203

Apply the finding below. After applying, append an `## Applied: F1` block under the finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

If the required change is ambiguous or unsafe to apply mechanically, append `## Blocked: F1` with a question instead.

Do NOT introduce new features, refactors, or stylistic changes beyond what the finding requires. Do NOT touch paper.tex, notes/, or any prose documents.

After editing, RUN the affected scripts (`python3 <path>` for SymPy, `math -script <path>` for Mathematica) and iterate until they exit 0 with all in-file checks passing.

## F1 — insufficient_verification

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage203_free_quintuple_scalar_closure_slice_and_crossing_theorem_sympy_audit.py:273-330` (Section VI)
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage203_free_quintuple_scalar_closure_slice_and_crossing_theorem_mathematica_audit.wl:263-306` (Section VI)

**Issue:**
Section VI is meant to exercise the stage's headline deliverable — the one-parameter crossing theorem on `\widehat\chi_Q(\mathbf y):=\chi_Q(\mathbf x_*^{\rm graph}(\mathbf y))` (paper Output line; appendix eq:app-part06-chi-hat-def, eq:app-part06-IVT-crossing-condition/result). As written it does NOT use the actual object. It takes the carried Stage 197 closure `\chi_Q=3(S\beta^5+9\Sigma_5)/(3S-\Sigma_0)`, hardcodes `\Sigma_0=0,\Sigma_5=0` to collapse it to `\beta^5` (sympy:289–291, wl:272), and moves only the `\gamma` coordinate through `beta_path` while the other four free coordinates (`lam_bar, cetaU_bar, KU_bar, KW_bar`, sympy:278–280 / wl:267) are inert. Consequently: (a) the Stage 202 graph map built in §IV (`T_graph, Keta_graph, mu_graph` → `Ctr, Cnt, epsEta`) plays no role in the slice scalar, so a wrong dependent-coordinate graph exponent would not be caught by §VI; and (b) the appendix's load-bearing structural claim — "the three target monomials are invariant for all `\tau` on this lift" (eq:app-part06-graph-lifted-ray), which is the *reason* closure reduces to one scalar — is never verified anywhere in the unit.

**Required change (REQUIREMENT + ACCEPTANCE CRITERIA — you design the route):**

Within Section VI of BOTH scripts, exercise the crossing/slice scalar on the genuine graph lift rather than a hand-zeroed `\beta^5`. The new content must satisfy:

1. **Monomial-invariance on the lift (new check, both engines).** Construct an explicit free-quintuple path `\mathbf y(\tau)` that moves AT LEAST TWO of the five free coordinates (e.g. `\gamma` and one of `c_{\eta U}, K_U, \lambda_W, K_W`) — not a `\gamma`-only path. Lift it through the Stage 202 graph map (reuse the §IV constructions of `T_graph, Keta_graph, mu_graph` / `logTGraph, logKetaGraph, logMuGraph`, now as functions of the path), form the three target monomials `Ctr, Cnt, epsEta` on the lift, and assert that each equals its target value identically in `\tau` (equivalently `q_tr=q_nt=q_eta=0` along the path). This must be a genuine `expect_zero` / `expectZero` that would FAIL if any graph exponent were wrong. ACCEPTANCE: a new zero-residual check whose residual is built from the §IV graph monomials, not from the error exponents `E_T,E_K,E_mu`.

2. **Crossing scalar sourced from the graph composition.** Drive the IVT demonstration from `\chi_Q` evaluated through the graph lift (so any `\Sigma_0,\Sigma_5` reduction is shown to FOLLOW from the lift's monomial invariance rather than imposed by `subs(Sigma0->0, Sigma5->0)` by fiat). If you retain `\Sigma_0=\Sigma_5=0`, add a one-line in-script justification tying it to the verified monomial invariance, AND keep the residual-bearing invariance check from (1) as the load-bearing assertion. ACCEPTANCE: `\widehat\Delta_Q` used in the sign/root checks traces to the graph composition, and the section no longer relies solely on a bare `\beta^5` with inert free coordinates.

3. **Preserve the working IVT checks.** The existing sign checks (`\widehat\Delta_Q(\tau_-)<0`, `\widehat\Delta_Q(\tau_+)>0`), the root at the interior `\tau`, the positive-denominator continuity check, and the unique-real-crossing check may remain as-is (they are correct); just ensure they apply to the graph-sourced scalar and to the multi-coordinate path.

**Anti-transliteration guard:** The Mathematica monomial-invariance residual must be formed in native log-space using `Log` / `PowerExpand` / `Simplify` (consistent with the existing §IV `qtr/qnt/qeta` log-space construction at wl:211–222), NOT by transcribing the SymPy power-space `sp.log(Ctr/Ctr_tgt)` construction. Both engines must derive the invariance independently in their existing representations (SymPy: multiplicative powers then `log`; Mathematica: additive logs).

**Claim manifest** (results the augmented §VI must independently verify):
- M1: Target-monomial invariance on the graph lift — for `\mathbf x_*^{\rm graph}(\mathbf y(\tau))` along a path moving ≥2 free coordinates, `q_{\rm tr}=q_{\rm nt}=q_\eta=0` for all `\tau`, i.e. `\mathfrak C_{{\rm tr},*}=\mathfrak C_{{\rm tr},*}^{\rm target}`, `\mathfrak C_{{\rm nt},*}=\mathfrak C_{{\rm nt},*}^{\rm target}`, `\epsilon_\eta=\epsilon_{\eta,*}^{\rm target}`.
- M2: One-parameter crossing on the graph-sourced scalar `\widehat\Delta_Q(\tau)=\widehat\chi_Q(\mathbf y(\tau))-1`: sign change at the endpoints and a zero `\tau_*` in the interior with `\widehat\Delta_Q(\tau_*)=0` (IVT crossing theorem).

**Verification command:**
After applying, the verifier runs `redteam exec-sympy 203` and `redteam exec-mathematica 203` and confirms (i) the new monomial-invariance check appears in §VI of both scripts and PASSES, (ii) the crossing scalar is sourced from the graph composition over a multi-coordinate path, and (iii) both scripts exit 0.

## Applied: F1

- files_changed:
  - `scripts/moving_throat_pde_stage203_free_quintuple_scalar_closure_slice_and_crossing_theorem_sympy_audit.py`
  - `mathematica/moving_throat_pde_stage203_free_quintuple_scalar_closure_slice_and_crossing_theorem_mathematica_audit.wl`
- summary: Section VI now lifts a two-coordinate free path through the Stage 202 graph, verifies the three target monomial residuals vanish, and sources the IVT scalar from the graph-composed beta coordinate.
- deviation: none

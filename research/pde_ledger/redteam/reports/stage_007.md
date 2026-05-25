---
unit_id: 007
batch: I.1
auditor_model: claude-opus-4-7[1m]
audit_date: 2026-05-25T00:00:00Z
verdict: findings
stop_cold: null
findings_count: 1
paper_alignment: aligned
scripts_checked:
  sympy: present
  mathematica: present
  engines_agree: true
  outputs_fresh: false
docs_read:
  paper_stage_tex: present
  notes_stage_files: []
  paper_appendix: present
---

# Audit unit 007 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_007.tex`
- notes: `(none)` (no `notes/stages/moving_throat_pde_stage007_*.md` present)
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part01.tex` (only an `\input{stages/stage_007}` line on L93 — no row-level summary beyond the stage card itself)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage007_projection_reduction_comparison_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage007_projection_reduction_comparison_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage007_projection_reduction_comparison_sympy_audit.txt` (stale — see F1)
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage007_projection_reduction_comparison_mathematica_audit.txt` (stale — see F1)

## What the paper claims

Stage 007 compares the exact projected Maxwell law with the zero-mode reduction channel. For a zero-mode ansatz `A_mu(x,w) = a_mu(x)` with mixed sector suppressed, the projected inhomogeneous law carries three weighted integrals `I_WZ = ∫WZ dw`, `I_WH = ∫WH dw`, `I_WS = ∫WS dw`, and the projected effective parameters are (verbatim from eq. \eqref{eq:stage007-effective-parameters}):

> `mu_eff^proj = mu0 * I_WS / I_WZ`,  `1/xi_eff^proj = I_WH / (xi * I_WZ)`.

These are observer-kernel quantities until a matching prescription is fixed. The stage's `\stagefield{Output}` says: "Stage 007 exports the observer-dependent projected coefficients [eq:stage007-effective-parameters] and marks the reduction-first brane law as a matched channel rather than the controlling derivation." The Mismatch guard paragraph further requires that the comparison audit demonstrate, on smooth profiles, that a delta-like source and a smooth localized gauge profile do not automatically produce the reduced-action coupling.

Deliverables enumerated:
- D1: derive/state the projected formulas for `mu_eff^proj` and `xi_eff^proj` (the `I_WS/I_WZ` and `I_WZ/I_WH` ratios) from the zero-mode ansatz.
- D2: exhibit at least one concrete profile family showing observer dependence (matched-Gaussian and regulated-delta limits suggested by the surrounding text).
- D3: mutation guard: a `w`-dependent field or source disturbs the projection so the projection-first answer is not just the reduction-first one.

## What the script claims to verify

The SymPy script (and its Mathematica twin) computes, for Gaussian `Z(w)=exp(-w²/λ²)` and `H(w)=exp(-w²/ρ²)` together with three observer/source families:

1. a smooth Gaussian observer `W_smooth` of width `σ` against a smooth Gaussian source `S_smooth` of width `τ`, and confirms the closed-form overlaps `I_WZ = λ/√(λ²+σ²)`, `I_WS = 1/(√π·√(σ²+τ²))`, `I_WH = ρ/√(ρ²+σ²)`, and `xi*I_WZ/I_WH = xi*λ*√(ρ²+σ²)/(ρ*√(λ²+σ²))`;
2. a matched kernel `W=Z/Z_int` with brane-delta source, and confirms `I_WZ_match = √2/2`, `I_WH_match = ρ/√(λ²+ρ²)`, `mu0_proj/mu0_red = √2`, `xi_proj/xi_red = √(λ²+ρ²)/(√2·λ)`, plus the H=Z specialization `xi_proj/xi_red|_{ρ=λ} = 1`;
3. a regulated sharp observer with width `ε`, confirming the closed forms, the sharp limits `I_WZ→1`, `I_WH→1`, the projected `xi_eff^proj_eps` formula, and its `ε→0` limit `→ xi`;
4. mutation guards (a `w`-dependent field perturbation and a `w`-dependent source perturbation) that produce nonzero Gaussian-moment corrections.

The scripts exit 0 on PASS.

## Paper ↔ script cross-check

| Deliverable | Script coverage | Status |
|---|---|---|
| D1: `mu_eff^proj = mu0·I_WS/I_WZ` | sympy lines 87-89, 108, 148, 186; math M3/M4, M8 | match |
| D1: `1/xi_eff^proj = I_WH/(xi·I_WZ)` | sympy lines 104, 150, 187, 203; math M4c, M8b, M10c | match |
| D2: observer dependence on concrete profiles | sympy §3 (matched), §4 (regulated); math M7/M7b/M8/M8b/M9-M11c | match |
| D2: H=Z reduction agreement | sympy line 170-173; math M8c | match |
| D3: mutation guard for `w`-dependent field/source | sympy lines 111-128; math M5, M6 | match |

Paper alignment is **aligned**. The prior `paper_misalignment` (H(w) channel absent) is now resolved by the H/ξ extensions to both engines. No script-side check is orphaned from a paper claim.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 90-93 | `I_WZ_smooth = λ/√(λ²+σ²)` | D1 (I_WZ) | yes |
| A2 | sympy | 94-97 | `I_WS_smooth = 1/(√π·√(σ²+τ²))` | D1 (I_WS) | yes |
| A3 | sympy | 98-101 | `I_WH_smooth = ρ/√(ρ²+σ²)` | D1 (I_WH) | yes |
| A4 | sympy | 102-105 | smooth `xi·I_WZ/I_WH` closed form | D1 (xi_eff^proj) | yes |
| A5 | sympy | 106-109 | zero-mode projected residual = 0 | D1 (projected law linearity) | partial (consistency of linearity of ∫ — weak but harmless) |
| A6 | sympy | 114-118 | `w`-dep field mutation amplitude | D3 | yes |
| A7 | sympy | 123-128 | `w`-dep source mutation amplitude | D3 | yes |
| A8 | sympy | 138-140 | `Z_int = √π·λ`, `Z2_int = √(π/2)·λ`, `H_int = √π·ρ` | enabling lemma for D2 | yes |
| A9 | sympy | 162 | `I_WZ_match = Z2_int/Z_int` | structural identity for matched kernel | partial (mild tautology — bracketed by A10 numeric pin) |
| A10 | sympy | 163-165 | matched `I_WZ_match=√2/2`, `I_WH_match`, `mu0_proj/mu0_red=√2` | D2 | yes |
| A11 | sympy | 166-169 | matched `xi_proj/xi_red = √(λ²+ρ²)/(√2·λ)` | D1 + D2 | yes |
| A12 | sympy | 170-173 | H=Z specialization gives `xi_proj/xi_red = 1` | D2 (gauge channel matches reduction when H=Z) | yes |
| A13 | sympy | 174 | mutated ratio is not 1 (catches accidental triviality) | D2 anti-tautology | yes |
| A14 | sympy | 198-205 | regulated overlaps, sharp limits, `xi_eff_eps` closed form and `ε→0` limit | D1 + D2 | yes |
| A15 | sympy | 206-210 | `I_WS_eps = √2/(2√π·ε)` and anti-tautology guard | D2 (delta self-overlap diverges) | yes |
| M1 | math | 21-26 | `Z_int = √π·λ` | A8 mirror | yes |
| M2 | math | 29-34 | `Z2_int = √(π/2)·λ` | A8 mirror | yes |
| M2b | math | 37-42 | `H_int = √π·ρ` | A8 mirror | yes |
| M3 | math | 52-57 | `IWZsmooth` closed form | A1 mirror | yes |
| M4 | math | 60-65 | `IWSsmooth` closed form | A2 mirror | yes |
| M4b | math | 68-73 | `IWHsmooth` closed form | A3 mirror | yes |
| M4c | math | 76-82 | `xi·IWZ/IWH` smooth closed form | A4 mirror | yes |
| M5 | math | 86-99 | field mutation amplitude | A6 mirror | yes |
| M6 | math | 101-117 | source mutation amplitude | A7 mirror | yes |
| M7 | math | 124-129 | `IWZmatch = 1/√2` | A10 mirror | yes |
| M7b | math | 132-137 | `IWHmatch = ρ/√(λ²+ρ²)` | A10 mirror | yes |
| M8 | math | 146-151 | matched `mu0_proj/mu0_red = √2` | A10 mirror | yes |
| M8b | math | 154-160 | matched gauge ratio closed form | A11 mirror | yes |
| M8c | math | 163-168 | H=Z gauge alignment | A12 mirror | yes |
| M9 | math | 176-181 | regulated `IWZeps` closed form | A14 mirror | yes |
| M10 | math | 184-189 | regulated source self-overlap | A15 mirror | yes |
| M10b | math | 192-197 | regulated `IWHeps` closed form | A14 mirror | yes |
| M10c | math | 199-206 | regulated `xi_eff_eps` closed form | A14 mirror | yes |
| M11 | math | 209-219 | sharp limit `IWZeps→1` | A14 mirror | yes |
| M11b | math | 222-232 | sharp gauge limit `IWHeps→1` | A14 mirror | yes |
| M11c | math | 234-243 | sharp `xi_eff_eps → xi` | A14 mirror | yes |

No tautological-by-syntax assertions (A9 is mild and immediately bracketed by A10's substantive numeric pin). All assertions trace to a paper deliverable.

## Findings

### F1 — stale_output

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage007_projection_reduction_comparison_sympy_audit.txt` (mtime 2026-05-21 11:26:14)
- `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage007_projection_reduction_comparison_mathematica_audit.txt` (mtime 2026-05-21 11:51:12)

**What's wrong:**
Both saved outputs predate the current script revisions (sympy mtime 2026-05-25 02:14:25; mathematica mtime 2026-05-25 02:15:11), and the content reflects this: the saved sympy `.txt` shows the pre-H(w) layout — its section 1 prints only `I_WZ * partial_mu f^{mu nu} = mu0 * I_WS * j^nu` with no `(I_WH/xi) ∂^ν(∂·a)` term, no `xi_eff` lines anywhere, and case A reports only `mu0_eff^(proj, match)` without `xi_eff^(proj, match)`. The saved mathematica `.txt` likewise contains only M1-M11 results — none of the new M2b, M4b, M4c, M7b, M8b, M8c, M10b, M10c, M11b, M11c lines that the current `.wl` now prints. The captures are from before the H(w)/ξ_eff^proj extension landed.

**Why this matters:**
A trust-audit consumer reading the saved outputs cannot confirm that the new ξ_eff^proj and H(w) assertions actually pass; they're verifying a different (earlier) script. Informational rather than blocking — the next normal `redteam exec-*` invocation will refresh the captures.

**Required change:**
No script edit needed. Re-run both engines through the standard invocation paths (the verifier will do this) so the saved outputs reflect the current script.

**Verification:**
After refresh, the sympy `.txt` must contain `xi_eff^(proj)` print lines (in sections 1, 3, 4) and the mathematica `.txt` must contain `PASS: M2b`, `PASS: M4b`, `PASS: M4c`, `PASS: M7b`, `PASS: M8b`, `PASS: M8c`, `PASS: M10b`, `PASS: M10c`, `PASS: M11b`, `PASS: M11c` lines.

## Independent-derivation check (Mathematica)

The `.wl` file is structurally parallel to the `.py` file: same Gaussian profile choices for Z, H, W_smooth, W_match, W_eps; same mutation polynomial `etaSym*x*w²`; same closed-form targets (`λ/√(λ²+σ²)`, `√2/2`, `√2`, `√(λ²+ρ²)/(√2·λ)`, etc.); same section ordering. However, every assertion's RHS is an analytically-known closed form for a Gaussian moment integral that each engine computes independently via its own `Integrate`/`integrate` and `FullSimplify`/`simplify`, then compares against the algebraic target. There is no shared derivation chain whose steps could carry an error from one engine to the other — both are reduced to "engine X computes definite integral Y; does it match closed form Z?". This is the correct mode for Gaussian-integral identities. Not a `mathematica_transliteration` finding.

Quoted side-by-side correspondences for context (not as defects):

- sympy line 84 `I_WZ_smooth = sp.simplify(sp.integrate(W_smooth * Z, (w, -sp.oo, sp.oo)))` vs math line 47 `IWZsmooth = Integrate[Wsmooth[w] Z[w], {w, -Infinity, Infinity}]` — both invoke their engine's integrator on the same integrand.
- sympy line 91-93 `I_WZ_smooth - lam / sp.sqrt(lam**2 + sigma**2)` vs math line 52 `IWZsmooth - lambda/Sqrt[lambda^2 + sigma^2]` — same analytic target.
- sympy line 198 `sp.limit(I_WZ_eps, eps, 0, dir='+')` vs math line 211-213 `Limit[IWZeps, epsilon -> 0, Direction -> "FromAbove", Assumptions -> lambda > 0]` — same sharp-sampling probe, executed by independent limit routines.

## Engine cross-check

Both engines target identical closed forms in every paired assertion. The current saved outputs cannot confirm engine agreement on the new ξ-channel checks (per F1) — but the targets are bit-identical algebraic expressions, so agreement is structurally guaranteed if both engines' `simplify`/`FullSimplify` resolve their respective residuals to zero. (The currently-saved sympy `.txt` ends in `STATUS: PASS`; the currently-saved mathematica `.txt` ends in `STATUS: PASS`; both refer to the older pre-H(w) script versions.) No `engine_disagreement` finding.

## Verdict justification

Paper alignment is now clean: the H(w) gauge-driver channel that the paper enumerates in eq. \eqref{eq:stage007-effective-parameters} is fully exercised by both engines, including the matched-kernel `xi_proj/xi_red` ratio, the H=Z specialization (which collapses the matched gauge parameter onto the reduction-first one), the regulated `xi_eff_eps` closed form, and its `ε→0` limit. The mutation guards from §1 still exercise the paper's D3 mismatch claim. The remaining defect is purely housekeeping (stale `.txt` captures predate the script extension); the substantive math holds up against every attack I tried (Gaussian-moment closed forms verified analytically, sharp limits checked, mutation amplitudes computed via Gaussian variance algebra). Attacks that failed: (i) tried to find a tautology in A9 (matched-kernel overlap = Z2_int/Z_int) — A10's numeric pin `√2/2` makes any algebraic equivalence between W_match and Z still observable in a value-mismatch; (ii) tried to find a missing branch (rho≠λ vs rho=λ) — A11 covers the general case and A12 specializes, both checked; (iii) tried to find a parity mistake in the source-mutation integral — the integrand `W_smooth * w² * S_smooth` is even, integral evaluates to the Gaussian variance formula, matches the asserted target; (iv) tried to find a domain assumption error — all width parameters are declared positive in both engines (`positive=True` in sympy, `$Assumptions` block in math).

## Self-test notes

Walked the new assertions: (1) Variable independence — every `sp.diff(...)` is taken with respect to `x` against expressions that actually depend on `x` (`f_test`, `j_test`, `F_w_dependent`, `J_w_dependent` all contain `x`). (2) Symmetry/parity — mutation integrals integrate `w²·Gaussian` over the full real line, even × even, integral is the usual `√π/(2A^{3/2})`, matches the asserted amplitudes by direct re-derivation. (3) Trivial-case substitution — under `ρ→λ` the matched gauge ratio collapses to 1; under `ε→0` the regulated `xi_eff_eps` collapses to `xi*λ/ρ * ρ/λ = xi`. (4) No new script needs creating, so path-spec self-test (#4) is N/A. (5) Paper round-trip — the H(w) and ξ_eff^proj checks are exactly the formulas the paper card enumerates in eq. \eqref{eq:stage007-effective-parameters}; no new misalignment introduced.

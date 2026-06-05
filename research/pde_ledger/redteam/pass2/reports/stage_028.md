---
unit_id: 028
batch: II.1
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-06-04T00:00:00Z
verdict: clean
stop_cold: null
findings_count: 0
paper_alignment: aligned
scripts_checked:
  sympy: present
  mathematica: present
  engines_agree: true
  outputs_fresh: false
docs_read:
  paper_stage_tex: present
  notes_stage_files: [/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage028_loaded_profile_selection.md]
  paper_appendix: present
---

# Audit unit 028 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_028.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage028_loaded_profile_selection.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part02.tex` (row at line 46)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage028_loaded_profile_selection_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage028_loaded_profile_selection_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage028_loaded_profile_selection_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage028_loaded_profile_selection_mathematica_audit.txt`

## What the paper claims

Stage 028 turns the wall-profile angle from a hand-chosen parameter into an output of the first loaded axial eigenproblem. In the two-mode N/N basis the bare wall matrix is diagonal `K_bare = diag(K0, K1)` with `K0 = K_eta + 6 T_Omega`, `K1 = K0 + Delta K_ax`, `Delta K_ax = pi^2 T_w / L^2`, and the support/mixed branch enters as a rank-one attractive load along the D/N overlap vector `v = (kappa0, kappa1)^T` with `kappa0 = 2 sqrt2/pi`, `kappa1 = -4/(3 pi)`. `\stagefield{Output}`: "Stage~028 outputs the rank-one wall matrix \eqref{eq:app-stage028-Keff}, exact eigenvalues \eqref{eq:app-stage028-eigenvalues}, the selected angle law \eqref{eq:app-stage028-angle-law}, and the softening threshold \eqref{eq:app-stage028-alpha-crit}." Concretely the four boxed deliverables are: (1) `K_eff = K_bare - alpha v v^T`, alpha>0; (2) `lambda_± = (1/2)[K0+K1-alpha(kappa0^2+kappa1^2) ± sqrt((Delta K_ax + alpha(kappa0^2-kappa1^2))^2 + 4 alpha^2 kappa0^2 kappa1^2)]`; (3) `tan(2 theta_-) = 2 alpha kappa0 kappa1 / (Delta K_ax + alpha(kappa0^2-kappa1^2))`; (4) `alpha_crit = K0 K1 / (K1 kappa0^2 + K0 kappa1^2)`. The card and notes also make the sign argument (`theta_- < 0`, away from the positive blind angle) and the strong-loading limit (`theta -> theta_max`, `tan(theta_max) = -sqrt2/3`). The appendix row (line 46) summarizes the stage as "Rank-one loaded wall eigenproblem, lower-mode angle law, and softening threshold."

## What the script claims to verify

Both engines build `K_bare` and `v` from the literal overlap constants, form `K_eff = K_bare - alpha v v^T`, and verify: (a) the sign ledger `kappa0^2-kappa1^2 = 56/(9 pi^2) > 0` and `2 kappa0 kappa1 = -16 sqrt2/(3 pi^2) < 0`; (b) trace and determinant of `K_eff` equal the paper's symbolic forms; (c) the characteristic factorization `(x-lambda_-)(x-lambda_+) = x^2 - tr x + det` reproduces the boxed eigenvalues; (d) the energy-stationarity derivative `dE/dtheta` equals the numerator form of the boxed `tan(2 theta)` law, and the weak-loading leading coefficient is `kappa0 kappa1 / Delta K_ax`; (e) the strong-loading limit `alpha->oo` of `tan(2 theta)` equals `tan(2 theta_max)` with `tan(theta_max) = -sqrt2/3`; (f) the `det(K_eff)=0` solution equals the boxed `alpha_crit`, with the determinant strictly positive on the `alpha<alpha_crit` stable side. The Mathematica script additionally cross-checks against the built-in `Eigenvalues[kEff]` (sum-vs-trace, product-vs-determinant) and against a finite-throat closed form `alpha_crit = 9 pi^2 K0 K1/(8(11 K0 + 9 Delta K_ax))`.

## Paper ↔ script cross-check

| Paper deliverable | Script check | Status |
|---|---|---|
| `v = (kappa0,kappa1)`, `kappa0=2√2/π`, `kappa1=-4/(3π)` | py L60-61, wl L33-34 (`kappa0`, `kappa1`) | match |
| `K_bare = diag(K0,K1)`, `K0=K_eta+6T_Ω`, `K1=K0+ΔK_ax`, `ΔK_ax=π²T_w/L²` | py L62-64, wl L35-37,40 | match |
| `K_eff = K_bare - α v vᵀ` (eq Keff) | py L76, wl L41 | match |
| `λ_±` exact eigenvalues (eq eigenvalues) | py L112-119 char-factorization; wl L60-76 char-factorization + `Eigenvalues[kEff]` sum/product | match |
| `tan(2θ_-)` angle law (eq angle-law) | py L142-146 (dE/dθ vs stationarity numerator); wl L84-93 | match |
| sign: θ_-<0 (κ0>0, κ1<0) | py L152-153, wl L97-98 sign checks; sign ledger py L87-88, wl L47-48 | match |
| strong-loading θ→θ_max, tan(θ_max)=-√2/3 | py L171-180, wl L104-113 | match |
| `α_crit = K0K1/(K1κ0²+K0κ1²)` (eq alpha-crit) | py L196-200, wl L115-120; stable-side det>0 py L204-209, wl L126 | match |
| (weak-loading leading term, notes §4.1) | py L156-158, wl L100-102 | match |

Dominant pattern: all paper deliverables faithfully exercised → `paper_alignment: aligned`.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored? |
|---|---|---|---|---|---|
| A1 | sympy | 109 | `trace - expected == 0` | trace of K_eff (eq eigenvalues prereq) | yes |
| A2 | sympy | 110 | `det - expected == 0` | det of K_eff (alpha_crit prereq) | yes |
| A3 | sympy | 119 | char-factorization `== 0` | eigenvalues (claim 2) | yes |
| A4 | sympy | 146 | `dE/dθ - stationarity/2 == 0` | angle law (claim 3) | yes |
| A5 | sympy | 152-153 | `κ0-|κ0|==0`, `κ1+|κ1|==0` | sign of θ_- | yes |
| A6 | sympy | 158 | `weak coeff - κ0κ1/ΔK == 0` | weak-loading limit | yes |
| A7 | sympy | 179-180 | `lim - tan(2θ_max)==0`, `tan(θ_max)+√2/3==0` | strong-loading limit | yes |
| A8 | sympy | 200 | `α_crit - expected == 0` | alpha_crit (claim 4) | yes |
| A9 | sympy | 207 | `det(α_crit)==0` | softening threshold | yes |
| A10 | mma | 57-58 | trace/det `== expected` | trace/det | yes |
| A11 | mma | 67-74 | `Eigenvalues[kEff]` sum/product vs λ_± | eigenvalues (independent) | yes |
| A12 | mma | 76 | char-factorization `== 0` | eigenvalues (claim 2) | yes |
| A13 | mma | 93 | `dE/dθ - stationarity/2 == 0` | angle law (claim 3) | yes |
| A14 | mma | 97-98 | sign checks | sign of θ_- | yes |
| A15 | mma | 102 | weak coeff `== 0` | weak-loading limit | yes |
| A16 | mma | 112-113 | strong limit + tan(θ_max) | strong-loading limit | yes |
| A17 | mma | 116-119 | `Solve[detEff==0]` vs ratio form | alpha_crit (claim 4) | yes |
| A18 | mma | 123 | `α_crit - finite-throat closed form == 0` | alpha_crit re-expression (extra) | yes |
| A19 | mma | 124 | `det(α_crit)==0` | softening threshold | yes |

No tautological rows: each assertion compares a *computed* quantity (from `K_eff`, `Eigenvalues`, `D[energy,theta]`, `Solve`) against an *independently written* target form, so a transcription error in either side would surface as a nonzero residual.

## Findings

None. (One non-finding informational note on output freshness is recorded under Engine cross-check; it does not rise to a `stale_output` finding because content does not disagree — see below.)

## Independent-derivation check (Mathematica)

The `.wl` is NOT a line-by-line transliteration. It shares the overall choreography (build `K_eff`, check trace/det/eigenvalues/angle/alpha_crit), which is unavoidable for the same physics, but it adds genuinely independent verification routes the SymPy script lacks:
- `eigvalsDirect = Eigenvalues[kEff]` (wl L66) then `Total[eigvalsDirect] - (lambdaMinus+lambdaPlus)` (L67-70) and `Times @@ eigvalsDirect - lambdaMinus*lambdaPlus` (L71-74) — Mathematica computes the eigenvalues with its own kernel solver and checks them against the hand-written boxed forms. SymPy only does the algebraic characteristic factorization (py L118).
- `alphaCritClosed = 9*Pi^2*K0*K1/(8*(11*K0 + 9*deltaK))` (wl L121) with check at L123 — an independent finite-throat re-expression of `alpha_crit` that SymPy does not carry.
- Threshold solved via `Solve[detEff == 0, alpha]` with `ConditionalExpression` stripping (wl L115) — SymPy uses `sp.solve(sp.Eq(det_eff,0), alpha)[0]` (py L196).
This satisfies the second-engine policy: the `.wl` derives the eigenvalues independently (built-in solver) and adds an extra closed-form, rather than echoing SymPy's algebra. No `mathematica_transliteration` finding.

## Engine cross-check

Both engines agree on every shared result (comparing the two saved outputs):
- `kappa0^2-kappa1^2 = 56/(9 pi^2)` (py L52, wl L8); `2 kappa0 kappa1 = -16 sqrt2/(3 pi^2)` (py L53, wl L9).
- trace `2*Keta+12*TOmega-88*alpha/(9*pi^2)+pi^2*Tw/L^2` (py L54, wl L10).
- `tan(2 theta) = -48 sqrt2 L^2 alpha/(56 L^2 alpha + 9 pi^4 Tw)` (py L76, wl L26).
- weak coeff `-8 sqrt2 L^2/(3 pi^4 Tw)` (py L80, wl L32).
- strong limit `-6 sqrt2/7`, `tan(theta_max) = -sqrt2/3` (py L82-83, wl L35-36).
- `det(alpha_crit*(1-eps))` factors to the same positive form (py L93-98, wl L48). All in-file checks report `0`/`PASS`.

**Output freshness (informational, not a finding):** the saved `.txt` files (mtime 2026-05-21 17:07) predate the scripts (mtime 2026-06-03 15:59). However the only post-output change to either script was the cosmetic stage-number relabel in commit `e2a4780` (banner `STAGE 011 → STAGE 028` in the `.wl`; docstring `Stage 11 → Stage 28` in the `.py`) — verified via `git show e2a4780`: 2 lines changed per file, "no equation, value, variable, \label, or code logic changed." Consequently the `.wl` output line 3 still prints the old banner `STAGE 011 — LOADED PROFILE SELECTION`, and the `.py` filename/docstring strings still read `stage11`/`Stage 11` (py L2-3 of source). No mathematical content of the captured output disagrees with what the current script would produce. Per the prompt's `stale_output` guidance ("informational, not blocking, unless the output's content also disagrees"), and since the content does not disagree, this is recorded as a note rather than a finding. The orchestrator's independent re-run will refresh the banner string.

## Verdict justification

`clean`. I read the paper card, the notes, and the appendix row before the scripts, and the four boxed deliverables (K_eff, eigenvalues, angle law, alpha_crit) each map to a non-tautological script check in both engines, with values agreeing across SymPy and Mathematica and against the docs. Attacks tried and failed: (1) zero-derivative trap on `D[energy, theta]` — `energy` genuinely depends on `theta` via `q=(cos θ, sin θ)`, so the derivative is non-trivial and the stationarity check can fail; (2) tautology hunt — every assertion compares a computed object against an independently typed target, not against itself; (3) symbol-domain attack — `alpha>0, Tw>0, L>0` are all justified by the physical setup (positive loading, positive tension/length) and are exactly what the paper assumes; (4) transliteration — the `.wl` adds an independent `Eigenvalues[kEff]` route and an extra closed form, so it is not an echo of the `.py`; (5) alpha_crit algebra — hand-verified `K1 kappa0^2 + K0 kappa1^2 = (8/(9 pi^2))(11 K0 + 9 ΔK_ax)`, confirming the finite-throat closed form is a correct re-expression. The only blemish is benign output staleness limited to a cosmetic banner string.

## Self-test notes

Checked: (1) variable independence — `diff(E,theta)` is non-trivial because E depends on theta through q; assertion can fail. (2) Parity/symmetry — no unbounded integrals in this stage; N/A. (3) Trivial-case pre-check — at `alpha=0`, `K_eff=K_bare`, `tan(2θ)=0` (numerator vanishes), `lambda_±={K0,K1}`, and `alpha_crit` reduces to `K0K1/(K1κ0²+K0κ1²)` finite — all consistent with the boxed forms. (4) Path specs — N/A (no missing-script finding). (5) Paper round-trip — no fix prescribed, so no new misalignment risk. Conclusion: no findings; output staleness is cosmetic-only.

## Value Reconciliation (pass-2 augmentation)

Outputs are stale by mtime but content-current except the cosmetic banner string (see Engine cross-check). Reconciliation is based on the script source plus the committed saved outputs, which agree on all emitted values.

| value | source (py/wl + output line) | .tex/.md location | status |
|---|---|---|---|
| `kappa0 = 2√2/π` | py L60; out py L22-23, wl out L5 | tex L47 (`\kappa_0=\frac{2\sqrt2}{\pi}`); md L78 | MATCH |
| `kappa1 = -4/(3π)` | py L61; out py L26-28, wl out L5 | tex L48 (`\kappa_1=-\frac{4}{3\pi}`); md L80 | MATCH |
| `ΔK_ax = π²T_w/L²` | py L62; out py L33-36, wl out L6 | tex L38; md L70 | MATCH |
| `K0 = K_eta+6T_Ω` | py L63; out py L30, wl out L6 | tex L32; md L66 | MATCH |
| `K1 = K0+ΔK_ax` | py L64; out py L31-36, wl out L6 | tex L33-34; md L68 | MATCH |
| `K_eff = K_bare-α vvᵀ` | py L76; out py L37-47, wl out L7 | tex L56 (eq Keff); md L100,104-105 | MATCH |
| `κ0²-κ1² = 56/(9π²)` | py L87; out py L52, wl out L8 | md L137 | MATCH |
| `2κ0κ1 = -16√2/(3π²)` | py L88; out py L53, wl out L9 | md L141 | MATCH |
| trace(K_eff) | py L107; out py L54, wl out L10 | md L117 | MATCH |
| det(K_eff) | py L108; out py L55, wl out L11 | md L119 | MATCH |
| `λ_±` (boxed eigenvalues) | py L113-114; out py L59-74, wl out L22-23 | tex L64-68 (eq eigenvalues); md L123-126 | MATCH |
| `tan(2θ) = -48√2L²α/(56L²α+9π⁴T_w)` (= boxed angle law) | py L142; out py L76, wl out L26 | tex L74-75 (eq angle-law, symbolic form); md L132-133 | MATCH |
| weak-loading leading term `-8√2L²α/(3π⁴T_w)` | py L156; out py L80, wl out L32 | md L163 | MATCH |
| strong limit `tan(2θ)→-6√2/7` | py L171; out py L82, wl out L35 | (intermediate; tan(θ_max) form carried) — md §4.2 L173,187 | MATCH (via tan(θ_max)) |
| `tan(θ_max) = -√2/3` | py L172,180; out py L83,86, wl out L36 | md L187 (`tan(theta_max) = -sqrt(2)/3`) | MATCH |
| `α_crit = K0K1/(K1κ0²+K0κ1²)` | py L197; out py L89-90, wl out L43-45 | tex L84-85 (eq alpha-crit); md L240 | MATCH |
| `det(α_crit*(1-eps))` stable-side form | py L208-209; out py L93-98, wl out L48 | md §6 L246 (qualitative "stable for α<α_crit") | MATCH (qualitative) |

INTERNAL (scaffolding / re-expression, no prose home expected, no finding):
- `alphaCritClosed = 9π²K0K1/(8(11K0+9ΔK))` (wl L121) — algebraic re-expression of `alpha_crit`, verified equal to the boxed ratio form; an extra second-engine cross-check, not a stated deliverable.
- pass/fail flags, `expectZero` residuals (all `= 0`), `eps` sign helper, banner/section strings.

reconciliation: complete; 17 deliverable values checked, 0 misaligned.

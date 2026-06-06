---
unit_id: 117
batch: IV.3
created_at: 2026-06-06T00:00:00Z
findings_count: 1
stop_cold: null
applied: true
applied_at: 2026-06-06T17:14:08Z
findings_applied: 1
findings_blocked: 0
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 117

Apply the finding below. After applying, append an `## Applied: F1` block under it with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

You DESIGN the route. This directive states the requirement and acceptance criteria only; choose the actual independent construction yourself. Do NOT introduce unrelated features, refactors, or stylistic changes. Edit only the `.wl`.

After editing, RUN the script (`timeout 600 math -script <path>`) and iterate until it exits 0 with all in-file checks passing. A timeout (exit 124) is a FAILURE — reformulate, never raise the cap. The orchestrator independently re-runs afterward.

Do NOT touch paper.tex, notes/, prose documents, or any numbering label. Do NOT edit the SymPy `.py` (it is the independent reference engine and stays unchanged).

## Context

The IV.2 transliteration-watch band (105-175) is being re-scrutinized for line-by-line `.wl`↔`.py` ports. The orchestrator's ground-truth read of stage 117 confirms the `.wl` is a **full line-by-line transliteration of the `.py` across ALL SIX sections** — not only §5. Every section reuses the *same* `Series[…,{z,0,5}]`/`Normal` black-box on the *same* generating function as the `.py`, the *same* `Coefficient[…,z,n]` extraction, the *same* `Solve` matching equations with the *same* RHS literals, and verbatim `expectZero` label strings; §5 additionally re-types the reduced shell/mixed core closed forms `rC/rhoC/sigmaC/kappaC/gammaC`. The math and every emitted value are CORRECT — the only defect is that the Mathematica engine re-runs the SymPy computation step-for-step, so a shared error (a mis-expanded series, a mis-typed reduced form, a wrong `deltaCoreExpected`) would pass both engines undetected.

The IV.2 precedent is direct and is the route-class to follow:
- **107 / 110** (`Series[ratio]` on a rational generating function) were re-authored to an **order-by-order undetermined-coefficient solve** of the *defining* relation `Λ · Y = numerator` (`CoefficientList` + `Solve` for the `Y` coefficients), so the canonical-even coefficients are *derived* in-engine rather than read off a `Series` black-box.
- **114** (re-typed 2×2 reduced Schur forms) was re-authored to an **explicit scalar `Solve`-elimination of the 2×2 core system** (no matrix `Inverse`, no re-typed reduced literal).

Stage 117 needs BOTH treatments: the §1–§4 class checks get the 107/110 treatment; §5 gets the 114 treatment.

## F1 — mathematica_transliteration (FULL re-author)

**Target:** `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage117_outlet_core_status_mathematica_audit.wl`

**Requirement (what must become true):** No section of the `.wl` may obtain the canonical-even (z²→1/9, z⁴→4/81) / odd-norm (z⁵→χ) coefficients by `Series`/`Normal`-expanding the *same* generating function the `.py` expands. Each section must extract those coefficients by a construction that is genuinely distinct from the `.py`'s `sp.series(...).coeff(z,n)` black-box, so that a series/expansion bug or a shared mis-statement of a closed form can no longer pass both engines. Concretely:

1. **§1 Harmless scale/argument, §2 Pure Robin, §3 Standalone mixed-pole, §4 Hybrid split** — replace the `Normal[Series[<ratio>, {z,0,5}]]` + `Coefficient` extraction with an **undetermined-coefficient derivation**: posit `Y(z) = Sum[c[k] z^k, {k,0,5}]` (or the appropriate odd/even ansatz), impose the *defining* relation `denominator · Y == numerator` and solve the resulting order-by-order linear system (`CoefficientList`/`Solve`) for the `c[k]`; then read the canonical-even / odd-norm conditions off the derived `c[k]` and `Solve` for the class parameter (`beta`, `rho`, `kappa`/`sigma`, `{rho,kappa}`) exactly as now. The matching targets (1/9, 4/81), the branch-root assertions, and the `chi_Q = 1` / `1 - 9 σγ` odd-norm checks are UNCHANGED — only the route to the coefficients changes.
2. **§5 Concrete core realization** — derive the concrete-core correction `δΛ_core(z)` by an INDEPENDENT elimination of the 2×2 core system rather than re-typing the reduced `sigmaC` literal: build `coreMatrix = {{ks, lam}, {lam, -kq dW}}` with `dW = 1 - kappa0 z^2 - I gamma0 z^5`, `coreSource = {gs, gq}`, and `Solve` (or explicit Cramer `adj/det`) the linear system `coreMatrix . {s, q} == coreSource` for `{s, q}`; form the mouth feedback `δΛ_core = gs s + gq q` (the Schur complement), `Together`/`FullSimplify`. `rhoC`, `sigmaC`, `kappaC`, `gammaC` must now be CONSEQUENCES of this elimination, not re-typed closed forms. Then apply the SAME §5 surface substitutions already present (`gq -> First[gqSolutions]`, `kappa0 -> kappa0FromTube`, `gamma0 -> (1 + rC)/9`) and assert `Normal[Series[deltaCore - deltaCoreExpected, {z,0,5}]] === 0`.

KEEP unchanged: the `lambdaOut` fingerprint definition (it is the canonical INPUT, identical by necessity); the §5 `lWForward`/`kappa0FromTube` forward tube-length construction; `sigmaStar`; the `core-balance surface branches` / `both branches give the same sigma_*` / `core-balance sigma_* value` checks (their `sigmaC` must now come from the §5 elimination); `deltaCoreExpected`; every `expectZero` label string; the §6 boolean wiring (it reads the §1–§5 residuals — leave the wiring, it inherits the new derivations); the §5 `Print` provenance lines for κ₀ (forward) and γ₀ (ansatz).

**Acceptance criteria:**
- No `Series[...]`/`Normal[Series[...]]` of the §1–§4 class ratios remains; the class coefficients come from a `Solve` of `denominator·Y==numerator`. (`Series` of the FINAL `deltaCore - deltaCoreExpected` residual in §5 is fine — that is a residual check, not the coefficient-extraction black-box.)
- §5 `sigmaC` (and the §5 collapse) is an in-engine consequence of a `Solve`/elimination on the 2×2 core matrix, not a re-typed literal `(ks gq - lam gs)^2/(ks^2 kq (1 + rC))`.
- Every `expectZero` still prints `= 0` / `PASS`; the printed deliverables (`scale/argument solutions`, the four class matches, `core-balance surface branches`, `sigma_*`, `kappa0FromTube`, the `concrete core collapses…` residual, the classification rows, the survivor set) are byte-identical to the current transcript.
- The script exits 0 under `timeout 600`.

**Claim manifest** (the re-authored `.wl` must independently verify):
- M1: in each of §1–§4, the canonical-even match (z²→1/9, z⁴→4/81) selects exactly the current branch root(s) (`beta∈{-1,1}`; `rho=0`; `kappa=-1/9`, `sigma=0`; `{rho=0,kappa=0}` cancellation and `{kappa=1/3}` compensated), derived from the undetermined-coefficient system — not a `Series` read-off.
- M2: the Schur complement of `{{K_s,λ},{λ,−K_q D_W}}` against source `{g_s,g_q}` at `u=1` equals `ρ_c − σ_c/(1−κ_c z²−iγ_c z⁵)+O(z⁶)` with `ρ_c=g_s²/K_s`, `σ_c=(K_s g_q−λ g_s)²/(K_s² K_q (1+r_c))`, `κ_c=κ₀/(1+r_c)`, `γ_c=γ₀/(1+r_c)`, `r_c=λ²/(K_s K_q)` — derived, not re-typed.
- M3: on the balance branch `g_q` (`ρ_c=4σ_c`), `σ_c = σ_* = g_s²/(4K_s)`; and with `κ₀=4 L_W²/(π²a²)` from `L_W=πa√((1+r_c)/3)/2` and `γ₀=(1+r_c)/9`, `δΛ_core = 4σ_* − σ_*/(1−z²/3−iz⁵/9) + O(z⁶)`.

**Verification command:** the orchestrator runs `redteam exec-mathematica 117` and the verifier confirms: (a) no `Series`/`Normal` extraction of the §1–§4 class ratios survives (only the final §5 residual `Series` may remain); (b) §5 contains a `Solve`/`adj-det` elimination of the core system in place of the re-typed `sigmaC` literal; (c) intermediate variable names / steps differ from the `.py`; (d) every `expectZero` prints `= 0` / `PASS`; (e) all printed deliverables are byte-identical to the current committed transcript; (f) the script exits 0.

## Applied: F1

- files_changed:
  - `mathematica/moving_throat_pde_stage117_outlet_core_status_mathematica_audit.wl`
- summary: Replaced the §1–§4 ratio-series coefficient extraction with undetermined-coefficient solves and derived the §5 core correction from an explicit 2x2 core-system elimination.
- deviation: none

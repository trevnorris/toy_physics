---
unit_id: 201
batch: VI.1
created_at: 2026-06-01T00:00:00-06:00
findings_count: 1
stop_cold: null
applied: true
applied_at: 2026-06-02T00:36:26-06:00
findings_applied: 1
findings_blocked: 0
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 201

Apply the finding below. After applying, append an `## Applied: F1` block under it with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

If the required change is ambiguous or unsafe to apply mechanically, append `## Blocked: F1` with a question instead.

Do NOT introduce new features, refactors, or stylistic changes beyond what is specified.

After editing, RUN the new Mathematica script (`math -script <path>`) and iterate until it exits 0 with all in-file checks passing. Getting the script to run cleanly is your job; the orchestrator independently re-runs afterward.

Do NOT touch paper.tex, notes/, or any prose documents. The red-team only modifies scripts.

## F1 — missing_verification_script (subtype: missing_mathematica)

**Target:** `mathematica/moving_throat_pde_stage201_explicit_realization_compiler_and_canonical_orbit_projection_mathematica_audit.wl`

(`.wl` files live in `mathematica/`. Do not place this in `scripts/`. The corresponding output will be written to `mathematica/output/` by the runner.)

**Issue:**
Stage 201 is non-status-only and not a checkpoint, but only the SymPy engine exists; the paper card itself records "Mathematica audit: none yet." The stage is finite symbolic linear algebra plus log/exp identities, which Mathematica can verify independently. Add a second-engine `.wl` that re-derives and checks the results below. This is a REQUIREMENT + ACCEPTANCE-CRITERIA directive: you design and write the script. State each result as a passing symbolic check (use a helper that prints the residual and calls `Exit[1]` on any nonzero entry, analogous to the existing `expectZero` idiom in this repo).

**Anti-transliteration guard (mandatory):**
The `.wl` must derive each claim independently using native Mathematica primitives — `Inverse`, `LinearSolve`, `Solve`/`Reduce`, `D[]`, `Series`+`Coefficient`, `Simplify`/`FullSimplify`, `PowerExpand`, matrix ops — via a DIFFERENT decomposition than the `.py`. A line-by-line port of the SymPy algebra (same variable choreography, same `Edep . Inverse[Pdep]` construction, same intermediate names) is rejected as `mathematica_transliteration`. In particular: obtain the dependent-triple solution (M6) with `LinearSolve` / `Solve` directly on the linear system, NOT by reconstructing the section matrix `S` and multiplying; and verify the log-ratio identities (M2, M5) by explicit exponent arithmetic with `PowerExpand`/`Simplify` rather than mirroring the `expand_log(..., force=True)` path. The script must build `M_*` from the physical column conventions, not copy the SymPy matrix literal verbatim as the sole source of truth (you may state the same matrix, but at least one load-bearing result must come from an independent route, e.g. solve the linear system rather than invert the pivot block).

**Setup the script must establish (symbols, all positive/real as appropriate):**
- Parameters `chi0_star, deltaU_star, E_star, F_star` (positive reals).
- The exact Stage-192 quotient map `M_*` (3×8) on the ordered microscopic drift basis `(Δλ, Δc, Δγ, ΔU, ΔKη, ΔW, Δμ, ΔT)`:
  - row tr: `(0, 1+deltaU_star, 1+deltaU_star, -(2+chi0_star+deltaU_star), 0, 0, 0, 1+chi0_star)`
  - row nt: `(2(1+E_star), 0, 2 E_star, F_star-E_star, -1, -(2+E_star), 1, -F_star)`
  - row eta: `(0, 2, 0, -1, -1, 0, 0, 0)`
- Quotient vector `q = (q_tr, q_nt, q_eta)` with `q_* = Log[R_*]`, `R_tr, R_nt, R_eta` positive.
- The dependent triple is `(T, Kη, μ)` = drift indices `(8, 5, 7)` in 1-based ordering (columns T, Kη, μ of `M_*`).

**Claim manifest** (each must be an independent passing check):

- **M1 — Right-section identity / cancellation.** There is a section `S` (support only on rows `T, Kη, μ`) with `M_* . S == I_3`. Equivalently, the canonical repair drift `Δx_rep` (defined below) satisfies
  `M_* . Δx_rep == -q`.
  Verify `M_* . Δx_rep + q == (0,0,0)`.

- **M2 — Mismatch chart.** With `q_* = Log[R_*]`, verify the three exact identities:
  - `m_T == R_tr^(1/(1+chi0_star))`, where `m_T := Exp[q_tr/(1+chi0_star)]`;
  - `m_K == R_eta^(-1)`, where `m_K := Exp[-q_eta]`;
  - `m_mu == R_nt * R_eta^(-1) * R_tr^(F_star/(1+chi0_star))`, where `m_mu := Exp[q_nt - q_eta + F_star q_tr/(1+chi0_star)]`.
  (Establish via `PowerExpand`/`Simplify` of explicit exponents, not by log-expansion mirroring the .py.)

- **M3 — Repair vector closed form and support.** The canonical repair drift vector (8 components, in the drift basis order above) is exactly
  `Δx_rep = (0, 0, 0, 0, q_eta, 0, -q_nt + q_eta - F_star q_tr/(1+chi0_star), -q_tr/(1+chi0_star))`.
  Verify it has nonzero entries only on `(Kη=row5, μ=row7, T=row8)` — i.e. rows 1,2,3,4,6 are identically zero — and that the three nonzero entries equal the stated expressions. Confirm the support claim by showing rows {1,2,3,4,6} are zero symbolically.

- **M4 — Repair vanishes iff on target orbit.** Verify `Δx_rep == 0` (all 8 entries) under the substitution `R_tr -> 1, R_nt -> 1, R_eta -> 1` (equivalently `q_tr=q_nt=q_eta=0`); and verify at least one entry is nonzero when exactly one ratio is perturbed (e.g. `R_eta -> r` with `r != 1` gives a nonzero `Kη` entry). This pins the iff direction.

- **M5 — Canonical orbit projection.** Define the projected state
  `x_proj = (λ, c, γ, K_U, Kη·R_eta, K_W, μ·R_nt^(-1)·R_eta·R_tr^(-F_star/(1+chi0_star)), T·R_tr^(-1/(1+chi0_star)))`
  over a positive free state `x = (λ, c, γ, K_U, Kη, K_W, μ, T)`. Verify componentwise that
  `Log[x_proj[[i]]/x[[i]]] == Δx_rep[[i]]` for all i = 1..8.
  Also verify the fixed-point reduction: `x_proj == x` under `R_tr->1, R_nt->1, R_eta->1`.

- **M6 — Same-free-quintuple uniqueness (INDEPENDENT SOLVE).** Set the free quintuple drifts to zero and let only `(dT, dKη, dμ)` be free: `Δx_dep = (0,0,0,0, dKη, 0, dμ, dT)`. Using `Solve` (or `LinearSolve`) on the linear system `M_* . Δx_dep == -q`, obtain the UNIQUE solution and verify it equals
  `dT = -q_tr/(1+chi0_star)`, `dKη = q_eta`, `dμ = -q_nt + q_eta - F_star q_tr/(1+chi0_star)`.
  Confirm the solution is unique (the 3×3 dependent block is invertible — `Det != 0`). This must be solved directly, not by reusing the section `S` from M1.

- **M7 — Intrinsic packet equals pairwise-witness packet.** With `Ctr_star, Cnt_star, epsEta_star` (target values) and `Ctr2, Cnt2, eps2` (candidate values), and a pairwise witness `(Ctr1,Cnt1,eps1)`, verify that under `Ctr1->Ctr_star, Cnt1->Cnt_star, eps1->epsEta_star` the pairwise packet
  `(chi_Q-1, Log[Ctr2/Ctr1], Log[Cnt2/Cnt1], Log[eps2/eps1])`
  equals the intrinsic packet
  `(chi_Q-1, Log[Ctr2/Ctr_star], Log[Cnt2/Cnt_star], Log[eps2/epsEta_star])`.

- **M8 — First-order linearized compiler.** With a general drift `Δx = (Δλ,Δc,Δγ,ΔU,ΔKη,ΔW,Δμ,ΔT)`, set `q^lin = M_* . Δx` and `Δx_rep^lin = -S . q^lin` (where `S` is the right-section from M1, or equivalently obtain `Δx_rep^lin` by solving the dependent system for `q^lin`). Verify
  `M_* . (Δx + Δx_rep^lin) == (0,0,0)`.
  Note for self-test: `q^lin` depends on all eight drift symbols, so this residual is not identically zero by construction — it vanishes precisely because the section is a right-inverse on the quotient. Do NOT collapse `Δx` to a single symbol.

**Verification command:**
After Codex applies, the verifier will run `redteam exec-mathematica 201` and confirm each of M1–M8 appears as an explicit residual check AND the script exits 0.

## Applied: F1

- files_changed:
  - `mathematica/moving_throat_pde_stage201_explicit_realization_compiler_and_canonical_orbit_projection_mathematica_audit.wl`
- summary: Added the missing Mathematica audit with independent constrained-solve, direct dependent-solve, exponent-arithmetic, packet, projection, and linearized compiler checks for M1-M8.
- deviation: none

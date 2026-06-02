---
unit_id: 217
batch: VI.1
created_at: 2026-06-01T00:00:00-06:00
findings_count: 2
stop_cold: null
applied: true
applied_at: 2026-06-02T12:14:18-06:00
findings_applied: 2
findings_blocked: 0
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 217

Apply EVERY finding below in order (F1, F2). After applying, append an `## Applied: F<n>` block under that finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

F1 was a `paper_misalignment` (value_mismatch) that the USER HAS NOW RESOLVED with direction (a): the script's `162` (= 3⁴·2) is correct and arithmetically forced; the published card/appendix `179` and notes `230` are typos. You are EXPLICITLY AUTHORIZED to apply the F1 paper/notes edits in F1's `## RESOLVED` block (card `stage_217.tex`, appendix `stage_appendix_part06.tex`, the stage-217 notes). Codex applies; Claude reviews. Do NOT change the SymPy script's `162` — the script is correct. (Stage 218 also touches `stage_appendix_part06.tex`; this directive OWNS the part06 appendix edits and runs before 218.)

If a non-paper_misalignment finding's required change is ambiguous or unsafe to apply mechanically, append `## Blocked: F<n>` with a question instead — skip that finding, continue with the rest.

Do NOT introduce new features, refactors, or stylistic changes. Edit exactly the targets named.

After editing, RUN the affected script (`math -script <path>` for the new Mathematica file) and iterate until it exits 0 with all in-file checks passing. Getting the script to run cleanly is your job; the orchestrator independently re-runs afterward.

Do NOT touch paper.tex, notes/, or any prose documents.

## F1 — paper_misalignment

**Subtype:** value_mismatch

**Paper side:**
- `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_217.tex:13` quote: "proves the preferred \(179\)-candidate bound per envelope"
- `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part06.tex:1200` quote: `3^4\cdot2=179`
- `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part06.tex:65` quote: "preferred bound \(179\) per envelope"
- `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage217_full_interior_five_coordinate_simplex_optimizer_and_finite_candidate_reduction.md:409` quote: `3\cdot 3\cdot 3\cdot 3\cdot 2 = 230` (also lines 5, 32, 616)

**Script side:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage217_full_interior_five_coordinate_simplex_optimizer_and_finite_candidate_reduction_sympy_audit.py:138` quote: `if bezout_lift != 162:`

## RESOLVED (user direction: (a) — correct card/appendix/notes to 162; AUTHORIZED)

The user confirmed direction **(a)**: the script's `162` (= 3⁴·2) is correct and arithmetically forced — and it is the only value consistent with the existing `2×162=324` budget and `1140+324=1464` total. The published card/appendix `179` and notes `230` are typos. Codex applies these prose edits; Claude reviews. No script change.

**Authorized paper/notes edits (content-anchor each; the line numbers are hints from the 2026-06-01 audit and may have drifted):**
- `paper/stages/stage_217.tex` (~:13): "proves the preferred \(179\)-candidate bound per envelope" → "...preferred \(162\)-candidate bound per envelope".
- `paper/appendices/stage_appendix_part06.tex` (~:1200): `3^4\cdot2=179` → `3^4\cdot2=162`.
- `paper/appendices/stage_appendix_part06.tex` (~:65): "preferred bound \(179\) per envelope" → "preferred bound \(162\) per envelope".
- `notes/stages/moving_throat_pde_stage217_full_interior_five_coordinate_simplex_optimizer_and_finite_candidate_reduction.md` (lines ~5, 32, 409, 616): every `3\cdot3\cdot3\cdot3\cdot2 = 230` (and any prose "230"-per-envelope) → `... = 162`. Correct ALL four occurrences; grep the notes for `230` to be sure none is missed.

Do NOT change `2\times162=324`, `1140+324=1464`, or any other figure — those are already correct (they use 162). After the edits, grep card+appendix+notes to confirm no stray `179`/`230` per-envelope figure remains. F2 (`.wl`) is independent and applied as written.

## Applied: F1

- files_changed:
  - `paper/stages/stage_217.tex`
  - `paper/appendices/stage_appendix_part06.tex`
  - `notes/stages/moving_throat_pde_stage217_full_interior_five_coordinate_simplex_optimizer_and_finite_candidate_reduction.md`
- summary: Corrected the stage-217 preferred lifted candidate bound from the paper/notes typos to the arithmetically forced value 162.
- deviation: none

## F2 — missing_verification_script

**Subtype:** missing_mathematica

**Target:** `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage217_full_interior_five_coordinate_simplex_optimizer_and_finite_candidate_reduction_mathematica_audit.wl`

**Issue:** Stage 217 is non-status-only and non-checkpoint, yet it has only a SymPy engine. Mathematica can independently verify every computable deliverable of this stage (stationary-numerator identities, lifted/projected polynomial degrees and Bézout products, and the two special reductions). The dual-engine contract requires a second engine wherever independent verification is possible; it is possible here. Design and write a new `.wl` at the Target path that independently establishes the claim manifest below.

**Required change:**
Create the Target `.wl`. It must independently set up the same physical objects from their definitions and assert M1-M6. State the bound as the literal arithmetic `3^4*2` evaluated by Mathematica (which yields 162) and the projected product as `5*5*5*6` (= 750) — do NOT hardcode the paper's 179/230; assert what the degree products actually evaluate to, matching the SymPy script. Use the symbol set: `H0`, ratio variables `r,s,t,u` and the lift variable `y` (all positive reals where the script uses positive reals), free slopes `kL,kc,kg,kU,kW` (positive), and the 15 discriminant coefficients `A..O` as free reals. Define
- `S = 1 + r^2 + s^2 + t^2 + u^2`,
- `Klin = kL + r kc + s kg + t kU + u kW`,
- `Δ = A + B r + C s + D t + E u + F r^2 + G r s + H r t + I r u + J s^2 + K s t + L s u + M t^2 + N t u + O u^2`,
- `Mr = S kc - r Klin`, `Ms = S kg - s Klin`, `Mt = S kU - t Klin`, `Mu = S kW - u Klin`,
- `Lr = S D[Δ,r] - 2 r Δ`, and analogously `Ls, Lt, Lu`,
- `Φ = (Klin + Sqrt[Δ])/Sqrt[S]`.

**Anti-transliteration guard:** The `.wl` must derive the results INDEPENDENTLY using native Mathematica primitives — `D[]`, `Simplify`/`FullSimplify`/`Together`, `Solve`/`Reduce`, `Series`+`Coefficient`, `Exponent`, `MonomialList`/`CoefficientRules`, `FindRoot`, matrix ops — via a DIFFERENT decomposition than the `.py`. A line-by-line port of the SymPy algebra (same `expect_zero` choreography, same intermediate `simplify(S*kc - r*Klin)` ordering, same `Poly(...).total_degree()` call rewritten as a one-liner) is rejected as transliteration. Concretely, the `.wl` must NOT mirror the SymPy script's exact verification expression for M1; it must reach the stationary-numerator identity by an independent route (e.g. clear the `Sqrt[Δ]` and `Sqrt[S]` denominators of `D[Φ,r]` symbolically with `Together`/`Simplify` and compare the resulting numerator to `2 Mr Sqrt[Δ] + Lr`, OR rationalize `y -> Sqrt[Δ]` and verify `Numerator[Together[D[Φ,r]]]` matches up to the known positive factor), and it must extract total degrees with a native primitive (e.g. `Max[Total /@ CoefficientRules[expr,{r,s,t,u,y}][[All,1]]]`) rather than transcribing the SymPy degree call.

**Claim manifest** (the new `.wl` must independently verify each):

- **M1 — exact stationary-numerator identities.** For each `v ∈ {r,s,t,u}`, the partial of the objective numerator vanishes iff the lifted numerator vanishes, i.e. the scaled stationary identity holds:
  \[ 2\sqrt{\Delta}\,S^{3/2}\,\partial_v\Phi_\star \;=\; 2 M_v \sqrt{\Delta} + L_v \quad\text{identically in all symbols.} \]
  Assert the residual `2 Sqrt[Δ] S^(3/2) D[Φ,v] - (2 Mv Sqrt[Δ] + Lv)` simplifies to 0 for `v = r,s,t,u`. (Reach this by an independent denominator-clearing route per the guard, not by copying the `.py` expression.)

- **M2 — face count.** The unique primitive five-coordinate simplex `{λ,c,γ,U,W}` has exactly five codimension-one quadruple faces (each face = the full set minus one axis). Assert `Length` of the generated face list `== 5` and that each face has cardinality 4.

- **M3 — lifted degree pattern and Bézout product.** With `Fv = 2 Mv y + Lv` (`v=r,s,t,u`) and `FΔ = y^2 - Δ`, the total degrees in `(r,s,t,u,y)` are `(3,3,3,3,2)` and the product is `3*3*3*3*2` which evaluates to **162**. Assert the degree tuple equals `{3,3,3,3,2}` and the product equals `3^4*2` (= 162). Do NOT assert 179 or 230.

- **M4 — projected square-root-free degree pattern and bound.** Define `Crs = Ms Lr - Mr Ls`, `Crt = Mt Lr - Mr Lt`, `Cru = Mu Lr - Mr Lu` (quintics) and `Sr = Lr^2 - 4 Mr^2 Δ` (sextic). Assert their total degrees in `(r,s,t,u)` are `(5,5,5,6)` and the product equals `5*5*5*6` (= 750).

- **M5 — diagonal-isotropic gradient-optimal reduction.** Under the diagonal-isotropic substitution `A=kL^2-2 H0 ν, F=kc^2-2 H0 ν, J=kg^2-2 H0 ν, M=kU^2-2 H0 ν, O=kW^2-2 H0 ν` and the off-diagonals `B=2 kL kc, C=2 kL kg, D=2 kL kU, E=2 kL kW, G=2 kc kg, H=2 kc kU, I=2 kc kW, K=2 kg kU, L=2 kg kW, N=2 kU kW`: (a) `Lv(diag) - 2 Klin Mv(diag)` simplifies to 0 for `v=r,s,t,u`; and (b) at the gradient-optimal ratios `r=kc/kL, s=kg/kL, t=kU/kL, u=kW/kL`, each `Mv` and each `Lv` simplifies to 0. Assert all eight residuals vanish.

- **M6 — fivefold-symmetric equal-mix reduction.** Under `kL=kc=kg=kU=kW=k`, diagonal coefficients `A=F=J=M=O=k^2-2 H0 ν_d`, off-diagonals `B=C=D=E=G=H=I=K=L=N=2 k^2-4 H0 ν_o`: at the equal-mix barycenter `r=s=t=u=1`, each `Mv` and each `Lv` simplifies to 0 (`v=r,s,t,u`). Assert all eight residuals vanish.

Each assertion must hard-fail the script (nonzero `Exit[1]`) on a nonzero residual or a wrong degree/product, mirroring the SymPy `raise`/`expect_zero` discipline but via native Mathematica `If[... =!= 0, Print[...]; Exit[1]]` (or equivalent) checks.

**Verification command:**
After Codex applies, the verifier will run `redteam exec-mathematica 217` and confirm the new `.wl` exists, asserts M1-M6, exits 0, and reports the lifted product as 162 (agreeing with the SymPy script) and the projected product as 750.

## Applied: F2

- files_changed:
  - `mathematica/moving_throat_pde_stage217_full_interior_five_coordinate_simplex_optimizer_and_finite_candidate_reduction_mathematica_audit.wl`
- summary: Added an independent Mathematica audit for M1-M6 covering stationary identities, face count, degree products, and the two special reductions.
- deviation: none

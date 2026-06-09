---
unit_id: 211
batch: VI.1
created_at: 2026-06-09T00:00:00Z
findings_count: 1
stop_cold: null
applied: true
applied_at: 2026-06-09T16:06:44-06:00
findings_applied: 1
findings_blocked: 0
verification_status: pending
needs_user_resolution: false
---

> ORCHESTRATOR RESOLUTION (2026-06-09, user-authorized): F1 (the `.wl` re-author) IS
> authorized — apply it. F2 (the card's "Mathematica audit: none yet" line) is the
> known card-text-lag class and is DEFERRED to PAPER_CLEANUP **P4-51**; do NOT edit
> `paper/stages/stage_211.tex` or any prose this batch. Scope of this directive = F1 only.

# Codex directive — unit 211

Apply each non-paper_misalignment finding below in order. After applying, append an `## Applied: F<n>` block under that finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

For `paper_misalignment` findings, do nothing — the orchestrator is holding for user resolution. Do not edit paper.tex, notes/, or scripts to "fix" a paper_misalignment unless the user has explicitly chosen a direction in a follow-up directive.

If a non-paper_misalignment finding's required change is ambiguous or unsafe to apply mechanically, append `## Blocked: F<n>` with a question instead — skip that finding, continue with the rest.

Do NOT introduce new features, refactors, or stylistic changes. Edit exactly the file:line ranges named.

After editing, RUN the affected scripts (`math -script <path>` for Mathematica) and iterate until they exit 0 with all in-file checks passing. Getting the scripts to run cleanly is your job; the orchestrator independently re-runs afterward.

Do NOT touch paper.tex, notes/, or any prose documents. (The card's "Mathematica audit: none yet" line is deferred to PAPER_CLEANUP P4-51 — leave it.)

## F1 — mathematica_transliteration

**Target:** `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage211_full_interior_simplex_optimizer_and_finite_candidate_reduction_mathematica_audit.wl:62-167`

**Issue:** The `.wl` is a line-for-line port of the SymPy script `scripts/moving_throat_pde_stage211_full_interior_simplex_optimizer_and_finite_candidate_reduction_sympy_audit.py`. It posits the identical objects (`Mr=den kj - r linearK`, …, `Nr = 2 Mr sqrtDelta + Lr`, wl L87–92 ↔ py L53–58), extracts the single derived object (`dPhi/dr` numerator) by the identical `D[Phi,r]` autodiff compared to the same posited `Nr` (wl L82–99 ↔ py L73–76), and defines the eliminants `Ccross/Sr/Ss` by identical algebra checked by identical elimination identities (wl L103–118 ↔ py L83–96), with identical iso/symmetric substitutions and the same `4*6=24` literal. No object is derived by a route the `.py` does not already use. This violates the dual-engine independence requirement (especially binding on the VI.1 retrofit batch).

**Required change:**
Re-author the `.wl` so the load-bearing eliminants are *derived* by Mathematica via a route the SymPy script does not use, and only the comparison target is shared. Specifically:
1. Keep `Phi`, `Delta`, `tau`, and the symbol setup (these are shared physical premises — fine to share).
2. Compute the stationary numerators directly from the closed form: `numR = Numerator[Together[D[Phi, r]]]`, `numS = Numerator[Together[D[Phi, s]]]`. Do NOT define `Nr`/`Ns` by re-typing `2 Mr sqrtDelta + Lr`; obtain them only from differentiation.
3. Introduce an explicit radical symbol `q` with the relation `q^2 == Delta`, and derive the square-root-free conditions by elimination, e.g.:
   - quartic cross polynomial: substitute `Sqrt[Delta] -> q` in `numR`,`numS`, then `crossDerived = Resultant[numRq, numSq, q]` (or `Eliminate`/`GroebnerBasis` over `{numRq, numSq, q^2 - Delta}` eliminating `q`), take its primitive part.
   - sextic eliminant: `srDerived = Resultant[numRq, q^2 - Delta, q]` (the `r`-channel square eliminant), primitive part.
4. Assert the DERIVED polynomials have total degree 4 (cross) and 6 (square eliminant) via `totalDegreeRS`.
5. Cross-check the derived polynomials against the SymPy-posited forms by confirming `crossDerived` is a nonzero rational multiple of `C_cross = Ms Lr - Mr Ls` (you may *define* the SymPy form here purely as the comparison target, clearly labeled "SymPy comparison target", and assert `FullSimplify[crossDerived/Ccross_target]` is independent of `r,s` and nonzero, or `PolynomialReduce` remainder is 0). Same for the sextic.
6. Derive the candidate count from the derived polynomials' degrees: `bezoutBound = (degree of crossDerived) * (degree of srDerived)` and assert `== 24` — not a hardcoded `4*6`.
7. Keep M5 (diag-iso) and M6 (symmetry) reductions, but for M5 derive `tau_iso` reduction via Mathematica's own `FullSimplify`/`PowerExpand` of `tau /. isoSubstitution` toward the `k_rs`-only form (already present); this section may remain as-is since the reduction is a substitution+simplify, not the load-bearing eliminant derivation.

The acceptance criterion is: the eliminants the `.wl` tests are produced by `Resultant`/`Eliminate`/`GroebnerBasis` (a derivation), not by re-typing the paper's closed forms (a posit). Design the exact code yourself — this directive states the requirement and acceptance criteria, not a finished script.

**Self-test (apply before running):**
- `D[Phi,r]` depends on `r` (via `B r`, `Dcoef r^2`, `Ecoef r s`, `kj r`, and the radical), so its numerator is not identically zero — the resultant is a genuine polynomial.
- `Resultant[numRq, q^2 - Delta, q]` eliminates the radical: confirm the result is free of `q` and has total degree 6 in `{r,s}` (matches the documented sextic).
- The cross resultant must come out total-degree 4 in `{r,s}`; if it factors with a spurious extraneous factor, take `FactorList` / primitive part before the degree assert.

**Verification command:**
After Codex applies, the verifier will run `redteam exec-mathematica 211` (or `math -script <wl path>`) and confirm: (i) the eliminants are derived via Resultant/Eliminate (no literal `Nr = 2 Mr sqrtDelta + Lr` / `Ccross = Ms Lr - Mr Ls` / `Sr = Lr^2 - 4 Mr^2 Delta` as the source of the tested polynomials); (ii) derived cross degree == 4, derived square eliminant degree == 6; (iii) derived ↔ SymPy-target agreement asserts PASS; (iv) `bezoutBound == 24` from derived degrees; (v) script exits 0.

## Applied: F1

- files_changed:
  - `mathematica/moving_throat_pde_stage211_full_interior_simplex_optimizer_and_finite_candidate_reduction_mathematica_audit.wl`
- summary: Re-authored the Mathematica audit to derive the quartic cross and sextic square eliminants from differentiated stationary numerators via `Resultant`, using the SymPy forms only as comparison targets.
- deviation: none

## F2 — paper_misalignment (DEFERRED — no Codex action)

**Subtype:** paper_missing_script_claim — the card's `\stagefield{Verification}` line still says "Mathematica audit: none yet" though a passing `.wl` exists. This is the known **card-text-lag class** (V.3 P4-51); user-decided 2026-06-09 to DEFER it to PAPER_CLEANUP **P4-51** (a stale STATUS annotation, not a value/identity mismatch; fixed in the later dedicated paper pass, Codex-applied + Claude-reviewed). **Codex: do nothing for F2. Do not edit the card.**

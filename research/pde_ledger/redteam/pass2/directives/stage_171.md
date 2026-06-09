---
unit_id: 171
batch: V.1
created_at: 2026-06-08T00:00:00Z
findings_count: 1
stop_cold: null
applied: true
applied_at: 2026-06-08T22:48:31Z
findings_applied: 1
findings_blocked: 0
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 171

Apply the finding below. After applying, append an `## Applied: F1` block under the finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

If the required change is ambiguous or unsafe to apply mechanically, append `## Blocked: F1` with a question instead.

Do NOT introduce new features, refactors, or stylistic changes beyond what is named. Do NOT touch paper.tex, notes/, or any prose documents — only the `.wl`.

After editing, RUN the Mathematica script (`math -script /var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage171_microscopic_grouped_obstructions_mathematica_audit.wl`) and iterate until it exits 0 with all in-file checks passing.

## F1 — mathematica_transliteration

**Target:** `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage171_microscopic_grouped_obstructions_mathematica_audit.wl`

**Issue:** The `.wl` is largely a line-by-line port of the SymPy `.py` (identical `D[expr,var]*dvar` summation choreography, identical typed `*Formula` targets, identical print strings). It already carries two genuinely independent re-derivations — the `Series`-in-`eps2` routes for the `Z` and `N` bundles (wl:127–130, 140–143) — but the BdG bundle and the weak-axisymmetric K/G checks remain pure ports of the `.py`'s `diff`-summation. Extend the independent-mechanism coverage to those two remaining sections, following the `Series`-route pattern already in the file.

**Required change:**
1. After the existing `BdG obstruction bundle` check (currently ending at wl:58), add a second route that builds the same first-order combination `(b2 + b0/9)` by Taylor expansion in a perturbation parameter, distinct from the `diff`-summation route already present:
   ```
   bCombSeries = Coefficient[Normal[Series[
     (b2 + b0/9) /. {c -> c + eps2*dc, w -> w + eps2*dw},
     {eps2, 0, 1}]], eps2];
   expectZero["BdG obstruction bundle (series route)", bCombSeries - bCombFormula];
   ```
   (Place it immediately after the existing `expectZero["BdG obstruction bundle", ...]` at wl:58. `b0`, `b2`, `bCombFormula` are already defined; `eps2` is a fresh local symbol — `Clear[eps2]` first if needed so a prior binding cannot leak.)
2. For the weak-axisymmetric section (wl:148–163), add — inside or after the `Do` loop — a route that derives `(𝔎₁,𝔊₁)` by substituting the Section-1 split structure rather than re-typing the collected `kMicro`/`gMicro` forms. Specifically, define once (outside the loop):
   ```
   kScalar = (k1 - b01 - z01)/9 + (-m1 - b21 - z21);
   gScalar = n01 - p0*(k1 - b01 - z01);
   ```
   and inside the loop add:
   ```
   expectZero[
     "weak-axisymmetric K scalar route lambda=" <> ToString[InputForm[lamVal]],
     eps*lamVal*kScalar - kRebuilt
   ];
   expectZero[
     "weak-axisymmetric G scalar route lambda=" <> ToString[InputForm[lamVal]],
     eps*lamVal*gScalar - gRebuilt
   ];
   ```
   This checks the rebuilt forms against the Section-1 split-derived scalar amplitudes, a route distinct from the literally-typed `kMicro`/`gMicro`. (`kRebuilt`/`gRebuilt` are already defined in the loop body.)
3. Do not remove any existing checks. The new checks are additive.

**Verification command:**
After Codex applies, the verifier will run `redteam exec-mathematica 171` and confirm the new PASS lines appear (`BdG obstruction bundle (series route)`, `weak-axisymmetric K scalar route lambda=...`, `weak-axisymmetric G scalar route lambda=...`) AND the script exits 0.

## Self-test notes

- `Series[..., {eps2,0,1}]` then `Coefficient[..., eps2]` reproduces the first variation by a mechanism distinct from `D[]`-summation; comparing against `bCombFormula` (the paper-typed target) is non-tautological — a wrong typed coefficient leaves a nonzero residual.
- `b2 + b0/9` depends on both `c` and `w`, so both substituted perturbations contribute; neither `δc` nor `δw` term is identically zero.
- The weak-axisymmetric `kScalar`/`gScalar` route compares the Section-1 split structure against the loop's `kRebuilt`/`gRebuilt`; since both are built from the same `(k1,m1,b01,b21,z01,z21,n01,p0)` symbols but via different groupings, the residual vanishes only if the split is algebraically correct.
- No paper-side constant is introduced; no `paper_misalignment` created.

## Applied: F1

- files_changed:
  - `mathematica/moving_throat_pde_stage171_microscopic_grouped_obstructions_mathematica_audit.wl`
- summary: Added the BdG series-route check and weak-axisymmetric scalar-route checks requested by the directive.
- deviation: none

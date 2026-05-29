---
unit_id: 171
batch: V.1
created_at: 2026-05-28T00:00:00-06:00
findings_count: 2
stop_cold: null
applied: false
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 171

Apply each finding below in order. After applying, append an `## Applied: F<n>` block under that finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

If a finding's required change is ambiguous or unsafe to apply mechanically, append `## Blocked: F<n>` with a question instead — skip that finding, continue with the rest.

Do NOT introduce new features, refactors, or stylistic changes. Edit exactly the file:line ranges named.

After editing, RUN the affected scripts (`python3 <path>` for SymPy, `math -script <path>` for Mathematica) and iterate until they exit 0 with all in-file checks passing. Getting the scripts to run cleanly is your job; the orchestrator independently re-runs afterward.

Do NOT touch paper.tex, notes/, or any prose documents.

## F1 — mathematica_transliteration

**Target:** `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage171_microscopic_grouped_obstructions_mathematica_audit.wl:113-126`

**Issue:** The `.wl` is a line-by-line transliteration of the SymPy script. The most factor-heavy comparison targets — `zCombFormula` (wl 113-119) and `nCombFormula` (wl 122-126) — are copied as literals identical to the SymPy file rather than independently reconstructed by Mathematica's own differentiation. A wrong hand-written target copied into both engines would pass in both, defeating the second-engine cross-check. The differentiation layer (`dZ2Exact`, `dN0Exact`, etc., wl 99-101) is already independent; only the bundle-target literals echo SymPy.

**Required change:**
Reconstruct the two bundle targets from the engine's own already-computed exact variations, so the second engine builds the bundle independently instead of restating the literal.

Before (wl 113-119):
```
zCombExact = Expand[dZ2Exact + dZ0Exact/9];
zCombFormula = Expand[
  (s/delta^2 + 1/(9*delta))*dQ +
  (q/delta^2)*dS -
  dG/delta +
  (gSym/delta^2 - q/(9*delta^2) - 2*q*s/delta^3)*dDelta
];
```
After:
```
zCombExact = Expand[dZ2Exact + dZ0Exact/9];
(* Independent reconstruction from the engine's own exact variations, *)
(* not a copy of the hand-written SymPy literal.                      *)
dZ2Indep = D[z2, q]*dQ + D[z2, s]*dS + D[z2, gSym]*dG + D[z2, delta]*dDelta;
dZ0Indep = D[z0, q]*dQ + D[z0, delta]*dDelta;
zCombFormula = Expand[dZ2Indep + dZ0Indep/9];
```

Before (wl 122-126):
```
nCombExact = Expand[dN0Exact + p0*dZ0Exact];
nCombFormula = Expand[
  (p0/delta)*dQ + (2*p/delta^2)*dP -
  (p0*q/delta^2 + 2*p^2/delta^3)*dDelta
];
```
After:
```
nCombExact = Expand[dN0Exact + p0*dZ0Exact];
(* Independent reconstruction from the engine's own exact variations. *)
dN0Indep = D[n0Expr, p]*dP + D[n0Expr, delta]*dDelta;
dZ0Indep2 = D[z0, q]*dQ + D[z0, delta]*dDelta;
nCombFormula = Expand[dN0Indep + p0*dZ0Indep2];
```

Note: the `expectZero["Z obstruction bundle", zCombExact - zCombFormula]` (wl 120) and `expectZero["N obstruction bundle", nCombExact - nCombFormula]` (wl 127) lines stay as-is. The difference `zCombExact - zCombFormula` now compares two engine-derived forms (one assembled before the `1/9` weighting, one after) and must still simplify to 0.

**Claim manifest:** n/a (both scripts already exist; this restores second-engine independence for two bundle targets).

**Verification command:**
After Codex applies, the verifier runs `redteam exec-mathematica 171` and confirms the `.wl` exits 0 and that `mathematica/output/...:39-42` still print `Z obstruction bundle = 0` / `PASS` and `N obstruction bundle = 0` / `PASS`.

## Applied: F1

files_changed: mathematica/moving_throat_pde_stage171_microscopic_grouped_obstructions_mathematica_audit.wl
summary: Replaced the two hand-copied SymPy bundle-target literals (zCombFormula, nCombFormula) with independent reconstructions built from the engine's own D[...] exact variations.
deviation: none

## F2 — insufficient_verification (mislabeled unit)

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage171_microscopic_grouped_obstructions_sympy_audit.py:27`
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage171_microscopic_grouped_obstructions_sympy_audit.py:5`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage171_microscopic_grouped_obstructions_mathematica_audit.wl:26`

**Issue:** Both scripts print the banner `"STAGE 154 — MICROSCOPIC GROUPED OUTLET OBSTRUCTIONS"` (sympy 27, wl 26) for what is Stage 171, and the SymPy docstring (line 5) says "from Stage 170." The wrong label is captured into both transcripts (`.txt:11`), breaking traceability between the saved output and the audited unit. No assertion changes.

**Required change:**
1. sympy line 27: change `banner("STAGE 154 — MICROSCOPIC GROUPED OUTLET OBSTRUCTIONS")` to `banner("STAGE 171 — MICROSCOPIC GROUPED OUTLET OBSTRUCTIONS")`.
2. wl line 26: change `banner["STAGE 154 — MICROSCOPIC GROUPED OUTLET OBSTRUCTIONS"];` to `banner["STAGE 171 — MICROSCOPIC GROUPED OUTLET OBSTRUCTIONS"];`.
3. sympy line 5: change `obstructions from Stage 170.` to `obstructions for Stage 171.` (drop the upstream-stage reference, keep the docstring otherwise unchanged).

**Verification command:**
After Codex applies, the verifier runs `redteam exec-sympy 171` and `redteam exec-mathematica 171`; both `output/...:11` lines must print `STAGE 171 — MICROSCOPIC GROUPED OUTLET OBSTRUCTIONS`, both scripts exit 0, and all `= 0` / `PASS` lines are unchanged.

## Applied: F2

files_changed: scripts/moving_throat_pde_stage171_microscopic_grouped_obstructions_sympy_audit.py, mathematica/moving_throat_pde_stage171_microscopic_grouped_obstructions_mathematica_audit.wl
summary: Corrected the mislabeled "STAGE 154" banners to "STAGE 171" in both scripts and changed the SymPy docstring's "from Stage 170." to "for Stage 171.".
deviation: none

## Orchestrator rework (F1 — directive route was tautological)
The directive's prescribed F1 fix rebuilt `zCombFormula`/`nCombFormula` from `D[z2,...]`/`D[z0,...]`/`D[n0Expr,...]` — the SAME calls that build `zCombExact`/`nCombExact` — making the bundle checks `X - X` (zero by construction, regardless of correctness). The first verification (needs_rework) caught this. Orchestrator replaced it with: (a) `zCombFormula`/`nCombFormula` = the paper's collected closed-form bundle target, native-typed, so `zCombExact` (the engine's own D[]-derived total differential) vs the collected form is load-bearing on the collected coefficients; PLUS (b) an independent second route per bundle — `zCombSeries`/`nCombSeries` extracted by `Series`-coefficient in a perturbation parameter (a distinct mechanism from the summed total differential), compared against the same target. Two new PASS lines per bundle. Re-run: 23 PASS, 0 FAIL, exit 0. F1 now resolved (load-bearing + genuinely independent, not tautological).

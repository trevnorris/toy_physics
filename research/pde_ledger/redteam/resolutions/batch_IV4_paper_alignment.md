---
batch: IV.4
created_at: 2026-05-27
clusters_resolved: 3
mechanical_groups: 6
status: directions_set
---

# Red-team batch IV.4 — paper-alignment resolution

12 audits landed (range 127-138). 22 findings raised. Three user-gated clusters were resolved on 2026-05-27 as `(Recommended)`; the remaining six mechanical groups apply without a gate.

## Cluster A — Mass renumbering (mechanical, 19 edits)

**Pattern:** Two uniform offsets surfaced across the batch.

- **Mathematica banners** in 9 `.wl` files have `STAGE N` where N = actual − 17:
  127→110, 129→112, 130→113, 131→114, 133→116, 134→117, 135→118, 137→120, 138→121.
- **Notes H1** in 10 files (127-136) have `# Moving-Throat PDE — Stage M` where M = actual + 102:
  127→229, 128→230, 129→231, 130→232, 131→233, 132→234, 133→235, 134→236, 135→237, 136→238.
  Files 137 and 138 have correct H1 (Stage 137, Stage 138).

Canonical numbering elsewhere (paper card filenames, appendix rows, MANIFEST, file paths) is 127-138.

**Direction:** mass-fix both groups in place.

- 9 `.wl` banner edits: `banner["STAGE <N-17> — …"]` → `banner["STAGE <N> — …"]`.
- 10 notes H1 edits: `# Moving-Throat PDE — Stage <N+102>: …` → `# Moving-Throat PDE — Stage <N>: …`.

No script content other than the banner string changes; no notes content other than the H1 line.

## Cluster B — Status-only carry-forward re-attribution (132, 136)

**Pattern:** Two status-only stages have notes that attribute load-bearing constants to *downstream* stages, breaking the carry-forward chain:

- Stage 132's notes cite "Stages 180-182" as upstream sources for `Π_*` and `g_Π`.
- Stage 136's notes cite "Stages 184-186" as upstream sources for `(M_s, M_q)`, `Σ_m^*`.

Both cited ranges lie numerically *after* their citing stage (132 < 180, 136 < 184), so they cannot be upstream.

The constants themselves actually originate within batch IV.4:
- `Π_* ≈ 1.50882951349316` is fixed at stage 131 (parent micro-threshold root).
- `g_Π` closed form is fixed at stage 130 (mouth-bias map).
- `g_-^{F1} ≈ 0.758035…` is fixed by upstream IV.3 stage 122.
- `(M_s, M_q)` reduction laws are fixed at stages 133-135.
- `Σ_m^* ≈ 0.451485277739090` is fixed at stage 135 (outlet-consistent closure).

**Direction:** edit notes for 132 and 136 to cite the actual IV.4 upstream stages.

- 132 notes: replace "Stages 180-182" attributions with "Stages 130-131" (for `g_Π` and `Π_*`).
- 136 notes: replace "Stages 184-186" attributions with "Stages 133-135" (for `(M_s, M_q)` and `Σ_m^*`).

The 132/136 paper cards themselves are clean (`Verification: SymPy audit: none yet. Mathematica audit: none yet.` is truthful).

## Cluster C — Stage 134 paper-card Checks downgrade

**Pattern:** Stage 134's paper card `\stagefield{Checks}` lists two items neither script exercises:
1. Outlet-consistency
2. Susceptibility closure

Outlet-consistency is the **subject** of stage 135 (the very next stage); susceptibility-closure runs through stages 137-138. Asking 134's scripts to exercise these would duplicate downstream work.

**Direction:** downgrade the two Checks at stage 134 to cite stages 135/137 as the carry-forward source, rather than requiring 134's own scripts to exercise them. Paper-card edit only; no script change for these Checks.

The stage 134 scripts still receive the F1+F2+F4 mechanical fixes (SymPy asserts + tautological replacement + independent `.wl` derivation) — only the paper_misalignment finding F3 is resolved by this paper-card edit.

## Mechanical groups (no user gate)

### M1 — `.wl` banner mass-fix
Same as Cluster A above. 9 edits.

### M2 — Notes H1 mass-fix
Same as Cluster A above. 10 edits.

### M3 — Add SymPy substance (130, 131, 134, 135, 137)
SymPy scripts have either zero meaningful assertions (131, 134, 135) or only check the tangential identity (130, 137). Add the following per-stage:

- **130**: assert sign of `d(g_Π)/dΠ` symbolically as `Cov(f, z)/L > 0` (replace numeric-only check); assert the closed form equals the integral evaluation.
- **131**: add asserts: `gPi(Pi_*) == g_minus`, `g'(Pi_*) > 0` (symbolic), and the parent-threshold identity at the root.
- **134**: assert `S_q(Pi_*) ≈ 0.658075937605428`, assert the fixed-point equation residual at `Pi_*` numerically, assert canonical-gain-line identity.
- **135**: replace tautological `Pi_* - sigma*(4 - S_*)` with five substantive asserts mirroring Mathematica parity: (i) `M_s = 4·Σ_m`, (ii) `M_q = -Σ_m`, (iii) `0 < S_q(Π_*) < 1`, (iv) `Σ_m^* ≈ 0.451485`, (v) mixed-lane correction.
- **137**: derive `rho_c`, `sigma_c` from variational principle (or cite Stage 97 Schur identity with a substantive equality check); add asserts that `M_s` and `M_q` closed forms match the paper card.

### M4 — Replace tautological wrappers (131, 134, 135, 137)
Each had at least one assertion whose RHS was derived from its LHS by construction. Replace with an independent re-derivation path, then assert equality between the two paths.

### M5 — Independent `.wl` derivation (127, 133, 134, 137)
`.wl` files mirrored `.py` line-by-line. Replace one section of each with an independent derivation:
- **127**: derive `g_slab(x)` and `g_exp(x)` from the source integral (`Integrate[sigma · cos(...), ...]`) symbolically in Mathematica, then `expectZero` the difference against the `.py`-computed closed form.
- **133**: use `DSolveValue` to derive `u(x)` from `u''-Π u' = 0` with BCs, replacing the hand-typed ansatz.
- **134**: use `Series` + `Limit` in Mathematica to extract `S_q` as the `kappa → π/2` boundary, rather than retyping the hand-derived expression.
- **137**: derive `rho_c`, `sigma_c` from free-energy variation in Mathematica (independent of `.py` formulas).

### M6 — Stale-output refresh
After all script edits land, run sympy then mathematica once per stage (single-seat, sequential) to refresh saved `.txt` outputs.

## Sequence

1. Apply M1+M2 (banner + H1 mass-fix), Cluster B notes re-attribution, Cluster C paper-card downgrade — all text edits, no script run needed.
2. Apply M3+M4+M5 per stage, in stage order (130, 131, 133, 134, 135, 137 for substance; 127 for transliteration).
3. Run each touched script via `python3` / `math -script` (single-seat); iterate until clean exit.
4. Render and dispatch verify prompts.
5. Close the batch.

The same orchestrator-direct math-authority pattern as IV.3 applies — Codex remains bypassed.

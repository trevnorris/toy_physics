# Grok review — force_visualizer physics fidelity (READ-ONLY)

You are an adversarial expert reviewer. READ-ONLY: do not modify files. Where you can, **COMPUTE-verify** the math with SymPy.

## What this is
A 2D phenomenology visualizer for a toy 4D-superfluid analog model with four force sectors (gravity, light, electric charge, magnetism). It implements **reduced/effective, calibrated** force laws (NOT the full PDE). Code: `software/force_visualizer/`. Spec: `software/force_visualizer/notes/build_spec.md`.

## The hard rule you are checking
Every physics law in `software/force_visualizer/physics/` must be **FAITHFUL to its cited authoritative doc** — correct functional FORM, sign, power, and any DERIVED constant. Magnitudes may be CALIBRATED (chosen) but must be **labeled** as such (never passed off as derived). Code must NOT invent physics beyond the docs and must NOT depend on `superfluid_lib`.

## Authoritative docs (source of truth) per sector
- Gravity: `software/stage1_solver/reports/pathA_29_brane_bulk_return.md` (+ `pathA_29_results.yaml`); PN two-body forms in `research/4d_1pn_full/paper/4d_1pn_full.tex`, `research/4d_2_5pn/paper/4d_2_5pn.tex`.
- Light: `software/stage1_solver/reports/pathA_36_c5_phase_potential.md` (+ results yaml); `research/pde_ledger_v2/notes/stages/ledger_stage003_transverse_photons_stray_longitudinal.md`, `ledger_stage005_sound_speed_light_ratio.md`; `research/1pn_optics/paper/1pn_optics.tex`.
- Charge: `software/stage1_solver/reports/pathA_38_throat_body_electric_localization.md` (+ `pathA_38_results.yaml`).
- Magnetism: `software/stage1_solver/reports/pathA_39_magnetic_force.md` (+ `pathA_39_magnetic_force_results.yaml`), `pathA_39_stage4_field_classification.md`, `pathA_39_scalar_admixture_screen.md`.
- Cross-sector calibration: `software/stage1_solver/reports/pathA_40_cone_lock.md`, `pathA_41_ng5_second_medium_drift.md`.

## Files to review
- `software/force_visualizer/physics/{gravity,light,charge,magnetism,departures,integrators}.py`
- `software/force_visualizer/params.py`
- `software/force_visualizer/report.py` + its output `software/force_visualizer/output/verification_report.txt`
- `software/force_visualizer/tests/*.py`

## Check and COMPUTE-verify
1. **Gravity:** Newtonian 1/r² + the 1PN two-body acceleration terms match the cited EIH/PN forms (cross-coefficients, signs). Verify the perihelion-precession formula the report checks against (expect 6πG(m1+m2)/(c²a(1−e²))). Confirm the monopole/dipole radiation-residual departure matches pathA_29's ε0-tied residual (form, not magnitude).
2. **Light:** transverse wave `c_γ²=μ_R/ρ_br`, dispersion ω=c_γk, exactly 2 transverse polarizations / no longitudinal; lensing deflection sign+scale; the stray-longitudinal departure matches pathA_36 (`FAIL_CAUCHY`). Check the longitudinal-speed value in the report.
3. **Charge:** 1/r² Coulomb; ±w sign (like-repel/unlike-attract); DERIVED constant `N0=8/(3ℓ)` (report: 3.333 at ℓ=0.8 — verify); Yukawa partner mass `√3/ℓ` (report: 2.165 — verify). Confirm the Yukawa *residue* is labeled calibrated (pathA_38 gives no universal residue).
4. **Magnetism:** current-current `1/R` potential → `1/R²` force; transverse channel like-currents-attract (correct sign); the unavoidable attractive longitudinal scalar admixture. **Verify the scalar/transverse ratio formula** — report shows ratio=49 at aL=0.35, aT=1, μ_R=100, B_eff=0.25 — and explicitly comment whether a scalar admixture DOMINATING the total force ~49:1 is doc-faithful or an artifact of the freely-chosen amplitudes aL/aT. Confirm aT, aL, c_E are labeled calibrated/sim-deferred.
5. **Provenance labels** (params.py / report): anything labeled DERIVED that is actually calibrated/imposed, or vice-versa? λγ and c_E MUST be CALIBRATED (pathA_40, the Δr=2 result).
6. **Tests:** are the golden tests genuinely able-to-fail, or tautological / pass-by-construction? Flag any test that cannot fail.
7. **Reconciliations Codex made — assess soundness:** (a) the 2.5PN Burke–Thorne normalization exposed as an optional calibrated benchmark, off by default (native normalization `GENUINE_BLOCKED`); (b) Yukawa residue labeled calibrated; (c) optics background speed identified with c_γ, default λγ=1 an explicit calibration.

## Output
For EACH sector: **FAITHFUL** or **DISCREPANCY** (with the doc reference + exactly what is wrong). Then: any rigged/tautological tests; any provenance mislabels; your magnetism scalar-dominance assessment; and an **overall verdict** (CLEAN / ISSUES with a prioritized list). Be specific — cite file:line and the doc where possible.

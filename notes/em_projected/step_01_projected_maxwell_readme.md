# Projected Maxwell SymPy bundle

This bundle contains three derivation note/script pairs plus one bundle-level
audit script.

## Files

- `step_02_projected_maxwell_covariant_notes.md`
- `step_02_projected_maxwell_covariant_sympy.py`
  - derives the exact projected inhomogeneous Maxwell law from the localized `4+1` bulk equation
  - keeps the projection kernel `W(w)` explicit
  - keeps the gauge term generic as `Gamma^nu`
  - derives the projected open-system charge continuity law

- `step_03_projected_maxwell_vector_notes.md`
- `step_03_projected_maxwell_vector_sympy.py`
  - rewrites the projected system in Maxwell-style `E/B` and `D/H` language
  - shows that the homogeneous equations project cleanly
  - shows that the inhomogeneous equations acquire leakage and gauge-driver terms
  - makes explicit that measured fields and source-coupled flux fields are not the same objects in general

- `step_04_projection_reduction_comparison_notes.md`
- `step_04_projection_reduction_comparison_sympy.py`
  - compares projection-first and reduction-first zero-mode couplings
  - derives the generic formula
    `mu0_eff^(proj) = mu0 * (∫ W S) / (∫ W Z)`
  - compares that with the reduction result
    `mu0_eff^(red) = mu0 / Z_int`
  - evaluates a Gaussian example

- `step_01_projected_maxwell_readme_sympy.py`
  - checks the bundle inventory and audits the main formulas summarized here

## Main takeaways

1. A projection-first Maxwell theory is naturally an **open-system electrodynamics**.
2. The homogeneous laws remain ordinary in terms of measured projected fields.
3. The inhomogeneous laws involve:
   - projected `Z F^{mu nu}` rather than only projected `F^{mu nu}`,
   - explicit transverse leakage terms from `Z F^{w nu}`,
   - projected gauge-driver terms.
4. In general, the brane-visible measured fields `(E_meas, B_meas)` are **not** the same as the source-coupled flux fields `(D_flux, H_flux)`.
5. Projection does **not** automatically reproduce the controlled reduction that yields standard brane Maxwell theory.

## Run

```bash
python step_01_projected_maxwell_readme_sympy.py
python step_02_projected_maxwell_covariant_sympy.py
python step_03_projected_maxwell_vector_sympy.py
python step_04_projection_reduction_comparison_sympy.py
```

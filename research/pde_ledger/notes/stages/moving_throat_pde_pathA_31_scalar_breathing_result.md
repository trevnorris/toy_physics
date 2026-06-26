# pathA_31 scalar breathing Gate 2 result

Verdict: `BREATHING_CALIBRATED`.

This records scaffold section 11.2. The scalar eta_00 wall channel was switched on, projected through the single ell=0 operator into `M_AB Qddot + K_AB Q = F_HF`, and only then compared to the legacy `(a,L)` closure.

Truncation at the predeclared Gate-1 anchor beta_L0=`1.85`: `o_1=0.993109102589`, `o_2=0.98776369936`, `min(omega_1^2,omega_2^2)=3.42251944599`, `gap=2.22787035351`, `epsilon_trunc=0.1`.

The result remains calibrated because `muEta`, `Tw`, and `K_eta` are frozen/calibrated wall inputs even though the two collective profiles are harmonic liftings of that operator.

Structure gate: `M_posdef=True`, `K_structure_ok=True`, `K_offdiag_negative=True`, computed from the derived `M_AB,K_AB` with probe `{'non_posdef_M_probe': {'mutation': 'M_aa -> -M_aa', 'M_posdef': False, 'leading_minor_positive': False, 'det_positive': False}, 'sign_flipped_K_probe': {'mutation': 'K_aL -> -K_aL', 'K_offdiag_negative': False, 'K_structure_ok': False}}`.

Artifacts:
- `software/stage1_solver/tools/pathA_31_scalar_breathing_sympy.py`
- `software/stage1_solver/tools/pathA_31_scalar_breathing.wl`
- `software/stage1_solver/reports/pathA_31_scalar_breathing.md`
- `software/stage1_solver/reports/pathA_31_results.yaml`

# S11c-b adjudication v2 builder report

- Command: `python research/pde_ledger_v3/scripts/S11c_b_adjudicated_comparison.py` (exit 0).
- Core multiset: `join=147 + py_only=36 + wl_only=48 = emitted=classified=231`; case-ID multiset equal.
- Routes: `MATCH=38`, `FLAG=12`, `REPRESENTATIONAL_DIVERGENCE=0`, `RESIDUAL_BULK=8`, `DIVERGENCE_INCOMPLETE=0`, `PROTECTED_UNREDUCED=32`, `STRUCTURE_INCOMPLETE=57`, `COVERAGE=84`; sum `231`.
- Namespace diagnostics outside that core multiset: `NAMESPACE_INCOMPLETE=12`.
- Bridge D is the imported 12-entry `PROFILE_GRADE_SUBS` object, including all four density bookkeepers and no `e_W_bg` (`scripts/S11c_b_brane_operator_sympy_audit.py:662-670,714`); anchored second-jet documentation is at lines 1850-1869.
- Fixtures: `a*φ_d1` -> `RESIDUAL_BULK`; `a_d1*φ+a*φ_d1` -> `REPRESENTATIONAL_DIVERGENCE`, `V=(a*φ,0,0)`, certificate `0`; the same residual on `SLAB_OPERATOR` -> exact `FLAG`.
- Anchored-profile fixture: `W_bg_d1*φ+W_bg*φ_d1` -> `REPRESENTATIONAL_DIVERGENCE`, `V=(W_bg*φ,0,0)`, Bridge-D certificate `0`; the printed Bridge-D/derivative noncommutator is nonzero.
- `FLAG`: `SLAB_OPERATOR_TERM_ORIGINS (OBJECT={ADVECTIVE,KINETIC}, BRANCH={LAB_HELD,MATERIAL_ADVECTED}, DENSITY={RHO4_CONSTANT,RHOBR_CONSTANT})`; `ADMISSIBILITY_OPERATOR_OPERAND (OBJECT=BODY_FORCE, DOF=THETA, BRANCH={LAB_HELD,MATERIAL_ADVECTED}, DENSITY={RHO4_CONSTANT,RHOBR_CONSTANT})`.
- `RESIDUAL_BULK`: `COUPLING_KERNEL (SECTOR=TRANSVERSE_TO_THICKNESS, BRANCH={LAB_HELD,MATERIAL_ADVECTED}, DENSITY={RHO4_CONSTANT,RHOBR_CONSTANT})`, both without `OBJECT` and with `OBJECT=ADJOINTNESS_OPERAND_FORWARD`.
- `PROTECTED_UNREDUCED`: `COUPLING_KERNEL (SECTOR=THICKNESS_TO_TRANSVERSE, OBJECT={absent,ADJOINTNESS_OPERAND_REVERSE}, BRANCH={LAB_HELD,MATERIAL_ADVECTED}, DENSITY={RHO4_CONSTANT,RHOBR_CONSTANT})`; `COUPLING_KERNEL (OBJECT=ADJOINTNESS_RESIDUAL, BRANCH={LAB_HELD,MATERIAL_ADVECTED}, DENSITY={RHO4_CONSTANT,RHOBR_CONSTANT})`.
- `PROTECTED_UNREDUCED`: `SLAB_OPERATOR (OBJECT={MASS_EVOLUTION_ROW,MU_THETA,THICKNESS_ROW,U_MOMENTUM_ROWS}, BRANCH={LAB_HELD,MATERIAL_ADVECTED}, DENSITY={RHO4_CONSTANT,RHOBR_CONSTANT})`; `MU_THETA_OPERATOR (BRANCH={LAB_HELD,MATERIAL_ADVECTED}, DENSITY={RHO4_CONSTANT,RHOBR_CONSTANT})`.
- Jet accounting: `JET_CONSERVED=231`, `JET_LOST=0`; full before/after occurrence multisets were printed per case.
- Runtime: `672.75058556 s` instrument time (`682.72 s` wall, peak RSS `567460 KiB`).
- New ablations: `--drop-bridge-d` touched 158 cases; `--drop-divergence` touched 8 and changed `RESIDUAL_BULK 8 -> 0`, `FLAG 12 -> 20`; both exited 0 with 231/231 accounting and printed before/after operands.
- Legacy ablations: `--drop-bridge-a` touched 52; `--drop-rename LWidth` touched 20 and changed `MATCH 38 -> 26`; `--collapse-jet w1_profile_d1=w1_profile` touched 88 and emitted `JET_LOST=56`; all exited 0 with 231/231 accounting.
- No residual target was supplied. Strong operators were compared exactly; HeldDiv nodes were formally expanded, never dropped.
- The divergence classifier ran only on source-proven scalar coupling-density paths; weak origin containers remained structural and `ADJOINTNESS_RELATION` remained ineligible.
- Protected sets were kept raw, and energy quotient representatives were not folded; only `ENERGY_BASIS_COUNT` entered exact comparison.

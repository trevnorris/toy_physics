# W0 round-2 repair — verification report

- Guard honored: I did not read `steps/S9_PILOT_ADJUDICATION.md`.
- E0 verified by execution: current S10-Python reached `emit_registry_comparison`, raised
  `TypeError: 'BoundDimensionLaw' object is not iterable`, and exited 1.
- E1 verified: the old builder required only the square-link residual; committed S9 MAIN and S10 MAIN
  D2–D5 outputs show the selected root quotient equals the emitted material recombination in both engines.
- E2 verified: selection was not bound to the predicate; S9 X8 emits the nonzero root `muG/rhoBr`
  (`mu_G/rho_br` in Python) with transverse nullity 0.
- E3 verified per engine/cell: S9 MAIN has predicate nullities `0,2`; S10 MAIN D2–D5 have respectively
  `0,1`, `0,2`, `0,3`, `0,4` in each engine, hence one predicate-true root in all ten cells.
- E4 verified from the old text: its sentinel entered `v_T^2`, and its selection mutation did not hold the
  solved spectrum or coefficients fixed.
- E5 verified algebraically from the emitted D-laws: solving for powers of the two named coefficients at
  dimension `[1,-1,0]` gives the unique exponents `{1/2,-1/2}` for symbolic `D`.
- E6 verified textually: the old boundary forbade all existing value changes while regression expressly
  allowed runtime movement and required inventory/count movement.
- E7 verified: the old manifest requirement did not require values and compared output only with a
  builder-authored manifest.
- E8 verified: the old scope/field text did not choose cell-versus-root placement or classify S9 direction
  specialisations; committed outputs contain both root-local and direction-specialised records.
- E9 verified: `L_T` is identically zero under its own square-root definition, while the repository's
  instantiated R4 residual can be nonzero for the current wrong bindings.
- Could not verify: post-repair engine emissions, manifests, mutations, or repointed coverage; no such
  builds exist, and W0 remains blocked until the W3 round-2 repair lands.
- Decision-list errors found: none.
- Changed only the two W0 statements and this report; no engine, `reduction/` file, committed `.out`, or
  forbidden adjudication file was modified, and no commit was made.

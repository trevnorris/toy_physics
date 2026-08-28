# Build brief — S11c-b hand-coded cross-engine reconcile layer (delegated; mirror the S11c-a pattern)

The T7 comparator (`scripts/S11c_b_cross_engine_comparator.py`, review-clean, committed) joins and prints
`operand_A`/`operand_B`/`A_minus_B` and DECIDES NOTHING. Its name-canonicalization inherits the S11c-a tables
but does NOT carry the many new S11c-b symbol spellings, so most core-family residuals print as
`Mismatch(kind='STRUCTURAL_ATOM_DISAGREE')` (a cross-engine SPELLING gap) rather than a real algebraic
difference. This layer is the reconciliation, written the S11c-a way: apply a HAND-VERIFIED, ENUMERATED
rename map and re-check the residual, and **FLAG everything the map does not cover — never massage.**

## Object
`research/pde_ledger_v3/scripts/S11c_b_handcoded_comparison.py`, mirroring the committed
`scripts/S11c_a_handcoded_comparison.py` (READ IT FIRST — same import-the-comparator + reconcile() +
residual_zero() + per-family verdict shape). It PRINTS a per-family, per-case verdict; it asserts nothing on
measured payloads and carries NO zero/nonzero target (rule 5). Exit 0 always except on operational failure.

## Scope — the SYMPY-PARSED CORE families ONLY
Reconcile exactly these families (the comparator parsed their WL operands into sympy — `operand_B` is
`Symbol(...)` srepr): `COUPLING_KERNEL`, `COUPLING_KERNEL_TERM_ORIGINS`, `SLAB_OPERATOR`,
`SLAB_OPERATOR_TERM_ORIGINS`, `MU_THETA_OPERATOR`, `ADMISSIBILITY_OPERATOR_OPERAND`, `ADMISSIBILITY_RESIDUAL`,
`ADMISSIBILITY_SUPPORT_OPERAND`, `ENERGY_BASIS_VARIABLE`, `ENERGY_BASIS_COUNT`, `ENERGY_BASIS_NEW_INVARIANTS`,
`ENERGY_BASIS_OMISSIONS`, `DIMENSIONS`, `HOMOGENEITY_BASE_OPERAND`, `HOMOGENEITY_CONTROL_OPERAND`,
`HOMOGENEITY_RESIDUAL`.
⛔ Do NOT attempt the control families `CONTROL_FORM_*`, `CONTROL_INDEPENDENCE_*`, `REP_INVARIANCE_*`,
`UNIFORM_LIMIT_*`: the comparator emitted their WL operands as RAW Wolfram InputForm (`<| … |>`
Associations it could not structurally parse), so they are namespace/structure-incomplete and OUT OF SCOPE
here. For each such family print one line `NAMESPACE_INCOMPLETE <family> (WL operand unparsed; cross-engine
control comparison owed; each engine's internal control verified in the build legs)` and move on. Do not
fabricate a comparison for them.

## The reconcile — mirror S11c-a exactly
1. Import the comparator: `import S11c_b_cross_engine_comparator as C`. Reuse `C.load_py`, `C.load_wl`,
   `C.extract_family`, `C.compare_family` (capture its `operand_A =`/`operand_B =` stdout lines the way
   `S11c_a_handcoded_comparison.py` does), `C.BOUND_BINDER`, and its integral combiner
   (`C.combine_bound_integrals` or the S11c-a equivalent it exposes) for the zero test.
2. Build `WL_TO_PY_RENAME` — an EXPLICIT dict, one entry per verified spelling, each with an inline comment
   citing WHERE in each engine source the two names denote the SAME physical quantity. Construct it by
   READING BOTH engine sources and the spec:
   - `scripts/S11c_b_brane_operator_sympy_audit.py` (PY snake_case emitted names),
   - `mathematica/S11c_b_brane_operator_mathematica_audit.wl` (WL camelCase emitted names),
   - `directives/S11c_b_SHARED_PHYSICS.md` (§1c energy, §1a/§3 — the canonical physics each name denotes).
   Confirmed starting spellings (VERIFY each against the sources; do not take them on faith, and ADD the rest
   you find — the background fields, the gamma-modulus/gamma-width couplings, the `A_T_s11cb_*` transverse
   amplitude jets, the `Lambda_*` face responses, `bRho`,`cCross`,`eWave`, etc.):
   `WZero→W_0`, `etaBg→eta_bg`, `sigmaW→sigma_W`, `kappaW→kappa_W`, `LWidth→L_W`,
   `w1ProfileZero→w1_profile`, `muR→mu_R_bg`. WL derivative jets collapse a base name (e.g.
   `widthProfileJet`, `m1ProfileZero`) while PY carries explicit `_dNdM` suffixes — reproduce the S11c-a
   `canon_wl_basic` PROFILE/jet decoding rather than hand-listing every jet, if that is how the WL source
   encodes them.
3. `reconcile(expr)`: apply ONLY the enumerated renames (and the S11c-a-style verified applied-vs-bare /
   jet-decode reproductions). ⛔ NO blanket `camelCase→snake_case` transform. ⛔ NO applied→bare collapse
   except where you have CITED that the engine never differentiates the object (S11c-a `INERT_APPLIED`).
4. `residual_zero(a,b)`: parse → reconcile both → `simplify(combine_bound_integrals(expand(a-b))) == 0`.

## ⛔ Two objects that must be FLAGGED, never folded (physics-bearing)
- **ENERGY_BASIS_* is a non-unique QUOTIENT** (modulo total in-plane divergences). ⛔ Do NOT name, assume, or
  fold any representative identity to force agreement — a variable-coefficient IBP generates first-background-
  jet terms that are PHYSICS. Compare raw after renames; a surviving difference is `FLAG`, surfaced for
  post-run adjudication.
- The **coupling-kernel ADJOINTNESS residual** is already reduced modulo compact-support IBP by the
  comparator. ⛔ Do NOT add any further collapse; if a residual survives it is the genuine non-self-adjoint
  `Λ_X` face response — `FLAG` it.

## Output (per family)
`MATCH` (all cases residual_zero True), `FLAG (n/m cases differ — <case keys>)` (print the surviving
`A_minus_B` residual for each differing case — that IS the candidate finding), `COVERAGE` (all cases
one-engine-only; from the comparator accounting), or `NAMESPACE_INCOMPLETE` (control families above). Also
emit a final `RENAME_MAP_SIZE=<n>` line and, for the one-sided ablation, a `--drop-rename <WLname>` flag that
removes one entry so a leg can confirm each rename is load-bearing (removing a genuine rename must turn a
MATCH into a FLAG, proving it was the same variable, not two different ones forced together).

## Definition of done (the build legs check these empirically)
Every in-scope case resolves to `MATCH` or `FLAG` (with the residual printed) or `COVERAGE`; control families
each print exactly one `NAMESPACE_INCOMPLETE`. The rename map is an explicit dict with per-entry source
citations. `--drop-rename` works. No case silently dropped; nothing asserted on measured payloads; no
`PASS`/`FAIL`/`VERDICT`/target. Runs against the committed transcripts by default (reuse `C.DEFAULT_PY`,
`C.DEFAULT_WL`).

## Builder report (≤25 lines)
Per-family verdict tally (MATCH/FLAG counts); the rename-map size and a few representative entries with their
source citations; which FLAGs survived and their family+case keys (NOT an interpretation of whether they are
findings — that is the orchestrator's post-run call); runtime. State that no residual target was given and
that ENERGY_BASIS/adjointness were compared raw (never folded).

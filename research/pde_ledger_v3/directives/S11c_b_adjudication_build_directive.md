# Build brief — S11c-b cross-engine ADJUDICATION layer (v4; folds decision-list rounds 1-3 + the comparator fix)

The T7 comparator (`scripts/S11c_b_cross_engine_comparator.py`, review-clean, committed) joins two engines'
emitted operands, losslessly decodes every WL `Derivative[...]` into a DISTINCT jet token (including time order,
after the `canon_jet_name` time-order fix `bb0bfc02`), maps applied base heads to bare PY symbols, and prints
`operand_A`/`operand_B`/`A_minus_B` per case while DECIDING NOTHING. The v1 reconcile
(`scripts/S11c_b_handcoded_comparison.py`, both build legs SOUND, committed) adds an 80-entry HAND-VERIFIED
spelling map (`H.WL_TO_PY_RENAME`) and re-checks zero. This layer adds ONE further source-verified algebraic
identity (Bridge A), classifies every remaining case by its operand SORT, compares what is genuinely
like-for-like, defers what is not, and accounts for every emitted case. It is the never-blanket-collapse danger
zone: **apply only the enumerated bridge and the enumerated source-cited container adapters; introduce NO
applied→bare or jet fold; account for every case; never massage.**

## The three script clauses (verbatim, non-negotiable)
1. PRINT computed objects; never state a conclusion. Every payload is a CAS object, never a prose verdict.
2. PRINT the residual; never assert it. This layer carries NO zero/nonzero target and NO expected value,
   classification, or ablation outcome for any case (rule 5). A surviving difference is a finding, not a build
   failure.
3. Interpretation belongs to the STEP RECORD, not the script.

## Object
`research/pde_ledger_v3/scripts/S11c_b_adjudicated_comparison.py`. IMPORTS the comparator
(`import S11c_b_cross_engine_comparator as C`) and v1 (`import S11c_b_handcoded_comparison as H`); REUSES
`H.WL_TO_PY_RENAME` UNCHANGED; relies on the comparator's committed jet decoding and adds no applied→bare or
jet-collapse logic of its own. Exit 0 except on operational failure. Default inputs `C.DEFAULT_PY`,
`C.DEFAULT_WL`. ⚠ The comparator's `--residual-leaf-budget` fallback makes residual serialization
timing-dependent on the heavy families; run/compute residuals with a GENEROUS or disabled leaf budget so the
output is reproducible run-to-run.

## Scope — the SYMPY-PARSED CORE families ONLY (identical to v1)
`COUPLING_KERNEL`, `COUPLING_KERNEL_TERM_ORIGINS`, `SLAB_OPERATOR`, `SLAB_OPERATOR_TERM_ORIGINS`,
`MU_THETA_OPERATOR`, `ADMISSIBILITY_OPERATOR_OPERAND`, `ADMISSIBILITY_RESIDUAL`,
`ADMISSIBILITY_SUPPORT_OPERAND`, `ENERGY_BASIS_VARIABLE`, `ENERGY_BASIS_COUNT`, `ENERGY_BASIS_NEW_INVARIANTS`,
`ENERGY_BASIS_OMISSIONS`, `DIMENSIONS`, `HOMOGENEITY_BASE_OPERAND`, `HOMOGENEITY_CONTROL_OPERAND`,
`HOMOGENEITY_RESIDUAL`. Control families (`CONTROL_FORM_*`, `CONTROL_INDEPENDENCE_*`, `REP_INVARIANCE_*`,
`UNIFORM_LIMIT_*`) each print exactly one `NAMESPACE_INCOMPLETE` line, as v1.

## LOCKED physics folds (verified refs supplied; a build leg re-verifies each)
- **Bridge A — the bRho / B_rho_3 normalization identity, as a SYMBOL substitution `bRho ↦ B_rho_3 / W_0`**
  (equivalently reduce modulo `W_0·bRho − B_rho_3 = 0`), applied to the renamed operands as an explicit,
  comment-cited substitution. Verified refs: WL θ² coefficient `bRho·anchoredWidth` with factor `1/2`,
  `anchoredWidth`=`W_bg` (`mathematica/S11c_b_brane_operator_mathematica_audit.wl:472`); homogeneous block
  `bRho·WZero` (WL:1621-1630); PY θ² coefficient `B_rho_3*W_bg/(2*W0)`
  (`scripts/S11c_b_brane_operator_sympy_audit.py:1130-1140`); spec `directives/S11c_b_SHARED_PHYSICS.md:102`
  (`B_ρ⁽³⁾ ≡ B_ρ W₀`). After it, no non-ablated ALGEBRAIC operand carries a residual `bRho` atom.
- **Bridge C — integral linearity** in the zero test (`H.COMBINE_BOUND_INTEGRALS`), unchanged from v1. No new
  fold.

## PROTECTED — never folded or adapted into agreement (verified refs; each carries physics)
PY quotient selected candidate set `{1,4,5,6,7,9,10,13}`, omitted `{2,3,8,11,12,14,15}`
(`scripts/S11c_b_brane_operator_sympy_audit.py` ~1273); WL's eight explicit invariants map to PY
`{1,4,5,6,8,9,11,13}` (WL:417-435). Therefore:
- `gamma{Width,Modulus}DivGrad{Theta,Ew}` are PY candidates **08/11** (WL:426,430) — omitted candidates, not
  the selected representatives v1 maps. Never map them.
- PY-selected representatives **07/10** have no WL counterpart — one-engine physics. Never fold or adapt them.
- **ENERGY_BASIS_\*** is a non-unique QUOTIENT (modulo total in-plane divergences); a variable-coefficient IBP
  generates first-background-jet terms that are physics. No substitution and no container adapter may pair or
  identify a representative.
- The coupling-kernel **ADJOINTNESS** residual is already reduced modulo compact-support IBP by the comparator
  (`*_DIVERGENCE_REDUCED` context); add no further collapse.

## Routing INVARIANT — classify every case by operand SORT (complete taxonomy, covers all cases)
For each joined case, determine the sort of the two operands (after `H.WL_TO_PY_RENAME` + Bridge A) and route to
EXACTLY one of:
- **ALGEBRAIC** — permitted ONLY when both operands are arithmetic, tested POSITIVELY as `isinstance(v, sp.Expr)`
  (a `Symbol` is arithmetic even though this SymPy also makes it a `Boolean` — test Expr, do not reject via the
  Boolean superclass), or same-shape matrices/tuples of arithmetic `sp.Expr`. Reconcile (renames + Bridge A +
  Bridge C) and re-check zero; `MATCH` if zero for all cases, else `FLAG` with the reduced `A_minus_B` residual
  printed per differing case. ⛔ A `Equivalent`/`Not`/relational Boolean, a `TextAtom`, a scalar-vs-relation
  sort mismatch, or a shape/arity mismatch is NOT arithmetic and MUST NOT enter this route.
- **CONTAINER** — like-typed containers (Association/tuple) whose leaves can be paired. May yield `MATCH` or
  `FLAG` ONLY through a source-cited **TOTAL BIJECTION** over ALL semantic leaves of BOTH sides (every leaf
  paired exactly once by a cited same-meaning label correspondence; NO unmatched, duplicate, or unlabeled leaf;
  NEVER pairing differently-labelled leaves or identifying a PROTECTED representative). Under such a bijection,
  compare each paired leaf ALGEBRAICALLY (as above) and print leaf residuals (`MATCH` if all zero, else `FLAG`).
  If no total-bijection adapter is enumerated and cited for a family, emit `STRUCTURE_INCOMPLETE <family> <key>
  (<computed shape diff>)`; a partial/subset leaf comparison may be printed as a diagnostic but NEVER yields a
  MATCH or FLAG verdict.
- **STRUCTURE_INCOMPLETE** — a JOINED case whose two operands are a genuine sort mismatch that no total
  bijection resolves: Boolean-vs-Expr, TextAtom-vs-Expr, relational-vs-scalar (e.g. the four
  `COUPLING_KERNEL/ADJOINTNESS_RELATION` cases: PY `NO_INDEPENDENT_SECOND_ROUTE` scalar `PY:1980` vs a WL
  relation `WL:1130`). Print the computed sort diff; do not fabricate an algebraic or container comparison.
- **COVERAGE** — a `<MISSING>` operand or a `TextAtom('UNDEFINED_UNJOINED')` join atom (one-engine-only, from
  the comparator accounting); never ALGEBRAIC.

**Enumerated container adapters are the builder's to find and cite; the build legs verify each is a faithful
total bijection.** The pattern: `SLAB_OPERATOR_TERM_ORIGINS/KINETIC` is a PY tuple (pos 0 = `U_MOMENTUM_ROWS`,
pos 1 = `THICKNESS_ROW`, `scripts/S11c_b_brane_operator_sympy_audit.py:1573`) versus a WL Association on those
same labels (WL:851) — a source-fixed total correspondence that MUST be adapted and compared, not deferred.
Read the sources and enumerate every such adapter with citations; where the correspondence is not total or not
source-backed, `STRUCTURE_INCOMPLETE`. Boolean admissibility operands (PY constructs logical `Equivalent`/`Not`
residuals, `scripts/S11c_b_brane_operator_sympy_audit.py:868-890`) are not arithmetic and are not adaptable by a
leaf bijection unless a separately enumerated, source-proved conversion exists.

## Accounting — EXACT case-ID multiset equality (no silent drop, no double-count)
The classification denominator is the comparator-emitted CORE case multiset = `join + py_only + wl_only` (NOT
`join` alone). Collect the multiset of case IDs the comparator emitted and the multiset of case IDs this layer
classified, and assert they are EQUAL including multiplicity; raise an operational error on any missing,
duplicate, or extra case ID. A sum of per-route counters is not sufficient (it cannot detect one-drop +
one-duplicate). Print the per-route counts AND the total, and confirm the total equals the emitted multiset
size.

## JET conservation — a COMPUTED invariant on the substitution map (the never-collapse guard)
Define a semantic jet ID `(canonical_base, derivative_multiindex)` where the base is taken AFTER applying
`H.WL_TO_PY_RENAME` (so a spelling rename is not a loss) and the decoder treats `theta_dN` and `grad_theta_N`
(and the general PY jet vocabulary — reuse `C.s11ca.canon_jet_name` and the `_dN`/`_tt` decoding) as the SAME
jet of order N. The layer's substitution map (renames + Bridge A) must be jet-order-preserving: every
substitution maps a token to one of the SAME derivative order (a rename permutes spelling: `theta_d1 →
grad_theta_1`, both order 1; it never lowers an order or deletes a jet without a same-order image). Take the
"before" jet-ID multiset from the comparator-captured operands (pre-reconcile) and the "after" from the
reconciled operands, and PRINT, per case, the computed `JET_CONSERVED`/`JET_LOST` from comparing those two
multisets — the label is derived from the computed comparison, not written in.

## Ablation hooks — each DECISIVE (each only changes the transformation and prints; states no outcome)
- `--collapse-jet <token>=<base>`: inject an order-reducing substitution (e.g. `w1_profile_d1=w1_profile` or
  `=0`) independent of all other logic.
- `--drop-bridge-a`: run without the bRho substitution.
- `--drop-rename <WLname>`: remove one spelling equivalence AND disable the comparator prepass for it (reproduce
  or delegate to `H._disable_imported_prepass_for_drop`, `scripts/S11c_b_handcoded_comparison.py:440`).
Each hook: its argument is drawn from the captured pre-transform inventory; an unknown or non-occurring argument
is an OPERATIONAL ERROR (nonzero exit), never a silent no-op; on success it reports the (nonempty) set of cases
it touched, shows a syntactically changed transformed operand, and prints recomputed before/after residuals and
jet multisets. A hook only CHANGES the transformation and PRINTS; it states no outcome.

## Definition of done (value-free structural conditions the build legs check empirically)
- Every emitted CORE case is classified into exactly one of `{ALGEBRAIC(MATCH|FLAG), CONTAINER(MATCH|FLAG via
  total bijection | STRUCTURE_INCOMPLETE), STRUCTURE_INCOMPLETE, COVERAGE}`; each control family prints exactly
  one `NAMESPACE_INCOMPLETE`. The classified case-ID multiset EQUALS the comparator-emitted `join+py_only+wl_only`
  multiset (asserted); no case dropped, duplicated, or invented.
- Bridge A is the symbol substitution `bRho ↦ B_rho_3/W_0`; after it, no non-ablated ALGEBRAIC operand carries a
  residual `bRho` atom; `--drop-bridge-a` removes it.
- The ALGEBRAIC route is entered only when both operands satisfy the positive `sp.Expr` arithmetic test (or
  same-shape arithmetic matrices/tuples); Boolean/TextAtom/relational/scalar-sort/shape mismatches never enter.
- Every CONTAINER verdict rests on a source-cited total bijection over ALL leaves; no subset comparison yields a
  verdict; no differently-labelled pairing; no PROTECTED-representative identification.
- Jet conservation is computed from before/after jet-ID multisets (base after rename; `theta_dN≡grad_theta_N`,
  time order via `_tt`); the substitution map is order-preserving; `--collapse-jet` forces the order-reducing
  case.
- `07/10`, every `gamma-DivGrad`, and ENERGY_BASIS representatives appear unfolded/unadapted wherever a case
  carries them.
- Each ablation hook is decisive: unknown/non-occurring argument → operational error; on success a nonempty
  touched-case set, a syntactically changed operand, and before/after residuals are printed.
- Nothing asserted on a measured payload; no `PASS`/`FAIL`/`VERDICT`/target/expected/predicted anywhere.

## Builder report (≤25 lines)
The routing accounting (per-route counts + total = emitted multiset size, asserted equal); the Bridge-A
substitution and its source citations; every enumerated container adapter with its citation; the computed
`JET_CONSERVED`/`JET_LOST` counts; which cases FLAG or are `STRUCTURE_INCOMPLETE` (family + keys only — NOT
interpreted); runtime. State that no residual target was given and that `07/10`, gamma-DivGrad, ENERGY_BASIS
reps, and adjointness were kept unfolded and unadapted.

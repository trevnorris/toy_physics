# Build brief — S11c-b cross-engine ADJUDICATION layer (v3; invariants + accounting, not a schema list)

The T7 comparator (`scripts/S11c_b_cross_engine_comparator.py`, review-clean, committed `5f01f9fa`) joins two
engines' emitted operands, losslessly decodes every WL `Derivative[...]` into a DISTINCT jet token
(`canon_wl_basic`, `FIELD`, `EXTRA_HEAD`; committed and reviewed), maps applied base heads to bare PY symbols,
and prints `operand_A`/`operand_B`/`A_minus_B` per case while DECIDING NOTHING. The v1 reconcile
(`scripts/S11c_b_handcoded_comparison.py`, both build legs SOUND, committed `82ec3b5f`) adds an 80-entry
HAND-VERIFIED spelling map (`H.WL_TO_PY_RENAME`) and re-checks zero. This layer adds ONE further source-verified
algebraic identity (Bridge A), then classifies every remaining case by an INVARIANT (its operand sort), compares
what is genuinely like-for-like, and defers what is not — with per-case ACCOUNTING so nothing is silently
dropped. It is the never-blanket-collapse danger zone: **apply only the enumerated bridge and the enumerated,
source-cited container adapters; introduce NO applied→bare or jet fold; account for every case; never massage.**

## The three script clauses (verbatim, non-negotiable)
1. PRINT computed objects; never state a conclusion. Every payload is a CAS object, never a prose verdict.
2. PRINT the residual; never assert it. This layer carries NO zero/nonzero target and NO expected value,
   classification, or ablation outcome for any case (rule 5). The diff is read on our side; a surviving
   difference is a finding, not a build failure.
3. Interpretation belongs to the STEP RECORD, not the script.

## Object
`research/pde_ledger_v3/scripts/S11c_b_adjudicated_comparison.py`. IMPORTS the comparator
(`import S11c_b_cross_engine_comparator as C`) and v1 (`import S11c_b_handcoded_comparison as H`); REUSES
`H.WL_TO_PY_RENAME` UNCHANGED; relies on the comparator's committed jet decoding and adds no applied→bare or
jet-collapse logic of its own. Exit 0 except on operational failure. Default inputs `C.DEFAULT_PY`, `C.DEFAULT_WL`.

## Scope — the SYMPY-PARSED CORE families ONLY (identical to v1)
`COUPLING_KERNEL`, `COUPLING_KERNEL_TERM_ORIGINS`, `SLAB_OPERATOR`, `SLAB_OPERATOR_TERM_ORIGINS`,
`MU_THETA_OPERATOR`, `ADMISSIBILITY_OPERATOR_OPERAND`, `ADMISSIBILITY_RESIDUAL`,
`ADMISSIBILITY_SUPPORT_OPERAND`, `ENERGY_BASIS_VARIABLE`, `ENERGY_BASIS_COUNT`, `ENERGY_BASIS_NEW_INVARIANTS`,
`ENERGY_BASIS_OMISSIONS`, `DIMENSIONS`, `HOMOGENEITY_BASE_OPERAND`, `HOMOGENEITY_CONTROL_OPERAND`,
`HOMOGENEITY_RESIDUAL`. Control families (`CONTROL_FORM_*`, `CONTROL_INDEPENDENCE_*`, `REP_INVARIANCE_*`,
`UNIFORM_LIMIT_*`) each print exactly one `NAMESPACE_INCOMPLETE` line, as v1.

## LOCKED physics folds (verified refs supplied; a build leg re-verifies each)
- **Bridge A — the bRho / B_rho_3 normalization identity, as a SYMBOL substitution `bRho ↦ B_rho_3 / W_0`**
  (equivalently reduce modulo `W_0·bRho − B_rho_3 = 0`). Verified refs: WL θ² coefficient `bRho·anchoredWidth`
  with factor `1/2`, `anchoredWidth`=`W_bg` (`mathematica/S11c_b_brane_operator_mathematica_audit.wl:472`);
  independent homogeneous block `bRho·WZero` (WL:1621–1630); PY θ² coefficient `B_rho_3*W_bg/(2*W0)`
  (`scripts/S11c_b_brane_operator_sympy_audit.py:1130-1140`); spec `directives/S11c_b_SHARED_PHYSICS.md:102`
  (`B_ρ⁽³⁾ ≡ B_ρ W₀`). Apply on the renamed operands as an explicit, comment-cited substitution.
- **Bridge C — integral linearity** in the zero test (`H.COMBINE_BOUND_INTEGRALS`), unchanged from v1. No new
  fold.

## PROTECTED — never folded or adapted into agreement (verified refs; each carries physics)
The PY quotient's selected candidate set is `{1,4,5,6,7,9,10,13}` and omitted `{2,3,8,11,12,14,15}`
(`scripts/S11c_b_brane_operator_sympy_audit.py` ~1273); WL's eight explicit invariants map to PY
`{1,4,5,6,8,9,11,13}` (WL:417–435). Therefore:
- `gamma{Width,Modulus}DivGrad{Theta,Ew}` are PY candidates **08/11** (WL:426,430) — omitted candidates, not the
  selected representatives (v1 maps the selected ones). Never map them.
- PY-selected representatives **07/10** have no WL counterpart — one-engine physics. Never fold or adapt them
  into a match; a case carrying them FLAGs.
- **ENERGY_BASIS_\*** is a non-unique QUOTIENT (modulo total in-plane divergences); a variable-coefficient IBP
  generates first-background-jet terms that are physics. No substitution and no container adapter may pair or
  identify a representative.
- The coupling-kernel **ADJOINTNESS** residual is already reduced modulo compact-support IBP by the comparator
  (`*_DIVERGENCE_REDUCED` context); add no further collapse.

## Routing INVARIANT — classify every case by operand SORT (this covers all cases uniformly)
For each joined case, determine the sort of the two operands (after `H.WL_TO_PY_RENAME` + Bridge A) and route:
- **ALGEBRAIC** — permitted ONLY when both operands are arithmetic `sp.Expr`, or same-shape matrices/tuples of
  arithmetic `sp.Expr`. Reconcile (renames + Bridge A + Bridge C) and re-check zero; print the reduced
  `A_minus_B` residual for each nonzero case. ⛔ A Boolean (`Equivalent`/`Not`/relational), a `TextAtom`, an
  atom-sort mismatch, or a shape/arity mismatch is NOT arithmetic and MUST NOT enter this route.
- **CONTAINER** — when the operands are like-typed containers (Association/tuple) whose leaves can be paired.
  A container comparison may yield a zero-residual verdict ONLY through a source-cited **TOTAL BIJECTION** over
  ALL semantic leaves of BOTH sides: every leaf paired exactly once, by a cited same-meaning label
  correspondence, with no unmatched leaf, no duplicate, no unlabeled leaf, and **never** pairing
  differently-labelled leaves or identifying a PROTECTED representative. Under such a bijection, compare the
  paired leaves ALGEBRAICALLY (as above) and print leaf residuals. If no total-bijection adapter is enumerated
  and cited for a family, emit `STRUCTURE_INCOMPLETE <family> <key> (<computed shape diff>)`. A partial/subset
  leaf comparison may be printed as a diagnostic but NEVER yields a zero verdict.
- **ONE-ENGINE / JOIN** — a `<MISSING>` operand or a `TextAtom('UNDEFINED_UNJOINED')` join atom routes to
  `COVERAGE` (from the comparator accounting), never ALGEBRAIC.

**Enumerated container adapters are the builder's to find and cite; the build legs verify each is a faithful
total bijection.** The pattern: `SLAB_OPERATOR_TERM_ORIGINS/KINETIC` is a PY tuple (pos 0 = `U_MOMENTUM_ROWS`,
pos 1 = `THICKNESS_ROW`, `scripts/S11c_b_brane_operator_sympy_audit.py:1573`) versus a WL Association on those
same labels (WL:851) — a source-fixed total correspondence that MUST be adapted and compared, not deferred.
Read the sources and enumerate every such adapter with citations; where the correspondence is not total or not
source-backed, defer (`STRUCTURE_INCOMPLETE`). Boolean admissibility operands (PY constructs logical
`Equivalent`/`Not` residuals, `scripts/S11c_b_brane_operator_sympy_audit.py:868-890`) are not arithmetic and are
not adaptable by a leaf bijection unless a separately enumerated, source-proved conversion exists.

## JET conservation — the never-collapse guard, made a COMPUTED invariant on the substitution map
Define a semantic jet ID `(canonical_base_after_rename, derivative_multiindex)` (reuse the comparator's jet-name
machinery, e.g. `C.s11ca.canon_jet_name` and the `_dN` decoding). The layer's entire substitution map (the
renames + Bridge A) must be **jet-order-preserving**: every substitution maps a token to one of the same
derivative order — a rename permutes a token's spelling (`theta_d1 → grad_theta_1`, same order 1); it never
lowers an order or deletes a jet without a same-order image. Take the "before" jet inventory from the
comparator-captured operands (pre-reconcile) and the "after" from the reconciled operands, and PRINT, per case,
`JET_CONSERVED`/`JET_LOST` computed from comparing those two jet-ID multisets — do not hardcode either label.

## Ablation hooks — each must be DECISIVE (a build leg breaks the reconcile with it)
- `--collapse-jet <token>=<base>`: inject an order-reducing substitution (e.g. `w1_profile_d1=w1_profile` or
  `=0`) independent of all other logic.
- `--drop-bridge-a`: run without the bRho substitution.
- `--drop-rename <WLname>`: remove one spelling equivalence AND disable the comparator prepass for it (reproduce
  or delegate to `H._disable_imported_prepass_for_drop`, `scripts/S11c_b_handcoded_comparison.py:440`).
Each hook: its argument is drawn from the captured pre-transform inventory; an unknown or non-occurring argument
is an OPERATIONAL ERROR (nonzero exit), never a silent no-op; on success it reports the (nonempty) set of cases
it touched, shows a syntactically changed transformed operand, and prints recomputed before/after residuals and
jet multisets. A hook only CHANGES the transformation and PRINTS; it predicts nothing.

## Definition of done (value-free; the build legs check these empirically)
- **Accounting:** every in-scope case is classified into exactly one of `{ALGEBRAIC(MATCH|FLAG), CONTAINER(MATCH
  via total bijection | STRUCTURE_INCOMPLETE), COVERAGE}`; the per-route counts are printed and sum to the joined
  case count; each control family prints exactly one `NAMESPACE_INCOMPLETE`. No case silently dropped.
- Bridge A is the symbol substitution `bRho ↦ B_rho_3/W_0`; after it, no non-ablated ALGEBRAIC operand contains a
  residual `bRho` atom; `--drop-bridge-a` restores it.
- The ALGEBRAIC route is entered only for arithmetic `sp.Expr` / same-shape arithmetic matrices/tuples; a
  Boolean, `TextAtom`, atom-sort, or shape mismatch never enters it.
- Every CONTAINER `MATCH` cites a total bijection over ALL leaves of both sides; no subset match, no
  differently-labelled pairing, no PROTECTED-representative identification. Each enumerated adapter carries a
  source citation.
- Jet conservation is computed from before/after jet-ID multisets; the substitution map is order-preserving;
  `--collapse-jet` forces the order-reducing case.
- `07/10`, every `gamma-DivGrad`, and ENERGY_BASIS representatives appear unfolded/unadapted wherever a case
  carries them.
- Each ablation hook is decisive: unknown/non-occurring argument → operational error; on success a nonempty
  touched-case set, a syntactically changed operand, and before/after residuals are printed.
- Nothing asserted on a measured payload; no `PASS`/`FAIL`/`VERDICT`/target/expected/predicted anywhere.

## Builder report (≤25 lines)
The routing accounting (per-route counts summing to the case count); the Bridge-A substitution and its source
citations; every enumerated container adapter with its citation; confirmation that every non-ablated case is
`JET_CONSERVED`; which cases FLAG or are `STRUCTURE_INCOMPLETE` (family + keys only — NOT interpreted); runtime.
State that no residual target was given and that `07/10`, gamma-DivGrad, ENERGY_BASIS reps, and adjointness were
kept unfolded and unadapted.

# Decision-list review — S11c-b requested-truncation row-level residual instrument (directive)

You are reviewing a **build directive** (an orchestrator-written decision list) BEFORE any builder runs. Your
job is to find every way the directive would make the builder compute the **wrong object**, leak an expected
answer, build a tautological check, or misread the spec. This is physics-bearing: an error here makes both a
later build and its two build legs agree on the same wrong thing.

Report a numbered list of findings (each: the defect, the spec/code evidence with file:line and quoted text,
and the correction). If the directive is sound on a point, do not pad — but a review that finds nothing is
weak evidence, so probe hard. **Change nothing on disk.**

## The artifact under review
`research/pde_ledger_v3/directives/S11c_b_row_residual_instrument_build_directive.md`

## What you are handed (read the sources of truth FIRST, form your own view, THEN read the directive)
- The spec: `research/pde_ledger_v3/directives/S11c_b_SHARED_PHYSICS.md` — §1d, §2a, §3a, §3b, §3c, §3d are
  the load-bearing sections. Read them and form your OWN statement of the requested truncation and of what
  object decides whether the two engines' operators agree.
- The two engine sources (what actually gets compared):
  `research/pde_ledger_v3/scripts/S11c_b_brane_operator_sympy_audit.py` (SymPy/PY, imports the prior ledger)
  and `research/pde_ledger_v3/mathematica/S11c_b_brane_operator_mathematica_audit.wl` (blind Wolfram/WL,
  re-derives). In particular the term-origin partition (PY L1574-1583; the strong-form `*_expanded` at
  L1470-1490) and the truncation primitives (PY `first_shape_series` L713-726; WL `truncateScalar` L99-116;
  jet grading PY L1865/L2117, WL L349).
- The reviewed alignment layer the instrument must reuse:
  `research/pde_ledger_v3/scripts/S11c_b_adjudicated_comparison.py` (rename map `WL_TO_PY_RENAME`,
  `transform`, `_bridge_d`, canonicalizer) and the step-1 multigrade instrument
  `research/pde_ledger_v3/scripts/S11c_b_background_multigrade.py`.
- Context (computed facts, not verdicts): the step-1 grade fingerprint at
  `~/.s11_build/S11c_b_grade_fingerprint.txt` shows the per-family per-bucket grade supports of the twenty
  differences. It states no verdict; treat it as measured context only.

## The questions you must answer with evidence (demand quotes and file:line; derive the truncation yourself)

1. **Truncation reading.** Independently, from §2a and §3a alone, state the requested truncation as a
   condition on the `(ε,η,σ_W)` exponents `(c,a,b)`. Does it match the directive's §1a formalization
   "retain iff `c≤1 ∧ a≤1 ∧ b≤1`, grading by amplitude-factor power not spatial-derivative order"? Quote the
   spec sentences you rely on. ⚠ Specifically: does §2a/§3a require **linearizing** a coefficient such as
   `W_bg²=W_0²(1+η w1)²` or the rational `1/(1+η w1)` in `η` (dropping `η²` and higher), or keeping the
   coefficient **in full**? If the spec is genuinely ambiguous on this, say so and name the ambiguous
   sentence — a spec ambiguity here is a SPEC defect to fix first, not a directive defect.

2. **Is the row-level object the right one?** The directive's object is the complete strong-form operator per
   DOF (the sum over ALL term-origin buckets, divergences carried out to strong form via
   `∇·(cF)=∇c·F+c∇·F`), compared at the requested truncation. Does this correctly decide whether the two
   engines' **operators** agree despite differing term-origin **bucketing** conventions? Or does it miss
   something — double-count, drop a face/boundary contribution that belongs to the strong-form operator, or
   conflate operands that live in different pairings?

3. **The coupling kernel (§3c) — weak vs strong.** §3c extracts the coupling as a **weak variational
   restriction** with the in-plane integration-by-parts boundary term fixed to zero (compact support). A weak
   operator is defined only **modulo total in-plane divergences**. Does the directive's residual, applied to
   the coupling kernel the same way as to the strong-form slab rows, correctly handle this? Or must the
   coupling residual be taken **modulo total in-plane divergence** (a variational-equivalence / weak-form
   residual) rather than as a strong-form difference? Quote §3c and §1d.

4. **Computability.** Is the object actually computable from the emitted operands? Check the engine code:
   does each engine emit the term-origin buckets and strong forms needed for the row assembly AND for the
   divergence expansion (do the operands carry the field spatial dependence so the product rule applies)?
   Name the file:lines. If PY's `*_expanded` already carries out the divergences but WL emits an unexpanded
   `Inactive` divergence form (or vice versa), does the directive's "carry out the divergence to strong form
   for each engine" step actually reconcile them, or is there a form the instrument cannot activate?

5. **Leak (rule 5).** Does the directive state, hint, or let the builder infer any expected residual value,
   any per-family verdict (representational vs genuine), or which engine is spec-correct for any row? Quote
   any leak. (The correct behaviour is compute-and-print; the diff is adjudicated off-instrument.)

6. **Tautological guards (rule 2 corollary 3).** Are `WINDOW_CLEAN` and `BUCKET_PARTITION_CHECK` genuine
   checks or zero-by-construction? A residual formed from operands produced by a single route is zero for any
   input and tests nothing. Is `BUCKET_PARTITION_CHECK` (bucket-sum minus the engine's own complete operator)
   an independent check, or does the assembly define the complete operator so the check is vacuous?

7. **Anything else physics-bearing** that would make the instrument answer the wrong adjudication question —
   a mis-specified pairing, a missing anchoring/density case, a truncation applied at the wrong stage (before
   vs after the row assembly), a canonicalization that would mask or manufacture a residual.

## Method
This is a document/directive review — no script to run. But your truncation reading (Q1) and your
computability claim (Q4) must be grounded in quoted spec text and cited engine file:lines, not asserted.
A prose claim with no quote/line behind it will be discarded. Physics filter: report a finding only if it
catches a way the instrument would compute the wrong object or leak the answer; do not report style.

# Decision-list review — S11c-b P2b §3a coefficient bridge (scale factors + certificate)

## Artifact
`research/pde_ledger_v3/directives/S11c_b_p2b_gamma_bridge_directive.md` — orchestrator-written DECISION LIST for a
PHYSICS-BEARING bridge (the manufactured-agreement danger zone). You are ONE of two independent decision legs (rule
7). The builder will TRUST this list. Find its defects now.

## Context
Both engines enumerate the same 15/source O(3)-Kronecker §3a family (#89/#89a), but with different coefficient
NAMES and — critically — different NORMALIZATION: a prior 2-leg review MEASURED `I_PY = W_0·I_WL` (W-family) and
`I_PY = μ_R·I_WL` (μ-family) — `EXACT_UNIQUE 0 / SCALED_UNIQUE 30`. WL carries the spurion normalization
`gradient[anchoredWidth]/WZero` (`…audit.wl:722`); PY uses raw jets (`sympy_audit.py:182`). The committed bridge
only expresses string→string renames (`handcoded_comparison.py:72,279`), so it CANNOT carry the scale. The directive
(B1-B6) is the fix. Your job: is it implementable and does it CATCH a wrong pairing / a folded normalization?

## Context you are handed (read the CODE, cite file:line)
- Bridge tables `scripts/S11c_b_handcoded_comparison.py` (`WL_TO_PY_RENAME` ~72-152, apply ~279-311; `BARE_APPLIED`
  ~159, arg-discard ~300) and `scripts/S11c_b_adjudicated_comparison.py` (`PROTECTED_ATOM_NAMES` ~54-65,
  `PROTECTED_UNREDUCED` ~953, ablation-touch exclusion ~1478).
- `scripts/S11c_b_cross_engine_comparator.py` `_extra_basic` (~379-429, the applied→bare collapse).
- PY invariants `scripts/S11c_b_brane_operator_sympy_audit.py`: `NEW_COEFFICIENTS` (~357-367), the invariant/coeff
  records + `ENERGY_BASIS_NEW_INVARIANTS` emit (~1819-1836, ~4101), the quotient selection (~1406-1439).
- WL invariants `research/pde_ledger_v3/mathematica/S11c_b_brane_operator_mathematica_audit.wl`: `newCoefficientSymbol`
  (~627), the new-invariant records incl. `COEFFICIENT_STANDARD_NAME` (~2136-2145), the spurion (~722).

## What to check (find defects; ground in code, cite file:line; run greps/probes and quote literal output)
1. **B1 scale representation:** is `gammaWJET* = W_0·gamma_s11cb_w_bg_NN` the correct direction (verify the algebra:
   energy term `γ·I` equal ⇒ which side carries `W_0`)? Can the bridge machinery express a scale factor, or does B1
   require a representation change the directive must name? Does the directive foreclose folding `W_0` OUT (the
   blanket-collapse trap) AND using `COEFFICIENT_STANDARD_NAME` as a no-op bridge?
2. **B2 pairing by invariant:** is matching the emitted invariants (up to the W_0/μ_R scale + existing renames)
   actually implementable — i.e. do PY and WL emit enough to compute the pairing WITHOUT a positional guess? Confirm
   the enumeration orders are independently computed (so positional matching is unsafe).
3. **B3 certificate:** does the certificate as specified actually FAIL a wrong pairing (the weak "residual moves"
   test does not)? Is anything missing from the certificate (coverage, uniqueness, both branches, the normalization
   factor, the invariant operands)?
4. **B4 argument-sensitive applied→bare:** is the fixture sufficient to catch the `widthBackground[base] −
   widthBackground[shifted] → 0` collapse? Are there other blanket-collapse sites the directive misses?
5. **B5 retire PROTECTED:** does dropping/re-adjudicating `07/10` + dead `gammaWidth*` risk regressing S11c-a or a
   legitimate quotient protection? Is the re-adjudication safe?
6. Any way the bridge MANUFACTURES agreement (positional guess, folded scale, no-op standard-name, protected atom
   hiding a residual) that B1-B6 do not catch. Any missing property.

## Required method
DOCUMENT review — do NOT modify the tree. Ground every claim in code (file:line); run greps/probes and quote LITERAL
output (e.g. compute one invariant pair's scale ratio; check the bridge apply path). Report defects as a numbered
list; name the directive item (B#) each fixes. A clean pass is weak evidence.

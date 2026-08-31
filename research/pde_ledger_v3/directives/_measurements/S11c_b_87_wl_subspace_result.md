# S11c-b #87 — WL's 8 hand-coded invariants ⊊ the correct 15 (RESULT)

## Headline (confirmed)
The Wolfram engine's 8 hand-coded spurion invariants (`newInvariantExpressions`,
`mathematica/S11c_b_brane_operator_mathematica_audit.wl:417-434`) span a **strict 8-dimensional
subspace of the correct 15/source** — genuinely undercomplete by exactly 7, with all 8 valid
(in-span), not 8 alternative-but-complete forms. This validates the #89 WL fix = "complete the
enumeration to 15" (nullity 0 ⇒ linear independence of all 15 gives 15).

## Computation (`~/.s11_build/S11c_b_87_probes/probe_wl_subspace.{py,out}`)
WL's 8 forms transliterated verbatim from `newInvariantExpressions` into the SymPy abstract space
(spurion→bg, u→bu, θ→btheta, e_W→be, ∇θ→bq, ∇e_W→br, divU→trace(bG), D[u_i,x_j]→bG[i][j]) and
compared to the complete enumeration `enumerate_new_candidates` (15/source).

    === RAW energy-candidate space ===
    rank(WL 8)            = 8
    rank(correct 15)      = 15
    rank(WL 8  UNION  15) = 15          # WL's 8 add NOTHING; all within span(15)
    each WL form in span(correct 15): True ×8

    === under the CORRECT (Hessian-retaining) quotient ===
    quotient-rank(WL 8)       = 8
    quotient-rank(correct 15) = 15
    quotient-rank(WL 8 U 15)  = 15

Both readings agree: WL's 8 are all in the span of the correct 15, span exactly 8 of the 15
directions, and none collapse under the quotient (raw rank = quotient rank = 8). The correct 15
are independent under the quotient (rank 15 = nullity 0, re-confirming #86).

## Rigor note / scope
This is an orchestrator transliteration-based confirmation of an already **code-resolved** fact
(#86: WL hand-codes 8, no divergence quotient — a DIFFERENT undercompleteness mechanism than PY's
frozen Hessian). The one hand step is the transliteration of WL's 8 forms from the `.wl` source;
it is re-verified independently when **#89's WL repair build runs the actual WL engine** (which
must emit the completed 15 and reduce to the committed 8 in the frozen limit). No engine changed.

## Next
#87 CONFIRMED. NEXT = **#89** both-engine §3a repair — fixes DIFFER per engine: **PY** correct the
quotient to retain the background Hessian (frozen 8 → 15/source); **WL** complete the enumeration
(8 → 15). Spec-confirm §1d/§3a intent → directive → 2 decision legs → build → 2 build legs each
(SERIALIZE Mathematica) → rebuild operator → re-run row-residual → **re-adjudicate KINETIC + θ**
(the #88 blast radius). Repaired engines must emit 40 (withheld from the builder). Then #90 PY §3c
content fix on the corrected basis.

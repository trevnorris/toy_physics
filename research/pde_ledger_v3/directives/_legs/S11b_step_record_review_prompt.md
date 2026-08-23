# Independent review — the S11b step record (physics-bearing prose)

## Artifact
`/var/projects/toy_physics/research/pde_ledger_v3/steps/S11b_interface_coupling_law.md`

This is the orchestrator-written step record for the unified S11b step ("the linear brane–bulk interface
coupling law"). It is the INTERPRETATION layer: it states conclusions the engines are forbidden to state, so
every conclusion must be faithful to what the engines actually COMPUTED. Your job is to catch any claim the
step record makes that the computed objects do NOT support — an overstatement, an understatement, a result
stated that no engine emitted, or a misrepresentation of the cross-engine comparison.

## What to check
Read the SOURCES OF TRUTH FIRST and form your own view of what the step establishes; only THEN read the step
record. For every physics claim in the record, find the emitted object that backs it and confirm they match.
Specifically verify against the `.out` transcripts (grep the tag; compare the parsed value):
1. **Energy basis count = 10** — `PY_S11B_ENERGY_BASIS_COUNT` (SymPy `.out`) and `WL_S11B_ENERGY_BASIS_COUNT`
   (WL `.out`). Confirm the record's X-1 narrative (11 pointwise → 10 modulo total divergence; the redundant
   invariant folded, not deleted; the representative split PY-`st²` vs WL-`(∇·u)²`) matches what the engines
   emit. ⭐ Derive the quotient dimension yourself if you doubt it (an O(3) enumeration + Euler–Lagrange
   nullity); save script + literal stdout.
2. **Transverse** — the record says the coupling is identically zero and the operator is
   `k²(μ_R + μ_S/2) − ω²ρ_br⁰`. Check `S11B_TRANSVERSE_DISPERSION` / `S11B_TRANSVERSE_COUPLING` in both `.out`s.
3. **The B2c/G13 slice map** — the record says both engines now emit static `Λ_p⁰ = −Λ_A⁰/ρ_m`. Check
   `ZPERM_SLICE_MAP` in both `.out`s.
4. **The passivity region and the reciprocity relation** — check the record's stated region
   (`Λ_A⁰ ≥ 0 and [Λ_V⁰=Λ_X⁰=0 or (Λ_X⁰=−Λ_V⁰, τ_V=τ_X=0)]`) and `Λ_X(ω)=−Λ_V(ω)` against the emitted
   admissibility / reciprocity objects. Is the "region, NOT a prohibition" framing faithful (the engines emit
   a classification, not an acceptance gate)?
5. **Breathing** — `K₀ = B_ρ⁽³⁾ − 2C W₀ + k_W W₀² > 0` on the `k=0` impermeable slice. Check the emitted
   stability object; confirm the record correctly limits it to the slice.
6. **The comparison adjudication** — the record claims ALL 108 disagreements in the re-run are
   representation/naming/convention, ZERO physics. Independently sample the comparison output
   `scripts/out/S11b_cross_engine_comparison.out`: are the 2 CONTENT items (`DEGENERATE_LOCI_EQUATIONS`,
   `ENERGY_BASIS_OMISSIONS`) really representation, or is either a genuine physics disagreement the record is
   waving away? Verify the DEGENERATE_LOCI "same equation under an i-factor and q_out≡q" claim yourself.
7. **Standing limits / owed items** — is anything misstated: the frozen-wall-width freeze
   (`ρ_br⁰ = rho_br`), the uncarried `O(v₀|q_n|/ω)` background-flow correction, the owed card items
   G12c/G12d, C's deferral and its bench-top-optics falsifiability?

## Fast path (the `.out` files are large — 18 MB SymPy)
⭐ A mechanical grep-digest of the load-bearing emitted tags (exactly the tags named below, from both
engines) is at `directives/_legs/S11b_step_record_emitted_digest.txt`. Use it as the quick reference; the
FULL `.out` files remain the authority — spot-check any tag there, and grep the full `.out` for any tag the
digest lacks. ⛔ Do not treat the digest as complete or as a substitute for your own verification.

## What you are handed
- The step record (path above).
- The emitted-tag digest `directives/_legs/S11b_step_record_emitted_digest.txt` (fast reference).
- The spec `directives/S11b_SHARED_PHYSICS.md` (the physics source of truth).
- The two engine transcripts: `scripts/out/S11b_interface_coupling_law_sympy_audit.out` and
  `mathematica/out/S11b_interface_coupling_law_mathematica_audit.out`.
- The frozen comparison `scripts/out/S11b_cross_engine_comparison.out`.
- `CLAUDE.md`. You have the full repo; ⛔ there is no do-not-read list.
⚠ You are deliberately NOT handed the build directives: a record can satisfy a directive and still
misrepresent what the engines computed, and that is exactly the case this review exists to catch.

## Required method — DOCUMENT review
Read the spec + the emitted objects first; form your own view; only then read the record. Where a claim is
checkable against a `.out` tag, VERIFY it (grep the tag, compare the parsed value) — quote both the record and
the emitted object. Where you derive anything (e.g. the basis quotient dimension, the DEGENERATE_LOCI
identity), save the script AND its literal stdout to named absolute paths and report them; a physics claim
without a computation behind it is discarded. ⛔ Do not modify the working tree.

## Physics filter
Report a finding only if it catches the record misrepresenting the physics or the computed result — an
overstated conclusion, a result no engine emitted, a mis-adjudicated disagreement, a dropped standing limit —
not prose style or wording preference.

## Output
Per item: verdict + the emitted object (tag + value) or your script path + literal stdout, and exact quotes
from the record and the source. End with a one-line bottom line: is the step record faithful to the computed
physics as written, or the specific corrections it needs.

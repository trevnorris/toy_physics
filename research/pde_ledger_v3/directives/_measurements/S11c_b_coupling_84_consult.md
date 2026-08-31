# S11c-b #84 — Codex+Grok consultation on the §3c coupling content (2 rounds; user-requested)

The user chose "pin physics intent first" and asked me to consult Codex and Grok on the two axes, then bring
back the best results. Two rounds: (1) independent, (2) adversarial cross-examination. This records the
converged verdict + the transcripts (rule 2: the claim carries its command).

⭐ UPDATE (post-consultation, current state): the §3c CONTENT verdict (INCLUDE/INCLUDE ⇒ WL correct) is a
two-consultant SPEC-READING that the USER SUBSEQUENTLY ENDORSED — it is NOT a from-scratch CAS which-engine
derivation, and it is no longer "awaiting adjudication" or conditional on user intent. Separately, the quantitative
bulk-basis probe flagged UNVERIFIED below was subsequently CAS-verified for PY (mechanism only) in
`S11c_b_coupling_84_basis_verification.md`; WL's freeze remains inferred (pending #87). Read the framing below as
the mid-consultation record.

## Prompts + transcripts (absolute)
- Round 1 prompt: `~/.s11_build/S11c_b_84_intent_consult.md`
  - Codex: `~/.s11_build/S11c_b_84_consult_codex.txt` (xhigh, 43k tok)
  - Grok:  `~/.s11_build/S11c_b_84_consult_grok.txt` (grok-4.5 high)
- Round 2 (adversarial) prompt: `~/.s11_build/S11c_b_84_round2.md`
  - Codex: `~/.s11_build/S11c_b_84_round2_codex.txt` (xhigh, 139k tok; read engine code + ran a live probe)
  - Grok:  `~/.s11_build/S11c_b_84_round2_grok.txt`

## Converged verdict (BOTH engines, BOTH rounds; survived adversarial refutation)
- **Q1 — reversible tilted-face geometry (`A_T`, no `Λ`): INCLUDE (WL spec-correct; PY under-extracts).**
- **Q2 — irreversible response coupling (`A_T·Λ(ω)`): INCLUDE (WL spec-correct; PY under-extracts).**
- **§3c is NOT genuinely ambiguous** (equal-support sense); the text determines INCLUDE/INCLUDE. But §0's
  "no permeability/memory kernel" wording is a genuine PROSE HAZARD (it uses the same English as `Λ(ω)`),
  inviting exactly PY's misread ⇒ a one-sentence §0 clarity pin is advised (not a prerequisite).

### The controlling argument (verified against the spec text myself)
"weak variational restriction under the stored-energy/kinetic pairing" names the EXTRACTION INSTRUMENT
(trial/test contraction, IBP, formal-adjoint comparison) — NOT a content filter. §3c: extract a block "of the
§3b operator itself", ⛔ "not from a parallel direct-variation route" (= PY's bulk-only recipe), ⛔ "do not
filter the kernel to a single channel"; §3b: consume the tilted-face substrate T-a..T-i "for every
boundary/face contribution"; §3c requires `TERM_ORIGINS` to classify bulk-energy/face-flux/advective (pointless
if face can't appear in the kernel); the §3c adjointness residual is expressly NOT ∂²U/∂u_T∂e_W (rules out an
energy-only Hessian object). VERIFIED VERBATIM: those §3c prohibitions exist (S11c_b_SHARED_PHYSICS.md §3c).

### Q2 decisive anchor (verified): the T-i seam
§1c SUPPLIES `Λ_I(ω)=Λ_I⁰/(1−iωτ_I)` as flat-face closure; S11c-a **T-i** shape-differentiates that closure and
states verbatim (S11c_a_SHARED_PHYSICS.md:~448) it is "**not** B0c's bulk-response assembly `δp=Z·v_bulk`; no
bulk DtN, impedance, or pressure-response solve belongs in T-i." §3b consumes T-a..T-i. ⇒ the `Λ` face closure
is inside the §3b operator and is distinct from the §0/S11c-c-excluded bulk solve. §0's ban targets the
curved-bulk DtN package, not the supplied `Λ`. VERIFIED (S11c-a T-i text quoted above exists).
⚠ Q2 PRECISION (Codex): keep `Λ(ω)` SYMBOLIC in the block; do NOT perform the bulk-elimination/DtN that
produces the effective curved-interface memory/self-energy (that IS S11c-c). WL keeps Λ symbolic — consistent.

### Suggested §0 clarity sentence (both gave near-identical; not a prerequisite to know the lean)
"The §0 exclusion of permeability/memory kernels applies only to kernels derived by the S11c-c curved-bulk
response solve (`δp=Z·v_bulk`); the supplied flat-face kernels `Λ_I(ω)` of §1c remain in every §3b face/flux
contribution and therefore in any §3c block that retains that contribution."

## The ~118-term bulk-core residual — a SEPARATE, still-open finding (NOT settled by Q1/Q2)
Both engines carry pure-bulk `γ·profile-jet·trial/test` terms and disagree on ~118 non-cancelling terms
(certified non-total-divergence). Round-2 assessment:
- Grok: likely PARTLY (maybe mostly) a basis-MISALIGNMENT artifact (equal counts 26=26 ≠ equal bases; §1d says
  identically-named `γ_i` can label different representatives differing by first-jet physics). Decisive check:
  compare invariant SPANS (the bilinear FORMS, not `γ` labels), fit a change-of-basis map, split residual into
  an alignment component vs an orthogonal physics remainder.
- Codex: ⚠⚠ UNVERIFIED-BY-ME CAS PROBE — claims a live Euler-signature probe gave PY selected span rank 8, WL
  selected span rank 8, UNION rank 10, neither containing the other, all 15 enumerated first-jet candidates
  independent once the generated background-Hessian jets are retained ⇒ "both 26-count bases may be
  UNDERCOMPLETE: their quotient tests effectively FROZE THE SPURION during IBP" (a shared §1d violation, a
  genuine construction DEFECT, not a renaming). Codex cites PY `..._sympy_audit.py:936` (15 contractions,
  selects {1,4,5,6,7,9,10,13}) and WL `..._mathematica_audit.wl:390` (DivGrad/ShearGrad reps — VERIFIED the
  file carries 8 `gammaWidth*`/`gammaModulus*` coeffs/source incl. DivGradTheta/DivGradEw).
  ⇒ ⛔ MUST verify Codex's rank-8/8/union-10 claim in CAS before acting (rule 13; "no rederivation trusted
  until in CAS"). If union rank > each ⇒ both §3a bases undercomplete (repair BOTH engines); if a change-of-
  basis kills the residual ⇒ it was alignment. This overlaps task #85 (energy-basis quotient).

## Net (for the user's INTENT decision — NOT yet acted on)
§3c-as-written mandates INCLUDE both axes ⇒ WL spec-correct, PY under-extracts (missing ALL face + response
coupling). IF the intended physics matches (coupling carries reversible face geometry + the supplied dissipative
face response), the path is: repair PY (+ §0 clarity pin), then separately verify/resolve the bulk-core basis
question. IF the intent is a conservative/bulk-only coupling, §3c-as-written contradicts it ⇒ a spec change.
The face/dissipation content decision is the user's; the bulk-core defect is a computation owed regardless.

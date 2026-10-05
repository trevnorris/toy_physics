# D4C.1–D4C.4 independent fidelity review

Both independent non-author reviews are **CLEAR**. No required finding remains.
The bounded full D4 constant-coefficient bulk contract is complete under
[FORMALIZATION_POLICY.md](../FORMALIZATION_POLICY.md). This closes D4C.1–D4C.4
only, not S11 as a whole or the separate CAS production work.

## Reviewed revision and reviewers

The user approved transfer of the fixed 52-file packet by answering
"Yes you can." to the exact-packet request. The authorization, run identities
and immutable snapshot are recorded in
`S11_lean_d4_bulk_review_state.json` (`archive/pre-cleanup-2026-10-04:research/pde_ledger_v3/_measurements/S11_lean_d4_bulk_review_state.json`).

- Packet: `S11 full D4 bulk D4C.1–D4C.4 fidelity contract v1`.
- Aggregate SHA256: `061d0d96f0c3b62439875877aa4a107b53b6870a206b308ee8a04056bd8e8a11`.
- Archive SHA256: `35bd9b379c9e93ba3b43fa5c6c6484f0931066b1462308133977fd1efe759eae`.
- Claude session: `c2238e91-147c-46a9-a659-e85c4f6ee195`, completed
  2026-09-17 17:34:59 UTC. Reported model usage includes
  `claude-opus-5[1m]` and `claude-haiku-4-5-20251001`.
- Grok session: `986e9e65-d39a-4dc5-92b7-8edfa5cfc61d`, completed
  2026-09-17 17:44:41 UTC. Reported model: `grok-4.6-build`.

These were separate, sequential source-only sessions with the same approved
packet. Neither session supplied the other review. Durable terminal reports:

- [Claude report](../../_measurements/S11_lean_d4_bulk_fidelity_claude_v1.md)
- [Grok report](../../_measurements/S11_lean_d4_bulk_fidelity_grok_v1.md)

Both read the contract, new statements and relevant definitions, native
identification, controls and recorded evidence. They independently checked the
mathematical interpretation and conventions; neither ran Lean, SymPy, Wolfram
or independent hash checks. Compilation, axiom and mutation results and file
preservation remain author-validated execution evidence, not a second build
performed by each reviewer.

The raw JSON, run records and stderr are retained alongside the reports.
Claude terminated with `success/end_turn/completed`, exit 0 and empty stderr.
Grok terminated with `end_turn` and exit 0 after startup warnings and a recovered
HTTP 503. Its final text contains the full substantive CLEAR report; earlier
commentary and its separate thought field are not the accepted report. The
extraction offsets and hashes are retained in the closure validation.

## Findings and dispositions

Both reviewers confirm the actual supplied action `L=-Q/2`, unrestricted
constant real coefficients, gradient and time-row conventions, actual integrated
first variation, exhaustive bulk equivalence, exact null family and dimensions,
both divergence identities, and the compact native sign/factor/basis map.
Neither requested a mathematical repair.

| Review note | Disposition |
|---|---|
| Claude: the new odd current normalization witness is at a value/jet, while the even witness is an actual smooth affine field. | Clarified in the live fidelity record. The frozen prompt's plural "actual smooth affine-current witnesses" is too broad for the new Controls module. The even affine-field divergence is 2; the odd current value/jet normalization is 1. The reviewed odd divergence theorem applies to every smooth field. No extra odd affine-field theorem is claimed or needed by D4C.3. |
| Claude: canonical definition-level mutations could strengthen evidence. | Not added. The sixteen paired controls mutate mathematical statements; two native tests mutate the original EL implementation. No Lean canonical source-replacement test is claimed. Existing normalization and passing controls meet the declared contract. |
| Claude: the direct bulk-equivalence pair tests an overly strong condition; another pair could test an overly weak one. | Retain the bounded suite. The universal iff is proved, while the separate longitudinal/transverse and null-sign controls discriminate the two response coefficients. No claim that the suite detects every wrong equivalence formula is made. |
| Claude: a field-level homogeneous identity would be stronger. | The agreed identification is the compact modal identity, including zero wavevector. No new field-level homogeneous theorem is claimed. |
| Claude: native null responses use mixed-derivative canonicalization. | Retain the explicit smooth-field interpretation and native normalization checks on nonzero density/momentum. Zero EL alone is not normalization evidence. |
| Grok: an explicit compact test witnessing nonzero variation could be constructed. | Not added. `exists_nonzero_firstVariation` proves actual existence from the exhaustive stationarity/null-family equivalence; no explicit test function is claimed. |
| Grok: broader null-Lagrangian, variable-coefficient, D5, interface or spectral results would be stronger. | These remain exclusions, not unfinished D4C obligations. |

The extracted reports are verbatim terminal text. Claude's description of
"third-derivative terms" in the even current calculation is a prose slip:
differentiating `u * partial u` produces second field derivatives, which cancel
by smooth mixed-partial commutation. The checked current identity and its
Lean statement are unchanged.

These dispositions clarify the evidence and retain the reviewed scope. They
change no definition, hypothesis, theorem, control or native computation. No
substantive re-review or build rerun is required for these documentation edits.

## Closure evidence and mathematical meaning

Author validation rechecked archive, snapshot, transport and all 52 live packet
files before closure edits; 33 current objects; 15 clean package pins and 12
direct Mathlib source/object pairs; 327 historical Lean source/evidence files;
three read-only native/generator inputs and 18 protected NP/T1/VC objects.
No shared-source drift was found. The full transitive Mathlib cache remains
the pinned baseline, not a fresh transitive rebuild.

Recorded verification remains 67 standard-axiom audits, sixteen mathematical
rejections and twenty positive executions (eighteen distinct source statements),
with 69 total records. Native evidence remains twelve exact identities, two
original-helper mutations and seven target checks. No instrument failure or
warning counts as a mathematical rejection.

For the classified four-coefficient family, bulk response depends exactly on
`c` and `a+b`. The response and null spaces each have dimension two; the full
null family is `(t,-t,0,beta)`, represented by the two explicit divergence
currents. This is exhaustive proof of the supplied family, not a new physical
discovery, pointwise vanishing of the null density, or absence of boundary
effects.

`S11_lean_d4_bulk_closure_validation.json` (`archive/pre-cleanup-2026-10-04:research/pde_ledger_v3/_measurements/S11_lean_d4_bulk_closure_validation.json`)
records the reviewed hashes, report extraction, final preservation checks and
documentation-only differences. The approved snapshot, archive and transport
are unchanged. No commit was requested or performed for this closure.

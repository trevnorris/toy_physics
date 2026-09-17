# D4B.1–D4B.4 independent fidelity review and closure

**The bounded D4 odd-density divergence and bulk-variation contract is
complete.** Claude and Grok independently returned CLEAR without required
corrections. Local verification and current source/object correspondence pass.
This closes [D4_ODD_COVERAGE.md](D4_ODD_COVERAGE.md) under
[FORMALIZATION_POLICY.md](../FORMALIZATION_POLICY.md).

## Fixed reviewed revision

The user approved this exact 38-file packet with “Yes you can”. Both reviewers
received the same isolated snapshot, sequentially in separate sessions with
read-only tools. Neither received the other's report or authoring duties.

- Revision: `S11 D4 odd divergence/variation D4B.1–D4B.4 fidelity contract v1`.
- Aggregate SHA256:
  `1c2fb3840e56f82e506adbddb13463a6fbc645193294a3173531026f0ae838d5`.
- Method: SHA256 of UTF-8 `json.dumps(files, sort_keys=True)`.
- Archive SHA256:
  `166f3de5984894b7b9383e87e03a7b1e5c731ea9a891ca1a724309836eb4b8cc`.
- Manifest SHA256:
  `970bb4ce1cc64423ca2da4c3c1fc86ad93fd515793c849f2fb56dea314ee7f14`.
- [Manifest](../../_measurements/S11_lean_d4_odd_review_packet.json),
  [archive](../../_measurements/S11_lean_d4_odd_review_packet_v1.tar.gz),
  [authorization and state](../../_measurements/S11_lean_d4_odd_review_state.json).

The repository manifest and the snapshot/transport `MANIFEST.json` are
byte-identical. Both terminal reports identify the exact aggregate above.

| Independent reviewer | Session | Result |
|---|---|---|
| Claude, `claude-opus-5[1m]` (run also records Haiku usage) | `f299d658-17d3-473c-90ac-442c76bd0aed` | CLEAR; terminal success, exit 0 |
| Grok, `grok-4.6-build` | `e77798a7-f754-4406-8520-377d3ee083da` | CLEAR; terminal `end_turn`, exit 0 |

Durable reports: [Claude](../../_measurements/S11_lean_d4_odd_fidelity_claude_v1.md)
and [Grok](../../_measurements/S11_lean_d4_odd_fidelity_grok_v1.md). Raw JSON,
stderr and run records remain unchanged. Grok's final report starts at character
offset 843 of `text`; preceding progress and the separate `thought` field are
not verdict evidence.

Both reviewers inspected the actual action, derivative-defined momenta,
divergence current, smoothness assumptions, finite relative action,
integration by parts, universal quantifiers, native conventions and mutation
diagnostics. They checked the algebra and counterexamples by hand. Neither
compiled Lean, executed instruments or recomputed hashes: these are independent
statement-fidelity reviews, not independent build attestations.

Author validation separately checked all 38 live files before closure edits,
archive/snapshot/transport, 43 check sources and logs, nineteen compiled objects,
the native/proof/dependency/generator/instrument hashes, fifteen clean tracked
package checkouts at their pins, and historical H/I/E/J/K/D4 proof bytes and
manifests. No shared-source change occurred during review.

Claude stderr is empty. Grok logged startup plugin-precedence, hook-loading
and `/tmp` repository-discovery warnings, then completed its substantive review
of the correct packet. Those warnings are neither mathematical evidence nor
a partial review. No cancelled, errored or plan-only response was accepted.

## Findings and dispositions

Neither reviewer requested a required fix. Optional suggestions were assessed
against the declared contract and its stopping rule.

| Observation | Disposition |
|---|---|
| Claude: package differentiability and zero derivative in a single `HasDerivAt` corollary. | Differentiability is already proved by `relativeAction_hasDerivAt`; `firstVariation_zero` is connected to that theorem through the actual first-variation identity. No use of an undefined derivative value supports the result. A wrapper is optional. |
| Both: add an explicit affine-field density witness alongside the jet witness. | The contract already proves the universal smooth-field identity and exhibits nonzero density and momenta. The affine field is admissible, and both reviewers checked its density by hand. A further convenience lemma is not a missing coverage obligation. |
| Claude: normalize the wrong-index diagnostic further; Grok: add the explicit field to the native wrong-index control. | The named false identity and independently checked field `u=(0,0,0,x1*x2)` establish the intended failure. Later rewrite/simp errors are explicitly excluded from evidence. No instrument failure is accepted as a rejection. |
| Grok: prove the entire finite relative action vanishes for compact variations. | The agreed contract claims its actual first variation, already proved for every smooth background and compact variation. Finite-action constancy is a stronger corollary and is not added under the stopping rule. |
| Claude: keep the limited Wolfram evidence explicit. | Retained: only the generator-transpose source-string check is claimed. No Wolfram execution or kernel certification of CAS code is implied. |

The current-factor mutation likewise rests on the false factor identity and
the field `u=(0,x1,0,x3)`, giving mutated divergence 2 versus density 1. Its
later witness-conversion error is not rejection evidence. Both reviewers
confirmed this distinction and the wrong-index counterexample.

Reviewer field subscripts occasionally use zero-based notation. Canonical
correspondence remains Lean component `j` to native `u_(j+1)`: in particular
the nonzero momentum witness is native matrix entry `[0,1]`, and the non-null
`G00²` control acts on the first component. No source convention changed.

No definition, theorem, assumption, instrument or native translation changed
after review. Re-review, proof rebuilding and mutation reruns are unnecessary
for these documentation-only dispositions.

## Closed result and limits

For the classified D4 odd density `P=F01 F23-F02 F13+F03 F12`, the supplied
`L=-beta P/2` has actual momenta `-beta M/2` and zero time row. The explicit
current `K_i=(1/2) sum_j u_j M_ij` satisfies `div K=P`. The actual local EL
and actual finite relative-action first variation vanish for **every constant
real beta, every smooth background and every smooth compact variation**.
There is no sign, nonzero-beta, time-independence or polarization premise.
The density and momenta can be nonzero; boundary effects are not excluded.

[D4_ODD_VERIFICATION.txt](D4_ODD_VERIFICATION.txt) records four new modules
plus root, fourteen guarded unchanged imports, 41 standard-axiom audits,
twelve mathematical rejections, twelve positives and compact native checks.
Run2 freshly built the new modules/root and ran every control; focused builds
were not accepted as recorded reuse evidence.

Only coverage, fidelity, verification and the README were updated after review;
this closure record and durable reports were added. The other 34 packet files,
all proof/instrument/native bytes, verification reports and nineteen objects
retain the reviewed correspondence. The approved snapshot/archive/transport
remain frozen. See [closure validation](../../_measurements/S11_lean_d4_odd_closure_validation.json).
The frozen local validation's pending-review status is historical; the review
state and this record establish completion.

This closes D4B.1–D4B.4 only. It is a checked constant-coefficient result, not a
new physical discovery or a full D4 bulk census. Variable beta, D5, general
null Lagrangians, spectra, interfaces and production/export work are excluded.
All S11c files and pinned exports remain untouched. The separately requested
[analytic assessment](S11C_ANALYTIC_ASSESSMENT.md) is a proposal, not another
completed proof contract. No commit was requested or made in this closure turn.

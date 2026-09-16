# D4.1–D4.4 independent fidelity review and closure

**The bounded D4 quadratic invariant classification D4.1–D4.4 is complete.**
Claude and Grok independently returned CLEAR with no blockers. Local checks
and current source/object correspondence pass. This closes
[D4_COVERAGE.md](D4_COVERAGE.md) under
[FORMALIZATION_POLICY.md](../FORMALIZATION_POLICY.md).

## Fixed reviewed revision

The user approved this exact 34-file packet with “Proceed”. Both reviewers
received the same isolated snapshot, sequentially, in separate sessions with
read-only tools. Neither received the other's report or authoring duties.

- Revision: `S11 D4 quadratic invariants D4.1–D4.4 fidelity contract v1`.
- Aggregate SHA256:
  `a91cfbf327a6d180e5f8d3e9e3dd9dc0037c426b71be41c0a55ec45d7a719762`.
- Method: SHA256 of UTF-8 `json.dumps(files, sort_keys=True)`.
- Archive SHA256:
  `cf69d3c8cc4b1b239052a8245d11f21bdfaa1e033e2fdd2806bb4e75c4c51229`.
- [Repository manifest](../../_measurements/S11_lean_d4_review_packet.json),
  [archive](../../_measurements/S11_lean_d4_review_packet_v1.tar.gz),
  [authorization and state](../../_measurements/S11_lean_d4_review_state.json).

The repository manifest is supplied to reviewers as `MANIFEST.json` at the
packet root. Those files are byte-identical, with SHA256
`b5e49158bb7ef2b3f81b7346ba449615e5d20327c6324d99993141e730b8a13b`.
This explicit alias resolves Claude's packet-name observation.

| Independent reviewer | Session | Result |
|---|---|---|
| Claude, `claude-opus-5[1m]` (run also records Haiku usage) | `bf27013a-bb48-4c14-84ad-71f23849893a` | CLEAR; terminal success, exit 0 |
| Grok, `grok-4.6-build` | `af92f927-4f14-4b00-81dd-116b01b7a034` | CLEAR; terminal `end_turn`, exit 0 |

Durable reports: [Claude](../../_measurements/S11_lean_d4_fidelity_claude_v1.md)
and [Grok](../../_measurements/S11_lean_d4_fidelity_grok_v1.md). Raw JSON, stderr
and run records are retained unchanged. Grok's final report begins at character
offset 1157 of its `text` field; preceding progress and its separate `thought`
field are not verdict evidence.

Both reviewers inspected the definitions, complete quadratic representation,
necessary constraints, reconstruction, full-group sufficiency, reflection
split, dimensions, normalization, native helper connection and control
diagnostics. They performed mathematical hand checks but did not compile
Lean, execute instruments or recompute hashes. Their contribution is independent
statement-fidelity review, not independent build attestation. The author
separately validated all 34 live packet files, archive/snapshot/transport, 44
check logs, fourteen compiled D4 objects, pins and fifteen clean tracked
dependency-package checkouts. No shared-source changes occurred during review.

Claude stderr was empty. Grok emitted startup plugin-precedence, hook-loading
and `/tmp` repository-discovery warnings, then completed its substantive
review of the correct packet. No partial, cancelled, errored or plan-only
response was accepted as clearance.

## Findings and dispositions

Neither reviewer requested a required fix. Optional observations were assessed
against the bounded contract and its stopping rule.

| Observation | Disposition |
|---|---|
| Claude A: repository and packet manifest names differ. | The byte-identical alias and exact aggregate are recorded above. No reviewed file needs correction. |
| Claude B: the determinant character on all orthogonal matrices is not packaged as a separate odd-space theorem. | The contract defines oddness using one reflection inside SO. Its complete classification and the all-matrix `orientation_conjugate` theorem are proved; no stronger wrapper is required. |
| Claude C: dimension and omission controls reuse proved lemmas instead of perturbing reconstruction certificates. | They test statement sensitivity. Two source mutations separately test form coefficients, and the kernel checks every reconstruction certificate. Existing controls satisfy D4.4; a further certificate mutation is optional. |
| Claude D: older step prose attributes the D4 extra to cross-pairings between isomorphic summands. | Recorded as a separate editorial follow-up. Claude notes that in D4 the self-dual and anti-self-dual SO(4) summands are non-isomorphic, and the odd invariant corresponds to the difference of their norms. The formal proof uses the exhaustive coordinate classification and does not depend on that prose. `steps/S11_stray_longitudinal.md` is preserved at its reviewed bytes. |
| Grok: a conceptual Pfaffian proof could replace the determinant expansion. | The existing polynomial identity is kernel-checked for every real matrix R. Replacing a valid proof would not close an open contract item. |
| Grok: native Q9 uses infinitesimal generators while Lean quantifies over the full group. | The compact check compares the actual native polynomial spaces to the completely classified Lean spaces. This tested translation boundary is explicit; no second formal proof of the native generator algorithm is required. |

Grok's probe summary is shorthand rather than a table of raw values. At
`G=diag(1,1,0,0)` the exact form value is `4a+2b+2c`; the canonical injectivity
proof and Claude's review use the full coordinates. No coefficient or theorem
change is indicated by that wording.

No substantive statement, assumption, normalization, instrument or native
translation changed after review. Re-review, proof rebuilding and mutation
reruns are therefore unnecessary for these closure-documentation changes.

## Closed result and limits

Every SO(4)-invariant real quadratic form on all real 4×4 gradients has unique
coefficients in `(tr G)²`, `tr(G²)`, `tr(GGᵀ)` and the displayed orientation
form P. Full O-invariance is exactly `d=0`; reflection-oddness inside SO is
exactly `a=b=c=0`. The actual invariant spaces have dimensions 4/3/1, and the
even/odd decomposition is complete with zero intersection. The native V1/V2/V6
spaces match exactly, including `P_D=P` and the fully summed epsilon
contraction `2P`. Wolfram evidence remains source inspection only.

[D4_VERIFICATION.txt](D4_VERIFICATION.txt) records the thirteen modules plus
root, 49 standard-axiom audits, fourteen mathematical rejections and sixteen
positive controls. Run5 reused run4's canonical builds under source/pin/generator
and input/output object guards and reran every control.

All 34 live packet files matched before closure documentation. Only D4's
coverage, fidelity, verification and README section changed afterward; this
closure record and durable reports were added. The other 30 packet files,
all proof bytes, native sources, instruments, local-check records and compiled
objects retain their reviewed hashes. The approved snapshot, archive and
transport remain unchanged. See
[closure validation](../../_measurements/S11_lean_d4_closure_validation.json).
The frozen local validation's pending-review label is historical; the review
state and this document record completion.

This closes D4.1–D4.4 only. The D4 odd term's divergence/zero bulk variation is
a separate next increment. D5, general null Lagrangians, variable coefficients,
dynamics/spectra/interfaces, production/export work and systematic CAS bridging
remain outside this closure. H/I/E/J/K evidence, all S11c files and pinned
exports are preserved. No commit was requested or made. The stopping rule
for this classification now applies.

# Bounded homogeneous S11 fidelity review and closure

Closed 2026-09-16 UTC. Author and adjudicator: Codex.
**H1–H4 are complete** under [FORMALIZATION_POLICY.md](../FORMALIZATION_POLICY.md)
and [COVERAGE.md](COVERAGE.md). Both independent non-author reviewers returned
**CLEAR**, with no mathematical/fidelity blockers and no required control gaps.
This closes the homogeneous Lean contract, not the entire S11 series.

## Fixed revision and review evidence

The user explicitly authorized this 47-file S11 packet for Claude and Grok.
That authorization is retained in the
[review state](../../_measurements/S11_lean_fidelity_review_state.json).
Both received the same bytes in fresh, sequential sessions. Neither was supplied
the other's findings. Transport was outside the repository to avoid discovering
parent project instructions as extra context.

- [Manifest](../../_measurements/S11_lean_fidelity_review_packet.json):
  `S11 homogeneous bounded fidelity contract v1`, aggregate SHA-256
  `4541cf8243b2058f813a09cf316b165a6d730ea53153f2284e27b39e5d4e51cf`.
  The aggregate hashes UTF-8 `json.dumps(files, sort_keys=True)`.
- [Exact archive](../../_measurements/S11_lean_fidelity_packet_v1.tar.gz): SHA-256
  `aea065d5dd6e5a4563f45ac5db33919b7ccdf9a13349e83b1c0a18091f6b2c15`.
  It contains the 47 files and manifest, including the pre-closure prose.
  Paths are relative to `research/pde_ledger_v3`, except `docs/model_map.md`,
  copied from the repository root under that packet path.

| Reviewer | Session and completion (UTC) | Verdict and retained evidence |
|---|---|---|
| Claude Opus 5 (`claude-opus-5[1m]` in metadata) | `c33ca478-e576-4edf-89b0-5c2c0ce3c815`; 2026-09-16 03:14:56 | **CLEAR**; [report](../../_measurements/S11_lean_fidelity_claude_v1.md), [raw JSON](../../_measurements/S11_lean_fidelity_claude_v1.json), [run](../../_measurements/S11_lean_fidelity_claude_v1_run.json) |
| Grok 4.6 (`grok-4.6-build`) | `d76418e2-171b-4df4-8522-7b1581655595`; 2026-09-16 03:24:02 | **CLEAR**; [report](../../_measurements/S11_lean_fidelity_grok_v1.md), [raw JSON](../../_measurements/S11_lean_fidelity_grok_v1.json), [run](../../_measurements/S11_lean_fidelity_grok_v1_run.json) |

Both processes exited 0 with `end_turn`; Claude additionally reports
`completed`, `success`, no error and no permission denials. These are completed
source reviews. Claude's Markdown preserves its decoded `result`; Grok's
preserves `text` from the final report heading onward. Its preceding progress
sentences remain in the raw JSON. Reports are retained without editorial changes.
The extraction offsets and artifact hashes are in the
[closure validation](../../_measurements/S11_lean_closure_validation.json).

Both reviewers report no Lean/CAS/script execution or hash recomputation. Each
states its inspected sources and limits. The author separately recomputed all
47 live, snapshot, transport and archived file hashes before disposition edits;
all matched. The local source/dependency pins, instrument, saved control sources
and outputs, and four compiled S11 output hashes also matched verification.

Claude was limited to Read/Glob/Grep in safe/restricted mode; Grok was configured
with Read only, no subagents and no web search. Grok's
[stderr](../../_measurements/S11_lean_fidelity_grok_v1.stderr) contains startup
warnings about plugin-name precedence, unrelated hook formats and git discovery
under `/tmp`. There is no tool-allowlist mapping failure in this run. These
warnings did not prevent the completed source review. Packet bytes remained
unchanged. The supervisor log is empty and the completion hook finished normally.

## Findings and dispositions

All submitted findings were nonblocking. No proof, action, parameter domain,
verification instrument or native source needed repair. Documentation changes
below clarify the already reviewed claim; optional additions are not promoted
to new completion requirements.

| Finding | Disposition |
|---|---|
| Claude R1: optional mutation of the branch in `phase_matching` | No added fixture. The branch identity is proved and audited, `positive_frequencies` gives nonvacuous longitudinal inputs, the longitudinal-coefficient control distinguishes the two branches, and the threshold controls pin `normalWaveSq`. The reviewer explicitly finds no required control gap. We do not claim a separate source mutation of this hypothesis. |
| Claude R2: make `modalAction = dot a (E a)/2` explicit | Documented its expansion from `lagrangian_split` and the existing two `modalAction_eq` theorems in [FIDELITY.md](FIDELITY.md). Both reviewers independently checked that algebra. No new theorem is required for the handwritten-reference connection; none is claimed. |
| Claude R3: optional Wolfram curl-normalization anchor | Retain the seven existing anchors. Both reviewers directly inspected `curlDensity`, including its double-sum half, at the frozen source revision. FIDELITY.md now states that this body is not separately guarded by an anchor. More drift protection does not close an outstanding present-revision fidelity gap. |
| Claude R4: name the determinant's matrix | Changed the wording to the determinant of the selected `M_B`, explaining the existing `/8` normalization. |
| Claude R5 / Grok recommendation 3: distinguish SymPy's phase rewrite from integration | FIDELITY.md now states that SymPy replaces sin² by 1/2 for quadratic MAIN, while Lean/Wolfram integrate. The checked agreement is exact symbolic agreement, not a numerical sample. |
| Claude R6: sharpen the longitudinal-gradient wording | Clarified that curl-freeness removes the shear contribution and compression acts through divergence. Retained “can contain”/“generally”: an unconditional claim that the traceless gradient is always nonzero would omit zero amplitude and phase nodes. See the precision note below. |
| Claude R7: status/sentinel flags and the route-A phase label | FIDELITY.md now explains that PASS and sentinel absence summarize successful assertions; they are not additional evidence. `route_A_stripped_phase` is the original routine's convention label, not a computed phase factor. The actual phase identity is proved in Lean and the native matrix checked separately. No instrument rewrite is needed. |
| Grok recommendation 1: calls the nonzero-wavevector guard mutant a dummy; proposes negating the zero-wavevector identity instead | Retain the pair with a precise interpretation. At the chosen parameters `B=1`, `rho=1`, `cs=2`, the allegedly unrelated equality `1=4` is exactly the coefficient-locus conclusion `B=rho*cs²`. The passing control establishes that `normalWaveSq=0` while that conclusion is false. Together they give a counterexample to dropping the `k≠0` guard; the rejecting conjunct alone would not. This is a concrete boundary-counterexample control, not a source-deletion mutation of the theorem. The proposed replacement tests the `k=0` identity itself and is not required for the existing guard counterexample. |
| Grok recommendation 2: optional static-transverse exclusion control at `B=0` | No added fixture. Exact action/operator/stationarity identities recover S10, whose `zero_frequency_iff` proves the entire longitudinal line. The current static longitudinal positive control is a witness to that established result, not evidence of basis completeness by itself. |

Two precision notes apply to the retained review prose. For longitudinal
`a=lambda*k` in D3, the traceless gradient is
`-lambda*sin(phase) * (k kᵀ - K I/3)`. The matrix in parentheses is nonzero when
`k≠0`, but the gradient vanishes when `lambda*sin(phase)=0`. Thus the
unqualified “nonzero” phrasing in both reports is too strong at phase nodes or
zero amplitude; the canonical proof makes no such claim. Also, the first row of
Claude's section 3 table is shorthand: `transverse_kernel_iff` assumes the
transverse root and classifies the whole kernel. It does not infer that root
from stationarity of a possibly zero amplitude. The actual declaration and
coverage table retain the correct quantifiers.

## Completion and stopping boundary

| Obligation | Evidence |
|---|---|
| H1: action and compact identification | Supplied curl-plus-compression density, integrated compact-test variation, local PDE and modal operator. Both native MAIN constructors identified, with exact selected SymPy checks and inspected Wolfram source. `M_A=-E`, `M_B=E/2`; no resolvent/residue normalization claim. |
| H2: exhaustive modes | Full transverse/longitudinal/coincident/off-root kernels and D3 dimensions 2/1/3/0, positive frequencies, exact `B=0` recovery and separate `k=0` boundary. No witness-only subspace claim. |
| H3: kinematic threshold | Longitudinal phase matching and exhaustive below/on/above-grazing sign equivalences under the stated positive-coefficient and nonzero-wavevector assumptions. |
| H4: controls and review | Three modules plus audit root compiled sequentially; 29 audits use only `propext`, `Classical.choice`, `Quot.sound`. All 15 mathematical mutants rejected, all 13 explicit positive controls passed. Both complete non-author reviews clear, with dispositions above. |

[VERIFICATION.txt](VERIFICATION.txt) retains the pre-review timestamp; its
“reviews pending” line describes that earlier state. This document closes that
obligation. Only three files from the packet changed after review:
`COVERAGE.md`, `FIDELITY.md` and `README.md`, for status and the clarifications
above. All proof sources, dependency pins, native sources, check scripts and
verification results retain the reviewed hashes. The closure validation records
the old/new documentation hashes. Under L3/L4 these changes require neither a
new build/mutation run nor re-review of a materially changed formal claim.

H1–H4 meet their finite stopping criteria. Invariant-count completeness, other
packages, S11b interfaces, S11c calculations/files, actual coupling, confinement,
bound poles and leakage remain excluded. This contract supplies no closure of
the broader CAS/comparator/export program. No commit was requested or made for
this increment.

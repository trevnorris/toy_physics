# S9 bounded fidelity review and closure

Closed 2026-09-15 local time (2026-09-16 UTC). Author and adjudicator: Codex.
**C1–C4 are complete** under [FORMALIZATION_POLICY.md](../FORMALIZATION_POLICY.md).
Both independent non-author reviewers returned **CLEAR**, with no mathematical
or statement-fidelity blockers. This closes the bounded Lean contract in
[COVERAGE.md](COVERAGE.md), not the entire S9 ledger step.

## Reviewed revision and independent evidence

Both reviewers received the same fixed 49-file packet, `S9 bounded fidelity
contract v1`, in separate sessions. Reviews ran sequentially; neither was given
the other's findings. The user's explicit authorization, “Yes you can send the
files for review”, is retained in the
[review state](../../_measurements/S9_lean_fidelity_review_state.json).

- [Manifest](../../_measurements/S9_lean_fidelity_review_packet.json), aggregate
  SHA-256: `6903a1fa2cedf302561deab7dd04c6e88f523125b372a416cd6448ed4325d0eb`.
  This hashes UTF-8 `json.dumps(files, sort_keys=True)`.
- [Exact reviewed archive](../../_measurements/S9_lean_fidelity_packet_v1.tar.gz),
  SHA-256: `10eb027fa73f93f6e469e0c8950abd5dc765bc7a984093adea7dd446beb1e5d4`.
  It preserves all 49 files and the manifest, including the pre-closure prose.
  A narrow ignore exception makes the archive eligible for a future checkpoint.
- Paths in the packet are relative to `research/pde_ledger_v3`, except
  `docs/model_map.md`, copied from the repository root under that packet path.

| Reviewer | Session and completion | Verdict and retained evidence |
|---|---|---|
| Claude Opus 5 (`claude-opus-5[1m]` in model metadata) | `a00d350a-61b5-4b28-87a7-d77ef4811ff9`; 2026-09-16 01:28:43 UTC | **CLEAR**; [report](../../_measurements/S9_lean_fidelity_claude_v1.md), [raw JSON](../../_measurements/S9_lean_fidelity_claude_v1.json), [run](../../_measurements/S9_lean_fidelity_claude_v1_run.json) |
| Grok 4.6 (`grok-4.6-build`) | `e1dc69d7-d887-44fc-89ca-aed0216d4503`; 2026-09-16 01:37:52 UTC | **CLEAR**; [report](../../_measurements/S9_lean_fidelity_grok_v1.md), [raw JSON](../../_measurements/S9_lean_fidelity_grok_v1.json), [run](../../_measurements/S9_lean_fidelity_grok_v1_run.json) |

Both CLIs exited 0 with `end_turn`. Claude additionally reports `completed`,
`success`, and no error. These are completed reviews, not plans or partial
responses. Claude's Markdown is its full decoded `result` field. Grok's is the
decoded `text` from the final report heading onward; the preceding progress
sentences remain in the raw JSON. Neither report has been editorially corrected.
Extraction details and artifact hashes are in the
[closure validation](../../_measurements/S9_lean_closure_validation.json).

Both reviewers inspected source and recorded evidence, and report that they
did not run Lean, CAS or hash computations. The author independently recomputed
all 49 live, packet and archived file hashes before the editorial closure edits;
all matched. Reviewer declarations about hashes are not a substitute for that
check. The supervisor log was empty and the completion hook finished normally.

The reviews were instructed to inspect without changing files. Claude was
limited to read/search tools. Grok's retained
[stderr](../../_measurements/S9_lean_fidelity_grok_v1.stderr) warns that two tool
allowlist names were unmappable and the full toolset remained available; a
technical read-only restriction therefore cannot be claimed for that session.
Its report states no execution, and all reviewed and live packet contents were
unchanged when validated. This launch warning does not alter its completed
statement-fidelity assessment.

## Findings and dispositions

All findings below were nonblocking in the submitted reviews. No mathematical
definition, proof or verification instrument required a change.

| Finding | Disposition |
|---|---|
| Claude F1: five Wolfram anchors do not separately cover the EL sign or route-B normalization | Clarified [FIDELITY.md](FIDELITY.md): only the five listed strings are literal-checked by the instrument. The whole source, including the other lines, is hash-pinned by the reviewed manifest, and both reviewers confirm the current mapping. A changed source requires a new fidelity assessment. Two extra anchors are optional drift protection, not a missing current-revision identification; no checker expansion was made. |
| Claude F2: off-cone zero kernel has no separate named corollary or concrete control | Added the explicit shared-control mapping in FIDELITY.md. The existing `propagating_mode_iff` proves the necessary cone condition for every nonzero amplitude, so off-cone exclusion is its direct contrapositive. Action/operator sign controls and specialized S10 coefficient controls test the same operator/cone; the exact determinant check also constrains the root set. No separate off-cone fixture is claimed. Claude marks the extra pair nonblocking; Grok explicitly finds it unnecessary for clearance. L3 permits shared controls for the same mathematical claim. |
| Claude F3: “amplitude range” could suggest two separately declared range theorems | Changed RESULT.md to “cosine-amplitude range”, matching `cosAmplitude_range`. Both quadrature membership theorems remain correctly described. No sine-range lemma was added. |
| Claude F4: the source instrument records Lean hashes but compares CAS objects to handwritten references | Made that boundary explicit. The instrument neither parses Lean nor checks its recorded hashes against an accepted revision. Independent source review and closure hash verification establish correspondence for this revision; a future instrument PASS alone does not renew it. |
| Claude F5: a generic rewrite failure need not show mathematical falsity | Retained the inspected diagnostic, not just its acceptance regex. For `phase_wrong_coordinate`, choose `hbar=mass=epsilon=B=1`, `A=0`, spacetime point `x=0`, `omega=0`, `k=(1,0,0)`, component `i=0`. The mutated time-derivative side is 0 while the claimed spatial side is 1. Both reviewers inspected the intended mathematical failure; the extra unused-variable warning is not the accepted evidence. |
| Grok: optional off-cone and sine-range lemmas or further Wolfram anchors | No expansion: these suggestions do not close an unresolved obligation. The scope policy's stopping rule applies. |

Two phrases in Grok's report need precision without changing the report itself:
at `hbar=0` the amplitude range is **{0}**, not an empty set; its nonzero part is
empty. Reversing the shear sign is a sign control within the curl-action family,
not a replacement of the stiffness form. The canonical definitions and author
records already use the correct meanings.

## Completion evidence and boundary

| Obligation | Evidence and result |
|---|---|
| C1: supplied action, complete modal classification, compact original-engine connection | Reused original integrated action/PDE and full-subspace theorems, plus exact S10 specialization. Three disjoint exhaustive frequency cases under positive coefficients and `k ≠ 0`. Exact native SymPy action/operator/determinant checks pass; Wolfram correspondence is source inspection at the pinned revision. Both reviewers clear the parameter map, sign, normalization and interpretation. |
| C2: scalar-phase velocity | `Madelung.lean` derives the actual spatial gradient and first perturbation derivative; both quadratures are longitudinal, with trivial transverse intersection for `k ≠ 0`, a nonzero witness and the `k=0` limit. Both reviewers clear this conditional kinematic claim. |
| C3: meaningful controls | Nine mathematical mutations rejected and five explicit positive controls compiled. Recorded diagnostics were inspected by both reviewers. Initial instrument attempts remain disclosed; they did not repair or alter the canonical mathematics. |
| C4: build, axiom audit and independent fidelity reviews | All nine S9 modules rebuilt sequentially and the audit root compiled; 32 selected axiom audits use only `propext`, `Classical.choice`, `Quot.sound`. No admissions or custom physics axioms. Both review legs complete and findings disposed above. |

The [verification record](CLOSURE_VERIFICATION.txt) and
[full control results](../../_measurements/S9_lean_contract_checks.json) retain
their pre-review timestamps. The verification record's “review pending” line
describes that earlier state; this document records its completion.

After review, only four files from the packet changed: `COVERAGE.md`,
`FIDELITY.md`, `README.md`, and `RESULT.md`. These are status, shared-control-map
and wording clarifications of the reviewed claim. Closure validation records
their old and new hashes. All packet Lean modules, dependency pins, CAS sources, check
scripts and verification results still match the reviewed revision. Per L3,
these documentation-only changes require no repeated build or mutation suite;
they introduce no materially changed physics-bearing claim needing re-review.

The conditional curl-action result is classical MacCullagh calibration; the
Madelung addition is the scalar phase-gradient constraint under its stated
ansatz. Neither derives the supplied physical premises. General PDE/Fourier
completeness, GNLS dynamics, the full spin-1/P2 assertion, ledger/export/comparator
reconciliation and S11 work remain outside this closure. C1–C4 now meet their
stopping criteria; no further Lean extension is required. This closure record
was prepared before the user-authorized checkpoint; consult git history for
the commit containing it.

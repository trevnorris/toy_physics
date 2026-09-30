# Fresh radiating-method review disposition — 2026-09-30

Both requested reports were delivered against the exact approved 37-file packet
`ad792d2e94fdbda24f17f633c22ae371c7c0bb78535a725513909a24edb6c4e1`.
The packet, archive, index, source hashes, independent session receipts and literal
final fields were verified. Claude stderr is empty; Grok's 2,863-byte stderr is
CLI plugin/hook diagnostics, with a complete final report delivered separately.
Raw reports are preserved in `S11c_d_numerical_radiating_replacement_review_claude.json`
and `_grok.json`; the full receipt/source disposition is in `_record.json`.

| Route | Claude literal verdict | Grok literal verdict |
| --- | --- | --- |
| Original tolerance | NEEDS REVISION | NEEDS REVISION |
| Prospective coarse branch | NEEDS REVISION | CLEAR SUBJECT TO THE STATED CURRENT-SIDE CONTROLS |

Both assess the finite-model transverse-current observable as coherent, subject
to the stated controls. Neither accepts a computed deficit or a physical leakage
claim. This local disposition is not independent clearance.

## Findings and actual current status

1. **Matrix quadrature path: missing packet evidence, not a demonstrated wrong
   rule.** The omitted parent source is
   `scripts/S11c_d_mixing_scattering_sympy_audit.py:1548–1570` under
   `research/pde_ledger_v3`. `BasisMomentum` inherits `ThreeMomentum.batches`
   unchanged. It calls `self.rule`, so the actual `Radiating.rule` override supplies
   the sin/cosh rule. Its centered inner-order selection matches `row_integral`.
   Three stdlib-only tests confirm the inheritance and exact traversal, order,
   point and weight sequences using stand-ins for one/two/three variables. No
   scientific imports, restoration or quadrature was run. A redundant batches
   rewrite is not justified. The requested actual Gaussian matrix-route comparison
   remains pending; a source-dispatch test does not replace that numerical check.
2. **Separate single-leg orders and monotone refinement: already implemented
   after the frozen packet.** The running `balance_v2` uses 48/base and
   64/refinements for single-variable groups only. Nested outer24/32 and inner8/12
   remain unchanged. The existing worker's tests pin this behavior. Its source
   and helper pins are left intact while it runs.
3. **Focused trial identities: actual run supports the intended sequence.** The
   saved four trials are (32,512), (48,512), (64,512), (64,768), all passing the
   original action tolerance; source-node movement also passes. The balance worker
   checks that set and every pass, rather than relying solely on the old support
   bit. A subsequent contained check must still assert actual row dimensionality
   and ordered settings from the saved operands.
4. **Momentum-zero and matrix-path Gaussian checks: pending.** Reuse saved source
   operands, row-integral returns and independent references. Check momentum zero
   at the adopted single-leg order, and compare actual single-leg and selected
   row70 Gaussian actions through the matrix route against their saved returns.
   Do not recompute independent references, accepted source construction, modes,
   paths, current maps or the previously completed row-integral calculations.
5. **Reporting corrections: accepted.** In the interpreted report, render legacy
   `NO_DEFICIT_RESOLVED` as `DEFICIT_NOT_NUMERICALLY_RESOLVED`; also report the
   separate `D < -E_num` flag, condition number and each setting's raw deficit.
   Preserve original result fields and every existing partition/sign/uniform gate.
   This is a reporting correction, not authorization to loosen any threshold.
6. **Coarse branch: unused.** Refinement passed the original tolerance. Claude and
   Grok disagree on handling a coarse uniform deficit between 1e-6 and 1e-5;
   their disagreement is retained. Neither alternative is adopted. The original
   hard 1e-6 uniform gate remains unchanged. No coarse policy or extra swap solves
   are needed for the present route.

## Execution decision and stop

At review inspection `central-balance-v2` was running in `base/group-2` with 22
complete operations and no completed finite-current solve. Leave that process
and its pinned sources intact: the parent source inspection found no wrong-rule
execution requiring an abort or restart. Completed work remains reusable.

Before interpreting a number or launching neighboring frequencies, complete the
narrow pending joins/checks above under the normal single-worker containment,
inspect the actual finite uniform/scaling/refinement/domain/regulator/current
partition/sign results, and reconcile all required findings. The completion
handler must read this disposition and the current execution plan. No extra
scientific job was launched during this review inspection, no current result is
accepted, and no reviewer rerun or external export is authorized by this record.

Remaining limits stay explicit: finite modal exterior and finite source interval;
shared inner rules in the middle comparison; no Gaussian coverage of all high-k
basis content; empirical current uncertainty rather than a physical bound;
retained-order physical-loss interpretation and analog-light calibration open.
Optional stronger controls are not turned into a new method campaign. The exact
omega1 track, full Green/FORM/A11/A12 and loss-mechanism separation remain parked.

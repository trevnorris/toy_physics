# Decision-list review — #90 PY §3c coupling content directive

You are one of two independent decision legs reviewing an ORCHESTRATOR-written builder decision list BEFORE any
engine is changed. Find every way it is wrong, incomplete, over-reaching, or leaks an expected value — the builder
trusts this list and nothing downstream re-checks it. Derive and quote; a leg that finds nothing is weak evidence.

## Artifact under review
`research/pde_ledger_v3/directives/S11c_b_90_coupling_content_directive.md`

## SETTLED — do NOT re-litigate
The §3c CONTENT verdict is INCLUDE/INCLUDE (WL spec-correct, PY under-extracts the reversible face `A_T` + the
irreversible response `A_T·Λ(ω)`) — adversarially confirmed (Codex+Grok ×2 rounds) and user-endorsed. You are NOT
asked whether to include face/response — that is closed. You ARE asked whether the directive correctly, completely,
and cleanly implements it in PY. Settled records (context, do not re-argue):
`research/pde_ledger_v3/directives/_measurements/S11c_b_coupling_84_{diagnosis,consult,basis_verification}.md`.

## What you are handed (read the sources, form your own view)
- The directive above.
- Spec: `research/pde_ledger_v3/directives/S11c_b_SHARED_PHYSICS.md` §0 (the permeability/memory-kernel exclusion),
  §1c (`Λ_I(ω)=Λ_I⁰/(1−iωτ_I)`; the balance-law method; "not by putting an irreversible response kernel in an
  ordinary action"), §3b (consume T-a..T-i for every face/boundary contribution), §3c (the weak-restriction
  instrument + its prohibitions: extract from the operator itself, no parallel direct-variation route, no
  single-channel filter, adjointness residual is NOT `∂²U`).
- `research/pde_ledger_v3/directives/S11c_a_SHARED_PHYSICS.md` T-i (`:~448`): the flat-face `Λ` is NOT B0c's
  bulk-response solve `δp=Z·v_bulk`.
- PY engine: `research/pde_ledger_v3/scripts/S11c_b_brane_operator_sympy_audit.py` — `task_coupling_kernel`
  (`:3874`), `build_kernel` (`:3254`), `weak_operator_blocks` (`:3404`, note it contracts only the operator's
  EXPANDED bulk-EL rows), `substrate_bundle`/`FACE_FLUX_BOUNDARY_OPERANDS` (`:2058`, the T-a..T-i face substrate PY
  already builds), the dead `bulk_kernel_from_density` (`:3177`)/`paired_kernel_from_density` (`:3227`), and the
  coupling depth `min(background_depth, STRONG_ROW_JET_DEPTH)` (`:3286`).

## Decide, WITH code and spec quoted for each claim
1. **Correct + complete + no over-reach.** Is directive change (2) the correct and COMPLETE realization of
   INCLUDE/INCLUDE — does applying the §3c weak restriction to the operator's FACE/FLUX contribution (`Λ` symbolic)
   yield the reversible face + irreversible response block §3c mandates, while RETAINING the prior bulk content and
   NOT over-reaching (no bulk elimination/DtN, no single-channel filter)? Is the root diagnosis correct (the emitted
   kernel is `build_kernel`→`weak_operator_blocks`, which ignores `FACE_FLUX_BOUNDARY_OPERANDS`)? Name any ripple the
   directive omits (the folded operator's new rows; the coupling depth cascade; the tower-depth control; the kernel
   `TERM_ORIGINS`; the §5c uniform limit that must still return the S11b decoupled zero; the dead functions).
2. **§0 pin accuracy.** Is the §0 clarity sentence accurate and consistent with §1c + the T-i seam (the supplied
   flat-face `Λ` stays in every face/flux contribution; only the S11c-c curved-bulk `δp=Z·v_bulk` solve is excluded)?
   Does it contradict §0/§1c anywhere?
3. **Rule-5 cleanliness.** Does ANY clause leak an expected value — a target kernel term, a term count, a specific
   `A_T`/`Λ`/`γ` structure, a channel, a sign, or a "match-WL" exit condition? Naming the SUPPLIED objects (`A_T`
   geometry, `Λ` kernels, bulk `γ·profile-jet`) is not a leak IF the acceptance does not gate on kernel content —
   judge whether it does. Name any leaking clause verbatim.
4. **Consistency.** Does anything contradict the #84 settled verdict, the T-i seam (Λ stays symbolic, no bulk
   solve), or §3c's prohibitions (operator-itself extraction, no channel filter, adjointness ≠ `∂²U`)?

## Output
Per question: CORRECT / DEFECT-with-fix / AMBIGUOUS, each with the governing spec/code line quoted. Give exact
replacement text for any change you would make. A prose claim without a quoted citation (or, for a spot-check, a
script + literal stdout) is discarded. Do not edit any file.

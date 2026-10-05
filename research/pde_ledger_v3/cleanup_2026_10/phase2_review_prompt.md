# Independent review — S11c cleanup, Phase 2 ledger records

## Artifact
Codex-written ledger records in the repository `/var/projects/toy_physics`, commits `a442df4e`, `b81d87a5` and
`8e0a8ee1`. View them with `git diff c8885c4f 8e0a8ee1`. The files:
- `research/pde_ledger_v3/steps/S11c_d_profile_conditioned_scattering.md` (new)
- `research/pde_ledger_v3/steps/S11c_PARTIAL_CLOSEOUT.md`
- `research/pde_ledger_v3/steps/S11c_b_variable_coefficient_operator.md`, its added "Changed after close" section
- `research/pde_ledger_v3/steps/S11c_c2_self_energy_fold.md`, the same
- `STATUS.md` (rewritten from 1,289 lines to 49)
- `research/pde_ledger_v3/V3_STEP_PLAN.md`, the S11/S11c sections
- `research/pde_ledger_v3/paper/parts/part01_light.tex`, the new S11c section
- `docs/s11_maccullagh_differentiation.md`, rewritten sentences
- `research/pde_ledger_v3/steps/S10_two_transverse_photons.md` and `research/pde_ledger_v3/lean/s11/README.md`
- `research/pde_ledger_v3/cleanup_2026_10/FINDINGS.md`, its new opening

## What to check
These records are physics-bearing. They state what the ledger now claims about light in this model. Report
anything that changes **what may be claimed**:
1. **Fidelity to sources.** Every result statement must match the source report it cites, including numbers,
   scope and qualifications: model point, frozen quantities, one engine vs two, and review status with literal
   verdicts. Flag any upgrade: "unresolved" becoming a result, "not detected" becoming "no leak", "scoped clear"
   becoming "cleared", or a CONDITIONAL result stated without its conditions. Flag any downgrade or omission that
   drops a still-valid result.
2. **The c2 mixed term.** The confirmed source fact is that `scripts/S11c_c2_selfenergy_fold_sympy_audit.py`
   (around lines 400–418) sets the direct three-leg `[0,2]` entry to zero while keeping the iterated product.
   Whether the complete closed response has a nonzero direct term is unresolved. Check that the records say
   exactly that, neither more nor less. Read the code and the joint disposition yourself.
3. **"Changed after close" sections (S11c-b, c2).** Do they list every post-close repair with its commit, its
   change, the scope of its review and its literal verdict? Do they make clear that the original verification
   text describes the pre-repair version?
4. **STATUS.md and the plan.** Is anything still open or owed from the old 1,289-line STATUS missing from the new
   front door? Compare it with `git show c8885c4f:STATUS.md`. Is there any stale "NEXT" or "in flight" left in the
   S11 sections of the plan? Are the plan's real downstream dependencies on S11c preserved?
5. **Plain-language openings** (FINDINGS, the S11c-d record, the closeout, STATUS). Are they accurate, and could a
   programmer who is not a physicist follow them? An opening that is readable but wrong is a finding.
6. **Paper.** Is the S11c section consistent with the records? Are its limits in the main text and not hidden in a
   macro field that the default build suppresses (`research/pde_ledger_v3/paper/macros.tex`)?

## What you are handed
The repository and its full history, including the tag `archive/pre-cleanup-2026-10-04`. The governing
instructions are `research/pde_ledger_v3/CLEANUP_2026-10_directive.md` and
`research/pde_ledger_v3/cleanup_2026_10/PHASE1_REVIEW.md`.

## Required method
This is a DOCUMENT review. For each record, read the sources it cites **first**, form your own view of what they
establish, and only then read the record. Quote both sides for every finding: the record's sentence and the
source's sentence, with paths and line numbers. For any claim about file content, show the command you used (for
example `grep -n` or `sed -n`) and its literal output. A finding that rests on paraphrase alone will be discarded.

## Physics filter
Report a finding only if it changes what the ledger claims, or would let a later step rely on something that is
not established. Do not report style, formatting or wording preferences.

## Sandbox
Read-only. ⛔ Never modify the working tree, commit, or run CAS engines or numerical workers.

## Bounds
Write your report and exit. ⛔ Do not spawn agents or build orchestration. Report format: a verdict line (CLEAR or
NEEDS REVISION), then the numbered findings with quotes and commands, then anything you checked and found sound.

# EM projected integration tracker

Created: 2026-05-11
Updated: 2026-05-11

## Scope decision

Use `notes/em_projected` as derivation source only through the derivation notes
ending at step 18.  The directory intentionally has no step 06 file.  Do
not import the later computational/runtime material into the PDE ledger paper:

- step 19 branch-export / Galerkin negative-control packet;
- step 20 and later reduced-family scans/frontiers;
- electron/SI falsification, dimensional runtime maps, CFD postprocessors,
  fail-fast classifiers, snapshot adapters, monopole JSONL screens, and related
  runtime tools;
- current simulator/CFD verdicts.

Those later files may remain useful as private diagnostics, but they are not
paper derivation support in this integration pass.

## Current ledger structure

The projected Maxwell derivation is now represented as ordinary global stages:

| Stage | Source note | Ledger role |
|---:|---|---|
| `004` | `step_01_projected_maxwell_readme.md` | Bundle index and derivation boundary. |
| `005` | `step_02_projected_maxwell_covariant_notes.md` | Covariant projected Maxwell law. |
| `006` | `step_03_projected_maxwell_vector_notes.md` | Vector projected Maxwell split. |
| `007` | `step_04_projection_reduction_comparison_notes.md` | Projection/reduction comparison. |
| `008` | `step_05_projected_maxwell_extension_notes.md` | Extension terms. |
| `009` | `step_07_projected_maxwell_near_throat_notes.md` | Near-throat packet. |
| `010` | `step_08_projected_maxwell_push_bundle_master_notes.md` | Push-bundle master. |
| `011` | `step_09_projected_maxwell_p2_bridge_notes.md` | Grouped `P_2` bridge. |
| `012` | `step_10_projected_maxwell_stage4_primitive_bridge_notes.md` | Primitive bridge. |
| `013` | `step_11_projected_maxwell_mouth_taylor_master_notes.md` | Mouth-Taylor master. |
| `014` | `step_12_projected_maxwell_mouth_taylor_gate_bridge_notes.md` | Mouth-Taylor gate bridge. |
| `015` | `step_13_parent_throat_action_master_notes.md` | Parent throat action master. |
| `016` | `step_14_parent_throat_action_candidate_notes.md` | Parent throat action candidate. |
| `017` | `step_15_parent_throat_action_weak_axisym_notes.md` | Weak-axisymmetric parent limit. |
| `018` | `step_16_parent_throat_action_bundle_master_notes.md` | Parent action bundle master. |
| `019` | `step_17_parent_throat_action_isotropic_bundle_notes.md` | Isotropic parent bundle. |
| `020` | `step_18_parent_throat_action_weak_axisym_packet_notes.md` | Weak-axisymmetric packet. |
| `021` | `notes/stages/moving_throat_pde_stage021_reduced_one_port_normal_form.md` | Reduced one-port normal form retained from the prior EM reduction stage. |

Former post-EM ledger stages start at Stage `022`.

## Verification policy

- SymPy audits are top-level files named
  `scripts/moving_throat_pde_stageNNN_*_sympy_audit.py` and should match the
  paper stage order.
- Projected Maxwell Stages `004--020` do not currently claim Mathematica
  mirrors.  Stage `021` retains the reduced one-port Mathematica audit.
- Mathematica audits must be regenerated one script at a time.
- Paper `\StageFile{...}` entries should resolve to real source files and to
  regenerated transcript files after the audit runs.

## Completed verification checks

- Full SymPy audit regeneration completed:
  `TOTAL: 241  PASS: 241  FAIL: 0  SKIPPED: 0`, including the repo-level
  master audit and `240` stage-level audits.
- Current Mathematica transcript set regenerated one script at a time:
  `TOTAL: 165  PASS: 165  FAIL: 0`.
- Numerical stress transcript summaries are current for both Python and
  Mathematica harnesses.
- Missing `\StageFile{...}` and stale old-numbering scans were run against
  paper-facing files and active metadata after the linear renumbering.
- `paper/pde_ledger.tex` was rebuilt in the post-renumbering verification pass.

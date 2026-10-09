# S9b spec amendment 1 and build directive amendment 4: review dispositions (orchestrator)

**Artifacts.** Two orchestrator-written amendments, reviewed together:
- amendment 1 to `directives/S9b_SHARED_PHYSICS.md`, on top of the Parts A–C acceptance at `97580010`;
- amendment 4 to `directives/S9b_repair_build_directive.md`, on top of `3699f4ef`.

**Why.** In the build's repair round 1, the SymPy builder stopped under directive item 12 after a second method
failure. For general radial profiles (directive item 13), it could not extract "the coefficient of `ln(1/b²)`" that
the spec named as the radar comparison object. Its report is `_scratch/s9b_build/s9b_py_amend3_stop_report.md`.
Its probe and literal stdout are in `_scratch/s9b_build/repair_build/py_amend3_stop/`. The builder changed no engine,
harness or delta (lookups). This is a gap in the spec, so it was repaired in the spec.

**The repair.**
- The spec now names the radar comparison object as the round trip's logarithmic slope in `b`, taken at fixed
  endpoints, for far endpoints: `𝒮_RT(b) ≡ lim_{Z_E,Z_R→∞} ∂Δt_RT/∂ln(1/b²)`.
- Where the round trip is a pure `ln(1/b²)` term plus terms that do not matter, the slope equals the old
  coefficient.
- The radar claim is limited to this slope.
- The directive's amendment 4 renames the object in items 4, 7, 10 and 14 and in Part 3.

This repair falls under the user's standing approval of technical repair rounds (2026-10-09). It changes what the
engines compute and what may be claimed, so it is reviewed until clear. The legs are Codex and Grok (G1).

Commands and literal output are in `S9b_SHARED_PHYSICS_amend1_review_lookups.md`. Generator:
`_scratch/s9b_build/gen/s9b_spec_amend1_r0_lookups.sh`.

## Round 0

Both legs used the identical prompt `_scratch/s9b_build/s9b_spec_amend1_review_prompt_r0.md`. They reviewed the
working-tree versions at the baseline `_scratch/s9b_build/s9b_spec_amend1_review_baseline_r0.sha256`: spec
`215ae8d7…`, directive `77869b47…`. Frozen copies are at `_scratch/s9b_build/S9b_SHARED_PHYSICS_amend1_reviewed_r0.md`
and `…/S9b_repair_build_directive_amend4_reviewed_r0.md`.
- **Codex** (gpt-6.1-sol, xhigh): CLEAR, no findings. Final report
  `_scratch/s9b_build/s9b_spec_amend1_review_r0_codex_final.txt`. Full report and evidence:
  `…_r0_codex_evidence/REVIEW.md`, with derivation scripts `derive_radar_v2.py` and `derive_radial_boundary.py`.
- **Grok** (grok-4.7): CLEAR, no findings. Report `_scratch/s9b_build/s9b_spec_amend1_review_r0_grok.txt`. Evidence:
  `…_r0_grok_evidence/radar_slope_derivation.py` and its stdout.

Both reported before I adjudicated. Neither document still names the old object (lookups).

**Recorded, not findings.**
- **The radar claim widens.** Both legs note that the slope also carries `b`-dependence from tails other than `1/r`,
  which the old logarithmic component dropped. This is the intended change, and the amended spec states it.
- **Radar and deflection may not be independent tests.** Codex's derivation finds that, at the compared order, the
  radar slope is fixed by the deflection, so the two effective `γ`s are redundant tests in this setting. The
  documents do not promise independent constraints. The engines compute both objects through separate
  constructions. The record must not present radar and deflection as independent constraints unless the engines'
  outputs support that.
- **A second-order piece.** Grok finds one piece of the slope that diverges, coming from a net kink in the ray. It
  is second order in the kink, so outside the retained grades.
- **K6's label.** K6 is still called a "coefficient knife". That is the pipeline's class for a knife no FORM knife
  can replace, not the old object. Its mutation now replaces the object whose slope is taken.
- **The builder's probe.** Both legs reproduced its two non-converging extractions. Neither extraction is `𝒮_RT`.

The legs' constructions, including their example profiles and the slope relation above, are reviewer measurements.
They stay out of builder packets and leg prompts.

**Outcome.** Nothing outstanding changes what is computed or what may be claimed. Both amendments are accepted at the
reviewed versions.

# O2 build r1: review dispositions (orchestrator)

**Artifacts:** the five O2 engine and harness files at `_scratch/s9b_build/o2_build_review_baseline_r1.sha256`,
the repair-round-1 output. They are preserved unaccepted as r1. The round-1 repairs answered every r0 finding
(`O2_build_r0_review_disposition.md`).

**Legs (Codex-written → fresh Claude + Grok, identical prompt `o2_build_review_prompt_r1.md`).** Both reported
before any adjudication.
- **Fresh Claude (opus): NEEDS REVISION, minor.** Seven findings
  (`_scratch/s9b_build/o2_build_review_r1_claude.md`). "None of it changes a computed object."
- **Grok: NEEDS REVISION.** One finding (`_scratch/s9b_build/o2_build_review_r1_grok.txt`).

**Agreement.** Both legs' own constructions match every explicit object in both engines with residual 0. The
Wolfram check is by PIT, at most 4.4e-16. The objects are the geometry, the material velocity, the mass law, the
in-plane carried momentum and the non-OPEN part of the hold. No OPEN operand is closed and no live quantity is
frozen. Every knife bites in both harnesses, no knife is all-zero, and K1, K3, K12 and K13 are selective as their
sites predict. Each leg also ran its own FORM ablations at sites off the knife list:
- Claude: an untilted graph normal, a non-radial graph, and an in-plane-only pairing velocity;
- Grok: an induced-measure divergence.

All of them bit, with diffs that equal the legs' independent expressions. The export has `IMPORT_KEYS = ()`, the
`o2_rho_br_live` rename, six roots at `F9A_ABSENT`, `ROUNDTRIP` all true, and Dummy-only drift between runs.

Each verification is a mechanical lookup. The commands and their literal output are in
`O2_build_r1_review_disposition_lookups.md`. I ran no CAS (E1).

| # | Finding (leg) | Disposition | What must be true after repair |
|---|---|---|---|
| C1 | The tag sets are not parallel: PY 32, WL 12. A hand-written key join is where a re-pointed name could manufacture agreement. (Claude) | **ROUTED to the comparator (sub-step 6), not an engine repair.** The r0 disposition stands. The blind WL builder cannot see PY's names, and making the sets parallel needs a shared vocabulary the directive did not supply. The comparator's reviewed join map must be **injective** and must pair objects by their computed content, never by name alone. Its review checks for vacuous joins. | (comparator) |
| C2 | WL represents each native face as a graph `qFace[s][x,t]` over the far-field coordinates. That excludes faces that are not single-valued over x, and the restriction is not stated. PY uses a general immersion. (Claude) | **ACCEPT, WL.** Lookup: WL lines 81–84 and 91. The spec leaves native geometry and measures to `𝒥_map` (§5). | Every restriction an engine's face representation places on the OPEN native face geometry is either removed or printed as a domain qualification on the objects it affects. |
| C3 | PY's `material_pairing = internal.dot(graph_velocity)` is only the graph-velocity contraction, and the normal generalized rates the graph does not fix stay OPEN. It could be read as the full material stress power, and WL emits no such contraction. (Claude) | **ACCEPT, PY (labelling).** Lookup: PY lines 258–261; WL lines 240–244. | The emitted pairing states what it contracts, and states that the unfixed normal generalized rates remain OPEN in the material-work action. |
| C4 | PY reduces the load with four independent OPEN operators (`FaceLoadReduction_0..3`) and the work with a fifth (`FaceWorkReduction`). The shared map exists only as an operand name. WL applies one `NativeToCoordinateDensity` to both. Spec §6: "with the same geometric map as its" force occurrence. (Claude) | **ACCEPT, PY.** Lookup: PY lines 225–227 and 255; WL lines 160–164 and 239; spec line 285. | The load and the face work are reduced by one and the same OPEN map, applied to their respective native integrands. |
| C5 | PY prints `OPEN_c_gamma(r)`, but spec §3.1 lists `c_γ` among the supplied optical-regime live identifications. (Claude) | **ACCEPT, PY.** Lookup: PY lines 286–287; spec line 111. | Supplied identifications are not labelled OPEN. |
| C6 | Both engines type the graph normal `(−∇ξ,1)/√(1+|∇ξ|²)` by hand, and WL also types its face normal. Both are correct, as verified by both legs, but a shared typing slip would hide behind agreement (structural rule: build-directive item 2). (Claude) | **ACCEPT, both (low).** Lookup: PY line 114; WL line 78. | Each normal is reached by computation from the computed tangents, not typed. |
| C7 | The K6 sites differ. PY replaces the profile in every momentum OPEN argument, and WL only in the differentiated section. Both bite. (Claude) | **ACCEPT, both (labelling).** Lookup: PY harness lines 37–42; WL harness lines 26–31. | Each harness's output records what its K6 knives replace and where, so the two engines' triples are not read as equivalent. |
| G1 | PY's `recorded_grades` binds δ, (∂ξ_w)² and V/c₀ to one bare `Symbol('epsilon')` (and its root), not to the computed `ε(r) = GM/(c₀²r)` printed beside them. That prints a stronger qualification than spec §7's three separate order statements. WL prints them as separate `OpticalGrade` objects. (Grok) | **ACCEPT, PY.** Lookup: PY lines 384–389; spec lines 323–325. | The model-point grades are printed as separate order statements in terms of the computed `ε(r)`. No two graded quantities are bound to one object. |

**Routing (G4).** No round-1 repair bred a defect: every r1 finding is in material the round-1 repair did not
introduce. So the same builders repair: WL gets C2, C6 and C7; PY gets C3–C7 and G1. The WL brief holds no PY
content.

**Round count.** This is review round 3 for the build. Its findings are labelling and representation, and no
computed object is in question. If round 3 does not clear, the orchestrator stops and brings the build to the
user. It does not fold a fourth time.

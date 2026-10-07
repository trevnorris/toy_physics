# O2 build r3: review dispositions and acceptance (orchestrator)

**Artifacts:** the five O2 engine and harness files at `_scratch/s9b_build/o2_build_review_baseline_r3.sha256`.
- The Wolfram engine and its harness are the output of repair round 3, which answered the four r2 findings
  (`O2_build_r2_review_disposition.md`). The user chose the same builder.
- The SymPy engine, its harness and `O2_exports.py` are unchanged from r2 (`42e85ccf`).

**Legs (Codex-written → fresh Claude + Grok, identical prompt `o2_build_review_prompt_r3.md`).** This is build
review round 4. Both legs reported before any adjudication.
- **Fresh Claude (opus): CLEAR.** No findings (`_scratch/s9b_build/o2_build_review_r3_claude.md`).
- **Grok: CLEAR.** No findings (`_scratch/s9b_build/o2_build_review_r3_grok.txt`).

**What the legs established.** Each leg built the objects independently, from the spec, before opening the
artifacts. Every closed §9 object in both engines matches that construction:
- the metric, inverse and `det g = 1+ξ′²`;
- the graph normal, computed from tangents;
- the material velocity, with `U^w = V_r ξ′` and a graph-normal component of 0;
- `(V·∇)U`;
- the coordinate-measure mass residual `∇·(ρV)+j_n`, emitted unsolved;
- the in-plane carried momentum `j_n V^i`.

The OPEN operands stay OPEN:
- the bulk carry, the momentum density and flux, and the internal force;
- the native face normal, tangents and measure, as `𝒥_map` actions in WL and from a general immersion in PY;
- the bulk-normal amplitude, now a native-point field;
- every §6 energy operand.

Every K1–K13 site occurs once, every knife bites, no knife is all-zero, and the selective knives (K1, K3, K8,
K13) move only the objects they should. Each leg ran its own FORM ablations at sites off the knife list, in both
engines. All of them bit:
- untilted graph normal;
- true-area element dropped;
- induced-measure mass law;
- frozen metric tangent;
- radial divergence without the spherical term;
- the metric installed as its own inverse.

Neither engine asserts an assembled balance, chooses a closure, or contains an expected value. `IMPORT_KEYS = ()`,
there are six `F9A_ABSENT` roots, and the upstream name matches are coordinates only.

The lookups confirm that the four r2 constructs are gone from the r3 Wolfram engine (counts 0): the face-label
amplitude, the `"Relation"` fields, and the height chart. K10's repair is attested by both legs' harness runs.
Commands and output are in `O2_build_r3_review_disposition_lookups.md`. I ran no CAS (E1).

## Notes from the legs (not findings), adjudicated

| Note | Disposition |
|---|---|
| Claude's seven named representational differences: face geometry (PY immersion vs WL OPEN `𝒥_map`), transport calculus, material power, where the graph normal enters, K2 on the face load, K9 site, and tag grain (PY 32 vs WL 11). | **ROUTED to the comparator (sub-step 6).** As with r1 C1, the comparator's join map must be injective and must pair objects by computed content, never by name. Ablation triples from different knife sites (K6, K7, K9) are not equivalent objects and must not be joined as agreement. |
| The WL K7 replacement is `orientation[s] unitNormal[…untilted tangents…]`, and r3 no longer defines `orientation` (count 0). | **Not outstanding.** The replacement is the untilted normal up to an unspecified sign, which is the directive's K7 form ("replaced by the untilted normal", directive line 57). The baseline normal is an OPEN `𝒥_map` action, so the bite shows the load does not use the untilted normal, whatever the sign. Recorded so the K7 triple is read with a symbolic sign factor. |
| `O2_exports.py` regenerates with different bytes, differing only in SymPy `Dummy` indices. | **Known (r1).** Dummy-only drift; the export round-trip is all true. |

## Outcome

**ACCEPTED: the O2 engines at r3** (baseline `o2_build_review_baseline_r3.sha256`). Both legs are clear, and
nothing outstanding changes what is computed or what may be claimed (G4).

Review history: round 1 (r0, 7 findings), round 2 (r1, 7 + 1), round 3 (r2, 4 Wolfram findings; the user chose a
round-3 repair by the same builder), round 4 (r3, clear).

**Next:**
1. Production runs: the `out/*.out` transcripts, via datalad.
2. The comparator (sub-step 6), carrying the routed items above.
3. The record (sub-step 7).

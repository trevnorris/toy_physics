# O2 build r0: review dispositions (orchestrator)

**Artifacts:** the five O2 engine/harness files at the sha256 baseline
`_scratch/s9b_build/o2_build_review_baseline.sha256`, preserved unaccepted as r0.
**Legs (Codex-written → fresh Claude + Grok, `CLAUDE.md` G1):** both reported before any adjudication. The prompt
was `_scratch/s9b_build/o2_build_review_prompt.md`, identical for both legs.
- Fresh Claude (opus): **NEEDS REVISION**, seven findings (`_scratch/s9b_build/o2_build_review_claude.md`).
- Grok: **CLEAR**, no findings (`_scratch/s9b_build/o2_build_review_grok.txt`).

**Agreement.** Both legs found that every explicit object matches their own construction with residual 0: the
graph geometry, the material velocity `(V^i, V·∇ξ_w)`, the coordinate-measure mass law and the in-plane carried
momentum `j_n V^i`. Both found that the bulk-direction carried component stays OPEN, that no body force is
present, and that every knife changes at least one tag in both harnesses. Grok's CLEAR brings no evidence
against findings 1–3 or 5–7: it did not examine them. On finding 4, both legs record the same fact, quoted below.

Each verification is a mechanical lookup (grep/sed). The commands and their literal output are in
`O2_build_r0_review_disposition_lookups.md`. The leg reports hold the CAS evidence; I ran none myself (E1).

| # | Finding (leg) | Disposition | What must be true after repair |
|---|---|---|---|
| 1 | Both engines leave a free native-face label `s` in the reduced load and the face work. PY: `s_face`. WL: `THold[s]`, `qFace[s]`, `orientation[s]`. No aggregation over faces. (Claude) | **ACCEPT**, both engines. Lookup: PY lines 158/171/200; WL lines 56/81/84. The faces that load the material element are not aggregated, so `ℬ_hold^live` and `ℬ_E^steady` describe the load from a single face. | The material element's mechanical load and face work, as they enter `ℬ_hold^live` and `ℬ_E^steady`, account for every native face that loads it. The face set stays OPEN (O6/`𝒥_map`): no face count, sheet, slab or face identification is chosen. No free face label remains in the assembled objects. |
| 2 | PY's premise-4 face normal is four independent OPEN functions, and its area factor is a separate OPEN function. The direction is named, not constructed. WL computes both from one face chart. (Claude) | **ACCEPT**, PY. Lookup: PY lines 160–164; WL lines 82–84. With an arbitrary OPEN 4-vector, PY admits tangential bulk traction, which premise 4 excludes (spec §3.3). | PY constructs the premise-4 traction: its direction is the unit normal of an OPEN native face geometry, and its area/measure factor comes from the same geometry. Its in-plane projection then arises from that face's tilt, and K7 acts on that tilt. |
| 3 | PY harness `difference()`: the matrix branch returns `corrupted − baseline` before the equality test, and does not canonicalize. Unchanged tags print a nonzero difference (14 false bites across 9 knives). (Claude) | **ACCEPT**, PY harness. Lookup: harness lines 60–78. The bite adjudication (pipeline §4) reads these triples. | A printed difference is zero exactly when the corrupted payload equals the baseline payload mathematically, for matrices as for every other payload type. |
| 4 | PY K6c replaces `xi_w(r)` by xreplace after the slope is already built as `Subs(Derivative(…))`, so the differentiated `ξ_w` in the momentum construction is untouched. (Claude) Grok: K6c "freezes the profile-tuple jacobian … and leaves already-built graph-slope atoms". | **ACCEPT**, PY harness. Both legs record that the slope atoms survive. The directive's site is "where it is differentiated" (item 8, K6). Lookup: PY lines 113/121/128/130 and harness lines 41–42. | K6c makes `ξ_w` constant at every site where it is differentiated in the material momentum storage and transport construction, as K6a and K6b do for `V_r` and `ρ_br`. |
| 5 | PY name `rho_br` (the live profile `rho_br(r)`) matches the upstream S11c-b row `rho_br`, a uniform KNOB, which is a different object. (Claude) | **ACCEPT**, PY. Lookup: S11c-b exports lines 9274–9281. Directive item 3 requires a rename. A name-keyed consumer would otherwise bind the uniform knob as the live profile, an M3 freeze. | No O2 name matches an upstream row that denotes a different object. |
| 6 | WL `"StressOccurrences" -> {fullStress}` and `"TransformationsUsedInMomentum" -> {}` are typed claims about the construction. (Claude) | **ACCEPT**, WL. Lookup: WL lines 150, 216–217. Line 217's qualification text also asserts "no such transfer used". `CLAUDE.md` E1/E2: a script prints computed objects and never states conclusions. | Every payload field that describes the construction is computed from it, or is absent. The recorded relative-`O(ε)` qualification is a supplied domain statement and stays. |
| 7 | PY carries `R_br` (O7) and `M_perp` (O1) only in the coupled-input register. WL carries them in every material action. (Claude) | **ACCEPT**, PY. Lookup: PY lines 96 and 123–125; WL lines 105–107; spec line 138 (`ℛ_br` "is input to the density in mass/momentum/energy accounting"). The two engines' objects will not join. | Each engine's material momentum, force and energy objects carry every §3.2 input the spec says enters them. |

**Note, not a repair item: tag grain.** WL emits 12 tags, about one per spec §9 row. PY emits 25 finer objects
(lookups, last section). Directive item 4 asked for parallel tag sets. Neither leg reports a physics effect. The
cross-engine comparator (sub-step 6) will declare an explicit, reviewed map from each engine's tags to the §9
rows. The engines are not re-cut for this.

**Repair routing (G4).** This is the first repair round, so the same builders repair, each in its own worktree.
The WL brief holds only findings 1 and 6 and no PY content, so WL stays blind. The PY brief holds findings 1–5
and 7. Each brief states the finding and what must be true, never the wording or form of the fix. The repaired
engines go back to a fresh Claude leg and a Grok leg under the same prompt, and the review continues until clear.

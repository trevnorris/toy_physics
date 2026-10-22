# Gravity prior-art section, review round 0: dispositions (orchestrator)

The mechanical lookups behind each disposition are in `gravity_prior_art_review_r0_lookups.md`, beside this file
(generator `_scratch/prior_art_gravity/gen_r0_lookups.sh`). "P-n" cites a section of that file. Every source text
quoted here is a file a leg fetched and saved.

## The round

- **Artifact.** The new section "ATTRIBUTION PASS — the gravity sector and S16 (the worldtube)" in
  `docs/medium_requirements_and_prior_art.md` (l.370–465), its cross-reference in the gravity row (l.246), and its
  sources file `docs/_measurements/gravity_prior_art_sources.md`. A Claude subagent of the orchestrator wrote them at
  the user's request (10-22). Both files matched the baseline after the round. The working tree was unchanged apart
  from the lookups file (P1).
- **Prompt and legs.** `_scratch/prior_art_gravity/review_prompt_r0.md` (P1). The orchestrator's side wrote the
  section, so the legs are Codex + Grok (G1). Both got the identical prompt, in parallel, with web access:
  - Codex (`gpt-6.1-sol`, effort xhigh, P1): `/var/projects/_review_scratch/prior-art-gravity-r0-codex/REVIEW.md`;
  - Grok (`grok-4.7`, effort high): `/var/projects/_review_scratch/prior-art-gravity-r0-grok/REVIEW.md`.
  - Both exited 0 (P1).
- **Verdicts.** Codex: four findings. Grok: five findings (P2). Two findings overlap, so there are seven.

## Findings

**R1. A flat, open-volume source for the centre relations exists. G3 and the NOT FOUND list say none was located.**
(Codex 1, Grok 2.) **CONFIRMED.**
- The section: G3 calls a flat, open-region source "not located", and NOT FOUND lists "the centre-of-content
  relations (`G3`): not located ⇒ unattributed" (P3).
- The source: Irschik–Holl (2004), p. 157, "we now introduce the center of mass of the non-material volume under
  consideration", with its definition (4.13), its rate (4.14)–(4.18) and its acceleration (4.20), with mass crossing
  the moving surface (P4). Irschik–Humer (2014) presents "Center-of-mass oriented formulations" (P4).
- Its conditions, by the legs: a mass density with momentum `∫ρv`, no interior source in that section, one control
  surface, and no split of the surface's motion into centre and radius parts.
- **What must be true.** G3 cites the located flat source, with the conditions that differ from S16's. The NOT FOUND
  list does not call these relations unlocated.

**R2. G5 and G6 are labelled KNOWN for an object that is not S16's. G6 omits the restriction of its result.**
(Grok 1; Codex 2 on G6.) **CONFIRMED.**
- The section: G5 "S16's architecture: what must be supplied about a body for its motion to follow from surfaces
  outside it" is **KNOWN** (GR). G6 is "the same, for a defect in a superfluid", **KNOWN**. The headline says the S16
  record should cite `G5`–`G7` (P3).
- The sources:
  - EIH: "the gravitational equations for empty space are in fact sufficient to determine the motion", with "surface
    integral conditions … which determine the motion" (P5).
  - Thouless–Ao–Niu: "the total nondissipative transverse force acting on a quantized vortex moving in a uniform
    background"; "the transverse force has the form ρs K × vV" (P5).
- The section's own NOT FOUND list gives Thouless–Ao–Niu as "a vortex, no drain" (P5).
- S16 supplies no field equation for the surroundings, and computes no force law and no passive factor (P6).
- So the question S16 asks is the problem of motion, which GR and the vortex literature asked. The results they
  publish answer it for other objects: motion from vacuum field equations, and the transverse force on a circulating
  defect.
- **What must be true.**
  - No row labels KNOWN a result whose object is not S16's object.
  - The section may record that S16's question was asked before, by whom, and what each answer rests on, as a
    statement about the question.
  - G6 names the transverse, nondissipative restriction of the force it cites.
  - The headline's citations for the S16 record follow the corrected labels.

**R3. G10 and the headline fuse two Riemann fragments, and date both with one fragment's footnote.** (Grok 3.)
**CONFIRMED.**
- The section: G10 cites "'Gravitation und Licht' and 'Neue mathematische Principien der Naturphilosophie' (footnote
  'Gefunden am 1. März 1853')". The headline gives one 1853 sketch with inflow into atoms ∝ mass and light as
  "another motion of the same substance" (P3).
- The 1876 edition: fragment 2, "Gravitation und Licht", begins at OCR line 35172. Its substance is conserved in a
  volume element: "bleibt daher constant", at line 35630. Fragment 3, "Neue mathematische Principien", begins at
  35687, and its footnote "Gefunden am 1. März 1853" is at 35706 (P8).
  - The drain into atoms is in fragment 3 (Grok 3's quotes).
  - The light motion superposed on the gravity motion is in fragment 2.
- Grok also quotes a footnote in fragment 3 that states a different inflow law (∝ the gravitational force), which
  the section does not mention.
- **What must be true.**
  - Each property is attributed to the fragment that states it, and the date only to the fragment it is printed on.
  - The two fragments' differences, the conserved substance and the drain, are stated where the section relies on
    both.
  - Both inflow statements printed in fragment 3 are given.
  - G10's label holds only for what one fragment, or the two fragments read together and said to be, establishes.

**R4. G4 quotes a Wikipedia sentence as the papers' text.** (Grok 4.) **CONFIRMED.**
- The quoted sentence "must therefore be supplemented by three additional constraints…" is in the Wikipedia wikitext
  (P9). The sources file records it as the Wikipedia quote, and records Mathisson's text as not opened (P9).
- **What must be true.** G4 marks that quotation secondary, by the section's own rule, or replaces it with text
  from a paper that was opened.

**R5. G7 identifies the Lagally and Wang force with S16's passive factor.** (Grok 5.) **CONFIRMED.**
- The section: "This is S16's **passive factor**" (P3).
- S16: it "computes no passive factor". A passive factor needs "a named field of the surroundings, a rule
  separating it from the body's own field, and its normalization". The force a drain exerts is S15's (P6).
- Writing the prior force as S16's object makes a prior result a premise for this program's object (M3).
- **What must be true.** G7 states the prior result under its conditions, and names S16's open passive factor as
  the requirement that the record will check against it. It does not identify the two.

**R6. G9's inverse square omits the velocity term in Wang's force between moving sinks.** (Codex 3.) **CONFIRMED.**
- Wang's Corollary 4, as Codex fetched it, adds `ρQ₂v₂` to the inverse-square term (P7).
- **What must be true.** G9 states the inverse square as holding for sinks at rest, or as the interaction term, and
  names the velocity term.

**R7. G8 omits a condition of Auton's result for the sphere.** (Codex 4.) **CONFIRMED.**
- Auton, Hunt and Prud'homme's abstract: "for the sphere it is also assumed that the changes in the ambient velocity
  field over the scale of the sphere are small compared with the velocity of the body relative to the flow" (P7).
- **What must be true.** G8 states that condition beside the one it gives.

## Not findings

- **Reading not completed.** Codex records that it could not open the full primary text of Reynolds, Lamb, Lagally,
  Papapetrou, Dixon or Taylor, and that those citations are not independently cleared by it. Grok opened Reynolds
  (1903) and found G1 matching, and found G8, G9 and G11–G13 matching the texts it opened. This limits what the legs
  cleared. It is not a defect in the section. It is carried: each citation shows what was opened for it.
- **Line 173.** Both legs found the light pass's NOT FOUND entry still accurate. Riemann's light is the resistance
  of a line element to change of length, not MacCullagh rotation on a sheet with throats.
- **No target.** Neither leg found a statement of what S16's engines should produce.

## Outcome

Seven findings, all confirmed. The section as reviewed is preserved with these dispositions. Its author repairs it to
what must be true, and the repair gets a full round 1 review by Codex and Grok.

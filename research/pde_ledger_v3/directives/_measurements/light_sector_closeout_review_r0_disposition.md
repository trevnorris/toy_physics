# Light-sector closeout, review round 0: dispositions (orchestrator)

The mechanical lookups behind each disposition are in `light_sector_closeout_review_r0_lookups.md`, beside
this file (generator `_scratch/light_closeout/gen/review_r0_lookups.sh`). "D-n" cites a section of that file.

## The round

- **Artifacts.** The packet the orchestrator wrote, uncommitted on top of `8d91f448`. It has seven files,
  hashed in `_scratch/light_closeout/packet_baseline_r0.sha256` (D1):
  - the new record `steps/S9b_flow_oracle_comparison.md`;
  - the edits to `SUBSTRATE_REQUIREMENTS.md`, `STATUS.md` and `DEFERRED_HEAVY_RUNS.md`;
  - three lookup files.

  Both legs left both repositories unchanged; the packet hashes still matched after the round.
- **Prompt and legs.** The prompt is `_scratch/light_closeout/review_prompt_r0.md`, sha256 `41388356…`
  (D1). Both legs received it identically.
  - **Codex:** `gpt-6.1-sol`, effort xhigh, sandbox workspace-write in its scratch directory. Report
    `/var/projects/_review_scratch/light-closeout-r0-codex/review_report.md`.
  - **Grok:** `grok-4.7`, effort high. Report `/var/projects/_review_scratch/light-closeout-r0-grok/grok_review.txt`.
- **Verdicts.** Codex: FINDINGS, 4. Grok: FINDINGS, 4 (D1). Several findings coincide, giving six distinct
  findings.

## Findings

**F1. The correspondence is stated with the mass identification alone.** (Grok 1, High; Codex 2, Medium.)
**CONFIRMED.**
- **What the packet says.** The record's "supports" paragraph, `STATUS.md`, the `R-S21-01` evidence line and
  the light-closing pass all say that under the `GM` identification "the kept flow-only profile is the
  river-model and Visser §7 profile".
- **What the check prints.**
  - The squared-speed zero test is `True` only with the common-parameter substitution. With a clock factor
    kept, it prints `GM == clockRate^2*gmH` (D2).
  - The signed comparison prints `False` for one orientation, even after the substitution (transcript lines
    724–729), and `True` for the other (line 741). In all, 40 signed tests print `False` (D2).
  - The map lists frame, free-fall versus lab time, and speed units as `UNMATCHED` (D2).
- **What the check supports.** Equality of the squared speeds, under the mass identification in common
  clock and speed units, with the signed match holding for one orientation per component.
- **Repair.** State that, everywhere the claim appears.

**F2. The exchange density law holds only where the velocity is nonzero.** (Grok 2, Medium; part of
Codex 2.) **CONFIRMED.**
- **What the packet says.** The `R-S21-01` evidence line says `jn` "vanishes only on `ρbr = C₁/r^(3/2)`".
- **What the check prints.** Of the exchange zero tests, 8 reduce to the density equation, all in cases
  1–2. The other 10 are `True`, in the zero-velocity cases, which impose no density law (D3). The record's
  own exchange paragraph already states this correctly; its "supports" paragraph does not.
- **Repair.** Restrict both statements to the nonzero-velocity components.

**F3. The falsifier in `R-S22-01`.** (Codex 1, High; Grok 3, High.) **CONFIRMED.** Neither the wording the
user approved nor the wording written is sufficient. The survey's conditions differ by measurement (D4):
- **E33** needs an additional ordinary thermal mode with the usual occupation, energy and dispersion, and
  a calibrated emission and absorption response.
- **E34** needs its power to enter the calibrated CMB radiance rather than cancel in the sky/reference
  comparison.
- **E38** needs extra gravitating relativistic radiation beyond the standard contents, and the survey says
  a weak population or earlier decoupling is not generically excluded.

Two further facts bear on any falsifier:
- In S11 the longitudinal speed is set by `B_comp/ρ_br` and the transverse speed by `μ_R/ρ_br` (D4). So the
  extra wave's speed generally differs from light's, and that changes both what it carries thermally and
  what it radiates.
- Both legs show that the population depends on more than the emitted share: on emission and uptake rates
  against each timescale, on escape into the bulk, and on initial population and expansion. The legs'
  scripts are cited as evidence of the argument only; nothing in the register quotes their numbers.

A "comparable share" is therefore neither necessary nor sufficient for a conflict, and "equilibrium with
matter" is not sufficient either.

**Repair.** Keep the object (the energy split), and add what applying any measurement needs: the spectrum
and rates of emission and uptake, the extra wave's own speed, and its escape into the bulk. Then state that
E33, E34 and E38 are each applied under their own conditions, and that neither the share alone nor
equilibrium alone decides a conflict. Report the change to the user.

**F4. The `c₀` sentence names the wrong entry.** (Grok 4, Medium.) **CONFIRMED.**
- **The problem.** `R-S1-02` is shear response as a function of the order parameter `χ_B` (D5). The printed
  relation `c0 == Sqrt[muPerp[r]/rhoBr[r]]` instead asks how `μ⊥` varies with the brane density.
- **Repair.** State what the printed relation requires within the flow-only restriction, and stop naming
  `R-S1-02`.

**F5. S22's test 2 is about the longitudinal wave escaping, not about incoming light.** (Codex 3, Medium.)
**CONFIRMED.**
- **The problem.** Plan line 1181 reads "does the longitudinal actually radiate into the bulk", and needs the
  interface law (D6). That escape rate enters the population of the extra wave (F3).
- **Repair.** Separate test 1 (incoming light converting at a defect) from test 2 (the longitudinal wave's
  escape), and make `R-S22-01`'s thermal use consume test 2's result.

**F6. The memory sentence about P2b overstates what was measured.** (Codex 4, Medium.) **CONFIRMED.**
- **What the record shows.** The round-3 runs were admitted at a 54 GiB cap. The guard's pool budget was
  66571993088 bytes, that is 62 GiB, and host memory available at admission was 60376612864 to 61462528000
  bytes (D7). So the runs failed at the cap tried. That does not show this box cannot finish them.
- **My earlier error.** I wrote "54 GiB = box maximum" in my report to the user and in session memory. That
  was wrong in the same way.
- **Repair.** "Full-family completion is unvalidated after OOM kills at 54 GiB; the memory needed to finish
  is unknown." Correct the user report and the memory.

## Notes that change no claim

- **Codex: the §7 conformal factor.** Visser says the §7 acoustic metric is "conformal to the
  Painlevé–Gullstrand form of the Schwarzschild geometry but not identical to it". The phonons are
  minimally coupled scalars, so the conformal factor affects the frequency-dependent greybody factors (D8).
  The packet already disclaims any line element, full wave equation or GR result. Adding this to the second
  route is faithful to the paper and changes no claim, so the repair adds it to the scalar-versus-transverse
  item.
- **Codex: the S9b record's own header.** It still says its revision awaits review. The record is outside
  this packet, and the acceptance commit `ff3bac63` settles its standing.
- **Grok: the reassessment lookup.** It lists the new record when re-run on the working tree. That is
  expected, because it was generated before the record existed.

## Outcome

Six findings, all confirmed. Nothing is accepted. The reviewed packet is preserved as is. The repair (r1)
then goes to a fresh round with the same authorship pairing: orchestrator-written, so Codex and Grok.

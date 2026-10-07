# Light-leakage inventory v3: final-pass dispositions (orchestrator, 2026-10-07)

**User rule (2026-10-07):** one more fold by a fresh Claude author, then a final Codex + Grok pass. Accept if
nothing outstanding changes what the bracket would compute or claim.

**Artifact:** `light_leakage_scoping.md` v3. `sha256sum -c scoping_v3.sha256` → `light_leakage_scoping.md: OK`.
**Legs (identical prompt `light_leakage_scoping_review_prompt_r3.md`):**
- Grok `scoping_review_r3_grok.txt`: **CLEAR**, no findings;
- Codex `scoping_review_r3_codex.txt`: **NEEDS REVISION**, two findings.

Both reported before adjudication.

## Verification (mechanical lookups; literal output)

**Finding 1 (Codex): spontaneous Brillouin is mapped to the background hold.**

```
$ nl -ba _scratch/light_leakage/light_leakage_scoping.md | sed -n '41p;92p' | cut -c1-…
41  … they hold the background profiles without a time dependence of their own (SPa:235–239). … pointing to C3
    (the nonlinear program, for stimulated gain) and O7 (a time-dependent background, for spontaneous initiation). …
92  … | Acoustic wave ↔ a time-dependent background profile, which the records do not supply
    (`directives/S11c_a_SHARED_PHYSICS.md:235–239`; `V3_STEP_PLAN.md:566–570`); `Ω` ↔ the frequency exchanged with
    that background; `q` ↔ its in-plane wave number. No coupling-strength identity. |
$ nl -ba research/pde_ledger_v3/directives/S11c_a_SHARED_PHYSICS.md | sed -n '194,195p;234,239p'
   194	Wave perturbations `u`, `ζ_c`, `δW`, `θ`, and the first variations of the face and bulk fields carry the
   195	independent amplitude bookkeeper `ε`. …
   234	Let `χ(x,t)` be the inverse material map, … For every background profile
   235	`Q_bg ∈ {W_bg, μ_R,bg, ρ_4D,bg⁰, ρ_br,bg⁰}`, the two supplied physical branches are
   238	LAB_HELD:          Q_bg^L(x,t) ≡ Q_bg(x) ,
   239	MATERIAL_ADVECTED: Q_bg^M(x,t) ≡ Q_bg(χ(x,t)) .
```

**Verified.** The records carry dynamical thickness and displacement waves (`δW`, `u`) as perturbations on the
held profiles. So an acoustic-type wave is not something "the records do not supply", and a light–phonon
interaction on a stationary background is a coupling among wave perturbations. That is the nonlinear program's
territory, not only a matter of releasing the profile hold. Row 4 already says "Branches exist, quadratic
cross-block is zero" (line 27). Lines 41 and 92 contradict it by sending spontaneous initiation only to O7 and
by mapping the acoustic wave to a background profile. **This changes a mapping and a destination**, so it is
inside the physics filter. The v2→v3 fold introduced this material in answer to r1 finding F3.

**Finding 2 (Codex): the Friedland–Giannotti Eq. (8) has no validity domain.**

```
$ grep -n -i "infinitely small\|unperturbed function\|energies" /tmp/llrev_r1leg.fgRA/fg.txt   # pdftotext of arXiv:0709.2164v1
172:bound : an infinitely small perturbation to this setup that lifts the zero mode, 0 → E ′ ≡
181:a(E ′ )e−i 2E |s| , while for s . |s0 | it can be approximated by the unperturbed function,
214:phenomenologically viable for energies ≪ k.
```

Lines 176–186 of the same file derive Eq. (8) from the turning point, the unperturbed solution and
`|a(E′)|² ∼ (E′)^{(n+1)/2}`.

**Verified.** Eq. (8) is a small-lift tunnelling estimate (`E′ = m²/2k²` small, i.e. `m/k ≪ 1`). Two caveats:
- Codex marks `m/k ≪ 1` as its own inference.
- Line 214's "energies ≪ k" concerns suppressing the continuum modes, not Eq. (8).

The inventory's validity cell (line 86) gives no domain. **This is a validity domain on an oracle**, so it is
inside the physics filter. It is small and localized. The r1 F4 repair added this row.

## Outcome under the user's rule

**Not accepted as-is.** Two verified findings remain, and each changes what the bracket would compute or claim:
1. the mapping and destination of spontaneous/thermal Brillouin initiation;
2. the validity domain of one oracle.

Both are in material that earlier folds added in answer to leg findings. The count has gone 19 → 16 → 8 → 2,
and Grok found nothing. Following the stop rule, the orchestrator does **not** fold again on its own. The
decision goes to the user.

What must be true after any repair:
- **F1.** Spontaneous and stimulated Brillouin (and Raman) both point to the generic nonlinear program, whatever
  their initiation or gain class. The missing pieces — the interaction and the fluctuation statistics of the
  receiving branch — are each named as a gap. O7 keeps only independently time-varying material profiles. The
  source's thermal/vacuum distinction is quoted, not asserted for v3. The rows stay NOT ADDRESSED with no
  invented record status.
- **F2.** Eq. (8)'s validity cell states the small-lift domain (`m/k ≪ 1`), marked as an inference from the
  quoted derivation. It keeps "no v3 correspondence" and "no inherited bound".

## User decision (2026-10-07)

"The items you listed under light-leakage table look fine. Let's go with your recommendation." That is:
1. A fresh Claude author repairs F1 and F2 only, in place. v3 is preserved as `light_leakage_scoping_v3.md`
   (`85f1d09c…`).
2. Codex and Grok check only those edits (Claude-written → Codex + Grok), with an identical prompt.
3. The table is accepted if nothing outstanding changes what the bracket would compute or claim.

**Scope of the repair.** Change only what F1 and F2 require, plus any cell elsewhere in the table that would
otherwise contradict the change (for example a status count, the O7 entry, the survival matrix row, or the WB
oracle row). Leave every other cell byte-identical.

The directive still binds:
- no calculations, no new physics, no verdicts;
- record disagreements quoted, not resolved;
- unopened sources marked UNVERIFIED;
- cite only records and published sources.

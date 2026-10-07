# Light-leakage inventory, repair 4: author report (fresh Claude author, 2026-10-07)

**Artifact:** `light_leakage_scoping.md`, edited in place from v3 (`85f1d09c…`, frozen copy `light_leakage_scoping_v3.md` untouched). No lines were added or removed. Eight lines changed (26, 27, 41, 71, 86, 92, 94, 132), and every table row keeps its column count. A cell-level comparison against v3 (split on unescaped `|`) shows that only the cells listed below differ.

Both findings **hold** against the records and sources, so both were applied.

## F1: Raman/Brillouin initiation and the O7 mapping

**Evidence that the finding holds:**
- `SPa:194–195`: "Wave perturbations `u`, `ζ_c`, `δW`, `θ`, and the first variations of the face and bulk fields carry the independent amplitude bookkeeper `ε`." `SPa:197`: "the requested truncation is first order in wave amplitude".
- `SPa:235–239`: the hold applies only to the background profiles `Q_bg ∈ {W_bg, μ_R,bg, ρ_4D,bg⁰, ρ_br,bg⁰}` (LAB_HELD / MATERIAL_ADVECTED). Holding these profiles does not remove the dynamical thickness or displacement waves.
- WB (published PDF, p.1246, §2.A; rendered and checked visually): "(In the context of Brillouin scattering, "spontaneous emission" of photons refers to any incoherent emission, be it true spontaneous emission due to phononic vacuum fluctuations or emission stimulated by pre-existing thermal phonons. Above a low threshold temperature T ≈ ℏΩ/k_B ≃ 1 K, the latter contribution dominates by far.)" Spontaneous initiation is therefore a property of the receiving branch's fluctuations, not of a time-dependent background.

**Changes:**

| Line | Cell | Change |
|---|---|---|
| 27 | Row 4, fiber mechanism (col 2) | Added WB §2.A's thermal/vacuum parenthetical, quoted and labelled "not asserted for v3". |
| 27 | Row 4, mapping (col 3) | Replaced "Thermal initiation needs a fluctuating background … **O7** … C3 … pointers not exclusive". Both forms now point to **C3**, whatever their initiation or WB gain class, because the records carry the acoustic-type fields as wave perturbations on held profiles (SPa:194–198,235–239; marked as this inventory's reading). Two gaps are named: the interaction (no vertex, partner unidentified) and the receiving branch's fluctuation statistics. WB's thermal/vacuum account is not asserted for v3. "Not O7". |
| 26 | Row 3, mapping (col 3) | Spontaneous and stimulated Raman both point to **C3** (S11:157–159; SC:50–51; SPa:194–198, marked as a reading). Three gaps are named: the receiving branch, the interaction and its fluctuation statistics. AB's "linear" is glossed as the first term of AB Eq. (4.9), an expansion in the optical field (AB §4.2: "The first term of the equation refers to the linear polarisation and depicts the spontaneous linear Raman scattering"). Without this gloss, column 2's "spontaneous linear Raman" quote would read as contradicting the nonlinear-program pointer. |
| 41 | Rows 3–4 taxonomy note | Rewritten to match the two rows: the SPa:194–198 wave-perturbation and truncation basis, both rows → C3, the two gaps, the WB thermal/vacuum account quoted and not asserted for v3, neither row → O7. Its closing no-verdict sentence is kept. |
| 71 | §2.1 O7, object (col 2) | O7 now covers only a background material profile that varies in time on its own, drifting or fluctuating (S11c's `Q_bg`, SPa:235, released from the SPa:237–239 hold), as distinct from the `ε` wave perturbations (SPa:194–195). "The analog of spontaneous Raman/Brillouin" is removed. Raman/Brillouin are stated to be C3. |
| 92 | WB oracle row, route (col 1) | `O7` → `C3 (Part 1 row 4)`. |
| 92 | WB oracle row, correspondence (col 4) | "Acoustic wave ↔ a time-dependent background profile, which the records do not supply" → "Acoustic wave ↔ a brane wave perturbation" (SPa:194–195,235–239). Which perturbation is the partner is a gap. The missing three-wave vertex (S11:157–159; SC:50–51) and the receiving branch's fluctuation statistics are gaps. "No coupling-strength identity" is kept. |

## F2: Friedland–Giannotti Eq. (8) validity domain

**Evidence that the finding holds** (local `arXiv:0709.2164v1` PDF, title page confirmed; pages 5–6 rendered and checked visually):
- p.5: "an infinitely small perturbation to this setup that lifts the zero mode, 0 → E′ ≡ m²/2k²";
- p.5: "(n + 1)(n + 3)/8(|s0| + 1)² = E′";
- p.5: "for s ≲ |s0| it can be approximated by the unperturbed function";
- p.5: "|a(E′)|² ∼ (E′)^{(n+1)/2}";
- p.6: "This suppresses the wave functions by the barrier penetration factor ∼ (m/k)^{(n+1)/2}, making the model phenomenologically viable for energies ≪ k". This sentence follows "the continuum modes residing in the bulk are suppressed on the brane", so it does not state Eq. (8)'s domain.

**Change, line 86, FG row, assumptions/validity (col 3):** the cell now states the small-lift domain `E′ = m²/2k² ≪ 1`, i.e. `m/k ≪ 1`. It is labelled "this inventory's inference from the quoted derivation (the source does not state a domain for Eq. (8) explicitly)", with the four p.5 quotes as basis. It also says that p.6's "energies ≪ k" concerns the continuum modes, not Eq. (8). It ends with: "The domain is in the source's variables: `m` and `k` have no record object (next column), and no bound is inherited (§2.2 preamble)." The existing no-correspondence column (col 4) is byte-identical.

## Other cells changed for consistency

- **Line 26, row 3, status (col 4), and line 27, row 4, status (col 4):** removed "O7 has no record status". The rows no longer point to O7, so the clause would imply that they still do. The "No recorded … status" wording and the C3 status text are unchanged, so no record status is invented.
- **Line 26, row 3, destination (col 6), and line 27, row 4, destination (col 6):** `(O7; C3; …)` → `(C3; …)`.
- **Line 94, §2.2 summary:** "O7 has Wolff et al. Eq. (1)" would now be false. It becomes "C3's Brillouin case (Part 1 row 4) has Wolff et al. Eq. (1), phase matching only, with no coupling strength". O7 joins the no-oracle **gap** list, and C3's no-oracle entry is qualified "beyond the Brillouin phase-matching condition below".
- **Line 132, matrix O7, frequency/polarization (col 5):** "oracle WB Eq. (1) (§2.2)" → "no oracle candidate (§2.2)".
- **Unchanged counts:** the mapping count (line 39: NOT ADDRESSED 6, rows 3–4 included) and the route count (line 75: 15 items) do not change.

## Mechanical lookup behind "no record line checked here states the receiving branch's fluctuation statistics" (line 92 col 4)

```
$ cd research/pde_ledger_v3 && grep -rn -i -E "thermal|temperature|fluctuat|zero-point|occupation|Bose|Boltzmann|phonon" steps/*.md directives/S11c_*.md SUBSTRATE_REQUIREMENTS.md V3_STEP_PLAN.md | cut -c1-220
steps/S11bB_interface_assembly.md:46:2. ⛔⛔ **Passivity and reciprocity are properties of fluctuations about an EQUILIBRIUM. This model's
steps/S11bB_interface_assembly.md:49:   the fluctuation subsystem need not be passive, because it can draw on the drive.
steps/S11b_interface_coupling_law.md:60:2. Passivity/reciprocity describe fluctuations about an **equilibrium**; the reference state carries `v₀ ≠ 0`
steps/S11b_HANDOFF.md:52:2. **Passivity and reciprocity describe fluctuations about an EQUILIBRIUM; our reference state carries
steps/S10_two_transverse_photons.md:871:transverse fluctuation” as KNOWN and attributes that object to Cembranos,
steps/S9_light_requires_shear.md:33:> one Goldstone phonon. **A scalar's spectrum cannot contain a spin-1 photon.** … This is
V3_STEP_PLAN.md:1090:**gravity-change/phonon speed** and `c_γ` is the **light-cone speed**
V3_STEP_PLAN.md:1093:> *"`c_s` is the phonon/gravity-change speed; `c_gamma` is the light-cone speed; `v_b` is the condensate
V3_STEP_PLAN.md:1105:> asserted `==5`) — **an untouched separate calibration** (light ↔ gravity-phonon), NOT the Part-V one."*
SUBSTRATE_REQUIREMENTS.md:198:  S11b-B states verbatim: *"⛔⛔ **Passivity and reciprocity are properties of fluctuations about an
directives/S11c_c2_export_repair_directive.md:13:2. Each exported row stores its value **fully `sp.expand`'d** (verbose srepr).
directives/S11c_b_wl_admissibility_repair_directive.md:49:that the mixed thickness/temperature gradient invariant `∇θ·∇e_W` is one of the independent §3a invariants.
```

None of these hits states a thermal or vacuum fluctuation population for any brane branch.

## Noticed but left alone

- **O7 wording elsewhere:** line 52 ("Time-dependent or fluctuating backgrounds are O7"), line 71 col 3 ("no record addresses a fluctuating background") and line 132's label ("O7 / time-dependent or fluctuating background") are consistent with the narrowed O7, a background *profile* that drifts or fluctuates. They are byte-identical.
- **Other O7 references:** FD (line 149) and TL (line 150) still list O7 as a frequency-exchange item. That remains true, so they are unchanged.
- **Matrix C3 rows:** line 124 does not list WB Eq. (1) or the fluctuation-statistics gap. Neither omission contradicts the change; no other matrix row lists oracles.
- **B:46–49:** "Passivity and reciprocity are properties of fluctuations about an EQUILIBRIUM … the fluctuation subsystem need not be passive". This concerns passivity and reciprocity, not a fluctuation population. It is not cited; the bracket author may want it when the fluctuation-statistics gap is taken up.
- **AB §4.3 (Boltzmann):** "The spontaneous Raman scattering generates a weak anti-Stokes signal, since the population of the excited vibration levels is low according to the Boltzmann distribution". Not added, because F1 asks only for WB's thermal/vacuum distinction.
- **WB §1:** WB separates Brillouin from "standard" opto-acoustic interactions (acousto-optic modulators), "in that the optical beams can themselves excite or influence the acoustic field". This bears on the C3/O7 boundary. Not added, because using it would create a new mapping.
- **WB row assumptions (line 92 col 3):** this cell still quotes §2.C, "initiated by thermal fluctuations". That is an accurate quote, and the vacuum caveat is now in row 4 col 2. Unchanged.
- **Truncation source:** the "first order in wave amplitude" truncation comes from S11c-a's shared physics (SPa:197). SC:48 ("In S11c (non-uniform, linear)") would also support it but is not cited.
- **Checksum file:** `scoping_v3.sha256` now reports FAILED for the edited file. This is expected; the frozen copy still matches `85f1d09c…`.

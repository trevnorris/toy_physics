# Light-leakage inventory v2: round-3 dispositions and repair brief (orchestrator, 2026-10-07)

**User decision (2026-10-07):** one more fold by a fresh Claude author, then a final Codex + Grok pass. The pass
is accepted if nothing outstanding changes what the bracket would compute or claim.

**Artifact:** `_scratch/light_leakage/light_leakage_scoping.md` (v2). Preserved as `light_leakage_scoping_v2.md`,
sha256 `7adcfdc1…` (`scoping_v2.sha256`).
**Round-3 legs:** Codex `scoping_review_r2_codex.txt` (C1–C2) and Grok `scoping_review_r2_grok.txt` (G1–G6). Both
said **NEEDS REVISION**. The orchestrator confirmed the record lines for G1, G5 and G6 by lookup:
- DEC:112–117 (a slit edge "not the small-`∇W₀` regime"; the O(1)-fraction argument is a "hint");
- B:200–202;
- d:64 (the frozen drain);
- SR:155–157 and S9:67–68 (shear-confinement stress at the throat);
- b:35 and b:44.

| # | Finding | What must be true in v3 |
|---|---|---|
| G1 | O6 ("conversion where the brane bends") hangs on S9:68 / SR:157. That sentence is about shear-confinement stress at the throat, which is O4's subject, not a bent light guide. | No route rests on a record line that does not concern it. The absence of a background centre-line in `Q_bg` (SPa:235) is stated as what it is: a freeze, or a gap. The throat stress stays with O4. Row 5's NOT ADDRESSED stays. |
| G2 | Row 1 (absorption) attaches B:74–78, which is about the uniform `e_W ↔ u_T` coupling, not about absorbers. | Each citation concerns the object of its row. |
| G3 | Rows 2 and 8 give c1's ESTABLISHED kernel status for an analog that also rests on the b and c2 objects. Those are marked CONDITIONAL (b:35, c2:42) and UNRESOLVED (b:44, c2:51, c2:55). | Every row that uses those objects carries each record's own status, side by side, as row 11 does. |
| G4 | E5 cites SC:17 for a radiation claim it does not make. | The radiation claim cites B:93–94. SC:16–17 is cited only for step ownership. |
| G5 | N11a, the S11c rest-frame limit, is attached to O2, the live-drain route. d:64 says O2's frozen-drain conclusion does not transfer to live order conversion. | N11a binds only the S11c-derived routes it governs. O2 states that the frozen-drain conclusion does not establish the live receiver. |
| G6 | The BO yardstick (grating efficiency) is assigned to the weak-gradient routes E1/E3 without DEC:112–117, which says a slit edge is order-unity, outside the weak-gradient regime, and that the O(1)-fraction argument is a hint. | Each yardstick names the bounded quantity as its record states it (B:200–202). Every assignment of the grating yardstick to a weak-gradient route carries DEC:112–117's regime limit. |
| C1 | Peyton Eq. (2): λ is "the S0 wavelength in this case", i.e. the converted wave, not the incident SH0 scale. | The parameter correspondence matches the source's own definition. The incident 6 mm SH0 wavelength is kept as the reported fit's input. |
| C2 | The Brillouin classification drops WB §2.C's steady-state qualifier ("after the system has achieved steady state"; backward SBS is "always a stimulated process, even if it is initiated spontaneously"). The C3 and O7 pointers read as exclusive. | Initiation (spontaneous or seeded) and classification (stimulated, by steady-state gain) are described as the source does. C3 and O7 are stated as non-exclusive, and the record-status mapping stays unchanged. |

The directive still binds: no calculations, no new physics, no verdicts; record disagreements quoted, not
resolved; unopened sources marked UNVERIFIED; only records and published sources cited.

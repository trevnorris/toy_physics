# Light-leakage inventory: PAUSED after review round 3, waiting on the user (orchestrator, 2026-10-07)

**State.** The current file is v2, written by a fresh Claude author after the Codex v1 repair bred defects. It is
preserved as `light_leakage_scoping_v2.md` (sha256 `7adcfdc1…`, in `scoping_v2.sha256`). Both round-3 legs say
**NEEDS REVISION**:
- Codex (6.1-sol): `scoping_review_r2_codex.txt`, 2 findings.
- Grok: `scoping_review_r2_grok.txt`, 6 findings.

Findings per round: r0 19 → r1 16 → r2 8. The rule on review rounds (more than 2 means stop and ask) applies, so
the loop stops here.

**Remaining findings** (my lookups confirm the record lines for G1, G5 and G6):

| # | Finding | Changes what the bracket would claim? |
|---|---|---|
| G1 | O6 ("conversion where the brane bends") rests on S9:68 / SR:157. That sentence is about shear-confinement stress at the throat, which is O4's subject, not a bent light guide. Row 5's NOT ADDRESSED is fine. O6 should go, or be restated as "no background centre-line profile exists (SPa:235)". | Yes: a spurious route. |
| G2 | Row 1 (absorption) cites B:74–78, the uniform-decoupling caveat, which is not about absorbers. | Minor. |
| G3 | Rows 2 and 8 give c1's ESTABLISHED kernel as their status, without b:35 / c2:42 (CONDITIONAL) and b:44 / c2:51 / c2:55 (UNRESOLVED). | Yes: overstates status. |
| G4 | E5 cites SC:17 for a radiation claim it does not make. B:93–94 is the right line. | Minor. |
| G5 | N11a (the rest-frame limit) is attached to O2, the live-drain route that d:64 says lies outside that limit. | Yes: a domain error. |
| G6 | The BO yardstick (grating efficiency) is assigned to the weak-gradient E1/E3 without DEC:112–117. Those lines say a slit edge is order-unity, "not the small-∇W₀ regime", and that the O(1)-fraction argument is a hint. | Yes: applies the yardstick outside its regime. |
| C1 | Peyton Eq. (2): λ is the converted S0 wavelength, not the incident SH0 scale. | Oracle correspondence only. |
| C2 | The Brillouin classification drops WB's steady-state qualifier. C3 and O7 overlap, so they should be described as non-exclusive. | Minor. |

**Options for the user:**
1. **One more fold (recommended).** The same Claude author, or a fresh one, applies these 8 findings. Then a
   final Codex + Grok pass, accepted if nothing changes what the bracket would claim. The findings are converging
   and are all record-fidelity issues, not a premise problem.
2. **Accept v2 as is.** Record these 8 findings as known defects in a header note, and let the bracket step
   re-check each route against its records when it uses it.
3. **Narrow the inventory before any fold.** For example, drop Part 1's fiber-textbook column structure and keep
   Part 2's route list and yardsticks, which the bracket actually consumes.

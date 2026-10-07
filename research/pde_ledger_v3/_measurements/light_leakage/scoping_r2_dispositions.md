# Light-leakage inventory r1: review dispositions and repair brief (orchestrator, 2026-10-07)

**Artifact:** `_scratch/light_leakage/light_leakage_scoping.md`, v1 (Codex repair 1). It is preserved as
`light_leakage_scoping_v1.md` (sha256 in `scoping_v1.sha256`). v0 is `light_leakage_scoping_v0.md`.
**Legs (Codex-written → fresh Claude + Grok, identical prompt):** both **NEEDS REVISION**.
- Fresh Claude: `scoping_review_r1_claude.md`, findings F1–F8 plus minors.
- Grok: `scoping_review_r1_grok.txt`, findings G1–G8.

**Author change (G4).** The r1 repair bred defects in the material it changed: the E4 rewrite (G1, G2, F8), the
new correlation P (G3), the new E5 route (G4), the rewritten Gu–Fuller cell (G8), the new PS row (F6), and review
reports cited as sources (F7). A fresh Claude author therefore writes v2, and its legs are Codex + Grok.
This is review round 3. If it does not clear, the orchestrator stops and brings the inventory to the user.

I checked by lookup the record lines behind the Claude-only substantive findings:
- `directives/S11c_decisions.md:130–133` (N13) and `:179–186` (N11a);
- `directives/S11c_a_SHARED_PHYSICS.md:235–239`;
- c1:181–182; S9:68; `V3_STEP_PLAN.md:566–573`;
- PC:7 and d:7 (both mark the uniform result **CONDITIONAL**);
- B:53–55.

Each says what the leg quotes, with one qualification at F3 (below).

## Dispositions (all ACCEPT; qualifications noted)

| # | Finding | What must be true in v2 |
|---|---|---|
| F1 | S11c's standing rest-frame decision (decisions N11a) and N13's definition of what does not survive are missing. | Every S11c-derived item quotes N11a as its drain freeze and validity domain. The survival convention cites both halves: PC:3 for what survives, N13 for what does not. |
| F2 | Bends are mapped to "curved faces", an analog the records lack. The background set has no centre-line, and the kernel's shape is the face-EVEN displacement. | Bends and face-odd roughness carry the mapping the records support, with record pointers (S9:68, SR:157 name brane bending). The flat background centre-line is named as a freeze. Conversion at a bend appears as a gap or open item. Marcuse's width-versus-straightness parity statement is an oracle candidate, never a bound. |
| F3 | The background's time dependence is never named as a freeze, so frequency exchange off a time-dependent background (the linear, spontaneous analog of Raman and Brillouin) is missing. **Qualification:** the records define two background branches, LAB_HELD `Q_bg(x)` and MATERIAL_ADVECTED `Q_bg(χ(x,t))` (S11c_a spec :235–239). Quote both, and the plan's hold (`V3_STEP_PLAN.md:566–570`) with its stated validity condition. | The background hold is named as a freeze, with both branches quoted. Frequency exchange off a time-dependent or fluctuating background is an inventory item, with its record status and its yardsticks. Rows 3–4 state the spontaneous (linear) and stimulated (nonlinear) conditions as their sources do. |
| F4 | Missing oracle candidates: Friedland–Giannotti's bound-photon tunnelling rate, Eqs. (7)–(9), for the bulk-shear leak; the odd/active-media prior art that B:53–55 names, for E4. | Both are listed as oracle candidates, with their conditions. Neither is a bound. |
| F5 | V0 says "no sourced observable", but `V3_STEP_PLAN.md:572–573` names the expansion rate (as a hook that is "not claimed, not derived"). | V0 quotes that line, with its caveat. |
| F6 | The PS row says the lifetime yardstick was "explicitly deferred". No record says that. | The row states what the records say: owner assignments only. Any yardstick choice is marked as the inventory's own gap. |
| F7 | Review reports and repair notes appear in the inventory as citations or residue. | The inventory cites only records and sources. No review report, finding number or repair note appears. Where a source's printed unit conflicts with its own definition (F–G's supernova criterion), quote both from the source. |
| F8 + G1 | E4's uniform-zero status is shown one-sided. B:191 says "unconditional on a uniform background"; PC:7 and d:7 say CONDITIONAL. U:77–80 adds that `Im ω = 0` holds only for `μ_⊥ ≥ 0`. | The coupling-zero statement and the stability condition are separate. The records' differing statuses are quoted side by side and left unresolved. |
| G2 | E4/P adopt one reading of the passive-region rule (B:34–38, "a classification, not an acceptance gate") without quoting U:70–72 ("the second law forbids … so the model refuses"). | Both are quoted and left unresolved. |
| G3 | Correlation P drops the supplied `τ_A ≥ 0` and the exclusion of negative relaxation times (U:53–55). | P quotes the record's full statement. |
| G4 | E5 writes `Re Z = ρ_m/a_q` and calls A:79 exact, while A:36/:41 give `Z = ωρ_m/q`. `a_q` and `b_q` are undefined in A, B and U. | Both impedance expressions are quoted. The undefined symbols are marked undefined, and neither expression is chosen. |
| G5 | E2 says "established … flux"; the established object is the relative flux `J_s`. The far-field/energy audit is UNDECIDED (c1:105–107, :112–115). | E2 names `J_s` relative flux and marks the energy/far-field audit with the record's own status. |
| G6 | E3 omits that F/G are withdrawn interpretations, not established physics (c2:53, :136–138). | E3 carries the withdrawal, so the kernel values cannot be read as an established loss. |
| G7 + minor | C1's token is not verbatim (S11:279–280 `INTERFACE_EQUATIONS_SUPPLIED={}`). | Record tokens are verbatim. |
| G8 | Gu–Fuller Eq. (15): `m` belongs inside the radical. | The equation cell matches the source. |
| Minors (Claude) | Peyton's experiments used aluminium (titanium was the simulation); More et al. also publish `nσ < 2×10⁻⁴ h Mpc⁻¹`; Marcuse 1969 Eqs. (75)–(76) assume `C₀=1, C₁=0` at `z=0`; SR:201–203 vs :207–212; the Corning row is a product specification, not an observation, and is unopened; row 13's ABSENT mapping is better supported by B:74–76, :191. | Each is corrected or qualified as its source states. A product specification is labelled as one, not as an observation. |

**Unchanged by this round:** the directive (`light_leakage_scoping_directive.md`) binds as before: no calculation,
no new physics, no verdicts; record disagreements quoted, not resolved; unopened sources marked UNVERIFIED.

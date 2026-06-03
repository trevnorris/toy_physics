---
unit_id: 221
batch: VII.1
verifier_model: claude-opus-4-8[1m]
verify_date: 2026-06-02T16:40:00-06:00
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 4
findings_total: 4
material_change: false
---

# Verification — unit 221

Checkpoint stage; verified at the higher bar (no rubber-stamping). All four findings
were applied by Codex and each holds up under independent scrutiny of the post-fix
scripts. The probe-only numeric slice (Aabs=7, gamma=2, r=3, eta=1/4, ...) is Print-only
and not load-bearing in either engine.

## Per-finding outcomes

### F1 — tautological_check

**Classification:** resolved

**What changed:**
`.wl` L173-180 (was L160-164). The old checks `survivalLeft - 2 UdispLowLossMax` and
`residueRequirement * eta/(1+eta^2) * Sfam^2 - 2 DeltaVreq` were replaced with:
- `((-Udisp) /. r -> 1/eta) - UdispLowLossMax`
- `(survivalLeft /. Aabs -> residueRequirement gamma) - 2 DeltaVreq`

**Assessment:**
Genuinely de-tautologized. `Udisp` (L148) is built from `reExpected` — the line shape
derived earlier (Section IV) — not from `UdispLowLossMax`, so check 1 fails if the line
shape is wrong; it requires `reExpected` at `r=1/eta` = `(Aabs/gamma)*eta/(1+eta^2)`, a
non-trivial reduction. Check 2 routes the saturation through `survivalLeft` (the physical
left side `Aabs/gamma * eta/(1+eta^2) * Sfam^2`) under `Aabs -> residueRequirement*gamma`,
so the equality to `2 DeltaVreq` holds only because `residueRequirement` has its specific
value — it is NOT a round-trip of `residueRequirement`'s own definition. No definition is
asserted against itself. Both PASS in the log (L74-77 of output).

### F2 — insufficient_verification

**Classification:** resolved

**What changed:**
`.py` L143-147: two new asserts after the survival-window definitions:
`assert simplify((-U_disp).subs(r, 1/eta) - Udisp_lowloss_max) == 0` and
`assert simplify((Aabs/gamma * eta/(1+eta**2) * Sfam**2).subs(Aabs, residue_requirement*gamma) - 2*DeltaV_req) == 0`.

**Assessment:**
Deliverable #9 (linear survival window) is now genuinely asserted in BOTH engines,
mirroring the F1 route. The first asserts the general line-shape `-U_disp` (L106, built
from `re_expected` at L75) evaluated at `r=1/eta` reproduces `Udisp_lowloss_max` — fails if
the line shape is wrong. The second saturates `2*DeltaV_req` through the physical left side.
Neither is print-only and neither is tautological (they reference `U_disp`/`re_expected`,
not a re-statement of `Udisp_lowloss_max`). SymPy assert count is now 14 (was 12), exit 0.
Coverage gap closed in both engines.

### F3 — mathematica_transliteration

**Classification:** resolved (re-authored to genuine independence)

**What changed:**
The `.wl` was re-authored across the three load-bearing sections via native primitives and
a different decomposition than the `.py`:
- **M1 (Section II, L96-102):** `NfunDerived = Together[D[QPi/DeltaPi, portPi]]` — native
  `D[]` of the raw `QPi/DeltaPi`, NOT the pre-written square. A separate "N is a perfect
  square" check (L98-101) then confirms `NfunDerived - (Afun GW + R GU)^2/DeltaPi^2 == 0`,
  and `NfunDerived` (not a hand-written form) feeds the `dD_Pi/dPi + N` check (L102).
  The `.py` (L50-52) hand-supplies the square as the target; the `.wl` derives it.
- **M2 (Section I, L78-81):** new `Residue[chiPassive, {delta, I gammaStar}] - Astar`
  check — a native `Residue` extraction the `.py` never performs.
- **M3 (Section IV, L121-131):** starts from the uncollapsed `chiBWgeneric = Aabs/(delta - I gamma)`,
  `ComplexExpand`s Re/Im in terms of `delta`, verifies the generic forms
  `Aabs delta/(delta^2+gamma^2)` and `Aabs gamma/(delta^2+gamma^2)` (L130-131), THEN
  substitutes `delta -> r gamma`. The `.py` (L72-73) starts from the pre-collapsed
  `Aabs/(r*gamma - I*gamma)` — a different decomposition.

**Assessment — independence:** INDEPENDENT, not back to transliteration. The previously
line-by-line-ported Sections I/II/IV now reach the load-bearing results through different
primitives (`D[]`, `Residue`, `ComplexExpand` of the un-collapsed pole) and different
intermediate quantities (`NfunDerived`, generic-detuning Re/Im) than the `.py`. The
equality cross-checks to the hand forms are retained as the new independent checks (the
native derivation must MATCH the hand form to PASS), which is the correct pattern. Sections
III/V/VI retain equality checks — acceptable per the directive (not the headline
independence concern). A transcription error in the `.py` target would now be CAUGHT, since
the `.wl` derives the form rather than copying it.

**Recorded deviation (sign) — CORRECT, not a defect:** Codex used
`NfunDerived = Together[D[QPi/DeltaPi, portPi]]` (no leading minus) rather than the
directive's `Together[-D[...]]`. This is reconciled correctly: since `D_Pi = K_B - Q_Pi/Delta_Pi`,
`∂_Π D_Π = -∂_Π(Q_Π/Δ_Π) = -NfunDerived`. The perfect-square check confirms
`NfunDerived = ∂_Π(Q_Π/Δ_Π) = (A G_W + R G_U)²/Δ²` (i.e. `NfunDerived = +N`), and the
L102 check `D[DPi, portPi] + NfunDerived` evaluates to `(-∂_Π(Q_Π/Δ_Π)) + ∂_Π(Q_Π/Δ_Π) = 0`,
i.e. it verifies the target identity `∂_Π D_Π = -N` with `N = (A G_W + R G_U)²/Δ²`. The
sign bookkeeping is internally consistent and the check genuinely compares the native
derivative to the target `N`. The directive's leading-minus form would have produced the
same final identity; Codex's form is mathematically equivalent and cleaner given the
`K_B - Q_Π/Δ_Π` definition. Deviation accepted as correct.

### F4 — paper_misalignment (notes_contradicts_script)

**Classification:** resolved (user-authorized notes-only edit)

**What changed:**
Per the directive `## RESOLVED` block (user-authorized 2026-06-02), Codex renumbered stale
"Stage 238" → "Stage 221" and "Stage 237" → "Stage 220" in the notes file only.

**Assessment:**
The captured diff patch (`stage_221_diff.patch`) touches ONLY the two script files
(`.wl`, `.py`) — no math/content change, no paper.tex, no appendix, no paper card. The
scripts themselves already used the canonical "Stage 220/221" numbering (e.g. `.py` L36
comment, `.wl` L87 banner), and they remain unchanged. The notes renumber is outside the
scripts-only verifier scope (I do not read notes/), but I confirm no stale numbering exists
in the scripts and the script + paper card are untouched. Math is identical either way;
this was a pure label fix. No script regression.

## Exec log assessment

**SymPy:** exit=0. Notable lines:
- `Verified exact Stage 220 derivative identity: dD_Pi/dPi = -N(omega)` (output L13-14)
- `|U_disp|_max in the low-loss window = Aabs*Sfam**2*eta/(2*gamma*(eta**2 + 1))` (L39)
- 14 asserts present (grep), all pass, exit 0. The two new survival asserts (L143-147)
  run silently (assert form) and the script completes — exit 0 confirms they held.

**Mathematica:** exit=0. 18 `PASS:` lines, 0 `FAIL` (grep-confirmed). Notable lines:
- `PASS: A_* is the residue at delta=I gamma_*` (output L19) — native `Residue` check (M2)
- `PASS: N is a perfect square` (L28) — native `D[]`-derived `N` (M1)
- `PASS: Re chi(delta) - generic expected` / `Im ...` (L44/L46) — uncollapsed line shape (M3)
- `PASS: low-loss |U_disp|_max recovered from line shape at r=1/eta` (L75) — de-taut survival
- `PASS: residue requirement saturates survival window via line shape` (L77) — de-taut survival

**Checkpoint PASS-count check:** 18 PASS in the `.wl` (12 original substantive + 4 new
genuine independence/de-tautology checks [A_* residue, N perfect square, Re generic, Im
generic] + 2 new de-tautologized survival checks). All 18 are substantive — none is a
round-trip of its own definition. SymPy has 14 asserts (12 + 2 survival). Claim count (9
deliverables) is fully exercised; deliverable #9 (the formerly-missing/tautological survival
window) now carries 2 genuine asserts per engine. No short or tautological PASS among the 18.
The probe-only numeric slice (Section VII, L185-206) is Print-only and not load-bearing.

**Output freshness:** confirmed fresh. Both committed `.txt` outputs are mtime 16:15:33
(2026-06-02), newer than the scripts (`.py` 14:27:23, `.wl` 14:29:26). Regenerated post-fix.

## Material-change assessment

`material_change`: **false**.

No derived constant, formula, or downstream-propagating result changed. All edits are
added-coverage (F2 SymPy asserts; F1/F3 new `.wl` checks) and independence re-authoring (F3)
that re-derive the SAME results by different routes, plus a notes-only label fix (F4). The
derived quantities (`A_*`, `gamma_*`, `gamma_wall`, the Re/Im line shape, `N = (A G_W + R G_U)²/Δ²`,
the survival-window bound) are unchanged. No downstream unit depending on stage 221 outputs
sees a different number or symbol. (The orchestrator will still mark units > 221 as
`upstream_stale` per policy, but there is no specific concern warranting a re-audit.)

## Side observations (non-blocking)

- The `.wl` top banner (L35) and output L8 still read "STAGE 204 — RESONANCE / LINEWIDTH
  TRADEOFF" rather than "STAGE 221". This is a cosmetic banner-label artifact in the script
  (a `Print` string only); it is NOT in any finding's scope, does not affect any check, and
  the per-section banners and all PASS labels correctly say "Stage 220"/"Stage 221". Flagging
  only — does not block verification. The auditor did not raise it and I do not introduce it
  as a new finding.

## Verdict justification

All four findings are `resolved`. F1's survival checks now route through the earlier-derived
line shape (`Udisp`/`survivalLeft`) and are non-tautological; F2 adds the matching genuine
asserts in SymPy so deliverable #9 is covered in both engines; F3's `.wl` is genuinely
re-authored to independence via native `D[]`, `Residue`, and `ComplexExpand` of the
un-collapsed Breit-Wigner through a different decomposition than the `.py` (not a port), and
the recorded sign deviation on M1 is mathematically correct (it verifies `∂_Π D_Π = -N` with
the perfect-square `N`); F4 is a user-authorized notes-only renumber with scripts/paper card
untouched. Both engines exit 0 with 18 substantive `.wl` PASS lines and 14 SymPy asserts,
outputs are fresh, and the diff shows no regressions. `material_change` is false. Verdict:
**verified**.

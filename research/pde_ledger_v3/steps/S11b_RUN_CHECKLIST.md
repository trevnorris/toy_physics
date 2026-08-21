# S11b run checklist — the full motion, mapped to the proven order

⭐ **Purpose:** make "we did not skip a step" auditable. Maps `.claude/skills/step-run/SKILL.md`'s 13-step
proven order (measured on S9) + the G-series (`directives/S11b_unified_decisions.md`) to S11b's concrete
obligations, with status and citation. ⛔ This is process tracking, not physics — it needs no legs.

Naming: **A+B unified → `S11b`** (G1); **C → `S11c`** (forward decision, no artifacts to migrate; record in
the plan when C opens). Slug = `S11b_interface_coupling_law`.

## The proven order (step-run table + runbook)

| # | step | status | S11b obligation + citation |
|---|---|---|---|
| 1 | author build directives | ✅ `9bd2f184` | both per-engine directives wrap spec §§0–13 |
| 2 | leak-gate directives | ✅ | SymPy gate EMPTY (2026-08-20); ⛔ re-gate the WL directive at WL launch |
| 3 | review directives Codex+Grok BEFORE build | ✅ `9bd2f184` | 2 legs each, killed the F9 paraphrase + T7 mis-attribution |
| 4 | fold once; stop build if question wrong | ✅ | folded once, not re-legged (rule 7) |
| 5 | build engine 1 (own call; prove prompt/log) | ⏳ SymPy running `b5itl5yqf` | prompt 8294 B proven; log ALIVE 237 KB @60s |
| 6 | baseline before any REPAIR overwrites it | ⬜ | ⚠ reconciles with rule 9 below: legs sandbox /tmp copies, so NO commit needed for review; commit AFTER fold; if a repair round, commit the reported baseline first |
| 7 | launch legs ON SIGHT (fresh agent + Grok) | ⬜ prompt staged | Codex-written ⇒ fresh Claude agent + Grok; SymPy is Python ⇒ parallel; FORM ablation MANDATORY (rule 14) |
| 8 | run the engine yourself; READ output | ⬜ | ⚠ capture stdout → committed `scripts/out/S11b_interface_coupling_law_sympy_audit.out` (S11b is "pre-chain" precisely for lacking a committed `.out`, STATUS:10) |
| 9 | build engine 2 (blind WL) barred from engine 1 | ⬜ WL directive done | ⚠ DEVIATION: we build `.py` first, not `.wl` first (see below) |
| 10 | cross-engine comparison + control/dim/parity | ⬜ | ⚠ `reduction/` harness DELETED — replacement mapping below |
| 11 | fold; repair; re-verify each fix myself (rule 13) | ⬜ | — |
| 12 | orchestrator writes the step record | ⬜ | 2 legs Codex+Grok (physics-bearing prose) |
| 13 | Codex writes the TeX card, its own 2 legs | ⬜ | re-point `paper/steps/S11b_interface_coupling_law.tex`; G12(c) owed fixes below |

## Export-chain integrity sub-checklist (the "import/export" the ask names)

- ⬜ **Carry-forward fidelity.** Every key in `S11_exports.LEDGER` (1663 rows) appears in `S11b_exports.LEDGER`
  with an identical value, EXCEPT any `F9c` collision (which prefixes `s11b_` and leaves the imported row
  untouched). ⛔ No S11 row silently dropped or mutated. Real-consumer check: `import S11b_exports` runs the
  module (`[[feedback-regression-bar-must-run-real-consumers]]`), ⛔ never read a cached `.out`.
- ⬜ **Bindings.** `c_s0`/`mu_R`/`rho_br` are the IMPORTED objects, ⛔ not re-declared; ⛔ no second inertia
  knob for `rho_br0`; `rho_m`/`v_dr` originate here with no upstream row (G7, directive §"Accumulated export").
- ⬜ **`v_0` false-merge guard.** The drain is `v_dr`/`v_bulk_normal_0`; ⛔ the key `v_0` never appears in the
  S11b rows (F9 object-compares `Symbol('v_0')` as EQUAL and silently merges, G3).
- ⬜ **Digests** pin exactly `{own source, S11_exports.py, S11b_SHARED_PHYSICS.md}`.
- ⬜ **D3** in-run round-trip present (G8b); **`_RELATIONALS`** reviver present (G8c); export FROZEN
  (`MappingProxyType` outer + per-record, `del _LEDGER`).
- ⬜ **F6 publish gate**: `S11b_exports.py` published only if every §9 primary task (all except B8 controls)
  completed.

## ⚠ Deviations & gaps from the literal proven order — disposition

1. **`.py` before `.wl` (step-run runbook 4 says blind `.wl` FIRST).** ⭐ DEFENSIBLE and adopted: the ".wl
   first" ordering is residue of the CUT quarantine era ("before any `.py` exists"). Under the current
   architecture the ONLY blindness control is that the WL engine imports NOTHING and re-derives (G7, G9,
   rule 12) — which is independent of whether `.py` exists. SAFEGUARD: launch the WL Codex build handed ONLY
   its self-contained blind directive; ⛔ never hand it the `.py`, `S11b_exports.py`, or the SymPy `.out`.
   Each engine is independently constructed and blind to the other's output — the "barred from engine 1"
   condition (step-run table row 9) holds regardless of order.
2. **`reduction/engine_output_checks.py` harness is DELETED (`fb29bba2`).** step-run 7c assumes a 4-layer
   harness (cross-engine agreement · dimensional homogeneity · control response · tag parity). Replacement:
   - cross-engine agreement + tag parity → the **T7 comparator** (G8a): join by name, residual paired
     payloads, ⛔ reject a native boolean as a residual operand, three-valued undecided, repoint ablation.
   - control response → the **FORM ablation** in the two script legs (rule 14) — each tag must move under a
     form change.
   - dimensional homogeneity → ⚠ THE GAP. A WRONG DIMENSION is the class caught ONLY by the second engine
     (step-run 7b, measured on S9), NOT by review or ablation. ⇒ the comparator MUST compare the per-row
     `dimension_key`/dimension tags across engines, and the WL engine must emit them. Verify dimension tags
     are in the comparator's join scope.
   - `derived_or_declared.py` triage is gone ⇒ lean on the FORM ablation as the defect-catcher.
3. **G13 acceptance diff is an ORCHESTRATOR step, not an automated gate.** After the build I compute the
   WITHHELD criteria myself and diff against the engine's emitted objects — the `mu_s=0` reduction relating
   `Lambda_p0` to `Lambda_A0`, and the A-slice reproduction of `Z_perm` (G13). A mismatch is a FINDING, ⛔
   not a build failure; the diff happens on our side (rule 5). ⛔ Do not skip — this is the physics
   acceptance.
4. **New-relation algebra (step-run 8 sub-bullet).** Each new relation S11b emits (admissibility region, any
   forced `Lambda_X`–`Lambda_V` reciprocity) needs its own algebra checked — substitute the derived root
   into the residual and confirm it vanishes; ablate to prove it can fail. Covered by legs + the G13 diff.
5. **Comparator freeze order (T7).** Author the comparator to the frozen T7 contract and COMMIT it BEFORE
   running the cross-engine comparison; its join/residual logic is the contract, ⛔ never tuned to the values
   ("frozen before it sees either output", `S11_C17_C18_spec_repair_decisions_v2.md:53-60`).
6. **Owed card items & standing scope limits — ⛔ do not silently drop.** G12(c): `Lambda_A0` used undefined
   + dropped `Lambda_p0=0` qualifier, fixed at card re-point. G12(d): B's uncarried background-flow
   convective correction `O(v0 |q_n|/omega)` recorded as a standing scope limit. The frozen-wall-width
   (`rho_br0 = rho_br`) is a recorded freeze, ⛔ not a fix (step-run static-or-instantaneous check) — name it
   in the step record.

## Commit / gate discipline (governing)

- ⛔ **No commit before BOTH legs report** (rule 9, `[[feedback-no-commit-before-legs-report]]`). My
  verification is not a leg. Legs review the working-tree artifact via /tmp copies (sandbox), so review needs
  no pre-commit; the working tree is the fixed target while legs run (I touch nothing).
- ⛔ **Fold ONCE, then go** (rule 7); ⛔ **verify each finding myself** (rule 13).
- Whatever writes does not review: engines/card (Codex) → fresh agent + Grok; directives/step record
  (orchestrator) → Codex + Grok.

## ✅ SymPy engine CLOSED (`864d6f41`, 2026-08-21)
Steps 1–9 of the proven order for the SymPy engine are done: directives (1–4) `9bd2f184`; build (5)
round-1; baseline commits (6) `7dd89076`/`6d57b27e`; script legs (7) each round; `.out` captured +
committed (8) `scripts/out/S11b_interface_coupling_law_sympy_audit.out`. Export-chain integrity
sub-checklist all ✅ (carry-forward fidelity: 1663 S11 rows, none dropped/mutated; bindings to imported
objects; drain `v_bulk_normal_0` no v_0 merge; digests(3); D3; _RELATIONALS; freeze; F6). G13 acceptance
diff ✅ (both legs independently reproduced `Λ_p⁰=−Λ_A⁰/ρ_m`). Two emission-fidelity repair rounds
(directives `0863645b`, `75560832`), HARD STOP, no round-3.
⬜ STILL OPEN for the STEP (not the SymPy engine): step 9 blind WL engine; step 10 T7 comparator (⚠ compare
dimension tags); step 11 fold; step 12 step record; step 13 card. The `.wl`-first deviation stands (blindness
by absence) — the WL Codex gets ONLY its self-contained directive.

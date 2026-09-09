# S11c-c2 N6 cross-engine reconcile — DISPOSITION = PATH B (adopted 2026-09-09; execute post-compact)

**Decision (user-adopted 2026-09-09):** STOP the collapse-instrument series at v3 NOT-SOUND. Surface the
surfaced cross-engine operand residuals as **representational-difference UNADJUDICATED**, and carry the
cross-engine operand corroboration as an explicit **open debt**. Reason on record = *an inadequately-specified
comparison + a scope decision*, ⛔ NOT computational intractability. This explicitly revises the earlier
A-scope adoption in `_measurements/S11c_c2_N6_reconcile_question.md:219`.

## Why B — the collapse bridge is the WRONG OBJECT (2 independent experts, both recommend B)
Strategy consult (advisory, full-context): astra `gpt-6-astra` high + Grok `grok-4.6` high, both **Path B**.
Reports `_legs/S11c_c2_N6_reconcile_strategy_consult_{astra,grok}.md`; prompt `_legs/…_strategy_consult_prompt.md`.
Two independent structural reasons the post-construction graded coefficient-table bridge cannot reconcile the
channel-(b) source operands:
1. **(astra) A density identity is not a source-operand rewrite.** The emitted source operands are
   `μ = EL(energy density)` (already Euler-Lagrange-differentiated). Map 4's 19-row table is correct at the
   DENSITY level, but it sends constant WL coefficients to POSITION-DEPENDENT expressions (`R_W = W_0/W_bg`), and
   `EL(T·L) ≠ T(EL·L)` for a spatially-varying `T` — applying the density-level table AFTER EL misses
   product-rule terms. Concrete E04/E14 counterexample (astra, verified with a small symbolic calc, γ_14 set to 0
   to isolate κ): with `W=W_bg, R=W_0/W, e=e_W, k=κ_θW`, `L_WL = a·θ'e' + b·eW'θ'`, `T(a)=kR, T(b)=−kR/W`;
   `EL_θ(TL) − T(EL_θ L) = kW_0/W²·W'e' − 2kW_0/W³·e(W')²`, and **the first term survives at σ_W^1**. This is
   BEFORE grading; no grading-stage repair fixes it. (Standard variational fact — EL does not commute with
   multiplication by a varying field; verifiable by reasoning, no CAS instrument needed.)
2. **(Grok) The ungraded `R_W` identity cannot be per-grade rewritten under the frozen contract.** `R_W·W_bg = W_0`
   is a product of two η-dependent series; both engines expand `WBg→W_0(1+η·w1)` BEFORE grade extraction (WL
   `profileRules`+`finish` `.wl:127-141`; SymPy `grades()` `xreplace(profiles)` `diagnostic:245-248`), so the
   graded leaves carry `W_0`, not live `W_bg`. Reproducing the ungraded identity coefficientwise requires
   cross-grade convolution `(TF)_1 = T_0F_1 + T_1F_0`, which the no-grade-mixing tripwire forbids. Every v4
   option (per-grade images / re-expand W_bg / apply-then-grade) reinserts a spurious `W_bg`, mixes grades, or
   leaves η^1 uncollapsed. (This is the v3 MUST I verified via G4 — `collapse_directive_gate.md`.)

⇒ The correct instrument would apply the thickness/basis transformation UPSTREAM of EL (a new emit path in BOTH
engines — cost ≥ Path A + 2 engine reviews), and even that would validate a REPLAY, ⛔ not retroactively
reconcile the already-emitted `.out` streams. A/C as scoped are not "one more fix"; A's real risk is a
**false-AGREE over-clear** (freeze `W_bg→W_0` at η^0, drop η^1 — same class as the L-CAS relapse). The recurring
map-4 difficulty (3 distinct findings across 3 rounds: v1 `7f1ed738`, v2 `adf178b9`, v3 `6cfa7148`) was the
architecture signalling the mis-decomposition.

## ⛔ How to write B HONESTLY (both experts; the disposition record + step record MUST state ALL of these)
1. **EARNED:** per-engine operator covariance (Reading B), BOTH engines, biting knives (SymPy δ≈2.6e-22; WL own
   primes/bound); cross-engine, every MATCHED `R_cov`/`R_cov_baseline`/`R_cov_control_delta`/`SOURCE_CONTROL_DELTA`
   and `CARRIER_BRIDGE_RESIDUAL` key vanished (160/0, 320/0) + matched STRUCTURAL support agrees (400/0). ⚠ These
   matched zeros are `(0)−(0)` — a dual-engine confirmation of the *vanishing statement*, ⛔ NOT operand AGREE.
   Reading B stands (`d21c8ff5`).
2. **NOT EARNED:** operand AGREE on `C_E`/`C_M` (40), `SOURCE_ACTUAL/PREDICTED` (76), Φ (18). Schema-unmatched
   `R_N6`/channels/RC-sources stay UNDECIDED as SCHEMA (different WL block/kernel/component vocab), ⛔ not
   heaviness (`reconcile_question.md` §1.iii).
3. **WHY THE COLLAPSE STOPPED:** the load-bearing representational identity (map 4: thickness `E=W_0·e_W/W_bg`
   vs `e_W`, `R_W`) is an ungraded DENSITY-level constructor fact that (a) does not commute with the EL
   differentiation that produced the source operands, and (b) cannot be reproduced coefficientwise under the
   retained-order no-grade-mixing contract. Continuing would reinsert `W_bg`, mix grades, or over-clear. Cite the
   Grok v3 MUST + my G4 verification (`collapse_directive_gate.md`) + the 2-expert consult.
4. **WHAT *WAS* CHECKED about constitutive physics:** the 19-row UNGRADED energy-basis table is independently
   re-derived and correct AS A CONSTRUCTOR IDENTITY (Opus v3 review). ⛔ That is NOT operand collapse — do NOT
   upgrade it to "the sources agree."
5. **⛔ THE B TRAP — do NOT let "representational-difference-unadjudicated" become "known to be just thickness."**
   The unresolved alternatives include an inconsistent constitutive coefficient-field CONVENTION *or* an actual
   IMPLEMENTATION error — two internally-covariant constructions can still produce different physical responses.
   Record that the leftover **SHAPE was NOT inspected** (a spelling-only pass would show whether the 76 source
   residuals ARE the predicted `R_W`/measure-expansion pattern or something else — not built). A later reader
   must ⛔ NOT upgrade UNDECIDED to "we know it's only thickness."
6. **⛔ Do NOT justify B as "c2 already has everything it needs."** S11c-d's downstream objects include
   gradient-driven mixing + leakage (`S11c_decisions.md:83`), so a retained gradient-term discrepancy is NOT
   dismissible just because covariance holds. Justify B as **an explicit decision to carry unresolved
   cross-engine corroboration as a DEBT while preserving the conditional per-engine result** (`RESOLVED.md:36`
   already says cross-engine agreement remains owed).

## CARRY OPEN (regardless of disposition — into the step record)
- The 3 per-engine premise caveats: Φ physical-correctness (derive from actual motion, not just faithfully
  implemented); V face-velocity transform (`V_E≡V_M` is builder agreement, not a derived transform);
  extracted-block omitted-block leakage.
- The cross-engine operand corroboration DEBT (this disposition) + the un-inspected leftover shape.
- Carried from earlier: the 2 S11c-b sign conventions; the 6 §3d questions; c1 ENERGY (UNDECIDED).

## Census (§3d) — do it as a FACT-LOOKUP in the step record, ⛔ no instrument
The 8 nonzero control/premise leaves: the frozen 6-row crosswalk (`kappa_a↔ADVECTION`, `kappa_j↔JUNK`,
`imported_theta_e_jets↔DOMAIN`, `coverage↔COVERAGE`, `uncovered↔UNCOVERED`, `max_present_jet_rank↔MAX_RANK`);
WL-only `THICKNESS`/`MATERIAL_NORMAL` stay one-sided; the `n6cov_J_mu`↔`junkMu` spelling. Surface, ⛔ don't
pre-adjudicate.

## Optional (⛔ default SKIP unless the user asks): a Φ-only spelling slice
Φ is excluded from map 4; a Φ-only pass (maps 1/3/5) would corroborate HALF of caveat 1 ("same declared map,
independently constructed") but NOT "is Φ physically correct." Both experts: don't spend a review cycle on it;
caveat 1 stays carried either way. astra's alternative if funded: a bounded upstream differential-compatibility
diagnostic (the E04/E14 sector, density + EL transport) to CONFIRM the bridge is ill-posed before any further
build — ⛔ not another grading-only v4.

## EXECUTE (post-compact, in order)
1. **Write the reconcile disposition record** — a short `steps/`-adjacent record (or fold into the step record)
   stating the 6 points above + carried-open. This is orchestrator prose ⇒ it gets 2 legs (O-written → Codex +
   Grok) review-until-clear before it's the governing disposition. ⛔ Do NOT over-claim; ⛔ do NOT call it "weak
   N6" — call it "covariance dual-engine confirmed; constitutive COMPONENTS not put in a common thickness
   coordinate at graded order; cross-engine operand agreement carried as a debt."
2. **Write the c2 step record** — `steps/…`; fix the misleading `I_{M→E}` "mapped-operand" terminology; preserve
   BOTH `R_N6`(=18/288, per-engine SymPy) and `R_cov`(=0, per-engine SymPy); carry the 3 premise caveats + the
   debt + the 2 S11c-b signs / 6 §3d / c1 ENERGY; census as fact-lookup; NO per-substep card (one S11c roll-up
   after S11c-e). Source-first fidelity; 2 legs (O-written → Codex + Grok) review-until-clear; both reports
   before commit.
3. ⛔ Do NOT rebuild/re-review the collapse directive; ⛔ do NOT re-litigate N6 covariance or c1.

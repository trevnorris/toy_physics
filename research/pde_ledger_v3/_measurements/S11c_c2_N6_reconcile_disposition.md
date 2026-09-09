# S11c-c2 N6 cross-engine reconcile — DISPOSITION (governing)

**Status:** orchestrator-written governing disposition of the SURFACED cross-engine residuals from the N6 T7
comparator run. Decision = **PATH B** (user-adopted 2026-09-09): STOP the collapse-instrument series at v3
NOT-SOUND; surface the surfaced cross-engine operand residuals as **representational-difference UNADJUDICATED**;
carry the cross-engine operand corroboration as an explicit **open debt**. Reason on record = *an
inadequately-specified comparison + a scope decision*, ⛔ NOT computational intractability (the run deferred 0,
peak ~330 MB). This supersedes the earlier A-scope adoption at
`_measurements/S11c_c2_N6_reconcile_question.md:219`.

⚠ This is the leg-reviewed governing disposition. The **decision + 2-expert rationale + write-up spec** live in
`_measurements/S11c_c2_N6_reconcile_disposition_PATH_B.md` (orchestrator planning doc). The **framing that
decides nothing** is `_measurements/S11c_c2_N6_reconcile_question.md`. This record states the disposition itself
and is what the c2 step record cites. It is orchestrator prose ⇒ **2 legs (Codex + Grok), review-until-clear,
both reports before commit** ([[feedback_specs_reviewed_until_clear]]); the c1 reconcile's correction-verify leg
scoped the orchestrator's verdict TWICE, so this disposition is NOT self-cleared.

**Grounding.** Every count below is from the committed mechanical tally
`_measurements/S11c_c2_N6_comparator_run_tally.txt` (run record `..._comparator_run_data.md`; comparator run
`a094b284`, joining the blind WL `.out` `ae73b884` against the 3 SymPy `.out` streams `7c0790ab`). ⚠ Read every
"0" as **"no nonzero found at the retained rectangle, under the adopted PIT-qualified disposition"**
(`S11c_c2_N6_covariance_sympy.py:8`, `_measurements/S11c_c2_N6_RESOLVED.md:23`), ⛔ NOT an exact symbolic zero.

---

## 1. EARNED — per-engine operator covariance, dual-engine confirmed (but a `(0)−(0)` confirmation)

- **Per-engine N6 = operator covariance (Reading B), in BOTH engines.** SymPy `R_cov` (source-naturality residual
  `ms − source_terms(μ_E.subs(Φ), V_E)`) gives no nonzero found at conditional δ≈2.6e-22 all 4 cases, knives bite
  (`d21c8ff5`, `_measurements/S11c_c2_N6_RESOLVED.md`); blind WL reproduces the carrier reconcile `C_E−C_M` and
  the source-naturality `R_cov` per-engine (`48a0b4e7`, own WL primes/degrees, import-free). Reading B stands
  (`d21c8ff5`, user-adopted).
- **Cross-engine, every MATCHED vanishing agreed.** SYMBOLIC channel: `N6COV_R_COV` /
  `R_COV_BASELINE` / `R_COV_CONTROL_DELTA` — 160 ZERO / 0 NONZERO each; `N6COV_SOURCE_CONTROL_DELTA` 160/0;
  `N6RC_CARRIER_BRIDGE_RESIDUAL` 320/0; `N6RC_ADVECTION_ABSENCE` 6/0; `N6RC_FROZEN_RELATIONS` 8/0. STRUCTURAL
  support: 400 ZERO / 0 NONZERO (every matched support key agrees; no support disagreement).
- ⚠⚠ **These matched zeros are `(0)−(0)`.** Each matched object (`R_cov`, the carrier bridge `C_E−C_M`) is
  *already 0 within each engine*, so its cross-engine difference is 0 **trivially** and says NOTHING about whether
  the underlying OPERANDS agree cross-engine (`reconcile_question.md` §2, grounded both engine sources). ⇒ this is
  a **dual-engine confirmation of the VANISHING statement (Reading B), ⛔ NOT operand AGREE.**

⇒ Honest headline: **operator covariance is confirmed by two independent engines; the constitutive COMPONENTS
were NOT put into a common thickness coordinate at graded order; cross-engine operand agreement is carried as a
DEBT.** ⛔ Do NOT phrase this as "weak N6."

## 2. NOT EARNED — operand agreement on the surfaced residuals

The substantive cross-engine signal is in the SURFACED operand residuals, which are UNADJUDICATED:
- **CARRIER operands** `N6RC_CARRIER_EULERIAN` / `CARRIER_MATERIAL` — 280 ZERO / **40 NONZERO** each (SymPy
  imported-slab carrier vs WL blind graph-geometry carrier). ⚠ The bridge `C_E=C_M` vanishing (§1) does ⛔ NOT
  dispose of these 40.
- **SOURCE + Φ** `N6COV_SOURCE_ACTUAL` / `SOURCE_BASELINE` / `SOURCE_PREDICTED` — 16 ZERO / **76 NONZERO** each
  (`SOURCE_BASELINE` is a nominal-control DUPLICATE of `SOURCE_ACTUAL` at shipped settings — ⛔ not a third
  discriminator); `N6COV_FROZEN_PHI` 42 ZERO / **18 NONZERO**. The constitutive source channel + the field map Φ.
- **`R_N6` itself, the channels, and the guards are ENTIRELY UNMATCHED** (`N6RC_R_N6`, `SPLIT_CHECK`, `SPLIT_SUM`,
  `*_OPERAND`, `*_CHANNEL`, `CROSS_CHANNEL`, `DIMENSIONS`, all guards — all UNDECIDED, no matched sibling, no
  structural support line). ⭐ This is a **SCHEMA non-join** — different WL block/kernel/component vocabularies,
  and a `(face,grade)` vs `{face,wave,grade}` key shape (`reconcile:105-111` vs `.wl:708-711`) — ⛔ NOT
  heaviness, ⛔ NOT a measured residual. ⇒ there is **NO direct cross-engine comparison of `R_N6` itself**; the N6
  cross-engine evidence rests on `R_cov` + the carrier bridge + the SOURCE/CARRIER support, ⛔ not on `R_N6`.

## 3. WHY THE COLLAPSE STOPPED — the load-bearing bridge is the WRONG OBJECT (structural, not cost)

The reconcile would collapse the surfaced operand residuals (§2) by applying a frozen justified-identity bridge
dictionary to both engines' operands and testing whether the residual vanishes coefficientwise per grade. Its
load-bearing entry (map 4) is the energy-basis coefficient table carrying the thickness-convention scale
`R_W = W_0/W_bg` (SymPy's rescaled local thickness `E = W_0·e_W/W_bg` vs WL's background thickness `W_bg`
directly). The collapse-instrument directive series reached **v3 NOT-SOUND** after 3 review rounds (v1
`7f1ed738`, v2 `adf178b9`, Codex-authored-v3 under rule 15 `6cfa7148`; v3 review legs SPLIT — fresh Opus SOUND /
Grok NOT-SOUND — and **my G4 confirmed GROK**, `_measurements/S11c_c2_N6_reconcile_collapse_directive_gate.md`).
A user-requested full-context strategy consult (astra `gpt-6-astra` + Grok, both **Path B**;
`_legs/S11c_c2_N6_reconcile_strategy_consult_{astra,grok}.md`) then gave **two independent structural reasons the
post-construction graded coefficient-table bridge cannot reconcile the channel-(b) source operands**:

1. **(astra) A density identity is not a source-operand rewrite.** The emitted source operands are
   `μ = EL(energy density)` — already Euler-Lagrange-differentiated. Map 4's table is correct at the DENSITY
   level but sends constant WL coefficients to POSITION-DEPENDENT expressions (`R_W = W_0/W_bg`), and
   `EL(T·L) ≠ T(EL·L)` for a spatially-varying `T` (product rule / IBP). Applying the density-level table AFTER
   EL misses derivative-of-`T` terms. E04/E14 counterexample (astra, small symbolic calc, γ_14→0 to isolate κ):
   `W=W_bg, R=W_0/W, e=e_W, k=κ_θW, L_WL = a·θ'e' + b·eW'θ'`, `T(a)=kR, T(b)=−kR/W` ⇒
   `EL_θ(TL) − T(EL_θ L) = kW_0/W²·W'e' − 2kW_0/W³·e(W')²`, and **the first term survives at σ_W^1**. This is
   BEFORE grading; no grading-stage repair reaches it. (Standard variational fact — EL does not commute with
   multiplication by a varying field; verifiable by reasoning, ⛔ no CAS instrument needed.)
2. **(Grok) The ungraded `R_W` identity cannot be per-grade rewritten under the frozen no-grade-mixing
   contract.** `R_W·W_bg = W_0` is a product of two η-dependent series; both engines expand
   `W_bg → W_0(1+η·w1)` BEFORE grade extraction (WL `profileRules`+`finish` `.wl:127-141`; SymPy `grades()`
   `xreplace(profiles)` `diagnostic:245-248`), so the graded leaves carry `W_0`, not live `W_bg`. Reproducing the
   ungraded identity coefficientwise requires cross-grade convolution `(TF)_1 = T_0F_1 + T_1F_0`, which the
   retained-order tripwire forbids. Every v4 option (per-grade images / re-expand `W_bg` / apply-then-grade)
   reinserts a spurious `W_bg`, mixes grades, or leaves η^1 uncollapsed. (This is the v3 MUST I verified via G4.)

⇒ The correct instrument applies the thickness/basis transformation **UPSTREAM of EL** — a new emit path in BOTH
engines (cost ≥ the doomed instrument + 2 engine reviews), and even that would validate a REPLAY, ⛔ NOT
retroactively reconcile the already-emitted `.out` streams. A/C as scoped are not "one more fix"; A's real risk
is a **false-AGREE over-clear** (freeze `W_bg→W_0` at η^0, drop η^1 — the same class as the L-CAS relapse). ⭐
The recurring map-4 difficulty (3 *distinct*, deepening findings across 3 rounds) was the architecture signalling
the mis-decomposition, ⛔ not a sequence of fixable nits. [[feedback_reconcile_bridge_must_precede_el]]
[[feedback_decompose_before_building_gates]]

## 4. WHAT *WAS* checked about constitutive physics (⛔ NOT operand collapse)

The 19-row UNGRADED energy-basis table is independently re-derived and correct **as a constructor identity**
(fresh-Opus v3 review). ⛔ That is NOT operand collapse — it does NOT establish that the surfaced source/carrier
operands agree cross-engine. Do NOT upgrade "the density table is a correct constructor fact" to "the sources
agree." The two are separated by exactly the density-vs-EL and graded-`R_W` obstruction of §3.

## 5. ⛔ THE B TRAP — "representational-difference-UNADJUDICATED" must NOT become "known to be just thickness"

The unresolved alternatives for the surfaced residuals include a consistent constitutive coefficient-field
CONVENTION difference **OR** an actual IMPLEMENTATION error — two internally-covariant constructions can still
produce different physical responses. ⭐ The leftover **SHAPE was NOT inspected**: a spelling-only pass would
show whether the 76 source residuals ARE the predicted `R_W`/measure-expansion pattern or something else — that
pass was NOT built. A later reader must ⛔ NOT upgrade UNDECIDED to "we know it's only thickness." The debt is a
genuine open item, ⛔ not a closed one dressed as open.

## 6. ⛔ Do NOT justify B as "c2 already has everything it needs"

S11c-d's downstream objects include gradient-driven mixing + leakage (`S11c_decisions.md:83`), so a retained
gradient-term discrepancy is NOT dismissible just because covariance holds. Justify B as **an explicit decision
to carry unresolved cross-engine corroboration as a DEBT while preserving the conditional per-engine result**
(`_measurements/S11c_c2_N6_RESOLVED.md:36` already records cross-engine agreement as owed).

## 7. CARRY OPEN — into the c2 step record

- **The cross-engine operand corroboration DEBT** (this disposition) + the **un-inspected leftover SHAPE** (§5).
- **The 3 per-engine N6 premise caveats** (Reading B does NOT close them —
  `_measurements/S11c_c2_N6_RESOLVED.md`): (1) is Φ itself physically correct (derived from actual material
  motion, not merely faithfully implemented — `R_cov` cannot exclude an error SHARED by the declared Φ and both
  of its own routes); (2) does face velocity V transform correctly (`V_E≡V_M` = builder agreement, prediction
  uses `V_E`, not `Φ(V_E)`); (3) extracted-block omitted-block leakage.
- **Carried from earlier** (surface, ⛔ do NOT pre-adjudicate): the 2 S11c-b sign conventions; the 6 §3d
  questions; c1 ENERGY (UNDECIDED).

## 8. Census (§3d) — a FACT-LOOKUP in the step record, ⛔ no instrument

The 8 nonzero control/premise leaves (`N6COV_ACTUAL_CONTROL_PARAMETERS` 4, `N6COV_PHI_DOMAIN_CENSUS` 4 [+44 BOOL
rejected]) are surfaced as **production-control equivalence** + **domain-coverage equivalence**, ⛔ not "do they
carry physics." The frozen 6-row crosswalk (`kappa_a↔ADVECTION`, `kappa_j↔JUNK`, `imported_theta_e_jets↔DOMAIN`,
`coverage↔COVERAGE`, `uncovered↔UNCOVERED`, `max_present_jet_rank↔MAX_RANK`); WL-only `THICKNESS` /
`MATERIAL_NORMAL` stay one-sided (PY `actual_amplitudes` retags only theta-advection + junk); the
`n6cov_J_mu`↔`junkMu` spelling. If a census leaf reflects differing μ-jet coverage it routes to §2 (source), ⛔
not to bookkeeping.

---

## Provenance
- Decision + 2-expert rationale + write-up spec: `_measurements/S11c_c2_N6_reconcile_disposition_PATH_B.md`.
- Framing (decides nothing): `_measurements/S11c_c2_N6_reconcile_question.md` + vet
  `_measurements/S11c_c2_N6_reconcile_question_vet.md`.
- Collapse-instrument closure: `_measurements/S11c_c2_N6_reconcile_collapse_directive_gate.md`; v3 review
  `_legs/S11c_c2_N6_reconcile_collapse_v3_review_{opus,grok}.md`; strategy consult
  `_legs/S11c_c2_N6_reconcile_strategy_consult_{astra,grok}.md`.
- Run data: `_measurements/S11c_c2_N6_comparator_run_{data,tally}.*` (reproducible —
  `python3 scripts/S11c_c2_N6_cross_engine_comparator.py`).
- Per-engine Reading B: `_measurements/S11c_c2_N6_RESOLVED.md` (`d21c8ff5`).

⛔ Do NOT rebuild/re-review the collapse directive; ⛔ do NOT author a v4 collapse instrument; ⛔ do NOT
re-litigate the N6 covariance resolution (Reading B, `R_cov=0`, `d21c8ff5`) or c1 (STANDS).

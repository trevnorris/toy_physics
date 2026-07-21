# Directive — magnetism as the moving throat, via boost-consistency (v1, DRAFT — pre-gauntlet)

**Status:** foundational derivation. **Target-blind — structurally.** Full gauntlet (Codex→Grok→Codex design-review BEFORE any build). Method = audit §F: **define the piece the model's OWN way, compute, ACCEPT whatever falls out — do NOT import Maxwell, do NOT import `j∝sV`, do NOT assume the magnetic sign or that the model is Lorentz-covariant.** Mirror the proven architecture of the gauntlet-clean electric directive (`directive_puncture_deflection_electric_sign.md`): every upstream question emits NEUTRAL FACTS; a single SEALED §4 owns all landing assignment; honest Tier-A/R1 tiering; atomic per-tooth acceptance; dual-engine.

**Why this exists.** Magnetism is the **twin** of the just-built electric sector: a committed mechanism (the *moving* `±w` throat — `pathA_39`'s `O(V)` operator on the throat-body) plus an **IMPORTED sign** (`j∝sV`, the audit's I1). The imported `j∝sV` is the *original anomaly* that triggered the whole EM reconsideration (the `force_visualizer` showed two parallel currents attracting on an imported current law). This directive tests, target-blind, whether the model's magnetism is **forced by its electricity under a boost** (emergent Lorentz structure — the wall the medium tradition had to clear) or an independent free knob / a characterized departure / a wrong-sign falsification.

**The core design — TWO independent derivations of the magnetic far-field; the COMPARISON is the result.**
- **Route A (boost):** apply a Lorentz boost to the *committed electric* interaction → the magnetic far-field the model MUST have if its EM is Lorentz-covariant (the Maxwell-consistent reference). Kinematic; robust to the electric sign's calibration (compare structure/velocity-dependence, not a calibrated number).
- **Route B (direct):** compute the magnetic far-field DIRECTLY from the moving-throat mechanism (the `O(V)` operator sourcing the transverse-vector on the brane), WITHOUT importing `j∝sV`.
- **Route B is computed BLIND to Route A** (no steering toward a match) — else the consistency test is rigged. §4 adjudicates whether they agree.

**Process (banked).** `codex exec -c model_reasoning_effort=xhigh … < /dev/null`; Mathematica-running builds add `--sandbox danger-full-access`; background/`setsid` for long runs (survives the harness kill); never kill mid-`WolframKernel` run (seat leak). Dual-engine (SymPy + independent Mathematica) wherever analytic; per-tooth ablation; dimensional check units-restored; NO fake/commentary scripts. **The BUILD session produces deliverables and EXITS — it does NOT run its own verification legs (self-verification = rig risk); independent 4-leg + Grok verify is arranged separately by the orchestrator.**

---

## 0. Read-first (build from the committed structure — do not paraphrase this directive)

- The gauntlet-clean electric build (reuse its machinery + its ARCHITECTURE): `software/em_charge_attribute/directive_puncture_deflection_electric_sign.md` (the neutral-facts→sealed-§4 discipline, honest tiering, atomic teeth) and `software/em_charge_attribute/puncture_deflection_electric_sign_{result.md,check.py}` (the committed electric far-field: field identity `ξ_w=ℓh`, the coupled `(u_L,h)` action, `m_gg=b z_g²/D`, the `1/R²` `h`-mediated interaction with its CALIBRATED sign). **The electric sign is calibrated (free `q/g`); Route A boosts this interaction structurally — the magnetism test must be robust to that calibration, not depend on a specific electric-sign value.**
- The committed magnetism / transverse sector: `software/stage1_solver/reports/pathA_39_*` (the `O(V)` moving-throat operator; note where it RESTS ON the imported `j∝sV` — that import is BARRED here) and `pathA_36_*` (the brane MacCullagh transverse shear = light; `c_γ²=μ_R/ρ_br`; magnetism propagates via THIS transverse-vector sector). The moving-throat = the `P_w`-ODD puncture orientation in motion; its `O(V)` operator is symmetry-allowed.
- The committed action + fields: `software/em_charge_attribute/g0_closure_card_v0.md` (fields, the transverse-vector DOF, the medium rest frame). The medium HAS a rest frame — so the "boost" that relates a static and a moving throat is the MEDIUM-frame boost; whether it reproduces the exact Lorentz boost is exactly what emergent-Lorentz-invariance means and is part of the test.
- The maintained map: `docs/model_definition_audit.md` §B (EM = characterized-departure; I1 `j∝sV` IMPORTED) + §G (the `~10⁴²` hierarchy; `F_e/F_g` = sharpest held-out both-sectors test, currently FIT), `docs/em_analog_next_phase_handoff.md`, `docs/two_throat_simulation_handoff_spec.md` (R1; the magnetic signs are `R_w`-ambient conditional; `pathA_39`'s kernel rests on imported `j∝sV`, barred from R2/R4 ancestry).

## 1. The object under test — state it with NEUTRAL labels (assume NOTHING about the magnetic sign or Lorentz covariance)

**Mechanism.** A charge is a static `±w` throat (its bend = the electric sector, built). A *moving* throat drags/shears the surrounding brane, sourcing a transverse-vector disturbance (the candidate magnetic field) via the same transverse sector that carries light (`pathA_36`). The magnetic FORCE appears between two throats when BOTH move (the velocity-dependent part of their interaction). **Do NOT name a magnetic sign, "attract/repel", "Maxwell", or "Lorentz-covariant" anywhere before §4.** Whether the model's magnetism is the boost of its electricity is UNKNOWN and is the object of the derivation.

**Q-CURRENT (define the moving-throat source the model's OWN way — do NOT import `j∝sV`).** Derive how a moving `±w` throat couples to the transverse-vector sector from the committed `O(V)` operator (`pathA_39`) + the throat's actual structure — i.e. derive the effective "current" (the source of the transverse-vector) as a computed object, NOT the asserted `j∝sV`. Preregister its field identity, units, and `R_w`/`P_w` parity BEFORE any force computation. Emit whether it reduces to `j∝sV`, or to something else, as a NEUTRAL FACT (not a landing).

## 2. Requirement — the questions (neutral throughout; two INDEPENDENT routes + the comparison)

**Q-BOOST (Route A — the Lorentz reference; kinematic, calibration-robust).** Apply a Lorentz boost (velocity `v`, to leading and next order in `v/c`) to the committed electric far-field interaction between two throats. Extract the velocity-dependent force (the candidate magnetic force) for parallel and antiparallel motion: its structure, falloff, sign relative to the electric term, and `v/c` order. **Emit this as the "Lorentz-consistent reference" neutral fact.** Robust-to-calibration: report the structure/sign *relative to the electric interaction* so the result does not depend on the specific calibrated electric-sign value. Do NOT compare to Maxwell here (that is §4).

**Q-DIRECT (Route B — the model's actual magnetism; computed BLIND to Q-BOOST).** From Q-CURRENT's derived moving-throat source + the transverse-vector sector, compute the DIRECT far-field magnetic interaction between two moving throats (parallel / antiparallel), in the MEDIUM frame: structure, falloff, sign, `v/c` order, coefficient. **This must NOT read, assume, or steer toward Q-BOOST's result** (the §6 anti-rig tooth enforces independence). Emit as neutral facts.

**Q-COMPARE (the neutral comparison — the crux; still no Maxwell).** Compare Route A and Route B as neutral facts: do their structure, falloff, sign (relative to electric), and leading `v/c` order AGREE? At what order do they first differ, if any? Emit `routes_agree(order)` / `routes_differ(order, characterized-delta)` / `route_B_R1(named gap)` — as facts. §4 owns the landing.

**Q-MAG (magnitude / free factors — mirror the electric Q-MAG).** Is the magnetic coefficient fixed by the throat geometry / the electric calibration with no residual free dimensionless factor, or does a free knob survive (a *second* calibration independent of the electric one)? Dimensionless-factor census + error bound. Emit `magnitude_forced_by_electric` / `magnitude_free_factor` / `R1(magnitude)` as a fact.

## 3. Method freedom + honest tiering (Codex designs the route)

Codex designs the computation. Tag every result with a neutral TIER fact (§4 maps tiers → landings):
- **Tier-A (`tier_A_conditional`, do now):** Route A (boost of the known electric far-field — analytic, calibration-robust); Route B's leading far-field IF the `O(V)`/transverse-vector operator is analytically tractable at leading `v/c`; the leading-order comparison. A Tier-A result carries `tier_A_conditional`; it does NOT itself reach a terminal verdict.
- **Tier-B (`tier_A_gaps_closed` once achieved — do NOT fake):** the full nonlinear moving-throat solve / higher `v/c` orders / the full transverse-vector dynamics that would DECISIVELY fix Route B and any characterized departure. If Route B's leading order is genuinely `R1`, say so and give Route A's prediction as the concrete R1 target Route B must match.
No clamping to Maxwell; derive the sector's own equations. `R1_REQUIRED` is a fully acceptable, first-class landing.

## 4. SEALED adjudication — Maxwell/magnetostatics target appears ONLY here — TOTAL precedence

Everything above is neutral. Only here map to the EM target (real magnetostatics: **parallel like-currents ATTRACT**, `1/R²`-class force, magnetism = the boost of electricity = `F_μν` one tensor). Build the full truth table over the neutral facts `{Q-CURRENT identity} × {routes_agree|routes_differ|route_B_R1} × {sign vs magnetostatics} × {Q-MAG fact} × {tier fact} × {internal_inconsistency?}`; first-match precedence:

1. **`internal_inconsistency(sector)`** (turning on the moving-throat coupling breaks a committed sector; instability; broken ledger row) → **`NO_GO(sector)`**.
2. **Unconditional-outcome predicate FIRST** (compute before R1-routing). The comparison is **unconditionally established** iff Routes A and B are both computed (not `route_B_R1`) at a common order AND `tier_A_gaps_closed`. If so, judge:
   - `routes_agree` at leading order AND the boosted-electric magnetic sign = magnetostatics (**parallel attract**) AND `magnitude_forced_by_electric` → **`MAGNETISM_LORENTZ_CONSISTENT`** (the win: magnetism forced by electricity, emergent Lorentz structure; one calibration → both sectors; triple-check).
   - `routes_agree` but `magnitude_free_factor` (a second independent knob) → **`MAGNETISM_CONSISTENT_FREE_MAGNITUDE(R1)`**.
   - `routes_differ` with the leading sign STILL correct (parallel attract) but a characterized `v/c` departure → **`MAGNETISM_DEPARTURE_CHARACTERIZED`** (a preferred-frame / non-Maxwell signature — a held-out prediction; quantify the delta).
   - leading magnetic sign WRONG (parallel like-currents REPEL), or Route A and Route B disagree in leading sign → **`MECHANISM_FALSIFIED(reason)`** (`wrong_magnetic_sign` / `routes_sign_conflict`).
3. **Else (predicate fails) — route to R1 by the blocking gap:** Route B leading order `route_B_R1` → **`R1_REQUIRED(direct_moving_throat)`** (state the concrete Route-A target it must match); Q-MAG unresolved → **`R1_REQUIRED(magnitude)`**; `tier_A_conditional` standing → **`R1_REQUIRED(consistency)`**.
4. Defensive backstop (should be unreachable) → `R1_REQUIRED(unclassified)` + diagnostic. No cell falls through; no code path branches on the EM target before this block.

## 5. Scope + prediction-hooks — HYPOTHESES TO CHARACTERIZE (allow YES/NO/UNDETERMINED)

- **Emergent Lorentz invariance (the deep one).** The Route-A-vs-Route-B delta (if any) at `O(v²/c²)` IS the preferred-frame / medium-rest-frame signature. Characterize it: is the model's effective EM Lorentz-invariant (delta → 0) or does the medium frame leak in (a held-out, in-principle-observable departure)? Report as a computed fact, not assumed.
- **Ledger consolidation (the endpoint).** Deliver the magnetic sector as a **term in the unified brane+bulk action** (the transverse-vector coupling to the moving-throat source), so the v2 ledger can absorb it — the four-sector action's completion.
- **Hierarchy capstone (flag only — do NOT compute here).** Once magnetism completes the four-sector action, `F_e/F_g` (§G) becomes the ratio of two coupling terms in ONE action → the sharpest held-out dimensionless test. Flag it as the downstream capstone; not in scope here.

## 6. Acceptance (able-to-fail, target-blind, ATOMIC per-tooth)

- **Faithfulness:** the `O(V)` operator + transverse-vector sector + electric far-field transcribed term-by-term from the committed sources; `j∝sV` NOT imported (Q-CURRENT derives the source; a tooth asserts the source is computed, not the asserted `j∝sV`).
- **Route independence (anti-rig, critical):** a computed guard proving Q-DIRECT (Route B) does not read/depend on Q-BOOST (Route A) — e.g. Route B's result is produced by a code path with no data dependence on Route A; a mutation that makes Route B copy Route A must be detected. The consistency claim is only meaningful if the two routes are causally independent.
- **Q-BOOST:** the boost is applied to the committed electric far-field, calibration-robust (result stated relative to the electric interaction); leading + next `v/c` order extracted symbolically; a wrong-boost mutation changes the answer.
- **Q-DIRECT:** the direct moving-throat far-field computed from Q-CURRENT's source; sign/falloff/`v/c` extracted symbolically; a wrong-source mutation changes the answer.
- **Q-COMPARE:** the agree/differ fact is a COMPUTED comparison of the two routes' symbolic results (not a literal); a mutation flipping either route's sign flips the fact.
- **§4 ownership + totality:** all landing assignment in §4; upstream emits neutral facts only; truth-table total, first-match, no fall-through; `NO_GO`=inconsistency-only vs `MECHANISM_FALSIFIED(reason)`=completed-nonmatch vs `R1_REQUIRED(named)`=unresolved; target-blindness guard (no magnetostatics/Maxwell/Lorentz token or desired-answer literal before §4); atomic tooth IDs, each with a production-path mutation failing only its own assert.
- **Dual-engine:** SymPy + independent Mathematica agree term-by-term where analytic; both exit 0.

## 7. Deliverables

- `magnetism_moving_throat_check.py` (SymPy) + `.wl` (independent Mathematica), both exit 0, emitting Q-CURRENT / Q-BOOST(Route A) / Q-DIRECT(Route B) / Q-COMPARE / Q-MAG / the neutral facts / §5 hooks + the §4 landing + all atomic teeth.
- `magnetism_moving_throat_result.md` — concise body (≤2pp) + appendix (exempt from the cap) with the two routes' full expressions, the comparison, the emergent-Lorentz delta, the ledger-ready action term, the §4 landing, and an explicit DECIDED-vs-R1 split. No softening.
- **After build (orchestrator arranges, NOT the build session):** fresh-agent 4-leg + Grok compute-verify.

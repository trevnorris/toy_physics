# Directive — the leftover-scalar (`u_L`) candidate for the electric SIGN

**Status:** foundational build (decides how far the leftover-scalar electric-sign candidate can be earned analytically). **Target-blind.** Full gauntlet. Method = audit §F: **define the piece the model's OWN way, compute, ACCEPT whatever falls out — do NOT import Maxwell, do NOT chase an exact match, do NOT assume the sign.**

**Scope note.** This is the **analytic** "solve given a supplied closure card" (sim spec R1) for the `u_L` mouth-ensemble closure, carried to the static two-body electric force — for what is analytically decidable, cleanly deferring the rest to R1. It does NOT re-run U2 (U2 proved the committed model does not *select* `𝔅` — 144/144 cells UNRESOLVED). What IS analytically decidable: (i) whether the committed action holds a nonrelaxable signed `u_L` datum **natively** (derive it; do not hardcode the landing); (ii) the **algebraic** static two-body sign/range for each precisely-specified held-boundary-datum, from the full coupled `(u_L,h)` functional. What is NOT (defer, do not fake): the "minimal" bolt-on (needs a supplied added action), and the assembled gravity/sleeve/momentum consistency (OPEN in G0 → R1).

⚠ **TARGET-BLIND DISCIPLINE (mandatory).** The build must compute Q1–Q3 from the committed ACTION only and emit NEUTRAL symbolic sign/range/admissibility results. It must NOT import any expected answer — not the sign, not a formula, not "repulsive/attractive," not "bolt-on expected." The prior writeup `docs/electric_sector_exclusions_and_leftover_scalar_hypothesis.md` records ONE un-verified clamp sketch (with a specific claimed formula); it is out-of-band SCOPE context only — do NOT open it during the derivation, and it must NOT enter the build's derivation ancestry or fixtures. The EM target (Coulomb like-repel, `1/R²`) may appear ONLY in the §5 verdict block.

**Process (banked).** `codex exec -c model_reasoning_effort=xhigh … < /dev/null`; Mathematica-running builds add `--sandbox danger-full-access`; background `run_in_background:true`, never wrap codex in shell `timeout`. Dual-engine REQUIRED (SymPy `.py` + independent Mathematica `.wl`), cross-reconciled term-by-term; anti-fake guards (couplings genuinely enter the brackets; controls share the code path; able to return ANY landing); dimensional-homogeneity check (units restored, able-to-fail) before any sign is trusted; NO `python3 -c` commentary scripts. **Codex→Grok→Codex design-review bookend BEFORE any build.**

---

## 0. Read-first (build from the REAL committed action — do not paraphrase from this directive)

- `docs/model_definition_audit.md` §B (electric row: falloff native-cond. on the scalar-`h` ansatz; SIGN ⚠ IMPORTED — like-repel only via a lock `G₀>0`; no-lock scalar limit *attracts*) + §F (I1/I2/I3; U1/U2/U3; the method).
- `docs/em_analog_next_phase_handoff.md` (maintained state) + `docs/two_throat_simulation_handoff_spec.md` (the closure-card; R1 "solve given a supplied closure card"; ⚠ **the maintained note that an ensemble/closure change need NOT flip the force sign** — do not assume it does).
- `software/em_charge_attribute/g0_closure_card_v0.md` — the committed EM-sector action + BCs: the coupled `(u_L,h)` block with `C_hu` mixing, `T_L=B_eff=ρ_B0²/χ_c>0`, `K_h`; the `h` mouth source `Σ_i η_i(k_mh − g_χ h s_i)`; `u_L` BC = natural continuity of `B_eff ∂_n u_L + C_hu ∂_n h` + `u_L→0` at IR; the declared-zero ledger (`r_B·u_L=0`, no direct wall–`u_L` coupling, `δE_g/δu_L=0`) — **every such entry is `[POSTULATE]`, not a native law**; the OPEN deferrals (§ dressed monopole / pair-force solve).
- `software/stage1_solver/reports/pathA_39_stage4_field_classification.md` — the `(u_T1,u_T2,u_L,h)` multiplet; `u_L` = the **brane longitudinal density mode** (NOT the bulk `c_s` gravity mode, NOT the committed `h` mediator); the conditional `q_L∝s`, `q_h`, `q_A^T` couplings; `q_M` tracked as a SEPARATE mass source, NOT a charge residue.
- `software/stage1_solver/reports/pathA_38*` — committed `THROAT_ELECTRIC_LOCALIZED_COULOMB` (the `h`-mediated `1/r²`, `M_h>0`); the clamp is a SECOND contribution → a double-count check is mandatory.
- `software/em_charge_attribute/directive_u2_boundary_adjudication.md` + report — the `𝔅` candidate universe (144/144 endpoint×ambient×stratum CELLS UNRESOLVED; U2 does not select); the clamp = the fixed-value/fixed-defect **mouth-ensemble** input keyed per `𝔅` branch (NOT one of the 144 cells).
- `software/stage1_solver/reports/pathA_24_T1_wall.md` — the `±w` **wall** result `T1_FAIL_NO_STABLE_WALL` (connected `S³`, zero barrier); relevant only as evidence about `±w`-orientation protection, **not** a direct proof about a `u_L` mouth datum.

## 1. The object under test (state it target-blind — no expected sign/formula)

**Candidate:** the electric charge's sign is carried by the leftover longitudinal scalar `u_L` when the throat imposes a signed boundary datum on `u_L` at the mouth (oriented by `±w`, `s=±1`). `u_L` is the brane longitudinal density mode, coupled to the committed `h` mediator via `C_hu`.

**⭐ Three DISTINCT closure choices — different held variable, reaction functional, and normalization, hence generally different two-body STRUCTURES. The SIGN may or may not differ between them — COMPUTE each with its correct (thermodynamically conjugate) reaction functional; do NOT assume they differ OR agree, and do NOT hardcode any sign:**
- **(V) fixed VALUE (Dirichlet):** `u_L|_mouth = s·u₀` held.
- **(M) fixed MONOPOLE / flux:** the far-field monopole strength (`∮` of the conormal) `= s·q` held.
- **(J) fixed SOURCE:** a `−J φ` Legendre term with `J=s·J₀` held.
For each, the held variable, the mouth profile, the conjugate reaction/work functional, and the normalization must be stated explicitly; the two-body energy is the reaction functional **thermodynamically conjugate to THAT held variable** (⚠ methodology, not the answer: e.g. a fixed-VALUE/Dirichlet ensemble uses the fixed-potential virtual-work functional = the Legendre dual of the bare field energy, NOT the bare field energy — using the wrong functional flips the sign). Do not assume the three functionals give equal energies.

## 2. Q1 — REALIZABILITY (native absence is decidable; the bolt-on is deferred)

**Q1a (DECIDABLE — derive it):** from the committed quadratic action + BCs + zero-ledger, does ANY committed ingredient hold a nonrelaxable signed `u_L` boundary datum against relaxation (a real energy barrier / protected constraint)? Enumerate the candidates the model actually has — the `u_L` self-stiffness, the `C_hu` mixing, the `h` mouth source + the natural `u_L` conormal-continuity BC, the IR `u_L→0` BC, the conditional stage-4 `q_L∝s` coupling (⚠ even if activated it is a SOURCE / J-type, NOT a nonrelaxable value/flux clamp — G0 has no direct `u_L` source), the wall/`χ_B` sector, the sleeve `±w` orientation, `S_hold`, geon couplings — and COMPUTE whether each can pin a fixed signed `u_L` datum. Emit `NATIVE_CLAMP_EXISTS(mechanism)` or `NO_NATIVE_CLAMP` (with the computed reason each candidate fails — e.g. `r_B·u_L=0`, `δE_g/δu_L=0`, `u_L→0` IR BC, connected-vacua from pathA_24 for the `±w` orientation).

**Q1b (DEFERRED — do NOT invent a "minimal" bolt-on):** there is no analytic minimality criterion without an exact added action + coefficient domain. Either (i) a SPECIFIC added-action term is supplied as input, and the build tests whether it holds the datum AND preserves Tier-A consistency (§4); or (ii) land `BOLT_ON_DEFERRED_TO_R1` and NAME what an added action would have to supply (a protected/disconnected `±w` sector + a direct `s↔u_L` datum coupling) — without asserting one exists or is minimal.

## 3. Q2 — the static two-body FORCE (full coupled functional; per held-variable; double-count guarded)

Use the **complete coupled `(u_L,h)` static functional** with the committed `C_hu` mixing and the existing `h` mouth source — NOT `u_L` in isolation (with `C_hu≠0` the relaxed static stiffness is `B_eff − C_hu²/K_h`, and the physical modes are mixed, not the bare `ω²=(B_eff/A_eff)k²` pole).

For EACH held-variable choice (V/M/J) compute, symbolically with `s₁,s₂` free:
- the static two-throat interaction energy + force (sign, far-field `1/R^n`, coefficient),
- via FOUR separated readouts under ONE fixed relative normalization: **`h`-only** (the committed pathA_38 contribution), **`u_L`-datum-only**, **`h`↔`u_L` interference**, **total** — so the candidate is neither double-counted against the committed `h`-Coulomb nor tuned against it.
Report all signs/ranges NEUTRALLY (the comparison to Coulomb is §5 only).

## 4. Q3 — CONSISTENCY (tiered; deferred checks block the strong verdict)

- **Tier-A (analytically testable — required):** scalar-sector Hessian positivity (`B_eff K_h − C_hu² > 0`, the committed stability bound); transverse (`u_T`) decoupling from the datum; G0 zero-ledger preservation (the closure introduces no forbidden coupling); `q_M`/mass-channel exclusion (the datum is charge-odd, NOT the mass leg — charge⊥mass guard); U2 non-selection (the closure is a supplied ensemble input, it does not claim to select `𝔅`).
- **Tier-B (DEFERRED to R1/sim — NOT analytically decidable; G0 defers these):** the assembled gravity/drain response under the datum, curved-sleeve stability, momentum/return closure, the dressed (non-perturbative) two-body response.
Any Tier-A failure ⇒ `NO_GO`. Tier-A all-pass with Tier-B deferred ⇒ at most `BOLT_ON_ALGEBRAICALLY_CONSISTENT / R1_R3_REQUIRED` — NEVER a full `BOLT_ON_CONSISTENT` from analysis alone.

## 5. Target-blind verdict block (the ONLY place the EM target may appear) — with explicit precedence

A **holder** = anything that pins a nonrelaxable signed `u_L` datum: a native ingredient (Q1a=`NATIVE_CLAMP_EXISTS`) OR a specifically-supplied bolt-on (Q1b) that the build VERIFIES actually holds a signed datum. Judge a holder by the force of the held-variable it ACTUALLY fixes (Q2) — never by whether some other hypothetical ensemble would have succeeded. Compare to the real-EM target (like charges repel; `1/R²`) ONLY here. Precedence, first match wins:
1. **Tier-A failure** (a committed invariant broken by the closure/bolt-on) → **`NO_GO(reason)`**.
2. **A HOLDER exists** (native, or a supplied bolt-on verified to hold a datum) — judge its ACTUAL held-variable force (Q2); Tier-A has passed (else item 1):
   - repel `1/R²`, **Tier-B passes** → **`NATIVE_DERIVED`** (native holder) / **`BOLT_ON_DERIVED`** (bolt-on holder) — the strong win; triple-check.
   - repel `1/R²`, **Tier-B deferred** → **`NATIVE_ALGEBRAIC_CONSISTENT / R1_R3_REQUIRED`** (native) / **`BOLT_ON_ALGEBRAIC_CONSISTENT / R1_R3_REQUIRED`** (bolt-on) — right algebraic sign; full consistency pending R1 (NEVER the strong landing).
   - wrong-sign / wrong-range / null → **`NO_GO(holder_held_variable_wrong_sign|wrong_range|null)`** — the holder exists but gives the wrong force (falsification, regardless of Tier-B or whether another hypothetical ensemble would work). (If both a native holder and a supplied bolt-on hold, the NATIVE holder is judged.)
3. **NO holder** — Q1a=`NO_NATIVE_CLAMP` AND (no bolt-on supplied, OR a supplied bolt-on the build finds does NOT hold a datum) → report the primary decidable result + per-held-variable (V/M/J) algebraic sign/range with the Coulomb comparison + the Q3 tier status. Landing: **`NO_NATIVE_CLAMP; ALGEBRAIC_SIGN={V:…,M:…,J:…}; BOLT_ON_DEFERRED_TO_R1`** (append **`BOLT_ON_INEFFECTIVE(reason)`** if a bolt-on was supplied and failed to hold the datum).
4. **`UNRESOLVED(named_datum)`** — first-class ONLY for a specifically-named analytic gap where no mutation/engine exhibits the derivation (anti-dodge).
This precedence is TOTAL over the decision domain (holder? native-vs-bolt-on; actual-holder Q2 sign/range; Tier-A pass/fail; Tier-B pass/defer): every point maps to exactly one landing. No code path may branch on the desired answer before this block. "Both signs reachable" is NOT required — do not mandate a flip.

## 6. Acceptance (exhaustive, per-tooth able-to-fail — one tooth per load-bearing assertion)

- **Faithfulness (transliteration audit):** the coupled `(u_L,h)` action, `C_hu`, `T_L=B_eff>0`, `K_h`, the three BCs (V/M/J), the `h` mouth source, and the zero-ledger `[POSTULATE]` entries match the committed sources term-by-term.
- **Q1a barrier tooth:** the realizability is a COMPUTED barrier/admissibility result; a mutation that INJECTS a holding coupling flips `NO_NATIVE_CLAMP→NATIVE_CLAMP_EXISTS` (proves it's not hardcoded).
- **Held-variable tooth:** each of V/M/J is computed from the reaction functional thermodynamically conjugate to that held variable; a mutation swapping in the WRONG (non-conjugate) functional changes the reported energy/force EXPRESSION (coefficient and/or sign) — proving the correct-functional choice is live. Do NOT require V vs M to differ in sign (the physical force is ensemble-independent; they may agree), and do NOT hardcode any expected V/M/J sign.
- **Double-count tooth:** `h`-only + `u_L`-only + interference = total (assert); a mutation dropping the interference term breaks it.
- **`C_hu` stability tooth:** violating `B_eff K_h − C_hu² > 0` fires `NO_GO(scalar_unstable)`.
- **Falloff + units teeth:** the `1/R^n` power and dimensional homogeneity each assert and each fails under a targeted mutation.
- **`q_M` guard tooth:** smuggling the mass leg into the charge datum fires the charge⊥mass guard.
- **Verdict-classification tooth:** a truth-table exercises EVERY §5 landing — `NO_GO(tier_A)` / `NO_GO(holder_held_variable_wrong_sign|wrong_range|null)` / `NATIVE_DERIVED` / `BOLT_ON_DERIVED` / `NATIVE_ALGEBRAIC_CONSISTENT_R1_R3_REQUIRED` / `BOLT_ON_ALGEBRAIC_CONSISTENT_R1_R3_REQUIRED` / `NO_NATIVE_CLAMP`+algebraic+`BOLT_ON_DEFERRED_TO_R1` / `BOLT_ON_INEFFECTIVE` / `UNRESOLVED` — spanning the full domain (holder? native-vs-bolt-on; actual-holder Q2 sign/range; Tier-A pass/fail; Tier-B pass/defer); the precedence is proven TOTAL (every domain point maps to exactly one landing) and enforced, not hardcoded.
- **Target-blindness tooth:** grep-independent guard that no §1–§4 quantity equals a hardcoded target sign/formula (the sign must be computed).
- Dual-engine SymPy + independent Mathematica agree term-by-term; both exit 0.

## 7. Deliverables

- `leftover_scalar_electric_sign_check.py` (SymPy) + `.wl` (independent Mathematica) — both exit 0, emitting Q1a / Q1b / Q2(V,M,J × 4 readouts) / Q3(tiered) + the §5 landing + every §6 tooth result.
- `leftover_scalar_electric_sign_result.md` — ≤2 pages: Q1a (native absence + computed reasons), Q1b (deferred / supplied-bolt-on result), Q2 (the V/M/J signs/ranges + the four readouts + the double-count check), Q3 (Tier-A pass/fail + Tier-B deferrals), the §5 landing, and an explicit honest statement of what is now decided vs what R1 must settle. No softening.

**After build:** fresh-agent 4-leg verification (arbiter re-run + fidelity + adversarial-recompute + read-only) + a Grok compute-verify pass, before the result is banked.

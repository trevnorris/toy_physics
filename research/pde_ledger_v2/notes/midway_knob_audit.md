# MIDWAY KNOB AUDIT — Parts I–II (a dry-run of the Part-VII codimension technique)

> **Status: COMPLETE (midway checkpoint). Layer-2 computed codimension check DONE + arbiter-verified; the Codex+Grok
> provenance re-audit DONE + reconciled; a Codex confirm-pass on the folded note landed `CONFIRM_ISSUES(1)` (a sub-split
> arithmetic fix, folded → now coherent). Corrected tallies + verdict in §2.6. Two items await the user: the C1/C2
> counting-convention decision (picks a point in ~34–43, post-Decision-16) + the R35-label refinement in the register.**
>
> **⭐ Decision-16 amendment (2026-07-21 — LANDED).** The brane polar field `P` is retired
> (`software/stage1_solver/decisions/16_retire_brane_polar_field.md`), removing `α_aniso`, `λ_Pu`, and the 3 `P`-machinery
> structural postulates (T0 P-reuse, massless spin-wave, parity-even P–u) = a **definite −5 route-less ACTION-input
> reduction**. All counts below are shifted accordingly: irreducible total ~39–48 → **~34–43**; route-less liability
> ~20–26 → **~15–21**; (c) postulated structure 19–25 → **14–20**; structural members 10 → **7**; continuous route-less
> params 10–16 → **8–14**. This −5 is the DEFINITE delta; the two parked counting-convention decisions still set the
> final point-in-range at Part VII. (Route A: the historical `SECOND_MEDIUM_DRIFT_AT_FREEZE(11)` / DOF=8 freeze record is
> RETAINED immutable as the freeze-as-run tier; these are the operative post-D16 numbers — `POST_D16_DRIFT(7)`, operative
> DOF=4, stage006 operative `DRIFT(5)`.)
> User-scheduled (2026-07-09) to run AFTER the Part-II gravity sector closes (done: count 29). This is a MIDWAY
> checkpoint over the BUILT stages only — NOT the final Part-VII count. Turns "are we accumulating too many knobs"
> into a number early enough to course-correct or hit a clean no-go (a bad number is a first-class result).
>
> **Approach (user-gated 2026-07-10):** hybrid — a reasoned tally over the register's edge classes (Layer 1, this
> chunk) + a block-decomposed COMPUTED RegionDimension/Gröbner dry-run (Layer 2, next chunk) as the anti-rig and
> Part-VII-technique validation. Highest-risk guard = a Codex + Grok provenance re-audit (an IMPOSED calibration
> mislabeled DERIVED would falsely shrink the irreducible count — a labeling error the dimension count cannot see).
>
> **Inputs:** `parameter_register.md` (the master table + R1–R48 edges + the §"two kinds of free" + the §"MIDWAY KNOB
> AUDIT" checkpoint spec), `notes/ledger_v2_blueprint.md` §8, `docs/conceptual_foundation.md` (design intent),
> and the reference method `software/stage1_solver/{tools,reports}/pathA_40_cone_lock*` (Δr=2, real-locus Krull 10→8).

---

## 0. Scope (what gets counted) + confirmed classification rules

**In scope — Parts I–II BUILT stages only:** `001, 002` (force primitives/assembly), `004–007` (I: medium +
brane freeze), `008–029` (II: gravity). **Out of scope:** `003` (light) is **Part III piloted**; Parts IV/V/VI/VII
are not yet built.

**Scope-boundary flags (on record):**
- The NG5 route-less trio `{ρ_B0, χ_c, C_hu}` — the register's *sharpest* health metric — lives in **Part III (003)
  + Part VI (041)**, so within a Parts-I–II-only audit it is a **forward-note**, NOT a counted item. Its only Part-I
  appearance is stage 006's **dead θ-branch** (`THETA_BRANCH_DEAD_NOT_ADMITTED`, R18 `CLOSED-NEG`).
- The Cluster-C SIM_DEFERRED symbols `{λ_c, I25, Ξ_Q, η_q, η_φ, rho_eff, K̄₀/₂/₄, …}` are Gate-6 branch data —
  **tracked-not-counted**, per the register. `rho_eff` ≠ stage005's `ρ0`.
- Control-construction scars `{k_warp, α, r_AL, I_wrong2, q_free, Xi_deferred, …}` — **tracked-not-counted**.

**Confirmed classification rules (user, 2026-07-10):**
- **Codimension / irreducible count** (tally #1) = nominal counted − DERIVED − CONV − discharged-reduction. **IMPOSED
  calibrations and PENDING debt STAY counted** (a tuning/debt is not a reduction). This is the predictive-surplus denominator.
- **(a)/(b)/(c) three-way split** of the irreducible set:
  - **(a) DECLARED universal constants** — fundamental-by-design, no route. Candidate: the medium primitives `{ħ, m, K, ρ0}`.
    ⚠ Even these carry open GAPs (`HBAR_FREE_SUBSTRATE_RELATION_MISSING`, `INFLOW_MASS_SOURCE_MISSING`); "declare
    fundamental vs keep reducing" is itself an open user checkpoint (register L46–47).
  - **(b) REDUCTION-DEBT** — `CALIB`/`FREE-UNREDUCED` with a **concrete** deferred computation as its route (a known
    address: the nonlinear-throat interior solve, the Gate-6 return closure, a named projection integral).
  - **(c) POSTULATED STRUCTURE** [NEW bucket] — free-postulated structure with **no concrete route** (route-less). The
    wall's route (c) `χ_B=|P_∥|²` is a high-risk future gate with falsified neighbors, so it does not make the wall (b).
- **Route-less = no CONCRETE equation/route** (a named *direction* without a concrete deferred computation is still
  route-less). Under this rule the throat packet (R30/R33/R36 → the one nonlinear-throat solve) and Route-A `{μ_R, ρ_br}`
  (R10 → same throat) are **route-ful (b)**; the NG5 trio (R12/R13, "no registered 4D→3D route yet") is **route-less**.

---

## 1. Counted inventory, by block (Layer-1 spine)

Blocks double as the Layer-2 computed-check decomposition (each ≤ ~10 params → RegionDimension tractable). "Class"
is the (a)/(b)/(c) bucket for **counted** items; DERIVED/CONV rows are shown as **collapsed** (removed from the count).

### Block M — medium primitives (I-1/I-2, stages 004/005)
| Param | Counted? | Class | Route / status |
|---|---|---|---|
| `ħ` | ✅ | **(a)** | GAP: `HBAR_FREE_SUBSTRATE_RELATION_MISSING` (declare-vs-reduce = open checkpoint) |
| `m_GNLS` | ✅ | **(a)** | GAP: `INFLOW_MASS_SOURCE_MISSING` (declare-vs-reduce = open checkpoint) |
| `K` (EOS) | ✅ | **(a)** | free input; **EOS exponent 5** is a separate `IMPOSED` structural choice → (c) |
| `ρ0` | ✅ | **(a)** | chosen asymptotic 4D-bulk number density |
| EOS exponent 5 | ✅ | **(c)** | `EOS_CLOSURE_IMPOSED` — an imposed structural exponent, route-less |
| `c_s0, ξ_h, h0, λγ` | ✗ collapsed | DERIVED | R1/R2/R3 → into `{K, ρ0, m}` (+ `c_γ`) |
| `a` (pin) | ✗ collapsed | CONV | R2 pin choice |

### Block B — brane / gauge moduli (I-3/I-4, stages 006/007)
| Param | Counted? | Class | Route / status |
|---|---|---|---|
| `μ_R` (brane shear) | ✅ | **(b)** | Route-A R10 PENDING (nonlinear-throat `μ_R`-as-bulk-defect integral — concrete deferred) |
| `ρ_br` (brane density) | ✅ | **(b)** | Route-A R10 PENDING (same throat) |
| `C_E, C_B` (bulk gauge metric) | ✅ | **(b)** | R6 brane-zero-mode reduction PENDING (`BRANE_ZERO_MODE_REDUCTION_UNDERIVED`) — ⚠ borderline-concrete |
| `μ_R⁽⁴⁾` (4D shear density) | ✅ | **(b)** | R17 projection `μ_R=∫χ_B μ_R⁽⁴⁾dw` PENDING (concrete integral); R22 firewall guards the ≠ |
| `W_slab` (slab width) | ✅ | **(b)** ⚠thin | R19 PENDING (kink admission ≠ slab stability; "requires dynamics", sim-deferred) — flagged un-earned by conceptual_foundation |
| `λ_Pu` (parity-repaired P–u) | ⊘ RETIRED | ~~(c)~~ | **RETIRED by Decision 16** (P–u twist coupling goes with the retired polar field `P`); was one of the G0 "11", drops from the operative count |
| `Ω_w` (u_w gap scale) | ✅ | **(c)** | one of the G0 "11"; **no reduction route named** → route-less |
| `g_ℓ(w; ℓ_g)` (freeze profile, 1 width knob) | ✅ | **(c)** | one of the G0 "11"; R21 scope-split; no route → route-less |
| `c_γ` (gauge cone) | (see M) | — | `= c_s·λγ`; the free thing is `c_γ`; cone-locks R7/R8 are **Part VI** (out of scope) |

### Block χ — the χ_B wall package (I-3, stage 006) — DRIFT(6) historical → **DRIFT(5) operative** (Decision 16 retires `α_aniso`)
| Param | Counted? | Class | Route / status |
|---|---|---|---|
| `χ_B` (order field) | ✅ | **(c)** | route (a) postulated; route (c) `χ_B=\|P_∥\|²` = high-risk future gate (neighbors falsified) → route-less |
| `a_B` (double-well) | ✅ | **(c)** | POSTULATED (parent `U(ρ)` single-well — wall cannot come from it) |
| `κ_B` (interface gradient) | ✅ | **(c)** | POSTULATED |
| `α_aniso` (P-orientation) | ⊘ RETIRED | ~~(c)~~ | **RETIRED by Decision 16** (P-orientation easy-plane term goes with the retired polar field `P`); drops the operative count `DRIFT(6)→DRIFT(5)` |
| `Γ_B` (conversion law) | ✅ | **(c)** | POSTULATED law; return/drain closure deferred |
| χ_B-gating structure | ✅ | **(c)** | structural 5th member of DRIFT(5) operative (was the 6th of DRIFT(6), `α_aniso` retired) |
| `δ, σ_wall` | ✗ collapsed | DERIVED | R20 → into `{a_B, κ_B}` (single-kink admission only) |

### Block W — the frozen-wall throat packet (II, stages 013/015/017)
| Param | Counted? | Class | Route / status |
|---|---|---|---|
| `μ_η` (wall inertia) | ✅ | **(b)** | R30 PENDING (nonlinear-throat solve) |
| `T_w` (wall tension) | ✅ | **(b)** | R30 PENDING |
| `β` (wall inverse-length) | ✅ | **(b)** | R30 PENDING |
| `Vp0/ℓ_c` (breathing drive) | ✅ | **(b)** | R33 PENDING (sibling of R30 — same throat) |
| `T_Ω` (ℓ=2 angular stiffness) | ✅ | **(b)** | R36 PENDING (Gate-1 `R0`-support-equation) |
| `β₂(w)` (ℓ=2 radial profile) | ✅ | **(b)** | R36 PENDING |
| `K_η` | ✗ collapsed | DERIVED | R29 `= T_w β²` |
| `M̃, K̃, T̃_Ω` (ℓ=2 radial scalars) | ✗ collapsed | DERIVED | R35 `= ∫density·β₂²` |
| `{B̃, Z̃}` (port-kernel support) | tracked | downstream-pinned | R17-adjacent; Part-VII adjudication (isotropy verdict value-independent of them) |

**⇒ The counted Part-II CALIB set = 6** `{μ_η, T_w, β, Vp0/ℓ_c, T_Ω, β₂}` — all (b), **all siblings of the ONE
deferred nonlinear-throat interior solve** (R30/R33/R36). This is block W's headline and the prime hidden-multiplicity
candidate for the Layer-2 computed check.

### Block R — return / radiation bookkeeping (II, stages 008/009/010/023)
| Param | Counted? | Class | Route / status |
|---|---|---|---|
| `d` (slab spacing) | ✅ | **(c)** | postulated executable-slab geometry (ACTION, not a medium constant); route-less |
| `ε0, ε1` (DC transmissions) | ✅ | **(b)** | R24-family / Gate-6 nonlinear return closure (concrete deferred) |
| `K0c` (ℓ=0 stiffness) | ✅ | **(b)** | R42-family scalar-reduction PENDING (⚠ dims ≠ 013/017 — convention trap; do NOT identify with raw densities) |
| `K1 = K_eta+2·T_Omega` (ℓ=1 sector) | ✅ | **(b)** | R42-family PENDING |
| `Z0_ret, Z1_ret` | ✗ | aliases | coordinate aliases of `ε0/ε1` (R42) — **zero new dofs** |
| `Z` | ✗ collapsed | DERIVED | R24 `= −M0·ε0/(1+ε0)` (accounting; `Z<0` premise stays a premise) |

### Calibrated magnitudes (tuned to benchmark, not postulates)
| Param | Counted? | Class | Route / status |
|---|---|---|---|
| force-magnitude norm (II, 002) | ✅ | CALIB (tuned) | matched to inter-defect force strength; route-less magnitude; form earned |
| `G` | ✗ | `GENUINE_BLOCKED` | the thing the PDE does NOT deliver — an external anchor, not a knob the model sets |

---

## 2. Provisional partition (pre-provenance-re-audit — ⚠ SUPERSEDED by §2.6)

> ⚠ **The buckets/counts below are the Chunk-1 provisional pass. The provenance re-audit (§2.6) corrected several:
> `C_E/C_B`→(c), `W_slab`→(c), `d`→(b), `K1`→2 dirs, +G0-structural-6→(c), ±`M̃/K̃/T̃_Ω`, ±`{B̃,Z̃}`. Cite §2.6 for the
> final range, not this section.**


> These are the Layer-1 hand counts. Layer 2 (computed RegionDimension per block) certifies the DERIVED collapses are
> genuine independent cuts and probes block W for hidden multiplicity; the Codex/Grok pass audits the provenance labels.
> **Final numbers land only after both.**

**(a) DECLARED universal constants (irreducible by design):** `{ħ, m_GNLS, K, ρ0}` = **4** (each with an open
declare-vs-reduce GAP flag).

**(b) REDUCTION-DEBT (route-ful — a concrete deferred computation):** the throat packet `{μ_η, T_w, β, Vp0/ℓ_c,
T_Ω, β₂}` (6) + `{μ_R, ρ_br, C_E, C_B}` (4) + `{μ_R⁽⁴⁾, W_slab}` (2) + `{ε0, ε1, K0c, K1}` (4) ≈ **16** — dominated by
the ONE nonlinear-throat solve (throat packet + Route-A) and the Gate-6 return closure. Honest, pending, "not-yet not never."

**(c) POSTULATED STRUCTURE (route-less liability):** the χ_B DRIFT(5) `{χ_B, a_B, κ_B, Γ_B, gating}` (5; `α_aniso`
retired by Decision 16) + `{Ω_w, g_ℓ}` (2; `λ_Pu` retired by Decision 16) + `{d, EOS-exp-5}` (2) ≈ **9** — the wall +
the imposed structural choices, no concrete route.

**Calibrated magnitudes:** force-magnitude norm (1); `G` = `GENUINE_BLOCKED` (not counted as a set-knob).

**Route-less irreducible subset (tally #2 denominator, Parts I–II):** the **(c)** set ≈ **9** (post-Decision-16) + the
force-magnitude tuning ≈ **10**. (The NG5 trio `{ρ_B0, χ_c, C_hu}` is the *canonical* route-less liability but is Part-III/VI — a
forward-note, adds ~3 when those Parts build.)

**Held-out (target-blind) predictions earned over Parts I–II** (tally #2 numerator, to be firmed up next chunk):
1. `p=2` → `1/r²` Newtonian gravity FORM survives the finite slab (pathA_29 / stage 010; counterfactual-guarded).
2. Attractive SIGN of the inter-defect force (stage 002).
3. Outgoing ℓ=2 DtN fingerprint `{u₂=1/9, u₄=4/81, v₅=1/27}` DERIVED (stage 018).
4. `χ_Q=+1` outgoing sign, from `j₂±i y₂` (stage 018).
5. The **27** in `54/5 = 2·27/5` (`derived_in_gate`; stage 020) — the earned part of the quadrupole magnitude.
6. Cross-ℓ fingerprints `{1, 1/2, 1/27}` at `{ω¹, ω³, ω⁵}` (stage 022).
7. SO(3) covariance `λ_m = ℓ(ℓ+1) = 6` (stage 016).
8. Squared-denominator prefactor factor-of-2 (`P₂` structure; stage 019).
9. Moment invariant `K̄₄ = 4K̄₂²/K̄₀` — the single open item `4d_2_5pn` left, reproduced (stage 028; CALIBRATED consistency).
10. **Falsifiable DEPARTURE**: bounded monopole/dipole `c_s`-radiation residual `∝ ε0` (GR forbids it; pathA_29).

⚠ **Honest caveat (falsification-first):** items 1–9 are **FORM / SIGN / rational-fingerprint** matches — the
**magnitudes** (`G`, the assembled `54/5`, the `2/5`) are **CALIBRATED / external_bridge_input**, NOT held-out. Item
10 is a prediction, not yet a confirmed match. So the *honest* held-out count is a handful of target-blind structural
matches + one departure — to be adjudicated against the ~10 route-less irreducible knobs (post-Decision-16) next chunk. If route-less ≳
held-out, that is the sober midway signal (reported flatly, per the model's own falsification-first stance).

**⭐ Design-intent frame (from `docs/conceptual_foundation.md`, per user steer):** the model does **not** claim
constant-minimization as the win ("postulate freely; calibration is fine; first-principles is not required; the test
is internal consistency"). So "too many knobs" is judged as **route-less liability vs held-out surplus**, with the
wall (route (a)/(c)) explicitly "the most imposed, least derived part." The doc names `W_slab`, the imposed wall
(`V_conf`/`Z/W/B_ℓ`/`k_w u_w²`), and throat `L/a` as the un-earned inputs to audit first — all present above.

---

## 2.5 Layer 2 — computed codimension dry-run (DONE, dual-engine, arbiter-verified)

Tool pair: `scripts/midway_knob_audit_codimension_sympy.py` (Gröbner-Krull dim of the initial-monomial ideal +
a positive-real-locus Jacobian-corank certificate) + `mathematica/midway_knob_audit_codimension_mathematica.wl`
(CAD `RegionDimension` over the real semialgebraic region `eqs==0 ∧ vars>0`). Materially independent routes; both
computed then **asserted** the payload (not typed), agreed on all 14 integers, and **re-ran to exit 0 from repo root
AND `/tmp`** under `timeout 600` (orchestrator arbiter re-run). Directive: `_scratch/midway_knob_audit_codimension_directive.md`.

| Block | Case | dim_before | dim_after | Δ | what it certifies |
|---|---|---:|---:|---:|---|
| **M** | baseline | 10 | **5** | 5 | the 5 collapses `{c_s0,a,ξ_h,h0,λγ}` are genuine independent cuts → residual free `{ħ,m,K,ρ0,c_γ}`. **λγ is NOT independent** (register headline). |
| M | C-M1 (vacuous R3) | 10 | 6 | 4 | replacing R3 with `λγ−λγ` leaves λγ free → R3 genuinely eliminates it (able-to-fail fires) |
| M | C-M2 (add entailed) | 10 | 5 | 5 | adding `ξ_h²−2a²` (entailed by R2a∧R2b) changes nothing → the count **cannot be inflated** by piling on dependent relations |
| M | C-M3 (inject `K=ρ0`) | 10 | 4 | 6 | a fake tie between primitives drops the dim → the machinery **would** catch real hidden multiplicity |
| **Wχ** | baseline | 10 | **7** | 3 | `K_η` collapses (R29 `K_η=T_wβ²`) + `δ,σ_wall` collapse (R20) → residual free `{μ_η,T_w,β,Vp0/ℓ_c,T_Ω,a_B,κ_B}`. Breathing set = `{μ_η,T_w,β}` = **3 not 4**. |
| Wχ | C-W1 (vacuous R29) | 10 | 8 | 2 | replacing R29 with `K_η−K_η` leaves `K_η` free → R29 genuinely eliminates it |
| Wχ | **C-W2 (inject `μ_η=T_w`)** | 10 | **6** | 4 | ⭐ a fake tie between two of the 6 CALIB knobs drops the dim → the check **is able to detect a collapse among the 6**; the real relation set finds NONE → **the 6 CALIB knobs are computationally certified independent** (answers audit item #4) |

**Findings.**
- **The DERIVED/CONV collapses used in §1–§2 are computationally certified genuine, independent cuts** — no
  vacuous relation is silently reducing the count (C-M1/C-W1), and the count cannot be inflated by entailed
  relations (C-M2). This is the falsification-first backing for tally #1 (truth in output, not prose).
- **No hidden multiplicity** among the 6 CALIB knobs or the medium primitives — the able-to-fail injections
  (C-M3, C-W2) prove the machinery *would* surface a collapse if one existed, and the real algebra shows none.
  (Contrast pathA_40, where the same technique *did* find `Δr=2` hidden multiplicity in the cone-lock.)
- **The block-decomposition works at this scale** (each block ≤10 vars → both engines trivial), validating the
  Part-VII scale-up strategy the pathA_40 scaling analysis flagged (full CAD over 15–30 vars would time out).
- **Honest limitation (stated in-engine):** `β₂(w)` is a profile and the R35 radial-scalar reductions
  `{M̃,K̃,T̃_Ω}=∫density·β₂²` are integral functionals — not polynomializable, so **out of the computed check**
  (reasoned-tally only). Block R (return admittances) and the (a)/(b)/(c) provenance LABELS are also out (the
  labels are the separate Codex/Grok re-audit). The computed check certifies ONLY the algebraic DERIVED collapses
  + the able-to-fail controls for blocks M and Wχ.

## 2.6 Provenance re-audit (Codex→Grok reconciliation) — corrected tallies + verdict

Two independent read-only reviews (Codex `gpt-5.6-sol` xhigh + Grok-4.5, each SymPy-compute-verifying the DERIVED
relations). **⭐ Highest-risk axis CLEAN:** R1, R2(×3), R3, R20 (incl. the wall-tension integral `σ_wall=√(2a_Bκ_B)/6`),
R29, R24, and the stage-024 two-port inverse **all compute-verified to 0** — no genuine algebraic IMPOSED-dressed-as-DERIVED
shrink. The catches are counting/classification, folded below. Reviews: `_scratch/midway_provenance_{codex,grok}.log`.

### Folded (uncontested — both engines, or register-backed)
- **G0 structural-3 → (c)** [Grok BLOCKING; **−3 by Decision 16**]: `SECOND_MEDIUM_DRIFT_AT_FREEZE(11)` originally included
  6 structural postulates (imposed ŵ axis + w=0 surface; free-slip uᵃ; ~~T0 P-reuse; massless spin-wave; parity-even P–u~~;
  no-C5-φ) — peers of the `χ_B`-gating structure already counted (c). I'd omitted them while counting their peers. **Decision 16
  retires the 3 `P`-machinery postulates** (T0 P-reuse, massless spin-wave, parity-even P–u), leaving 3. **(c) +3 route-less**
  (operative; was +6 historical/freeze-as-run).
- **K1 → two directions** [Codex C3, register row 166]: `K_eta`/`T_Omega` are separate generator dirs; `K1` is their
  manifestation → count 2 not 1. **(b) +1.**
- **`C_E, C_B` → (c)** [Codex C4 / Grok #3]: R6 (brane-zero-mode) names a direction with no concrete equation →
  route-less under the confirmed rule. **(b) −2, (c) +2.**
- **`W_slab` → (c)** [Codex C5]: R19 "requires dynamics", no concrete selection equation → route-less (Grok read (b)
  defensible; the strict rule + conceptual_foundation's "un-earned input" tip it to (c)). **(b) −1, (c) +1.**
- **`d` → (b)** [Codex C6]: the register's `d` row assigns the real geometry to the Gate-6 nonlinear closure (a concrete
  route). **(c) −1, (b) +1.** (`W_slab`/`d` roughly cancel in the route-less tally.)
- **Held-out list 10 → 8** [Codex C7/C8, Grok #6]: remove item 9 (`K̄₄=4K̄₂²/K̄₀` = internal calibrated consistency, not
  target-blind) and merge item 5 ("the 27" = `1/v₅`, already inside item 3's fingerprint).
- Wording: R36 ≠ the nonlinear throat (a separate Gate-1 `R0`-support equation) [C9/Grok #2]; `c_γ` collapsed by R4 into
  `{μ_R,ρ_br}` — a block-local Layer-2 residual, not double-counted [C10/Grok #7]; force-mag made explicit in the partition
  [Grok #5].

### Disputed — a Part-VII counting-convention bracket (surfaced, not force-resolved)
Two count-changing items where **Codex counts** but the **register + Grok defer/collapse** — both hinge on ONE
methodology question: *does an un-executed / downstream-deferred reduction count as a discharged collapse, or as counted debt?*
- **`M̃, K̃, T̃_Ω` (±3)** [Codex C1 — the register's-own-mislabel claim]: register calls R35 (`M̃=∫μ_η β₂²`)
  DERIVED-not-counted (the ℓ=2 analogue of `K_η=T_wβ²`). **VERIFIED at source:** the stage017 script *freezes* them as
  `calibration_inputs` (L41–42, 434, 441–442); the moment-integral is **never evaluated to bind them**. So R35 is a genuine
  moment-*definition* (functional of the counted `{μ_η,T_Ω,β₂}`) but **un-executed**. Register-methodology (reducible-in-principle
  → not counted, like K_η) vs strict (un-executed + frozen → counted debt). **±3 to (b).**
- **`{B̃, Z̃}` ℓ=2 port-kernel support scalars (±6)** [Codex C2]: register (row 181) explicitly DEFERS — "irreducibility a
  Part-VII adjudication … downstream-pinned, NOT counted in the clean headline." Codex: PENDING doesn't discharge, no
  concrete equation → count 6. **±6.**

⚠ **Register-refinement flag (for the user):** Codex C1 is a legitimate concern that R35 reads DERIVED while un-executed.
I did **not** unilaterally alter the register ([[feedback-never-alter-calibrated-process]]) — recommend a Codex-agreed
clarification of R35's status (DERIVED-in-form vs PENDING-debt) as a follow-up.

### Corrected tallies (an honest range)
| Bucket | Register-methodology | Strict-count (Codex) |
|---|---:|---:|
| (a) declared-universal | 4 | 4 |
| (b) reduction-debt (concrete route) | 15 | 18 (+`M̃/K̃/T̃_Ω`) |
| (c) postulated structure (route-less) | 14 | 20 (+`{B̃,Z̃}`) |
| force-mag (route-less CALIB) | 1 | 1 |
| **Irreducible total** | **~34** | **~43** |
| **Route-less liability (tally #2)** | **~15** | **~21** |

**⭐ The route-less set (register-methodology = 15 = (c) 14 + force-mag 1, post-Decision-16) is NOT one kind of thing — sub-split it 7 + 8:**
- **(c-struct) = 7 discrete STRUCTURAL postulates** — the χ_B wall's field/law/gating (`χ_B`, `Γ_B`, gating structure) +
  the **G0 structural-3** (Decision 16 retired the 3 `P`-machinery postulates: T0 P-reuse, massless spin-wave, parity-even
  P–u) + the imposed EOS-exponent-5. These are *modeling choices*, not continuous tuning knobs — they do
  **not** live in a continuous codimension variety (the RegionDimension technique does not apply to them).
- **(c-param) = 8 continuous route-less parameters** — `{a_B, κ_B, Ω_w, g_ℓ-width, C_E, C_B, W_slab}` (7,
  inside (c); `α_aniso`+`λ_Pu` retired by Decision 16) **+ force-mag** (1, tallied separately) — the genuine continuous
  predictive-surplus-eaters. The strict reading adds `{B̃,Z̃}` (+6) here → **~14**, taking route-less to 21.
- **Reconciles:** register 7 (c-struct) + 8 (c-param) = **15**; strict 7 + 14 = **21** — matching the tally table.

### ⭐ Verdict — the plain-number answer
**Held-out target-blind surplus (tally #2 numerator):** ~6–7 structural matches (`1/r²` exponent, attractive sign, the
ℓ=2 + cross-ℓ DtN rational fingerprints incl. the earned `27`, `χ_Q=+1`, SO(3) `λ_m=6`, the squared-denominator
factor-of-2) + **1 falsifiable departure** (monopole/dipole residual `∝ε0`). The magnitudes (`G`, assembled `54/5`,
`2/5`, `Γ̄₅`) are calibrated/external — NOT held-out (both engines confirm the note does not over-claim these).

**Are we accumulating too many knobs at the Parts-I–II midway?**
- **Route-less liability (~15–21, post-Decision-16) ≫ held-out surplus (~6–7)** — robust across both readings; the
  un-earned load still dominates the earned surplus (though the margin narrows after the `P` retirement).
- **The composition is the real story — route-less splits ~7/8:** **7 discrete STRUCTURAL postulates** (the postulated
  χ_B wall's field/law/gating + the G0 shear-surface freeze, now `P`-retired — exactly what `conceptual_foundation.md`
  names "the most imposed, least derived part") and **8 continuous route-less parameters** (~14 strict, with `{B̃,Z̃}`).
- **⭐ Even discounting the discrete postulates** (arguably not codimension-countable), the **continuous route-less tunings
  alone (~8, or ~14 strict) already exceed the held-out surplus (~6–7)** — so it is not merely a "structural postulate"
  artifact: the genuine continuous-knob load out-runs the earned predictive surplus at this midway.
- **Most of the (b) debt (~15–18) is honest, route-ful** — the throat-interior solve + the Gate-6 return closure,
  "not-yet not never"; it converts CALIB→DERIVED in one stroke if the deferred solves land. ⚠ **Post-Decision-16 the
  route-less liability (~15–21) no longer clearly EXCEEDS the route-ful debt (~15–18) — the two are now comparable**
  (before the `P` retirement, route-less ~20–26 clearly dominated it).

**⇒ The midway number is sobering but diagnostic, NOT a no-go.** The liability concentrates in (1) the *postulated wall +
freeze* — where route (c) `χ_B=|P_∥|²` is the named high-risk gate that would collapse ~half of (c-struct) — and (2) the
*un-reduced throat packet* (one deferred solve). The predictive-surplus denominator is dominated by **structural
postulates, not continuous tunings**, which sharpens where to push: **earn the wall (route c) and execute the
throat/moment reductions (R30/R33/R35/R36), or the analog stays "four sectors from one medium + an imposed wall."**

## 3. NEXT — remaining chunks (resume point)

- ✅ **Layer 2 — computed codimension dry-run: DONE + arbiter-verified** (§2.5). The DERIVED collapses are certified
  genuine independent cuts; no hidden multiplicity among the 6 CALIB knobs or the medium primitives;
  block-decomposition validated for Part VII.
- ✅ **Provenance re-audit (Codex + Grok) — DONE + reconciled** (§2.6). Highest-risk axis CLEAN (every DERIVED relation
  SymPy-verified to 0); uncontested catches folded; the two count-changing disputes (C1/C2) surfaced as a Part-VII
  counting-convention bracket. Corrected tallies + verdict landed.
- ⏭ **Two open items (both for the user):** (1) the **C1/C2 methodology decision** — do un-executed (`M̃/K̃/T̃_Ω`, R35) /
  downstream-deferred (`{B̃,Z̃}`) reductions count as discharged collapses or as counted debt? This picks a point in the
  ~34–43 range (post-Decision-16) and is a genuine Part-VII convention. (2) the **R35-label refinement** — whether to reclassify R35 in the
  register from DERIVED to DERIVED-in-form/PENDING-debt (a Codex-agreed edit, not done unilaterally).
- ⏭ **Optional closing bookend leg:** a Codex confirm-pass on the folded note (the folds are orchestrator prose over two
  surfaced-not-hidden disputes, so this is a lower-risk final check).

> **⏸ PAUSE POINT (2026-07-11):** Both audit layers DONE + verified. Layer 1 (inventory) user-approved; Layer 2 (computed
> codimension check) dual-engine arbiter-verified; the provenance re-audit reconciled (Codex+Grok) + a Codex confirm-pass
> (one sub-split arithmetic fix folded) → corrected range + verdict in §2.6. **Headline (post-Decision-16 −5): route-less
> liability ~15–21 ≫ held-out surplus ~6–7 (sobering but diagnostic; now comparable to the route-ful debt ~15–18) — split
> ~7 discrete structural postulates (the postulated wall + G0 freeze, `P`-retired) vs ~8 continuous route-less params; the
> continuous half alone (~8–14) already out-runs the held-out surplus.**
> Resume = the two user items above.

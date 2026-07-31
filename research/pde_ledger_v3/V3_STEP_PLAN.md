# v3 STEP PLAN — gravity, light, gravitomagnetism

**Read `CHARTER.md` first** (scope, what it excludes, the two hard constraints). This file is the
step-by-step route. **Method**: `docs/derivation_walkthrough_plan.md` — ⛔ not restated here.

---

## 0. ⭐⭐ How a step is run — SIDE BY SIDE (user decision, 2026-07-31)

⛔⛔ **This supersedes the v2 working style, and the reason is measured.** On v2 the user delegated
heavily to go fast. The cost: the `a`-pin sat at the foundation for eleven months, the lepton tower's
falsification never reached the front-door docs, and the Spin Problem sat in an orphan note. ⇒ **Speed
bought nothing and hid three structural findings.**

⭐ **The user works through every step with me. Not a gate at the end — a participant in the
derivation.** The purpose is that *they understand everything that is going on*, because they are the
only physics reviewer this project has.

### What "side by side" means operationally

1. ⛔ **Do NOT pre-derive a step and present conclusions.** Derive it in the open, one move at a time,
   with the reasoning stated **before** the result so it can be caught in flight.
2. **Open each step with the setup, not the answer:** what we are deriving · what we already have and
   from which step · what is missing · where I expect trouble.
3. **Stop at every substantive move**, not only at step boundaries. A step is not an atomic unit of
   conversation; a *derivation move* is.
4. ⭐⭐ **Flag every identification before making it.** *"I am about to treat X and Y as the same
   quantity — here is why, and here is what would make that wrong."* Ten pin-shaped identifications are
   already on record; this is the exact move that produces them, and it is the one place a second pair
   of eyes in real time is worth more than any downstream review.
5. **Ask rather than assume** on physical-picture calls. The maths is checkable; the picture is the
   user's, and this session showed their instinct catching things the corpus had wrong.
6. **Sub-agents are for LOOKUP, not judgement** — "find me this locus", "what does file X say about Y".
   ⛔ Not "decide whether this is right".
7. **Codex and the external reviewers come AFTER we have walked it** — to check our work, not to do it.
8. **Size steps to fit one sitting.** ⛔ If a step below is too big to walk in one go, split it and say
   so; the 23 steps are a route, not a contract.

### Then, and only then, the mechanical part

9. **Write the step record** — the eight fields the method doc specifies. The record is the deliverable,
   not a report about it.
10. **Add quantities and equations to `reduction/`** and run the dimensional gate. Must pass before the
    step is banked.
11. **Confront the register rows this step touches.** Each step below names them. ⛔ A step that walks
    past a row without resolving *or explicitly deferring* it is not done.
12. **Then the review legs** (below), then commit.

⚠ **This is slower, deliberately.** ⛔ Do not optimise it back. If a session covers two steps properly,
that is the plan working. The failure mode this replaces is *twenty-three steps banked and nobody able
to check them.*

⚠ ⛔ **Do not read this as licence to build apparatus instead.** Slower means *more explanation per
derivation*, not more machinery around it.

**Review load, scaled to what the step is:**

| step contains | legs |
|---|---|
| bookkeeping / re-banking / prose | one fresh reviewer |
| ⭐ **new physics, or a number that meets data** | **two independent fresh legs, both blocking** |
| a claim that something is *derived* which was previously *postulated* | two legs + explicit ablation |

⛔ **The builder never reviews its own work.** Codex builds; Grok and/or a fresh agent review; the
orchestrator adjudicates and does not type.

**Definition of done for a step:** record banked · gate green · register rows confronted · reviews clean
· committed.

---

## PHASE 0 — the substrate (4 steps)

Re-banked from v2, forward, with the corrections this session established. ⭐ Cheap — but it is what
makes v3 readable without v2, which is acceptance criterion (3).

### S1 · The medium, and the brane as its ordered state
One GNLS condensate in a 4D bulk; `ρ = |ψ|²`; `P = Kρ⁵`; the brane is the **ordered state** of the same
medium (`χ_B = 1`), not a second object.
**Expected new:** 4 scalars `{ħ, m_GNLS, K, ρ0}` · 3 discrete `{D=4, n=5, the two-phase split}` · 3 fields
· 1 BC.
**Register:** none. **Carry forward:** the two-phase split is `postulated` — the largest tier-1 item.
⚠ **Correct v2's stated reason.** `model_map.md:26` justifies the postulate by *"`U(ρ)` is single-well"*.
That reasoning is **wrong** — the brane comes from `V_χ(r_B)`, a potential on a *different* field. The
conclusion survives, the argument does not. → S5.

### S2 · The sound speed
`c_s² = (1/m)dP/dρ = nKρ^(n-1)/m`. **Expected new: nothing** — pure consequence.
**Finding to re-bank:** `[K] = M L^(4n−2) T⁻²`, so `n` and `[K]` are one structural choice.
**Register:** none. **Open:** O-02 (is `K`+`n` one entry or two — user call, still open).

### S3 · The dimensional foundation — three pins, not four
`{c_s0, ħ, m_GNLS}` over `{L,T,M}`: `det = 1`, rank 3, **nullity 0**. A complete unit system.
⭐ **This is the step that must state A1's lesson**, not just its fix: a fourth pin of dimension `L`
*must* produce `a·m·c_s0/ħ` — forced arithmetic. The **form** was never in doubt; the error was reading
it as a defining equation for a **defect** quantity.
**Register:** **A1** (retired — record why), and the ⭐ **kind-test** that A7's investigation produced:
> ⛔ Ask *"is this the same KIND of thing?"* — one number for the whole medium, or one per particle —
> **not** *"is this value right?"*. The latter passes the pin.

### S4 · Healing length and enthalpy scale
`h0 = m c_s0²/4`; `ξ_h = √(ħ²/(2 m h0)) = √2 ħ/(m c_s0)` — derived from core balance.
**Expected new: nothing.**
**Register:** **A7** — record both conventions and *why* the no-√2 form is forced by the BdG ratio
`(kξ/2)²`. ⛔ Do not "unify" them; they are different lengths with different jobs.

---

## PHASE 1 — the brane (4 steps) ⚠ the honest core

### S5 · The order field and the wall potential
`χ_B = r_B e^{iθ_B}`, `r_B ∈ [0,1]`; `S_χ = ∫ ½Z_χ(∂_t r_B)² − ½κ_χ|∇₄r_B|² − (λ_χ/4)r_B²(1−r_B)²`.
**Expected new:** the double-well **form** (a function, not a scalar) + `{Z_χ, κ_χ, λ_χ}` or the
`{a_B, κ_B}` pair.
**Register:** **C8** (the wall rests on two *postulated* constants — derived-from-postulates is not
derived-from-primitives). **B2** — ⛔ **state the no-go here**: the medium program that would have
grounded this returned `FAIL_COUPLE_STRESS_NOGO` and is superseded.

### S6 · The kink — wall thickness and surface tension
Solve the EL equation: `δ = √(κ_B/(2a_B))`, `σ_wall = √(2a_Bκ_B)/6`. **Genuinely derived** — from S5's
postulates.
**Register:** **A2** (`ℓ = δ` is an *identification*, not a derivation — confront it here or nowhere).

### S7 · ⛔ The slab width is NOT selected
`W_slab` is `FREE-UNREDUCED`; *"double-well selects NO width"*. The kink gives one interface; the brane
is a finite slab.
⭐ **A step whose result is "this is not determined" is a real step.** Bank it as such.
**Register:** **C7**, **A9** (`W_slab` merged with the `L/a` debt without an equation).

### S8 · The quadratic brane Lagrangian
Transverse and longitudinal sectors; the reduced `h`/`u_L` operator.
⭐⭐ **The step that defines v3's ceiling.** State plainly: **quadratic = linear response about an
assumed equilibrium.** Everything downstream in this ledger is small-oscillation physics on a brane
that S5–S7 postulated. It is why gravity and light are tractable, and why the defect is not.
**Register:** **C6** (no closed parent action — the coupling to sleeve/geon/drain does not exist).

---

## PHASE 2 — light (3 steps)

### S9 · The transverse shear speed
`c_γ² = μ_R/ρ_br`. **Expected new:** `c_γ` — or rather `{μ_R, ρ_br}`.
⭐⭐ **This closes O-01, the "universe hole."** v2's registry seed required `c_γ` as a supplied input
while **no step ever introduced it**. It enters *here*, in the excitations phase, with provenance —
⛔ never back-filled into the medium block for bookkeeping convenience.
**Register:** the `{μ_R, ρ_br}` R10 debt begins here → S22.

### S10 · Two transverse photons
The earned target-blind result: brane shear gives exactly two transverse polarisations.
**Expected new: nothing.** ⭐ A genuinely earned step — say so.

### S11 · The stray-longitudinal departure
Characterized, first-class. ⛔ Not a defect to fix — a recorded departure from the reference theory.

---

## PHASE 3 — the drain (3 steps) ⚠ where the interior is deferred

### S12 · The throat as a boundary condition
`∂_t ρ + ∇·(ρv) = S_drain + S_leakage`, and the momentum/energy partners.
⛔ **All six source terms are `[OPEN]`.** The throat enters this ledger **only** as a boundary condition
of given strength — that is the whole content of the scope boundary.
**Register:** **C6**; **C1** (⚠ the drain law is *parameterized, not derived* — and `g_phys` was never
mapped; carry the caveat, not the top line).

### S13 · The drain rate `J`, and what it is not
`J` is the invariant label. ⛔ **`m_defect` is a GAP** — `INFLOW_MASS_SOURCE_MISSING`; only a dimensional
bridge `α_J ħJ/c_γ²` exists, and pathA_21's verdict is `MASS_BRIDGE_FORM_NOT_DERIVED`.
⭐ **This is the single most load-bearing gap in the gravity sector**: the source of gravity is not
derivably connected to mass. Bank it as an unresolved step.
**Register:** **C3**.

### S14 · The far field
`1/r²` on the brane, `1/R³` in the bulk. The exponent comes from the slab DC-sink completions and the
zero mode — ⛔ **not** from the throat.

---

## PHASE 4 — gravity and gravitomagnetism (6 steps)

### S15 · Two drains attract
Noether stress at infinity, drain strengths `Q_i` **given**: `F_12 = −(m N_∞,3 Q1Q2/4πr²) r̂`.
Form and sign **earned target-blind**; magnitude falls to the interior.

### S16 · ⭐⭐ The worldtube reduction — the theorem the whole scope boundary rests on

⛔ **The load-bearing step of this entire ledger.** It is what makes the scope boundary legitimate rather
than convenient — and equally, the reason a clean v3-gravity is **not** evidence the model is sound.

⚠ **Read the actual statement, not the slogan** (`research/4d_2_5pn/paper/4d_2_5pn.tex:605-614`):

> *"For a **compact defect of size `a ≪ r`** in a smooth external field, the center-of-mass reduction
> gives the usual point-particle force law … `+ O(a²/r² · Φ_ext/r)` … so that the universal
> point-particle source **at leading order** is controlled by the worldline/worldtube multipoles rather
> than by arbitrary internal details of the defect."*

⇒ **The independence is (i) leading-order only and (ii) conditional on compactness.** The first
correction is `O(a²/r²)` — controlled by the throat radius `a`, which is **undetermined** (**C2**,
**D1**). So "gravity does not depend on the interior" is precise only to leading order; beyond it,
gravity depends on exactly the quantity this project cannot compute.

⭐⭐ **And it collides with C10.** The worldtube reduction's premise is that the defect is **compact**.
The Spin Problem's verdict is that *"you cannot get frame dragging from a compact 4D bubble; you need
the tail."* ⛔ **If the tail is real, S16's premise fails** — and with it the justification for scoping
gravity apart from the interior at all.

⇒ **S16 must be run BEFORE S19 is trusted, and its result feeds back into the charter.** Two blocking
legs. ⚠ Possible outcome: the scope boundary is narrower than the charter claims, in which case say so
rather than preserving the boundary.

### S17 · The PN ladder — cite-only
Seven separately-published papers, Zenodo-DOI'd (`ledger_stage029_pn_corpus_doi_cite.md:32-40`).
⛔ **Do not re-derive.** Record what each establishes and its inputs.

⚠ **Two precisions the word "audited" hides — checked 2026-07-31:**
- **`4d_2_5pn` and `4d_4pn` are *conditional* derivations** (their own titles say so). ⛔ Do not carry
  them forward as unconditional; the condition is part of the imported result.
- **The DOIs are README-authoritative, not source-cross-validated.** *"none of the seven papers declares
  its own Zenodo DOI in source"* (`:51`). ⛔ Repeat that caveat rather than dressing them as
  source-verified.
- stage029 itself is **cite-only**: no scripts, no dual-engine, no re-audit. *"029 does NOT re-audit the
  PN physics."* ⇒ The audit lives in the papers, not in the ledger.

### S18 · 2.5PN Burke–Thorne matchback
⚠ **Record honestly**: `INV2` is `a`-free — `K̄₀ = (54/5)Gc_s⁵/(a⁵c⁵)` times `a⁵/(27c_s⁵)` is
`2G/(5c⁵)` identically. It is the **rational identity** `54/(5·27) = 2/5`, i.e. calibrated closure
consistency, ⛔ **not** a first-principles derivation of `G`.
**Register:** **E1** — v2's count claims here rest on tagging a physical `a` as `CONV`.

### S19 · Gravitomagnetism ⛔ NOT a cite-and-record step — it has a known open problem
The frame-dragging / velocity-dependent PN terms. ⛔ Not the EM sector's moving throat.

⚠ **Checked 2026-07-31, and the plan's original assumption was wrong.** This is *not* already settled
inside S17's ladder. `conceptual_foundation.md:589` still lists spin as *"**not yet placed in the
picture**"*, while an **orphan** note (`research/4d_1pn_bridge/notes/tadpole.md`, referenced by nothing
in the corpus) claims a *"Standard Model of the Defect"* that solves it.

⭐⭐ **And the note identifies a structural tension inside the gravity sector itself:**
- the **inertial** sector (1PN precession) forces a **compact, "stubby"** throat;
- the **spin** sector needs a **dipole** gravitomagnetic falloff, and a compact source gives a
  *"gravitomagnetic monopole"* — *"physically inadmissible"*;
- verdict: *"You cannot get frame dragging from a compact 4D bubble; you **need** the tail."*

⇒ **Inertia wants compact, spin wants extended.** If that holds, it is a conflict **within the sector
v3 calls solid** — exactly the kind of thing this ledger exists to surface.

⚠ The proposed fix (a composite "stubby head + infinite vortex-filament tail") re-introduces
**quantized circulation**, which the charge sector explicitly disclaims. ⛔ Do not adopt it; **audit it**.

**Register:** **C10**. ⭐ **Two independent legs, both blocking** — this is new physics, not bookkeeping.
⚠ **Possible outcome: a NO-GO for the sector.** That is a first-class result, not a failure.

### S20 · The falsifier
stage009 `RETURN_RESIDUAL_PREDICTION`: `1−T_ℓ(0) = ε_ℓ/(1+ε_ℓ)` at ℓ=0/1, orders `p_res = 1,3` — a
departure **GR forbids**. ⭐ The sector's one live able-to-fail prediction. Bank it prominently.

---

## PHASE 5 — the knit and the deliverable (3 steps)

### S21 · What light and gravity share
Both ride the brane; both consume `{ρ_br, μ_R}`. ⛔ Per the method doc's retraction, integration may
introduce no new action/constitutive/source/BC input **only** given a frozen complete parent action —
and **C6** says there isn't one. So this step's honest form is: *name what integration needs that phase
0–1 did not supply.*

### S22 · ⭐⭐ THE DELIVERABLE — the interior debt list
One table. Every quantity the gravity and light sectors **assume from the throat interior**, with locus
and what would discharge it:

| debt | what it blocks | discharge |
|---|---|---|
| `{μ_R, ρ_br}` (R10) | `c_γ` and every cone lock | the nonlinear throat solve |
| `{μ_η, T_w, β}` (R30) | frozen-wall response | same solve |
| `{Vp0/ℓ_c}` (R33) | breathing drive | same solve |
| `{T_Ω, β₂}` (R36) | ℓ=2 support | ℓ=2 support equations |
| `m_defect` ↔ `J` (**C3**) | ⭐ **the source of gravity itself** | unbuilt — no named route |
| density-port magnitudes | Gate-6 numbers | SIM-deferred |

⇒ **This is what v3-gravity is for.** Today it is scattered across six documents; here it is one page.

### S23 · The count — sim-input set, partitioned
⛔ **Never one integer.** Scalars · profiles/functions · BCs and domains · discrete choices.
⛔ **No number leaves without the §7a closing certification** (admission gate · transitive leaf closure ·
block rank · top-down reconciliation · the sim-input-vs-residual diff). The walkthrough produces an
*inventory*; only §7a certifies it.

---

## 1. Order of attack, and the honest expectation

**S1–S4** are cheap re-banking; they exist to make v3 self-contained.
**S5–S8** are the honest core — three of the four bank a *postulate* or a *gap*.
**S9–S11, S15–S20** are where the earned results live.
**S22 is the deliverable.**

⚠ **Expect roughly this shape:** ~6 steps introducing genuinely new inputs, ~8 pure-consequence steps,
~4 steps whose entire result is *"this is not determined"*, and 1 (S22) that is the point of the
exercise. ⛔ If S5–S8 start producing derivations instead of postulates, something is being smuggled —
that is the phase where the model is known to be weakest.

## 2. ⛔ Failure modes to watch, specific to this plan

1. **Re-deriving the PN ladder.** S17 is cite-only. It is audited and DOI'd; touching it is pure cost.
2. **Letting S22 become a census.** It is one table with loci. It is not a provenance system.
3. **Treating a clean close as validation.** See `CHARTER.md` §1.1. A tidy v3-gravity means the far field
   does not depend on the interior — already known.
4. **Back-filling `c_γ` into phase 0** to make the registry balance at S2. That is the exact error O-01
   recorded.
5. **Building registry machinery as a precondition** for banking a step. Add one quantity and one
   relation when a step needs them; nothing more.

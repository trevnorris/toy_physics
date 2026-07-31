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

## S0.5 · ⛔ PRE-S1 REGISTRY SURGERY — required, and the plan originally missed it

**Found by review.** *"Reuse `reduction/`"* (charter) and *"`c_γ` enters at S9, never back-filled into the
medium block"* (S9) **collide**, because the live seed **already** carries `Q.medium.c_gamma` and
`Q.medium.lambda_gamma` as *medium* quantities with `R3: λγ = c_γ/c_s0` marked `DERIVED-EXECUTED`.

⇒ Executing S1–S3 against the seed as-is either **accepts the O-01 universe hole** or requires surgery
that no step names. Do the surgery **first**:

- `c_γ` and `λγ` leave the medium counting contract until S9 / the cone-lock step introduces them with
  provenance;
- recompute acceptance (⛔ never preserve — same rule as step ①);
- ⛔ **do not "resolve" it by leaving them where they are.** That is the back-fill S9 forbids.

⛔⛔ **This is NOT a one-line scope edit — round 2 confirmed it will break code.** `active_variables`
is computed from `counting_axis`, and `C-M1`, `able_to_fail.py` and `test_registry.py` are all built
around `R3`. Removing `c_γ`/`λγ` from the medium contract therefore requires, in one pass: the
counting-axis change · a replacement for the `R3`-dependent mutation case · the harness and tests
updated · acceptance recomputed and **independently re-derived**, never preserved. ⇒ Treat S0.5 as a
step in its own right with its own review legs, ⛔ not a preamble to S1.

⛔ **"Change the counting axis" is NOT an executable prescription.** Moving `c_γ`/`λγ` to
`convention-orbit` or `discrete-structural` **misclassifies** them, and changing `scope` alone does not
affect `active_variables`. ⇒ Specify the lifecycle exactly: **retire the two quantities and deactivate
/remove `R3`**, then re-introduce `c_γ` at S9 and `λγ`/`R3` at S20a.

⛔ **Enumerate every touchpoint before starting** — at minimum `registry_read.py`'s runnable propagation
smoke test (`:1014`) and `README.md`'s canonical example (`:174`), **both hardwired to `lambda_gamma`**,
plus `acceptance_check.py`, `able_to_fail.py` and `test_registry.py`.
⭐ A mutation case built on `R3` needs a replacement built on a surviving relation.
⛔ Derive the new payload independently; ⛔ never preserve, and ⛔ never read it off a prior document.

---

## PHASE 0 — the substrate (4 steps)

Re-banked from v2, forward, with the corrections this session established. ⭐ Cheap — but it is what
makes v3 readable without v2, which is acceptance criterion (3).

### S1 · The medium, and the brane as its ordered state
One GNLS condensate in a 4D bulk; `ρ = |ψ|²`; `P = Kρ⁵`; the brane is the **ordered state** of the same
medium (`χ_B = 1`), not a second object.
**Expected new:** 4 scalars `{ħ, m_GNLS, K, ρ0}` · 3 discrete `{D=4, n=5, the two-phase split}` · 3 fields
· 1 BC.
**Register:** ⛔ **A13 — GATE. Resolve before this step banks.** `χ_B` is used here as a **real material
fraction** (so `χ_B·n` is real) while S5 writes it **complex**. Choose the branch — real/dissipative or
complex/inertial — propagate the choice into the field count, the action and the conversion balance, and
**remove the unchosen inputs**. ⛔ S1 cannot bank on "Register: none".
**Carry forward:** the two-phase split is `postulated` — the largest tier-1 item.
⚠ **Correct v2's stated reason.** `model_map.md:26` justifies the postulate by *"`U(ρ)` is single-well"*.
That reasoning is **wrong** — the brane comes from `V_χ(r_B)`, a potential on a *different* field. The
conclusion survives, the argument does not. → S5.

### S1.5 · The substrate action and the Madelung balances ⭐ ADDED — was owed and unnumbered
The GNLS action, the quantum-gradient term, and the Madelung mass/momentum/energy balances, with the
pressure identity that fixes `U(ρ) = Kρ⁵/4`.
⛔ **Without this step the chain is broken**: S4 invokes a "core balance" and S12 needs the
momentum/energy partners, and neither exists in S1's output. *"Derive forward from what previous steps
banked"* is false at S4 until this is in.
⛔ `U(ρ)` is **not** a third field or an independent free function — it is fixed by the EOS.
**Expected new:** nothing beyond S1's primitives, if the action is faithfully transcribed. ⚠ If it turns
out to introduce something, that is a finding.

⛔⛔ **Order: S1 THEN S1.5.** This step consumes S1's primitives, so it cannot precede it. ⚠ Earlier
handoff text said `S0.5 → S1.5 → S1`; that was **wrong** — the correct order is
`S0.5 → S1 → S1.5 → S2`.

⛔ **`U(ρ)` is NOT uniquely fixed by the pressure identity.** `P = ρU′ − U = Kρ⁵` gives
`U = Kρ⁵/4 + Cρ` — the linear term is unconstrained. `C` is a **chemical-potential / phase-reference
choice**, and it matters to the energy and source ledger. ⇒ **State the reference that sets `C = 0`
explicitly**; do not present `Kρ⁵/4` as forced.

⚠ Also **fix the field count**: S1 lists three fields including `U(ρ)`, while this step says `U` is not a
field. Name the two actual fields and correct S1.
⇒ Bank an explicit, consistent `Π_Q`, `ε_Q`, `j_ε^Q` improvement convention — the three form one package.
⚠ This step supplies the **conservative left-hand sides only**; S12's non-variational source partners
stay open.

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
⚠ **The equation below is a TEMPLATE, not a pre-decided branch** — S1's A13 gate chooses first, and this step rewrites to match.
`χ_B = r_B e^{iθ_B}`, `r_B ∈ [0,1]`; `S_χ = ∫ ½Z_χ(∂_t r_B)² − ½κ_χ|∇₄r_B|² − (λ_χ/4)r_B²(1−r_B)²`.
**Expected new:** the double-well **form** (a function, not a scalar) + `{Z_χ, κ_χ, λ_χ}` or the
`{a_B, κ_B}` pair.
**Register:** ⛔ **A13 (gate — must already be resolved at S1; this step must use the chosen branch, and
⛔ `{a_B, κ_B}` cannot parameterize an action containing `Z_χ`)** · **C8** (the wall rests on two *postulated* constants — derived-from-postulates is not
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

### S8 · The quadratic brane Lagrangian ⭐ and `{ρ_br, μ_R}` enter HERE
Transverse and longitudinal sectors; the reduced `h`/`u_L` operator.
⛔ **Introduce and classify `{ρ_br, μ_R}` in this step** — the transverse action contains both *before*
`c_γ² = μ_R/ρ_br` can be derived, so introducing them at S9 would invert provenance.
⇒ **The R10 debt starts here**, not at S9.
⭐⭐ **The step that defines v3's ceiling.** State plainly: **quadratic = linear response about an
assumed equilibrium.** Everything downstream in this ledger is small-oscillation physics on a brane
that S5–S7 postulated. It is why gravity and light are tractable, and why the defect is not.
**Register:** **R10 starts here** (`{ρ_br, μ_R}`) · **C6** (no closed parent action — the coupling to sleeve/geon/drain does not exist).

---

## PHASE 2 — light (3 steps)

### S9 · The transverse shear speed
⛔ **Provenance order corrected (round 3): `{ρ_br, μ_R}` must be introduced and classified in S8, not
here.** S8 cannot write the transverse quadratic brane Lagrangian without them, so introducing them
afterwards violates forward provenance and puts the apparatus before its own coefficients.
⛔ **COMMITTED LAYOUT (round 4 — the earlier "pick one" left the ambiguity live):**
`{ρ_br, μ_R}` are **introduced and classified in S8**, and the **R10 debt starts in S8**.
⇒ S9 is **pure consequence**: `c_γ² = μ_R/ρ_br`. **Expected new: nothing.**
⭐⭐ **This closes O-01, the "universe hole."** v2's registry seed required `c_γ` as a supplied input
while **no step ever introduced it**. It enters *here*, in the excitations phase, with provenance —
⛔ never back-filled into the medium block for bookkeeping convenience.
**Register:** none — R10 began at S8. → S22.

### S10 · Two transverse photons
The earned target-blind result: brane shear gives exactly two transverse polarisations.
**Expected new: nothing.** ⭐ A genuinely earned step — say so.

### S11 · The stray-longitudinal departure
Characterized, first-class. ⛔ Not a defect to fix — a recorded departure from the reference theory.

---

## PHASE 3 — the drain (3 steps) ⚠ where the interior is deferred

### S12 · The throat as a boundary condition
⛔⛔ **CORRECTED 2026-07-31 — this step originally wrote a drain law the user had already RULED OUT.**

The committed non-variational drain is the **dynamical `Γ_B` order-conversion** law
(`stage045_nonvariational_block_prep.md:23`):

```
∂_t(χ_B n) + ∇₄·(χ_B n u + J_χ) = n Γ_B ,    Γ_B = Γ_return − Γ_drain ,    [Γ_B] = T⁻¹
```

An **internal order-conversion at the throat** — the dynamical wall converting ordered↔disordered
material. Phase-conversion, **not suction**.

⛔ **RULED OUT** (`:19`, user, 2026-07-24): the G0-card `S_drain` **ρ-mass-sink + remote return**, i.e.
`∂_t ρ + ∇·(ρv) = S_drain + S_leakage`. It is **frozen-wall-premised**, and the standing principle is
*"the drain cannot be a frozen wall — several instances where we tried to freeze the wall and it screwed
up all the calculations."* ⇒ If a static limit is ever wanted it must be **derived** as the `ω→0` limit
of the dynamical law, ⛔ never posed as an independent frozen solve.

⚠ The G0 card's *structure* is still reusable — return controllers, `f`/`w` balances, `𝔅`/mouth/collar/IR
BCs, the `F_var+F_flux+F_𝔅+F_rad` partition — but it must be **re-hosted** on the dynamical law.

⭐ **Note how this error happened**: I had the "never freeze the wall" rule in memory and wrote the frozen
form anyway, because I transcribed the equation from a spec section rather than deriving it forward. That
is precisely what §0's side-by-side rule exists to prevent.
⛔ **All six source terms are `[OPEN]`.** The throat enters this ledger **only** as a boundary condition
of given strength — that is the whole content of the scope boundary.
**Register:** ⛔ **A13 (the conversion balance must use the branch chosen at S1)** · **C6**; **C1** (⚠ the drain law is *parameterized, not derived* — and `g_phys` was never
mapped; carry the caveat, not the top line).

### S13 · The drain rate `J`, and what it is not
`J` is the invariant label. ⛔ **`m_defect` is a GAP** — `INFLOW_MASS_SOURCE_MISSING`; only a dimensional
bridge `α_J ħJ/c_γ²` exists, and pathA_21's verdict is `MASS_BRIDGE_FORM_NOT_DERIVED`.
⭐ **This is the single most load-bearing gap in the gravity sector**: the source of gravity is not
derivably connected to mass. Bank it as an unresolved step.
**Register:** **C3**.

### S14a · ⭐ RUN AFTER S13, BEFORE S14 · ⛔⛔ THE DRAIN BRIDGE — BLOCKING, and created by fixing S12

**Round 3 found that correcting S12 severed the chain to the imported gravity results.**

The committed dynamical `Γ_B` law **conserves total material and converts only order**
(`ledger_stage006_two_phase_chiB_ontology.md:81`). But the imported force derivation assumes a **Gauss
drain removing number flux** (`ledger_stage002_matter_stress_force_assembly.md:53,:111`), and S14/S20
use frozen DC sink/return completions. `Q = Θ_Q J / N_∞,3` *labels* the far-field flux; it does **not**
derive it from `Γ_B`.

⇒ **Derive, from the dynamical system:** the projected order-loss source · the definition of `J` · the
profile-dependent `J → Q` map · the controlled `ω → 0` return law.

⛔ **Until this succeeds, S14, S15 and S20 are CONDITIONAL, not carried results.** ⚠ Mark them so in
their records. ⭐ This is a real physics obligation produced by removing a frozen-wall shortcut — exactly
the kind of debt this ledger exists to surface, not a bookkeeping chore.

### S14 · The far field
⛔ **CONDITIONAL ON S14a.** Not a carried result until the dynamical-`Γ_B` bridge succeeds.
`1/r²` on the brane, `1/R³` in the bulk. The exponent comes from the slab DC-sink completions and the
zero mode — ⛔ **not** from the throat.

---

## PHASE 4 — gravity and gravitomagnetism (6 steps)

### S15 · Two drains attract ⚠ RUN S16 FIRST, or fold it in as the premise
⛔ **CONDITIONAL ON S14a.** Not a carried result until the dynamical-`Γ_B` bridge succeeds.
⚠ **Corrected twice.** Round 1: run S16 first, since the corpus establishes the worldtube closure before
the particle law. ⛔ **Round 2 corrects that correction: S16 does NOT license S15.** S16 is
*response-side* (how a body moves in a field); S15 is *source-side* (what force two drains exert) and
carries **its own** Noether / control-surface assumptions. ⇒ Run S16 first for regime context, but
⛔ derive S15's premises independently — do not inherit them.

Noether stress at infinity, drain strengths `Q_i` **given**: `F_12 = −(m N_∞,3 Q1Q2/4πr²) r̂`.

⛔ **Two corrections to the original wording:**
- **The source map is missing.** S13 introduces `J`; this step uses `Q_i` — the corpus relation
  `Q_i = Θ_Qi J_i / N_∞,3` must be an explicit step, not an unremarked substitution. Add
  `{J, Θ_Q, N_∞,3, Q}` to S22.
- **The sign claim is overstated.** This is the **leading matter-stress** contribution only; the source
  says the full sign retains quantum, confinement, Maxwell and profile residuals, and the normalization
  is **calibrated**. ⇒ Bank it as *"leading matter-stress sign"*, ⛔ not *"the force sign, earned"*.

### S16 · ⭐⭐ The worldtube reduction — the theorem the whole scope boundary rests on

⛔ **The load-bearing step of this entire ledger** — but ⛔ **not** because it "legitimises the scope".
Round 2: it does not establish interior-independence at all, so do not carry that framing forward. What
it gives is a **conditional response-side approximation**, and the honest use of it is to enumerate what
must be *supplied* for the far field to be computable without the interior.

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

⛔⛔ **And it is weaker still — an external review found the theorem is on the wrong SIDE.** What the
corpus establishes is a **response-side** result: how a compact body *moves in* a given external field
(`M_A ẌA = −M_A∇Φ_ext`). It is **not** a source-side statement about what field a defect *produces*.
`4d_1pn_full.tex:657,:673,:695,:734,:825` additionally assume compact support, smooth external fields,
a vanishing dipole, a quasi-static conservative regime, **a calibrated potential `−Gm/r`**, and
**supplied mass and multipole moments**, with boundary flux discarded by assumption.

⇒ **Rewrite S16 as a conditional response-side monopole approximation** — ⛔ not a source-generation or
interior-independence theorem — and list `M`, `a`, the multipoles, compactness, source normalization and
boundary-flux suppression as **supplied inputs** that belong in S22's debt table.

⇒ **S16 must be run BEFORE S19 is trusted, and its result feeds back into the charter.** Two blocking
legs. ⚠ Likely outcome: **the scope boundary is materially narrower than the charter claims.** Say so
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
⛔ **CONDITIONAL ON S14a.** Not a carried result until the dynamical-`Γ_B` bridge succeeds.
stage009 `RETURN_RESIDUAL_PREDICTION`: `1−T_ℓ(0) = ε_ℓ/(1+ε_ℓ)` at ℓ=0/1, orders `p_res = 1,3` — a
departure **GR forbids**. ⭐ The sector's one live able-to-fail prediction. Bank it prominently.

---

## PHASE 5 — the knit and the deliverable (4 steps)

⚠ `S20a` is numbered from phase 4 but belongs here: the cone lock is a **knit** question, not a
gravity-sector derivation.

### S20a · ⭐⭐ THE CONE LOCK — the actual cross-sector question, and it had no step
**Found by review, and it is the load-bearing knit for *this* scope.** For a gravity + light ledger the
cross-sector identification that matters is **not** "both ride `{μ_R, ρ_br}`" — it is **whether the light
speed and the gravity-signal speed are forced equal**: `λγ = c_γ/c_s`.

`conceptual_foundation.md:595` still lists **"Cone lock `λγ`"** as **open**; the corpus records
`CONE_LOCK_CALIBRATED` with `λγ = 1` as a **calibration, not a derivation**; and `R3` already encodes the
ratio in the registry seed.

⛔ **Two different objects — round 2. Do not merge them:** the **ratio** `λγ = c_γ/c_s` is a *derived
definition* (registry `R3`); the **equality `λγ = 1`** is a separate *calibrated / uncommitted* cone
lock. Classifying the ratio says nothing about the lock.

⇒ Introduce `λγ` with provenance, **classify it** (`calibrated` / `derived` / `debt`), confront
`CONE_LOCK_CALIBRATED`, and ⭐ **apply the kind-test**: is this a medium-wide ratio, or a relation between
two sector speeds that happen to share dimensions? ⚠ That is exactly the shape of the ten identifications
in register §A.

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

⛔ **The table above is INCOMPLETE — both reviews said so independently.** It was built from the v2
R-number habit rather than from what phases 3–4 actually consume. ⇒ **Build it from the imported PN and
worldtube input lists**, and add at minimum:

| also owed | why |
|---|---|
| **`a`, the throat radius** | controls the first worldtube correction `O(a²/r²)`, and the finite-size/radiation channels |
| `M` (worldtube mass) and the supplied multipoles | supplied inputs of the S16 reduction, not outputs |
| `G` / `N_∞,3` normalization | the force magnitude is calibrated through it |
| `Θ_Q` and the `J → Q` map | S15 silently substitutes one for the other |
| `κ_add`, `κ_PV`, `κ_ρ`, optical `n = 5` | the PN response packet S17 imports |
| the **cone lock** `λγ` | S20a |
| `I_F`, `W_eff`, profile closure, source-vs-inertial mass, EP, Newton-`G` normalization | named in the force-derivation source |

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
**S9–S11, S16–S19** are where the earned results live. ⛔ **S14, S15 and S20 are CONDITIONAL ON S14a**, not earned.
**S0.5 and the substrate-action step are prerequisites** — ⛔ nothing is banked before them.
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

# v3 STEP PLAN — gravity, light, gravitomagnetism, charge, magnetism

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
   quantity — here is why, and here is what would make that wrong."* Thirteen pin-shaped identifications
   are already on record; this is the exact move that produces them, and it is the one place a second pair
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

⚠ **`R`-number namespace hazard — `R⟨n⟩` and `R1`-the-rung are DIFFERENT namespaces.** `R⟨n⟩` (`R1`,
`R3`, `R10`, `R60`–`R73`) are **parameter-register edges** (`parameter_register.md:268` —
`| R1 | c_s0 = √(5K ρ0⁴/m_GNLS) | DERIVED |`); `R1`-the-**rung** is the **one nonlinear throat solve**
(`parameter_register.md:329` — `R1_REQUIRED(bc_selection)`; `docs/model_map.md#shared-r1-throat-solve`).
⛔ `R1_REQUIRED(x)` means **blocked on the throat solve**, ⛔ **not** on the sound-speed edge `R1`. {#s05-r-namespace-hazard}

⇒ Executing S1–S3 against the seed as-is either **accepts the O-01 universe hole** or requires surgery
that no step names. Do the surgery **first**:

- `c_γ` and `λγ` are **deleted** from the registry until S9 / the cone-lock step introduces them with
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
affect `active_variables`. ⇒ Specify the lifecycle exactly: **delete `Q.medium.c_gamma`,
`Q.medium.lambda_gamma` and `R3` from the registry outright** (user decision, 2026-07-31), then
introduce `c_γ` at **S9** and `λγ`/`R3` at **S20a**, each with its own provenance.

⭐⭐ **The right word is not "retired" but PREMATURE.** `c_γ` is not dead — it was placed in the medium
block **before the walkthrough reached the step that introduces it**. ⇒ It is removed so that **S9** can
introduce it properly with its own provenance, and `λγ`/`R3` likewise at **S20a**. ⛔ "Retire" implies
dead; this is **not yet born**.

⭐ **Nothing dead goes into the registry, the scripts, or the eventual `.tex` that composes the ledger.**
A retired-but-present row is dead weight a reader must filter.

⭐ **The precedent is deletion.** Step ① deleted `Q.medium.a_pin` outright — *"`Q.medium.a_pin` are
**deleted**, and no class needs picking"* (`docs/derivation_walkthrough_plan.md#classification-test`).

⚠ This also **moots** any question about whether removed quantities keep their `source_loci` — there are
no removed-but-present rows.

⛔ **Enumerate every touchpoint before starting** — at minimum `registry_read.py`'s runnable propagation
smoke test (`:1014-1018` — the `evaluate_output("lambda_gamma", …)` call plus its `EVALUATION:` print) and
`README.md`'s canonical example (`:179-181` — the `assert … == 1` inside the block that opens at `:174`),
**both hardwired to `lambda_gamma`** *and* to `c_gamma` in the input dict,
plus `acceptance_check.py`, `able_to_fail.py` and `test_registry.py`.
⭐ A mutation case built on `R3` needs a replacement built on a surviving relation.

- ⭐⭐ **`C-M1` dies with `R3`** — discover this now, not mid-step. Its entire content is *"drop `R3`"*:
  `acceptance_check.py:54` builds it as `without_r3 + (lambda_gamma - lambda_gamma,)`, where
  `without_r3` (`:29`) is the baseline minus `R3`'s residual, and `constraint_dimension`
  (`registry_read.py:638-644`) discards zero-valued constraints — `if sp.simplify(expression) != 0` — so
  the appended term is **inert**. ⇒ Its replacement must be a **plausible physics error on a surviving
  relation** — `R1` (the sound speed) or `R2.xi_h`/`R2.h0` — ⛔ not a relabelling and ⛔ not another
  inert term. ⭐ **Choosing which mutation is a user call**, not the builder's.

  > ⚠ **Measured — and stated in the harness's own terms.** ⛔ **Do not restate this as coverage of the
  > method doc's five named error classes:** that mapping was attempted twice and was wrong both times, in
  > opposite directions. What was actually measured, on the post-S0.5 block:
  > - The **count** (`constraint_dimension`) moves only when the number of **independent relations**
  >   changes — dropping a whole relation gives `7 → 5`, adding a genuinely new one `7 → 3`.
  > - ⛔ It does **not** move for a wrong **coefficient**, a wrong **sign**, a wrong **branch**, a term
  >   **omitted from inside** a relation, or even a **changed symbol set** — every one leaves `7 → 4`.
  > - ⭐ The **certifier** (`certify_positive_real_dimension`, `acceptance_check.py:66`, run on every case)
  >   **raises** on a wrong **sign** and a wrong **branch** while the count sits at 4.
  >
  > ⇒ **The mutations reachable for `C-M1`:** drop a whole relation · flip a sign · take the wrong branch ·
  > substitute a **circular** relation the survivors already imply. ⛔ **Do not assert that any error class
  > is or is not covered without measuring it** — this passage was wrong twice for exactly that reason.

⭐ **What the count currently means.** The medium block's free count is exactly the number of declared
primitives — the three relations determine three intermediates and nothing else. ⛔ The counter has
**nothing to say yet**; it starts earning its keep when a step introduces quantities with no route back
to primitives (`{a_B, κ_B}` at S5/S6, `{ρ_br, μ_R}` at S8).
⛔ **And it could not have caught the `a`-pin** — a quantity introduced with its own defining relation and
consumed by nothing is **count-neutral**. → `DEFECT_REGISTER.md` row **F6**.

⭐ **And this whole step is bookkeeping repair, ⛔ not physics.** Keep the `C-M1` replacement to the
mechanical minimum: the mutation slot guards little today (the paragraph above) and is ⛔ not worth
optimising.

⛔ Derive the new payload independently; ⛔ never preserve, and ⛔ never read it off a prior document.

---

## PHASE 0 — the substrate (4 steps)

Re-banked from v2, forward, with the corrections this session established. ⭐ Cheap — but it is what
makes v3 readable without v2, which is acceptance criterion (3).

### S1 · The medium, and the brane as its ordered state
One GNLS condensate in a 4D bulk; `ρ = |ψ|²`; `P = Kρ⁵`; the brane is the **ordered state** of the same
medium (`χ_B = 1`), not a second object.
**Expected new:** 4 scalars `{ħ, m_GNLS, K, ρ0}` · 3 discrete `{D=4, n=5, the two-phase split}` ·
**2 fields** · 1 BC.
⛔ **Two, not three.** `ψ` plus the order field in whichever form S1's **A13 gate** selects — state the
real-vs-complex degree count explicitly. `U(ρ)` is an **EOS energy density, not a field**
(`stage004:56-69`); it was miscounted in v2.
**Register:** ⛔ **A13 — GATE. Resolve before this step banks.** `χ_B` is used here as a **real material
fraction** (so `χ_B·n` is real) while S5 writes it **complex**. Choose the branch — real/dissipative or
complex/inertial — propagate the choice into the field count, the action and the conversion balance, and
**remove the unchosen inputs**. ⛔ S1 cannot bank on "Register: none".
**Carry forward:** the two-phase split is `postulated` — the largest tier-1 item.
⚠ **Correct v2's stated reason.** `model_map.md#two-material-states` justifies the postulate by *"`U(ρ)` is single-well"*.
That reasoning is **wrong** — the brane comes from `V_χ(r_B)`, a potential on a *different* field. The
conclusion survives, the argument does not. → S5.

### S1.5 · The substrate action and the Madelung balances  ⚠ runs AFTER S1
The GNLS action, the quantum-gradient term, and the Madelung mass/momentum/energy balances.

⛔ **The pressure identity does NOT determine `U(ρ)`.** `P = ρU′ − U = Kρ⁵` has general solution
```
U(ρ) = Kρ⁵/4 + Cρ
```
The linear term is unconstrained by the EOS. `C` is a **chemical-potential / energy-reference choice**
and it enters the energy and source ledger. ⇒ **Choose it explicitly and classify the choice**; setting
`C = 0` is a decision to be recorded, ⛔ not a consequence to be asserted.

**Why this step exists:** S4 invokes a "core balance" and S12 needs momentum/energy partners, and
neither has an antecedent in S1's output. Without this step *"derive forward from what previous steps
banked"* is false at S4.

**Expected new:** the reference choice `C` (one decision), and an explicit `Π_Q` / `ε_Q` / `j_ε^Q`
improvement convention — ⭐ **the three are one package**, and picking them separately is how a
convention becomes an unnoticed input.
⚠ This step supplies the **conservative left-hand sides only**. S12's non-variational source partners
stay open.
⚠ Confirm S1's **two-field** inventory (`ψ` + the A13-selected order field); `U` is not a field.


### S2 · The sound speed
⛔ **Derive FORWARD from S1/S1.5. Do NOT re-bank v2's `stage005`** — it is wholly historical: it consumes
the retired pin relation and classifies `λγ` as a single underived calibration input, which would undo
S0.5 by reintroducing `λγ`. ⚠ Cite it for provenance only.

`c_s² = (1/m)dP/dρ = nKρ^(n-1)/m`. **Expected new: nothing** — pure consequence.
**Finding to re-bank:** `[K] = M L^(4n−2) T⁻²`, so `n` and `[K]` are one structural choice.
**Register:** none.

⭐⭐ **DECIDED 2026-08-01 — `K` and `n_eos` STAY in the medium block.** The question raised was whether
they should leave the way `c_γ` did at **S0.5**. They should not, and the reasoning matters more than
the verdict:

- ⛔ **S0.5's criterion was TWO criteria, and only one applies here.** **(a) Category** — `c_γ² = μ_R/ρ_br`
  is built from two **brane** properties, and the shear-free bulk cannot produce a shear modulus at all,
  so it was in the wrong block. **(b) Prematurity** — it sat in the block before the walkthrough reached
  the step that introduces it. ⭐ `c_γ` failed **both**; `K`/`n_eos` fail **only (b)**, and under
  requirements-first (b) indicts **all eight** medium rows equally, so it cannot discriminate.
- **Measured on the live registry** (⛔ not reasoned about): `K` appears in **exactly one** relation,
  `R1: c_s0 = √(5Kρ0⁴/m)`, whose designated output is `c_s0` and nothing else. `R2.h0` and `R2.xi_h`
  both reach it **through `c_s0`**. `R4` does not touch it. ⇒ **Every half-one consumer of `K` reaches it
  through `c_s0`.** And `n_eos` is `counting_axis: discrete-structural` — ⛔ **not** in
  `active_variables`, so dropping it moves **no number**; its only live job is `R1`'s
  `literal_consistency` guard tying the literal `5` and exponent `4` to `n=5`, `n−1=4`.
  ⇒ ⭐ **`K`, `n_eos` and `R1` stand or fall together.**
- ⭐⭐ **The deciding argument is HALF TWO.** Half one is linear response about `ρ0`, which only ever sees
  `dP/dρ` **at `ρ0`** — i.e. `c_s0` — never the **shape** of `P(ρ)`. ⛔ But half two needs a **parent
  action**, and the parent action contains `U(ρ) = Kρ⁵/4`. ⇒ The EOS **is requested — by half two.**
  ⚠ **Precision, so this is not overclaimed:** light's nonlinearity is a nonlinear **shear** term
  (`μ_R`'s extension), ⛔ not the bulk EOS directly; they connect through **`R17`,
  `μ_R = ∫ χ_B μ_R⁽⁴⁾ dw`** (`parameter_register.md:284`), which is **`PENDING`, "dim-consistency
  asserted only."** ⇒ The request is real but **routes through an open debt** — say so.
- ⚠ **The swap was NOT free, and the cost is S0.5-shaped:** `acceptance_check`'s **C-M3 mutation is
  `K − ρ0`** and dies with `K`; `able_to_fail`'s **`provenance` tooth names `R1`**. Same shape as
  *"C-M1 dies with `R3`"*. ⛔ The free **count** does not change either way (`{ħ,m,K,ρ0}` →
  `{ħ,m,c_s0,ρ0}` is still four with one dimensionless combination) — ⇒ **this was never a count
  question.**

⭐ **What the record must now say:** the medium block is a **v2 inheritance carried provisionally**;
the **knit** re-derives it. ⛔ No medium row reads as earned before then.

**Open:** **O-02** — is `K` + the EOS exponent **one entry or two**? ⚠ **Still open, and it is a
DIFFERENT question** from the one decided above; ⛔ do not read the decision as closing it.

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
derived-from-primitives). **B2** — ⛔ **state the no-go at its ACTUAL scope**: the GNLS-polar-smectic
gate program that would have grounded this wall is *"SUPERSEDED at the brane-existence level"*, which is
why the wall's constants stay **postulated** — ⛔ but `FAIL_COUPLE_STRESS_NOGO` itself is **narrower**
than that packaging: it is a route-failure for *deriving* `μ_R` from a polar substructure `P`, and
*"light stands on the bare postulated modulus"* (`DEFECT_REGISTER.md#B2`). ⇒ That half of **B2** is
`μ_R` provenance and belongs at **S8**, ⛔ not on the wall step.

### S6 · The kink — wall thickness and surface tension
Solve the EL equation: `δ = √(κ_B/(2a_B))`, `σ_wall = √(2a_Bκ_B)/6`. **Genuinely derived** — from S5's
postulates.
**Register:** **A2** (`ℓ = δ` is an *identification*, not a derivation — confront it here or nowhere).

⭐⭐ **INHERITED OBLIGATION FROM S11 — `B_comp` (`Q.brane.B_comp`) MUST BE RETIRED OR RE-AFFIRMED HERE.**
{#s6-b-comp-callback}
⚠ **Written at S6 on purpose. It was decided at S11, and a note living only in S11's record is a note
nobody reads when they arrive here.**

S11 entered the brane's **compression modulus `B_comp`** as a **postulated** knob, on the user's
explicit call (2026-08-02): *"pose B as a postulate knob, then we can claim it fully derived once we
reach that point in the program, and show it's no longer postulated. That way the ledger stays linear
and clear in what it claims."* ⇒ The knob count is an **upper bound that can only improve**, and the
retirement must be **visible** when it happens.

⭐ **What S11 predicts this step should find** — and it is able to fail:
in-plane compression is relieved through **two channels in series**, the medium's EOS and **widening of
the wall**, so compliances add: `1/B_comp = 1/B_EOS + 1/B_wall`. ⇒ `B_comp` is **softer than either
channel alone**, and its wall part is governed by **`σ_wall`, which THIS STEP derives**.
⛔ **If `B_comp` comes out of `σ_wall` plus the medium EOS with no new constant, the knob RETIRES** —
record the count moving, do not leave it postulated out of habit. ⛔ If it needs a genuinely new
constant, say so plainly; that is a result too.

⚠ **And S11's own cone is conditional on this step.** S11 computed `ω² = (B_comp/ρ_br)k²` with the wall
width **FROZEN**. Move 5 predicts that unfreezing it **softens the longitudinal at long wavelength**
(`ω ∝ k²`, flexural rather than acoustic, from the `σ_wall|∇W|²` cost of a modulated width).
⛔ **If the dispersion stays a clean cone once the width is dynamical, S11's move 5 was WRONG** — say so
rather than reconciling it. ⇒ `steps/S11_stray_longitudinal.md`.

### S7 · ⛔ The slab width is NOT selected
`W_slab` is `FREE-UNREDUCED`; *"double-well selects NO width"*. The kink gives one interface; the brane
is a finite slab.
⭐ **A step whose result is "this is not determined" is a real step.** Bank it as such.
**Register:** **C7**, **A9** (`W_slab` merged with the `L/a` debt without an equation).

⭐⭐ **S11 DEPENDS ON THIS STEP'S ANSWER, and reads it as a LIVE PHYSICAL QUESTION, ⛔ not bookkeeping.**
{#s7-b-comp-callback}
If the slab width is a genuine **flat direction**, then the wall offers **no resistance to thickening**,
`B_wall = 0`, and by the series law `B_comp = 0` — ⇒ the longitudinal zero of **S10 is never lifted** and
S11's propagating mode **does not exist**. ⚠ S11's answer is that the flatness is **lifted by
gradients**: a wave modulates the width, tilting the interfaces and stretching them at a cost
`∝ σ_wall|∇W|²`. ⇒ **flat at `k = 0`, stiff as `k²`.**
⛔ **So "not selected" is not the end of it — this step must say whether the flat direction survives
when the width is made position- and time-dependent.** A uniform-width statement does **not** answer
S11's question. ⇒ `steps/S11_stray_longitudinal.md`, moves 4–5.

### S8 · The quadratic brane Lagrangian ⭐ and `{ρ_br, μ_R}` enter HERE
Transverse and longitudinal sectors; the reduced `h`/`u_L` operator.
⛔ **Introduce and classify `{ρ_br, μ_R}` in this step** — the transverse action contains both *before*
`c_γ² = μ_R/ρ_br` can be derived, so introducing them at S9 would invert provenance.
⇒ **The R10 debt starts here**, not at S9.
⭐⭐ **The step that defines the LINEAR/NONLINEAR SEAM.** State plainly: **quadratic = linear response
about an assumed equilibrium.** Everything in **half one** is small-oscillation physics on a brane that
S5–S7 postulated. It is why gravity and light are tractable, and why the defect is not.

⛔⛔ **A SEAM, ⛔ NOT A CEILING** (user, 2026-08-01 — `CHARTER.md#two-halves`). ⛔ Do **not** write that
this ledger "stops" at linear response, and ⛔ do **not** record missing nonlinearity as a *blocker*.
Linear is **half one's scope**; nonlinear is **half two's subject** — the simulation-setup half. ⚠ The
wrong framing already cost a session: it recorded the absent nonlinear shear action as *"⛔ the blocker"*
on the photon-simulation track, when that action **is** what half two exists to specify.
**Defect register:** **C6** (no closed parent action — the coupling to sleeve/geon/drain does not exist) ·
**B2** (⛔ the *only* thing `FAIL_COUPLE_STRESS_NOGO` closes is deriving `μ_R` from a polar `P` — so
`μ_R` enters here **postulated** and stays postulated; ⛔ that is a missing reduction, **not** a
falsified foundation for light — `DEFECT_REGISTER.md#B2`).
**Parameter-register edges:** **R10 starts here** (`{ρ_br, μ_R}`).

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
**Defect register:** none.
**Parameter-register edges:** none — R10 began at S8. → S22.

### S10 · Two transverse photons
⛔⛔ **REWRITTEN 2026-08-09 — ⚠ the previous two lines here asserted what the closed record REFUTES**, and
they were still live while the paper card cited this very range. ⇒ ⭐ `steps/S10_two_transverse_photons.md`
is the authority; ⛔ this entry may not outrun it.

⭐ **What S10 measured:** the conditional map `D ↦ D − 1` transverse null directions **for the cases it
swept** (`D = 2,3,4,5`), so the `D = 3` member has two — ⛔ **conditional on the supplied action, the
supplied `[u]`, and BOTH structural premises.**

⛔ **What this entry used to claim, and why each is wrong:**
1. ⛔ *"brane shear gives exactly two transverse polarisations."* ⚠ **The form control refutes the
   attribution**: `FULLGRAD` returns the **same root and the same `D − 1` count**. ⭐ What curl-only
   determines is that the **longitudinal does not propagate** — a different claim about a different root.
   ⚠ And the count needs a **second** premise the slogan never mentions: **isotropic inertia**. `ANISO`
   keeps `S_curl`, satisfies every other condition, and drops the **exactly-transverse** count to `D − 2`.
2. ⛔ *"A genuinely earned step — say so."* ⚠ The record says of itself that it does **NOT** establish a
   viable physical light sector — only a conditional mode count — with **seven OPEN substrate obligations**
   ⇒ `SUBSTRATE_REQUIREMENTS.md`. ⚠ And the algebra is **standard**: MacCullagh 1839
   ⇒ `docs/medium_requirements_and_prior_art.md`. ⛔ **Expected new: nothing** was right; ⛔ *"earned"* was
   the wrong word for why.

⚠ **`D = 3` is NOT selected here** ⇒ `R-S1-01`, target **S6**, status OPEN.

### S11 · The stray-longitudinal departure
Characterized, first-class. ⛔ Not a defect to fix — a recorded departure from the reference theory.

⛔⛔ **The token `FAIL_CAUCHY_STRAY_LONGITUDINAL` reads as a defect. It is not.** ⚠ Record the **naming
defect** explicitly: the token is spelled `FAIL_…`, and that spelling already caused a reader this
session to treat a **required feature** as a defect. ⛔ **Do not rename the token** — it is committed in
engines and reports — ⭐ but flag at **every point of use** that the "failure" names a **characterized
departure the ledger keeps as first-class** — ⛔ not a bug to erase. ⚠ ⛔ Do **not** upgrade that into
*"the charge sector requires it"*: it does not (see the DIFFERENT OBJECTS block below).

⭐ **The computed result: one medium carries both mode families consistently**
(`software/stage1_solver/decisions/15_em_medium_native_physical_picture.md:266-272`):

> *"Spectrum of the symmetric law is `(λ−μ_br k²)²(λ−(K_br+4μ_br/3)k²)` … **μ_br>0** (ordinary/Cauchy
> symmetric shear …) → two transverse modes `μ_br k²` **but also** a longitudinal mode
> `(K_br+4μ_br/3)k²` = a stray **"second photon"** → `FAIL_CAUCHY_STRAY_LONGITUDINAL` at Stage 4."*

⇒ ⛔ **The obstruction is to SUPPRESSING the longitudinal mode** — i.e. to getting **transverse-only**
Maxwell — ⛔ **never** to carrying both. Carrying two transverse modes **and** a longitudinal one
simultaneously is exactly what the `μ_br > 0` branch *does*.

⭐⭐ **The mode is REAL and IRREMOVABLE, and its role was CATALYTIC — ⛔ not constitutive** (user's own
history, recalled 2026-07-31). ⛔⛔ **It is not charge; it does not become charge; it supplies no part of
charge.** ⚠ An earlier version of this step said it "becomes the **±w** displacement" — ⛔ that was a
false identification, corrected here. The history, in the order it happened:

1. A **gauge field** had been used to obtain exact Maxwell, and had been in use for **months**. The user
   pressed on **what it meant physically**, and there was no compelling answer.
2. The then-current view of charge was a **4D interaction seen as a 3D shadow**. ⛔ Those calculations
   gave **`1/r³`** — the wrong power (falsified-route box below).
3. The calculations were redone **from the core concept of the model**, to see whether the two transverse
   modes would come out. ⭐ **They did — but an extra scalar came with them, and could not be removed.**
   ⇒ Exact Maxwell was **not going to happen**.
4. ⭐⭐ **Because the scalar was irremovable, it became clear the charge picture itself was wrong.** The
   replacement: charge as the **tugging of the brane into the bulk (`w`) direction** — which gives
   **`1/r²`**, sits **on the brane**, and does **not disturb the physics already built**.
5. ⇒ *"The fact that it popped out unlocked a view of the conceptual model that didn't exist before."*

⭐⭐⭐ **AND HERE IS THE VIEW IT UNLOCKED — the user's own words, 2026-08-02. ⛔ Record it at every point
of use; it has had to be explained THREE times:**

> **The brane is like a DRUM HEAD. The defect is like your FINGER pushing into it.** The two different
> directions you can push are **the two charges**. Pushing creates an **energy potential** in that
> direction. **Same charge → energy INCREASES → they push away. Opposite → brings them together.**

⇒ ⭐ **Coulomb's SIGN is ENERGY GEOMETRY** — two same-`w` deflections superpose and raise the elastic
energy, so the configuration relaxes by separating; opposite-`w` deflections cancel, so it relaxes by
approaching. ⛔ **No gauge field anywhere in that account.**

⚠⚠ **AND THE GAUGE FIELD IS THE POINT OF THE WHOLE STORY.** The model carried *"exact Maxwell"* for
**months, across multiple papers**, obtained by **INSERTING a gauge field** — ⚠ uncaught at the time. The
question that broke it open was ⭐ ***"what physically motivates this gauge field?"*** Answer: nothing. It
was inserted **in order to get Maxwell**. ⇒ ⭐⭐ **Generalise the diagnostic: of any object a derivation
leans on, ask whether it is a THING or MACHINERY placed there to reach a known answer.**

⇒ ⭐⭐ **What the mode did was falsify the TARGET.** Exact Maxwell gives two transverse photons **and
nothing else**, with charge entering as an **external source put in by hand** — so matching it exactly
would have **inherited that gap**. ⛔ The point is that **the target was wrong**, ⛔ *not* that the
longitudinal mode is the thing that carries charge. Its irremovability **closed the exact-Maxwell route
and forced the rethink**; the rethink produced the `±w` mechanism (PHASE 4b, `Q2`).

⛔⛔ **DIFFERENT OBJECTS — ⛔ do not weld them, at any strength.** `u_L` is an **in-plane longitudinal
brane displacement**; **`±w`** is the **throat's normal/orientation direction**, and the charge the
rethink produced is carried by the **`h`-branon**. The committed sources keep them apart:

- `research/pde_ledger_v2/notes/stages/ledger_stage031_puncture_deflection_field_identity_source.md#h-a-distinct-from-u-l`
  — *"Define the held mouth datum `h_A = ξ_w|_A/ℓ = P₀H|_A`, `[h_A]=1` **(distinct from `u_L`)**"*.
- `docs/em_analog_next_phase_handoff.md#u-l-clamp-not-the-charge-scalar` — of the `u_L`-clamp candidate:
  *"It is NOT the committed charge scalar (that's the `h`-branon, `h≠u_L`; `h` remains the committed
  mediator for the conditional `1/R²` falloff); `u_L` is a separate charge-odd density mode BC'd to
  relax to 0."*
- `docs/model_map.md#superseded-charge-routes` — the *"leftover-scalar `u_L`-clamp"* is listed among the
  **superseded** charge routes.

⇒ ⛔ **Never write that the stray longitudinal *is*, *becomes*, or *supplies* the `±w` displacement.**

⚠⚠ **FALSIFIED ROUTE — the wrong-power signature** (user, recalling the history 2026-07-31).
Charge-as-a-4D-interaction-seen-as-a-3D-shadow gives **`1/r³`**; the brane-tug (`±w` puncture) view gives
**`1/r²`**. ⛔⛔ **This is RECALLED HISTORY, ⛔ not a located computation.** The *power-counting statement*
is in the corpus — `docs/conceptual_foundation.md:236` (*"a source acting through the **4D bulk** naively
gives the **wrong** falloff (4-space → `1/R²` potential → `1/R³` force)"*) and
`docs/em_sector_reconsideration.md:62` (*"**CORRECTED (was backwards).** A finite bulk source gives
`1/r³` at large `r`, not `1/r²`"*) — but the **calculation that landed `1/r³` for charge is
`UNLOCATED IN CORPUS`**; ⛔ do not cite one, and ⛔ do not manufacture a locus for it.

⭐⭐ **A SECOND DIAGNOSTIC QUESTION (register idiom).** A **mathematical device can survive for months
without a physical referent** — the gauge field did, until it was asked what it *meant* physically
(step 1 above). ⚠ This is a **different question** from A7's KIND test
(`DEFECT_REGISTER.md#a7-kind-test`): the KIND test asks ***"is this the same KIND of thing?"*** and
**presumes both are things**; this one asks ***"is this a thing at all, or only machinery?"***
⇒ Ask it of every object a derivation leans on — ⛔ not only of the two being compared.

⚠ **And the Maxwell locus IS reachable** — `pathA_36` emits `C5_RESOLVED_MAXWELL_BY_TUNING`
(`research/pde_ledger_v2/paper/stages/stage_003.tex:25-28`):

> *"the tuned Maxwell locus (`K_θ=C_J²/ρ_br`, `B_eff=0`, `m_θ²=0`) is reachable and emits
> `C5_RESOLVED_MAXWELL_BY_TUNING` — so the FAIL is not hardcoded — but `BY_TUNING`, not
> `WITH_PROVENANCE`."*

⇒ The departure is a **choice the model makes**, ⛔ not an inability.

### S11b · ⭐⭐ The brane–bulk interface coupling — LINEAR, and it closes S11 {#s11b}

**Added 2026-08-02, after S11 walked. ⭐ It is not scope creep, and the test is one line: it is LINEAR,
and by this ledger's own definition linear is HALF ONE** (`CHARTER.md#two-halves`). It was deferred by
**choice**, ⛔ not by difficulty.

⭐⭐ **Why it is worth its own step: S11 ends with three apparently separate open questions that are ONE
question**, and this is that question:

| S11 leaves open | reduces to |
|---|---|
| does the longitudinal radiate into the bulk, or stay bound? | **the coupling law** |
| is light's confinement **unconditional**, or does it rest on a polarization-overlap argument? | **the coupling law** |
| does a second characteristic speed break Lorentz invariance **for us**? | **the coupling law** |

⇒ ⭐ **The second mode's entire physical status — radiative, observable, Lorentz-breaking, or none of
the above — is decided here.**

**What it needs:** interface conditions at the wall (continuity of displacement and of stress), the
scalar bulk sound mode already in the registry, and the brane sector from S8–S11. ⛔ It does **not**
require the wall solved — standard elastic/acoustic matching is enough for the conditions themselves.

⛔ **What S11 established that this step must NOT re-assume:** phase matching is **kinematic only**. It
settles whether a propagating bulk channel **exists**; it settles ⛔ **neither** that an existing channel
is used, ⛔ **nor** that an evanescent solution forms a bound eigenmode. Both are this step's job.

⚠ **Expect the transverse answer to be strong and the longitudinal answer to be weaker** — the bulk has
no transverse mode to match at all, whereas the parallel pairing is like-for-like in polarization.
⛔ But S11's own engines refused to conclude even the transverse case without this law; ⛔ do not import
that refusal as an answer in either direction.

⭐ **What it unlocks beyond the ledger** (→ **S22**): three of the five near-term *linear* simulations —
longitudinal radiation, the width mode's inertia and stiffness, and the flexural crossover of move 5.
⇒ `steps/S11_stray_longitudinal.md`.

**Expected new:** interface/matching conditions; ⛔ no new medium constant is anticipated — if one is
required, that is a result, record it.
**Register:** **C13** is adjacent but ⛔ **not** this step's job (a gravitational wave is not a brane
mode).

#### ⭐⭐ S11b IS SPLIT INTO THREE — ⛔ AND S11b IS **NOT CLOSED** UNTIL C IS {#s11b-split}

> ⭐⭐ **UNIFIED REWRITE UNDERWAY (2026-08-19) — this section's A/B/C framing is SUPERSEDED.** A+B are being
> rebuilt as **ONE unified export-chain step** (decision list `directives/S11b_unified_decisions.md`
> `ddd0ae4c`; unified spec `directives/S11b_SHARED_PHYSICS.md` `1a2395a3`, two legs Codex+Grok folded once).
> ⭐ **The immediate NEXT is the two per-engine build directives** (each reviewed before its build); the two
> blind engines, comparator, step record and card follow. ⚠ The **"✅ CLOSED"** tags below are the **OLD
> pre-chain** sense — the physics stands and is being re-packaged, but the ledger rows **do not exist yet**.
> ⛔ **SUPERSEDED below:** *"C runs immediately after B's rebuild"*, *"C is not a separate step"*, and the
> table's **"S11b-C ▶ NEXT"** — per G1, **C is a separate LATER step** (a later build **after** the unified
> A+B build closes; ⛔ not immediately after B, and ⛔ not the immediate NEXT); no renumbering (still
> S11b-C). ⇒ live NEXT in `STATUS.md`'s top block and memory `project-s11b-interface-law-result`.

⛔⛔ **CORRECTED 2026-08-05 (user).** ⚠ Several docs — including a memory titled *"S11b CLOSED"* — read as
though S11b were finished. ⛔ **It is not.** The split below was for **specification tractability**: three
attempts to spec the whole interface in one pass were rejected. ⇒ ⭐ **A, B and C are ONE STEP**, and
⛔ **S11b closes only when C closes.** ⚠ A ledger that calls a step closed when a third of it was never
built overstates what the light sector establishes.
⇒ ⭐ **C runs immediately after B's rebuild** — ⛔ it is a FIRST BUILD, not a rebuild, since C never existed.
⭐ ⛔ **No renumbering:** C is genuinely part of S11b, ⛔ not a separate step.



⚠ **2026-08-03.** Two attempts to specify the whole interface in one pass were rejected before any build,
and a third was rejected for mandating a non-uniform background while fixing plane waves. ⇒ three steps:

| | scope | state |
|---|---|---|
| **S11b-A** | the bulk's response to moving faces + the projection identity | ✅ **CLOSED** — `steps/S11bA_interface_response.md` |
| **S11b-B** | the **homogeneous** assembly ⇒ does the longitudinal radiate or stay bound | ✅ **CLOSED** — `steps/S11bB_interface_assembly.md`. ⚠ **A's velocity leak lies OUTSIDE the passive region ⇒ it costs a named reservoir; ⛔ it is NOT forbidden.** Passivity + Onsager are classifications, ⛔ not gates, and both assume an equilibrium reference state this model does not have (`v₀ ≠ 0`) |
| **S11b-C** | the **non-uniform** transverse coupling ⇒ is light's confinement unconditional | ▶ **NEXT.** ⚠ B proved the uniform coupling is **identically zero**, so a uniform-limit control is now **known-vacuous** |

⭐ **The seam is real:** the longitudinal mode's fate needs no gradients; light's confinement needs them.
⇒ Full state, traps, and C's requirements: **`steps/S11b_HANDOFF.md`**.

#### ⭐ `S_leak` EXISTS ALREADY, and it is an IDENTITY — ⛔ not an ansatz {#s11b-sleak}

⚠ **Found by corpus search during the S11b walk, 2026-08-03**, and it is the mass-channel expression of
this very step's coupling. Projecting 4D continuity onto the slab against a window `W(w)` leaves an
integration-by-parts remainder — `research/pde_audit/scripts/stage_v2_06_07_poisson_newtonian_sympy_audit.py:53-62`.
⭐ **Nobody set out to derive a leak; it is what is left over because the window has edges.**

⛔⛔ **DO NOT IMPORT IT.** ⚠ **User, 2026-08-03:** *"We did a bunch of work on pde_audit. It was completed
at the time, but I'm unsure how solid it is since this model is always evolving in minor ways."*
⇒ It is a **lead, not an authority**. Re-derive it in **both** engines; it is cheap, and blind dual-engine
is exactly what settles whether it still holds. ⭐ That also means S11b carries **no dependency** on that
tree's solidity.

⚠ **It is adjacent to a RULED-OUT object and the distinction is load-bearing.** S12's `S_drain` mass-sink
was ruled out as **frozen-wall-premised**. `S_leak` is an *identity*, so the ruling does not obviously
reach it — ⛔ but that has **not** been verified, and it must be, not assumed.

#### ⭐⭐ THE USER'S DARK-ENERGY POSTULATE — banked 2026-08-03, ⛔ do not lose it {#dark-energy-postulate}

> *"this mechanism I am postulating is the reason for dark energy. As the bulk reorders onto the brane,
> the brane expands to accommodate. Causing expansion. But that's for later."*

⇒ `S_leak` carries a **steady DC part**: the bulk continuously re-orders onto the brane, and the brane
expands to accommodate it. ⭐ **Re-ordering, ⛔ not material moving between two substances** — brane and
bulk are one medium in two phases, so the leak is the order parameter converting at the wall
(user-confirmed 2026-08-03; likely the same physics as S12's `Γ_B`, ⚠ unverified).

⚠ **It has an IMMEDIATE technical consequence for S11b, so it is not purely "for later":** the wave is a
perturbation on a background that is **not static**. ⭐ Separable, and generously — a leak timescale of
order the Hubble time (`~4×10¹⁷ s`) against an optical period (`~2×10⁻¹⁵ s`) is **~32 orders**. ⇒ hold the
background fixed over a wave period and take `δS_leak` as the coupling, ⭐ **but record
`wave period ≪ leak timescale` as a stated validity condition** rather than assuming it silently.

⭐ **Why it is worth keeping past this step:** it ties the DC leak rate to an **observable** (the expansion
rate), which is a future calibration hook — ⛔ not claimed, not derived, and out of scope here.

---

## PHASE 3 — the drain (3 steps) ⚠ where the interior is deferred

### S12 · The throat's non-variational SOURCE, plus its boundary data
⛔ **Not merely a boundary condition.** `n*Gamma_B` is the RHS of a **local order-balance PDE** — a *source*. The return controllers and the mouth/collar/IR conditions are **separate** objects.
⇒ Inventory the source functions and the boundary data **separately**; conflating them hides which of the two is open.
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
⛔ **Two inventories, kept apart:**
- **non-variational source / controller functions** — `Γ_return`, `Γ_drain`, `J_χ`, and the return
  controllers. ⛔ All `[OPEN]` in form.
- **boundary / domain data** — the mouth, collar and IR conditions, and the `𝔅` branch selection.

⇒ The throat enters this ledger with its **source form open and its strength supplied**. ⛔ Not "only a
boundary condition" — that framing hides which of the two is unresolved.
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

### S15 · Two drains attract ⚠ RUN S16 FIRST for regime context — ⛔ NOT as S15's premise
⛔ **CONDITIONAL ON S14a.** Not a carried result until the dynamical-`Γ_B` bridge succeeds.
⚠ **Corrected twice.** Round 1: run S16 first, since the corpus establishes the worldtube closure before
the particle law. ⛔ **Round 2 corrects that correction: S16 does NOT license S15.** S16 is
*response-side* (how a body moves in a field); S15 is *source-side* (what force two drains exert) and
carries **its own** Noether / control-surface assumptions. ⇒ Run S16 first for regime context, but
⛔ derive S15's premises independently — do not inherit them.

Noether stress at infinity, drain strengths `Q_i` **given**: `F_12 = −I_F,12 · (m N_∞,3 Q1Q2 / 4πr²) r̂`  ⛔ **`I_F,12` is not 1** — it is listed as debt in S22; ⛔ do not display the coefficient-one form.

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

⚠ **Read the actual statement, not the slogan** — ⛔ and note the transcription caveat below (`research/4d_2_5pn/paper/4d_2_5pn.tex:605-614`):

> *"For a **compact defect of size `a ≪ r`** in a smooth external field, the center-of-mass reduction
> gives the usual point-particle force law … `+ O((a²/r²)·Φ_ext/r)` … so that the universal
> point-particle source **at leading order** is controlled by the worldline/worldtube multipoles rather
> than by arbitrary internal details of the defect."*

⇒ **The independence is (i) leading-order only and (ii) conditional on compactness.** The first
correction is `O(a²/r²)` — controlled by `a_WT`, ⛔ **not** the mouth radius (**A11**); `a_WT` is itself undetermined. So "gravity does not depend on the interior" is precise only to leading order; beyond it,
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

⛔ **QUOTATION INTEGRITY — the block above is the source AS WRITTEN.** `4d_2_5pn.tex:605-614` writes the remainder **without** the mass factor,
which makes it dimensionally inhomogeneous (force on the left, acceleration in the remainder). ⛔ **Do
not present the repaired formula as a quotation.** Quote the source as written, record the missing-mass
**typo**, and cite the dimensionally correct form from `4d_1pn_full.tex:886`.

⛔ **A11 GATE — this step must not substitute one length for another.** The source's `a` is a
**Gaussian/profile support width** controlling `Q_ij = M a² δ_ij / 2`; the model's **mouth radius** is a
different object, and the bridge report says outright it *"is not an invariant reduction width"*.
⇒ **Use distinct symbols — `a_WT` (worldtube profile width) and `a_mouth`** — and list the **missing
bridge between them** as debt in S22. ⛔ Never substitute one for the other. Same dimension is not same
object.

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

## PHASE 4b — charge and magnetism (7 steps) {#phase4b-split-register-fields}

⚠ **Dependencies here are PER STEP**, ⛔ **not phase-wide** — it is **not** true that all of Q1–Q7 wait on
PHASE 3's drain:

- **Q1–Q5** do **not** require the drain bridge. Stage 030 folds **ONLY** the charge-new electric-scalar
  subset and leaves the drain predicates (`PASS_DRAIN_KERNEL_NORMALIZATIONS`, the drain mass / momentum /
  energy controllers) explicitly **out of scope**
  (`ledger_stage030_electric_scalar_localized_h_closure.md:242`).
- **Q6** needs **one** earlier thing: PHASE 2's **transverse-vector row** (imported and cited, ⛔ not
  re-earned — `parameter_register.md:233`).
  ⛔ **`τ_d` is NOT a PHASE 3 inheritance.** The active-drain time arrow is **Q6/stage-034's own
  POSTULATED input**: 034's magnetism-NEW content is *"exactly `{moving coupling, q_T, τ_d}`"*
  (`ledger_stage034_transverse_move_action_row.md:38`, `:189`), and the register carries `τ_d` as
  *"structural / postulated (**NOT** a reducible knob)"* with ⛔ **no reduction route**
  (`parameter_register.md:232`, edge **R68**). ⚠ The `O(V)` moving row **requires** it — *"a passive
  T-even throat would not supply it"* (`parameter_register.md:232`) — ⛔ but **requiring is not
  deriving**, and PHASE 3 supplies no time arrow.
- **Q7** checks coexistence against the **assembled** substrate's row list (gravity's scalar / drain /
  return / wall rows among them, byte-for-byte unchanged — `parameter_register.md:234`), so it runs last
  in the phase. ⚠ That is the row **inventory**, ⛔ not the S14a drain **bridge**.

⛔ **The whole phase is independent of the PN ladder.**

⛔ **It does NOT sit next to light in PHASE 2.** Light is brane shear and is genuinely insulated from the
interior; **charge is throat-coupled**. ⇒ Placing it beside light would import an insulation this sector
does not have.

⭐ **Its far-field FORM is earned target-blind, while its SIGN and MAGNITUDE route through the interior
solve.** Hold those two apart in every step below; they have different classes and different blockers.

⚠ **The `Q` namespace is deliberate.** `S14a` already runs *before* `S14`, and `S20a` sits in a
*different phase* from `S20`; more `S`-suffixes would compound a confusing convention.

### ⛔⛔ READ BEFORE STARTING THIS PHASE — ⚠ prior art says LINEAR CANNOT CARRY CHARGE

⚠ **Raised 2026-08-09. ⛔ Not a blocker, ⛔ not a verdict — an OPEN QUESTION this phase must answer early,
because the answer decides whether the charter's linear scope survives contact with charge.**

⭐ **Unzicker**, *"Nonlinear continuum mechanics with defects resembles electrodynamics — A comeback of the
aether?"*, *ZAMM* `10.1002/zamm.202100280`. ⇒ detail and citation corrections in
`docs/medium_requirements_and_prior_art.md`.

⛔⛔ **HIS STATED REASON FOR LEAVING LINEAR ELASTICITY:** *"MacCullagh's theory based on linear elasticity
**cannot describe charges**"* — so he goes to a **twist disclination**, a topological defect *"causing
large deformations requiring nonlinear description."*

⚠ **Why this reaches us specifically:**
- ⭐ **Our charge is a throat / puncture — a LARGE deformation of the brane**, ⛔ not a small oscillation
  about a flat state ⇒ [[project-puncture-deflection-charge-mechanism]].
- ⛔ **And our brane–bulk interface law (S11b-A/B) is LINEAR**, and it is what the charge story runs on.
- ⇒ ⭐ the **object** is plausibly nonlinear while the **instrument** we built for it is linear.

⭐⭐ **THE DECISIVE TEST IS CHEAP AND ALREADY NAMED: ⛔ S11b-C, non-uniform coupling, WAS NEVER BUILT** —
and it is the **MacCullagh differentiator**. ⇒ ⭐ **Run it before committing this phase to anything.**
- ⭐ If S11b-C anchors charge linearly, we have what he did not, and half one's linear scope holds.
- ⛔ If it cannot, that is **independent convergence on his conclusion**, and **C6** (the absent nonlinear
  brane-shear action) moves onto the critical path **with a measured reason**, ⛔ not an assumption.

⛔⛔ **AND "GO NONLINEAR" IS NOT A MATTER OF KEEPING MORE TERMS.** ⚠ **C6**: every form in the corpus is
quadratic and the quadratic Lagrangian was written **directly**, ⛔ not expanded from a nonlinear parent
⇒ it ⛔ **cannot be reached by un-linearising.** Authoring a nonlinear parent is a **modelling decision**,
⛔ not a derivation ⇒ see the half-two section below.

⛔ **TWO THINGS THIS NOTE DOES NOT SAY:**
1. ⛔⛔ **Prior art is an ORACLE, ⛔ NEVER a premise** (`CLAUDE.md` rule 16). ⚠ His conditions are **not**
   ours: ⭐ **his medium is a 3D incompressible elastic continuum with ⛔ NO bulk, ⛔ no brane, ⛔ no
   codimension.** ⇒ the confinement architecture is genuinely not his; ⛔ only the charge-from-defect move
   is. ⚠ **The paper is PAYWALLED and has NOT been read** — this is abstract-level.
2. ⛔⛔ **Nonlinearity does NOT discharge `R-S8-04`.** ⚠ The couple-stress / internal-angular-momentum
   objection is a **separate axis**: Cosserat buys **microrotation**, ⛔ not amplitude. ⇒ a fully nonlinear
   theory can still be inadmissible for the same 19th-century reason ⇒ `SUBSTRATE_REQUIREMENTS.md`.

### Q1 · The electric-scalar substrate
The **static electric scalar**, closed by a **localized-`H` / PT** construction. ⚠ This is the phase's
entry point onto PHASE 3's throat, ⛔ not onto PHASE 2's brane-shear apparatus.

**Locus:** `research/pde_ledger_v2/notes/stages/ledger_stage030_electric_scalar_localized_h_closure.md:15`.

The reduction

```
M_h = N₀M₄ ,    K_h = N₀K₄ = M_h c_E²
```

is **EARNED — *given the postulated action***.

**Expected new:** the localized-`H`/PT closure, and the **throat Green speed `c_E`**.
**Class:** `M_h`, `K_h` — **derived, conditional on the action**; `c_E` — an **interior** quantity
(`parameter_register.md:135`).
**Regime:** static.
**Carry forward:** ⛔ **the action is postulated — a tier-1 item.** ⇒ "EARNED" here means earned *from a
postulate*, not from primitives. ⛔ Do not promote it in transcription. → S22.
**Defect register:** **C6** (*"No closed parent action."*) — ⛔ Q1 does **not** resolve it: `M_h`/`K_h`
are earned *given* a postulated action, which is exactly the closure gap C6 names. ⇒ **Explicitly
deferred**, not closed. → S22.
**Parameter-register edges:** none.

### Q2 · The puncture deflection and its source ⛔ the holder is a DEBT, not a result
A **±w puncture geometrically bends the brane into ±w**: the field identity `ξ_w = ℓh` and the
orientation-odd mouth source. Token: **`THROAT_H_SOURCE_1_OVER_R2`**.

**Loci:** `ledger_stage031_puncture_deflection_field_identity_source.md:12`, `:18`, `:25`.

**Expected new:** the **source identity** and the **far-field FORM** — ⛔ **not a holder.**
**Class:** the **`1/R²` falloff** and the **`s₁s₂` product** are **target-blind EARNED** (`stage031:25`) —
⭐ this is the **FORM** the phase preamble promises, and it is the part that does *not* wait on the
interior. ⚠ EARNED ***within Q1's postulated G0 closure*** (`stage031:18`), ⛔ not from primitives. {#q2-earned-within-g0}
**Regime:** static, one puncture against another.

**Open — ⛔ the core holder is carried as an R1 DEBT:**
⛔ **`NONZERO_HA_REQUIRES_CORE_HOLDER`** — a nonzero signed `h_A` **is not a stationary point of the
core-less exterior energy**, and ⛔ **031 does NOT resolve it**: the guaranteed-nonzero amplitude needs the
sim-deferred core holder / nonlinear throat solve (`stage031:170`, `:173`; `parameter_register.md:327`
R61 — a sibling of R10/R30/R33).
⛔ **`h_A` is a held mouth *datum*** (`parameter_register.md:204`), and ⛔ **a held boundary datum is not a
physical holder** — ⛔ do not write "the deflection is held by the puncture" and count a holder as earned.
⚠ `z_g > 0` is **POSTULATED** (`stage031:249`).

**⭐ Candidate holder mechanisms — RECORDED, ⛔ NOT SETTLED.** ⛔ Do not pick a winner in this step; which
agent holds the throat open is exactly what the deferred nonlinear throat solve is for.
1. **Trapped standing-wave pressure.** The equilibrium throat radius `R*` is set by **outward**
   trapped-wave pressure against **inward** brane TENSION plus ground-state superfluid BACKPRESSURE
   (`docs/conceptual_foundation.md:483`); the trapped mode is a **brane-shear standing wave** with
   Wheeler-**geon** ancestry
   (`software/stage1_solver/decisions/15_em_medium_native_physical_picture.md:304-306`), whose *"only job
   is **structural** — it holds the throat **open**"*
   (`software/stage1_solver/directives/pathA_27_drain_sector.md:34`).
2. ⚠ **CONTESTED — the drain-flow reading.** A later document says instead that *"what holds the throat
   open is the drain-flow/amplitude balance — a sim question, nothing static"*
   (`research/pde_ledger_v2/notes/stage044_v2_unfreeze_prep.md:27-28`), and the directive leaves the agent
   **disjunctive and unsolved**: *"the **holding-open** agent (geon standing-wave pressure **and/or**
   drain flow)"* (`software/stage1_solver/directives/pathA_24_brane_existence_defect_structure.md:307`).
3. ⚠ **The one COMPUTED result does not settle it either.**
   `FLUID_ONLY_COLLAPSE_NO_INTERIOR_STATIONARY` — no wave term ⇒ no stationary throat
   (`software/stage1_solver/reports/pathA_26_derrick.md:49`) — was run with a **declared scalar
   surrogate**, ⛔ **not** a brane-shear mode (`docs/conceptual_history.md:132`).

**Defect register:** **C4** and **C1** — ⛔ both **explicitly DEFERRED**, ⛔ neither resolved. Candidate
mechanism 1 invokes the **geon**, whose *"profile is a **declared OPEN input**"* (**C4**); mechanism 3
rests on `pathA_26_derrick.md`, the same report **C1** records as `THROAT_DRAIN_DESTABILIZED`. ⇒ Q2
**names** them; ⛔ it discharges neither. → S22.
**Parameter-register edges:** **R61**.

### Q3 · The BC ensemble and the electric sign ⛔ the result is *"this is not determined"*
Four coefficients `A_V / A_J / A_M / A_MIXED`, **DERIVED decided-conditional**
(`parameter_register.md:227`, edge **R62** `:328`).

⛔⛔ **The sign is NEITHER earned NOR calibrated.** It is **`outcome_not_invariant`** across the four
branches and carries **`R1_REQUIRED(bc_selection)`** (`parameter_register.md:329` R63). The selection
data `{φ, q, j, λ}` is **R1-deferred** — *"tracked, **not counted-clean**"* (`:228`).

**Expected new:** the four-branch BC ensemble; ⛔ **no sign.**
**Class:** the coefficients — `derived` **decided-conditional**; the **sign** — `R1_REQUIRED`.
⭐⭐ **This step's entire result is *"this is not determined"* — ⛔ that is the content, not a failure.**
Bank it as such.
**Defect register:** none.
**Parameter-register edges:** **R62**, **R63**. → S22.

### Q4 · Charge magnitude ⛔ the result is *"this is not determined"*
`Q_E` is **`R1_REQUIRED(magnitude)`**, and the core normalizations `c_a`, `c_ξ` are *"unbounded at
tier-A"* (`parameter_register.md:143`, edge **R65** `:331`;
`software/em_charge_attribute/puncture_deflection_electric_sign_result.md:51`).

**What it waits on — ⛔ THREE DIFFERENT CLASSES, not one bucket.** ⛔ Give each item the method doc's
class (`docs/derivation_walkthrough_plan.md#classification-test`); ⛔ the shared `FREE-UNREDUCED` register label is **not**
a class, and one of these is not an interior debt at all:
- **`a`** (throat mouth radius) — **`postulated`**: `FREE-UNREDUCED`, ⛔ **no defining equation**, value
  undetermined, and ⛔ **no named reduction route** (`parameter_register.md:132`). ⚠ **Count `a` only** —
  `c_a` is then just `a/r_e`, ⛔ not a second independent supplied quantity.
- **`Q_E`** (charge magnitude) — **`debt`**: the *same* `FREE-UNREDUCED` label as `a`, ⛔ **but a
  different class** — `Q_E` *has* a named route, `R1_REQUIRED(magnitude)` under the deferred throat
  solve (`:143`, edge **R65**). ⚠ `c_a`, `c_ξ` ride **inside** it, *"unbounded at tier-A"* (`:143`) —
  ⛔ not a class of their own.
- **`ℓ`** (embedding / PT healing length) — **`calibrated`**: ⛔ **`IMPOSED`/`CALIB`**, ⛔ **not R65** and
  ⛔ **not a registered throat-solve output**: the RATIO `ℓ/a = 1/20` is a **frozen handoff scale — a
  tuning, not derived** (`parameter_register.md:200`). ⇒ A **calibrated input** to this step, with
  ⛔ **no reduction route registered**.

⚠ **A LIVE TENSION on `a` — ⛔ state it, do not hide it.** `parameter_register.md:132` says `a` has **no
defining equation**, while `software/em_charge_attribute/puncture_deflection_electric_sign_result.md:51`
writes **`a = c_a r_e`**. ⇒ That is ⛔ **not** a defining equation under the method doc's test 1: `r_e =
k_e e²/(m_e c²)` is an **external empirical** scale, ⛔ not a model quantity, and *"neither `c_a` nor
`c_ξ` is fixed by the reduced action"*. ⇒ It **parameterizes** `a` against an outside length; `a` stays
**postulated**. ⛔ Do not read `a = c_a r_e` as a reduction.

⭐ **Entire result: *"this is not determined."*** ⛔ Do not close it with a normalization choice.
**Defect register:** **A8** — ⛔ **explicitly DEFERRED**, and ⛔ the earlier *"none"* was **FALSE**. A8 is
exactly this relation: `r_e` is called *"the throat-**body** size"*, from the balance `ke²/a = m_ec²`,
while `a` is the **mouth** radius (`conceptual_foundation.md:449-455`; `DEFECT_REGISTER.md` A8, **OPEN**).
⇒ Q4 **uses** the relation A8 flags and ⛔ resolves none of it. → S22.
**Parameter-register edges:** **R65**. → S22.

### Q5 · Charge universality ⛔ the result is *"this is not determined"*
**`SIM_GATED_REQUIRED_UNIVERSALITY`** — `q_h/Q_E` universality is **required, not earned** from the
current `b/ℓ` ledger (`ledger_stage042_charge_coupled_scalar.md:29`).

⭐ **Entire result: *"this is not determined."*** ⛔ *Required* is an assumption the ledger carries, not a
result it produced.
**Defect register:** **C5.**
**Parameter-register edges:** none.

### Q6 · Magnetism — the moving throat
The electric twin: **the same ±w throat, moving**.

```
J_{T,i} = q_T s_i η_a V_i
```

⛔ **FOUR different objects on two different routes — do not merge them:**
1. **The current FORM** `J_{T,i} = q_T s_i η_a V_i` — **DERIVED** from **defect continuity** (signed-dent
   continuity fixes the unique isotropic flux coefficient `α = 1`; `parameter_register.md:235`, R69
   `:335`). ⚠ Its **magnitude still rides the free `q_T`** ⇒ the compact-limit verdict is
   `CONVECTION_LIKE_CONDITIONAL`, ⛔ not a value.
2. **The Maxwell–Darwin kernel** `I_ij = (δ_ij + n_i n_j)/8πR` — a **Route-A REFERENCE**, `DERIVED` by
   **boosting the electric interaction** and ⛔ **explicitly tier-A conditional**: its overall coefficient
   rides the electric `A_E` (R63, `R1_REQUIRED(bc_selection)`) (`parameter_register.md:237`, R70 `:336`).
   ⛔ **Route A never touches `J_T`** (`ledger_stage036_route_a_maxwell_darwin_reference.md:214`) — it is a
   reference object, ⛔ not the moving-throat result.
3. **The direct moving-throat route (Route B)** — ⛔ **R1**: `R1_REQUIRED(direct_moving_throat)`, its
   magnitude carried by `q_T` (`parameter_register.md:238`, R71 `:337`; `docs/model_map.md#q-direct-route-b`).
4. **The structural COMPARISON of the two routes** — `BOOST_STRUCTURAL_RELATION_HOLDS`, **EARNED
   target-blind**: tensor structure, `R⁻²` falloff and `O(V₁V₂)` order agree, with the prefactor-stripped
   kernels symbolically equal (`parameter_register.md:238`; `docs/model_map.md#magnetism-part-v`). ⛔ **Emergent Lorentz
   is NOT claimed** — that needs `δ_BA = 0` **and** `r_cone = 1` **and** a closed `O(v⁴/c_γ⁴)`.

**Expected new:** the throat current **form**, the Route-A **reference** kernel, and the **structural**
route comparison. ⛔ **No magnitude and no sign.**
**Class:** the `J_{T,i}` form — `derived`, magnitude conditional on `q_T`; the Route-A kernel — `derived`
**tier-A conditional**; Route B — **`R1`**; the route comparison — `earned`; `q_T` — **FREE-UNREDUCED**,
`R1(throat)` (`:231`).
**Regime:** the **moving** throat. ⛔ Not gravitomagnetism — S19 is the PN frame-dragging step and is a
different object.
**Departure:** **`B_TIME_REVERSAL_EVEN`** — **DERIVED** (`:240`, R73 `:339`).

**Open:**
- ⛔ **The magnetic sign inherits the electric R1** — it is ***"doubly-R1"***, with co-blockers
  `direct_moving_throat` / `magnitude` / `consistency` (`model_map.md#q-compare-crux`, R72 `:338`).
- ⚠ `τ_d` is *"structural / postulated (**NOT** a reducible knob)"* (`:232`, R68 `:334`).

**Defect register:** none.
**Parameter-register edges:** **R68**, **R69**, **R70**, **R71**, **R72**, **R73**. → S22.

### Q7 · Coexistence, and the departure
**Computed:** **`internal_inconsistency = none`** — gravity and charge **coexist in one action**
(`parameter_register.md:229`, R64 `:330`). The **neutral far-field shell** `A = m_gg·C` is **DERIVED**
(`:224`, R60 `:326`).

**Departure: `NATIVE_P_NO_EMERGENT_GAUSS`** — the native `P^a` is **second-class at quadratic order**
and develops **NO emergent first-class U(1) Gauss constraint** (`parameter_register.md:230`, R66 `:332`;
decided 2026-07-12, commits `fd38efe8`, `d010e977`).

⇒ ⛔ **The target is the characterized-departure Maxwell analog, NOT exact Maxwell.** ⚠ Carry that in the
record's wording; a step that writes "Maxwell" unqualified has already lost the departure.
**Defect register:** none.
**Parameter-register edges:** **R60**, **R64**, **R66**.

---

## PHASE 5 — the knit and the deliverable (4 steps)

⚠ `S20a` is numbered from phase 4 but belongs here: the cone lock is a **knit** question, not a
gravity-sector derivation.

⭐⭐ **THIS PHASE CARRIES BOTH HALVES** (`CHARTER.md#two-halves`, user decision 2026-08-01):

| | half | what it is | steps |
|---|---|---|---|
| **1** | ⭐ **one medium supports the LINEAR part of every force** — all far-field effects | a **VERDICT**, able to fail: a no-go between two sectors' requirements **IS** the falsification | **S20a**, **S21** |
| **2** | ⭐ **what is left, that only a SIMULATION can settle** | an **INVENTORY** that sets up future work | **S22**, **S23** |

⛔⛔ **Nothing requiring a simulation is DONE here.** Half two specifies the work; it does not attempt it.

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

⇒ ⛔ **Both classifications are fixed, not open:** the **ratio** `λγ = c_γ/c_s` is **`derived`** (it is `R3`); the **equality `λγ = 1`** is **`calibrated`/uncommitted**. ⛔ Do not re-open the ratio's class. Introduce `λγ` with provenance, confront
`CONE_LOCK_CALIBRATED`, and ⭐ **apply the kind-test**: is this a medium-wide ratio, or a relation between
two sector speeds that happen to share dimensions? ⚠ That is exactly the shape of the thirteen
identifications in register §A.

⭐⭐ **`λγ = 1` is not a free landing — observation already constrains it.** In this model `c_s` is the
**gravity-change/phonon speed** and `c_γ` is the **light-cone speed**
(`research/pde_ledger_v2/notes/stages/ledger_stage005_sound_speed_light_ratio.md:75-77`):

> *"`c_s` is the phonon/gravity-change speed; `c_gamma` is the light-cone speed; `v_b` is the condensate
> flow."*

Gravitational-wave and electromagnetic arrival-time observations (**GW170817 / GRB 170817A**) constrain
those two speeds to agree to roughly **1 part in 10¹⁵**. ⚠ **The observational value IS known** — `λγ = 1`
to that precision is not in doubt as a *target*.

⛔⛔ **But `λγ` is currently a CALIBRATION, not a prediction.** `pathA_40` landed `CONE_LOCK_CALIBRATED`,
and **Lock A** (`λγ = 1`) is recorded as *"an untouched separate calibration"*
(`research/pde_ledger_v2/paper/stages/stage_040.tex:80-85`):

> *"**Lock A** (`λ_γ=1 ⟺ L_A`) is `WITNESSED` … at the witness, `L_A = m μ_R − 5Kρ⁴ρ_br = 5` (COMPUTED,
> asserted `==5`) — **an untouched separate calibration** (light ↔ gravity-phonon), NOT the Part-V one."*

⇒ ⭐ **A calibrated value cannot falsify — it was fitted.** At S20a the model **spends a knob** here and
**declines to predict** `λγ`. ⛔ Only a **DERIVED** value materially different from `1` would be a
falsification; landing on `1` by calibration tests **nothing**.

⚠ **That is the survey's standing verdict, not a new concession**
(`docs/medium_requirements_and_prior_art.md:181-182`):

> *"**Gate K — cone-lock `c_γ≈c_s` (B11).** In-plane rotational modulus vs bulk compressibility → two
> speeds, **equal nowhere automatically** (survey). Almost certainly a **calibration gap**, not a
> derivation — consistent with `λγ`'s current status. **Flag, don't bank**."*

⚠ **The specific worry to check:** `c_s` is **longitudinal compression** of the bulk while `c_γ` is
**transverse shear** of the brane. In an ordinary elastic solid those two differ by an O(1) factor
(`c_T/c_L ≈ 0.5–0.6`). ⛔ **That is a generic-elasticity expectation, NOT a model calculation** — carry it
as the risk to check, ⛔ never as a model prediction. If the model's shear modulus turns out to relate to
the bulk EOS the ordinary way, it lands near `0.55` — ⛔ which today's calibration would **absorb**, not
report as a failure.

⭐ **And the knob is LOAD-BEARING, not cosmetic** — `λγ³` enters the gravity normalization
(`ledger_stage005_sound_speed_light_ratio.md:80-82`):

> *"`lambda_gamma` enters the gravity-sector normalization to the fifth power (`(c/c_s)^3 =
> lambda_gamma^3`, and further factors downstream), so mis-stating it as `1` would silently smuggle a
> calibration."*

⇒ Spending the knob at `λγ = 1` propagates into the gravity-sector normalization. ⛔ Do not carry it as a
cosmetic convention.

⚠ **And `stage005` already REJECTS inferring `c_γ = c_s` from shared dimensions**
(`ledger_stage005_sound_speed_light_ratio.md:198`, `:204-205`):

> *"That `[c_s] = [c_gamma] = L T^-1` is explicitly **NON-EVIDENTIARY for equality**."*
>
> *"The parent action does not pin it; `c_gamma = c_s` from shared dimensions or from legacy weak-field
> prose is **REJECTED**."*

⇒ Observation requiring `λγ ≈ 1` and the ledger forbidding a dimensional shortcut to it are **compatible,
and both must hold** — ⛔ but **the timing is the whole point, so state it**: **today** `λγ` is a **spent
knob**, and all S20a does is *record that spend*. **Earning** the number that experiment already knows is
a **future route**, ⛔ not something S20a attempts or can be graded on. ⇒ That route is **Route A** —
registry edge **R10**, *"derive `λγ` from a nonlinear-throat `μ_R`-as-bulk-defect integral"*, status
`PENDING`, and it *"needs the deferred nonlinear throat"*
(`research/pde_ledger_v2/notes/parameter_register.md#R10`; the alternative Route B is `CLOSED-NEG` at
**R11**). ⛔ Nothing here licenses reading S20a as able to fail by landing off
`1` — a spent knob absorbs whatever it lands on.

⭐ **The SI anchor is calibration, not prediction.** Fixing model units by the measured `c` costs nothing
predictive — absolute constants are convention here, and only dimensionless ratios test the model.

### S21 · What light and gravity share
Both ride the brane; both consume `{ρ_br, μ_R}`. ⛔ **Use the method doc's CORRECTED rule**
(`docs/derivation_walkthrough_plan.md#knit-rule`), ⛔ **not** the earlier blanket ban on anything new:

> ⛔ the knit may not introduce a new **input** — action, constitutive, source or BC — **that revises what
> an earlier sector already derived** (such an input falsifies **completeness of the proposed substrate**,
> ⛔ not the one-medium hypothesis by itself); ⭐ it **is expected to produce new *consequences***,
> including constraints on the throat interior — that is the knit's **purpose**, not a violation.

⚠ **C6 still bites:** there is **no frozen complete parent action** to test "revises an earlier sector"
against. ⇒ This step's honest form is: *name what integration needs that phase 0–1 did not supply*, and
⛔ sort each item into a **revising input** (a substrate-completeness failure) or a **new consequence**
(the expected product) — ⛔ never bank an item unsorted.

### S22 · ⭐⭐ HALF TWO — what only a simulation can settle

#### ⭐⭐ THE NEAR-TERM LINEAR SIM INVENTORY — banked 2026-08-02, ⛔ do not re-derive it

⚠ **These need NO nonlinearity.** They are blocked by **one LINEAR object** — the brane–bulk interface
coupling law (**S11b**) — ⛔ **not** by the missing nonlinear shear action. ⇒ Much cheaper than the geon.

| # | what it measures | needs | blocked by |
|---|---|---|---|
| 1 | ⭐⭐ **transverse → longitudinal conversion at a defect** — how much light converts to the stray mode passing a particle. ⚠ The families do **not** mix in a *homogeneous* brane (coupling `∝ k·a`); `∇μ_R ≠ 0` mixes them | a defect profile | **form** available now with a generic profile; **magnitude** needs the interior (`R1`) |
| 2 | ⭐⭐ **does the longitudinal actually radiate into the bulk** | the interface law | **S11b** |
| 3 | ⭐ **the flexural crossover** — `ω ∝ k²` at long wavelength once the wall width is dynamical. ⛔ **Able to fail:** a clean cone means move 5 was wrong | width mode inertia + stiffness | **S11b** + **S5–S7** |
| 4 | **birefringence near a defect** — the two polarisations are degenerate *by symmetry* in a homogeneous brane; a defect splits them | a defect profile | as (1) |
| 5 | **Test 5** — set `χ_B = 0`, confirm shear propagation dies (`notes/brane_bulk_handoff.md:880`) | nothing new | ⭐ **runnable now**, designed, never run |

⛔⛔ **TAUTOLOGY GUARD:** a sim fed `μ_R` that reports `c_γ = √(μ_R/ρ_br)` has measured **nothing** — a
linear integrator already does exactly that. ⇒ Postulated parameters are fine to *run* with; they bound
what the run may **claim**. Conversion efficiency, mode structure, stability and confinement are genuine
outputs; ⛔ **the speeds are not.**

#### ⭐⭐ "STAYS QUANTIZED IN A PACKET" IS THREE DIFFERENT QUESTIONS — ⛔ do not merge them

⚠ The user's ambition is *"a photon travelling through space, self-sustaining."* Three readings, **three
different blockers**, and conflating them caused a false alarm that `ħ` was on the critical path:

| reading | needs | status |
|---|---|---|
| **a bound state** — a localized mode with a **discrete** frequency spectrum, trapped by a defect | ⭐ a linear wave equation in a **varying** medium (Sturm–Liouville) | ⭐⭐ **reachable, LINEAR.** No `ħ`, no nonlinearity |
| **a soliton** — one allowed lump, energy fixed by a size–amplitude relation | the nonlinear shear action | ⛔ **C6** — the real wall, and what the user actually wants |
| **`ħω` quanta** — energy in discrete lumps | field quantization | ⛔ **not a classical PDE sim at all**, at any effort |

⭐⭐ **⇒ Shelving `ħ_model` does NOT block the simulation.** A classical field has continuous energy no
matter how good the code is; row three was never reachable this way. ⇒ `ħ_model` connects the model to
quantum mechanics — a **separate project** from watching the brane move. ⛔ Do not re-open it here.


⛔⛔ **BROADENED 2026-08-01 (user decision: broad, ⛔ not photon-specific).** This step was *"the interior
debt list."* It is now **every nonlinear gap** — the throat interior, the geon, the drain law, the
**nonlinear brane-shear action**, and the **packet/soliton**.

⭐ **Why broad and not two separate parts:** S22's organising spine is already *"one nonlinear throat
solve"*, and the **nonlinear shear is the SAME missing object as the geon** (**C6** — no closed parent
action; `em_gravity_mined_verdicts.md:38` calls the trapped-shear geon *"intrinsically NONLINEAR"*).
⇒ Splitting them would file **one** missing parent action in **two** places.

⛔ **The light-sector rows this broadening ADDS** — ⛔ they were previously absent because the list was
interior-only:

| debt | what it blocks |
|---|---|
| the **nonlinear brane-shear action** (**C6**) | ⛔ every form in the corpus is **quadratic**, and the quadratic Lagrangian was written **directly**, ⛔ not expanded from a nonlinear parent ⇒ it ⛔ **cannot** be reached by un-linearising |
| the **packet / soliton** | linear `ω² = c_γ²k²` is exactly non-dispersive, so any shape translates unchanged — ⛔ but with **no preferred size or amplitude**. A packet with a characteristic scale **requires** a nonlinear term |
| **`R17`**, `μ_R = ∫ χ_B μ_R⁽⁴⁾ dw` (`parameter_register.md:284`) | the bulk→brane projection that would ground `μ_R` — **`PENDING`, "dim-consistency asserted only"** ⇒ the route by which the EOS reaches the brane sector |
| **`ħ_model`** (`DEFECT_REGISTER.md#C11`) | ⭐ **`postulated`, and a half-two row by user decision (2026-08-02).** ⛔ Not identified with `ħ_physical` — nothing forces it, and doing so later is a **separate calibration**. Either **derive** it (needs the substructure, which is EMPTY) or **define it operationally well enough for the sim to run**. ⛔ Do not spend derivation time on it |

⛔ **GATED on A4, C4 and C9** — this step classifies exactly the quantities those rows dispute, so it
cannot bank without confronting them. ⚠ **C9 is a live source disagreement** (`FREE-UNREDUCED` vs
*"likely DERIVED"* for the same three symbols): preserve both readings as unresolved; ⛔ do not pick one
silently.
⭐⭐ **The organizing spine is ONE nonlinear throat solve** — ⛔ not a per-sector debt list
(`docs/model_map.md#shared-r1-throat-solve`):

> *"One nonlinear throat solve is the shared R1 for gravity {μ_R, ρ_br} … electric `bc_selection`, and
> magnetism `q_T` — one interior solve collapses several knobs at once"*

⇒ `Q_E` and `bc_selection` are **siblings of R10/R30** (`parameter_register.md:143`, `:329`), ⛔ **not a
separate charge-sector programme.**

⭐ **Section the table by what each sector is waiting on** — gravity-side · charge-side · magnetism-side ·
the shared throat packet — and then keep the two non-throat loci separate: brane constitutive inputs ·
external/calibrated PN inputs. A flat list hides which debts one solve would actually discharge.
⛔ **`J`, `m_defect` and the geon are THREE distinct unbridged objects**, not one debt: the geon profile
is undeclared (**C4**), `m_defect` has only a dimensional bridge (**C3**), and nothing links either to
`J`.

**(i-a) Throat interior — gravity-side:**

| debt | what it blocks |
|---|---|
| `{mu_R, rho_br}` (R10, opened at S8) | `c_gamma` and every cone lock |
| `{mu_eta, T_w, beta}` (R30) | frozen-wall response. GATED ON A4 (`beta*a = 1` is an unearned lock) |
| `{Vp0/ell_c}` (R33) | breathing drive |
| `{T_Omega, beta_2}` (R36) | l=2 support. C9 IS LIVE: the same symbols are `FREE-UNREDUCED` in one document and "likely DERIVED" in another. Preserve both readings; do not pick one silently |
| **`J`** the drain rate | gravity's source term |
| **`m_defect`** | a SEPARATE object: only a dimensional bridge `alpha_J hbar J/c_gamma^2` exists (C3) |
| `INFLOW_MASS_SOURCE_MISSING` | the S13 gap itself: the source of gravity is not derivably connected to mass |
| **the geon** | a THIRD separate object, profile UNDECLARED (C4) |
| the `J` / `m_defect` / geon BRIDGES | none exist. Three unbridged objects, not one debt |
| density-port magnitudes | Gate-6 numbers, SIM-deferred |

**(i-b) Throat interior — charge-side** (PHASE 4b; siblings of R10/R30,
`parameter_register.md:143`, `:329`). ⛔ **One heading is not one class.** The class column is the method
doc's — `debt` · `postulated` · `calibrated` (`docs/derivation_walkthrough_plan.md#classification-test`) — and ⛔ a
**blocked OUTPUT is not an input debt at all**:

| quantity | class | what it blocks / why listed |
|---|---|---|
| `bc_selection` (R63) | `debt` | the electric SIGN — and ⛔ that sign is a BLOCKED OUTPUT, not an input. `outcome_not_invariant` across the four BC branches (Q3); the sign is NEITHER earned NOR calibrated |
| `Q_E` (R65) | `debt` | charge MAGNITUDE. `R1_REQUIRED(magnitude)` — a named route under the deferred throat solve |
| `c_a`, `c_xi` (R65) | `debt`, ⛔ INSIDE `Q_E` | core normalizations, "unbounded at tier-A"; ⛔ not a separate debt — they ride `Q_E`'s magnitude |
| **`a`** (throat mouth radius) | `postulated` | Q4's magnitude waits on it. `FREE-UNREDUCED`, no defining equation, ⛔ NO named route (`parameter_register.md:132`); `a = c_a r_e` parameterizes it against the EXTERNAL `r_e` and does not reduce it. GATED ON **A8** (mouth radius vs throat-body size) |
| `ell` (embedding / PT healing length) | `calibrated` | ⛔ NOT an interior debt and NOT a throat-solve output: `IMPOSED`/`CALIB`, the ratio `ell/a = 1/20` is a frozen handoff scale, a tuning (`parameter_register.md:200`). Listed only because Q4 waits on it; ⛔ NO reduction route is registered |
| the `𝔅` / mouth boundary class | `debt` | the branch selection S12 left open on the boundary-data side |

**(i-c) Throat interior — magnetism-side** — ⛔ same class discipline as (i-b):

| quantity | class | what it blocks / why listed |
|---|---|---|
| `q_T` (R1(throat), `parameter_register.md:231`) | `debt` | FREE-UNREDUCED; the moving-throat current strength |
| `tau_d` (R68) | `postulated` | ⛔ NOT a debt and ⛔ NOT a Phase-3 inheritance: "structural / postulated (NOT a reducible knob)" with NO reduction route (`parameter_register.md:232`); Q6/stage-034 supplies it itself |
| the magnetic SIGN | ⛔ BLOCKED OUTPUT, not an input | inherits the electric R1 — "doubly-R1", co-blockers `direct_moving_throat` / `magnitude` / `consistency` (R72) |

**(i-d) The shared throat packet** — `{mu_eta, T_w, beta, Vp0/ell_c, T_Omega, beta_2}`: all six are
recorded as **siblings of the ONE deferred throat-interior solve**, via R30/R33/R36
(`parameter_register.md:49`), which is why they re-appear sector-by-sector above.
⛔ **That sibling relation is WEAKER than a six-way collapse, and the map says so in as many words:**
gravity's shared R1 is `{mu_R, rho_br}` via **R10 + R30 + R33 — *"not all six"***
(`docs/model_map.md#shared-r1-throat-solve`). ⇒ Claim the siblinghood; ⛔ do **not** claim one solve is known to discharge
all six.

**(ii) Brane constitutive** — would survive a throat solve untouched:

| debt | note |
|---|---|
| `W_slab` | `FREE-UNREDUCED`; the double-well selects no width (C7) |
| `{a_B, kappa_B}` | the wall's own postulated constants (C8) |
| the cone lock `lambda_gamma = 1` | calibrated/uncommitted (S20a) |

**(iii) External / calibrated PN inputs** — not interior debts at all:

| debt | note |
|---|---|
| `a_WT` (worldtube profile width) | controls the `O(a_WT^2/r^2)` correction |
| the `a_WT` / `a_mouth` BRIDGE | DOES NOT EXIST (A11). Never close it by substitution |
| `M` and the supplied multipoles | inputs of the S16 reduction, not outputs |
| `G` / `N_inf,3` normalization, `I_F,12` | the force magnitude is calibrated through them |
| `Theta_Q` and the `J -> Q` map | S15 must not substitute silently |
| `kappa_add`, `kappa_PV`, `kappa_rho`, optical `n = 5` | the PN response packet S17 imports |

⇒ **This is what v3 is for.** Today it is scattered across six documents on the gravity side alone, and
across the charge/magnetism register rows besides; here it is one page.

### S23 · The count — sim-input set, partitioned
⛔ **Never one integer.** Scalars · profiles/functions · BCs and domains · discrete choices.
⛔ **No number leaves without the §7a closing certification** (admission gate · transitive leaf closure ·
block rank · top-down reconciliation · the sim-input-vs-residual diff). The walkthrough produces an
*inventory*; only §7a certifies it.

---

## 1. Order of attack, and the honest expectation

**S1–S4** are cheap re-banking; they exist to make v3 self-contained.
**S5–S8** are the honest core — three of the four bank a *postulate* or a *gap*.
**S9–S11, S16–S19** are where the earned results live — and so are **Q2** (the `1/R²` falloff and the
`s₁s₂` product, target-blind) and **Q6** (the boost structural relation).
⛔ **S14, S15 and S20 are CONDITIONAL ON S14a**, not earned.
**Q1–Q7 (PHASE 4b)** carry ⛔ **per-step** dependencies, ⛔ **not a phase-wide one**: **Q1–Q5** do **not**
need the drain bridge; **Q6** needs PHASE 2's transverse row, and its active-drain arrow `τ_d` is ⛔ **its
own POSTULATED input, not a PHASE 3 inheritance** (`parameter_register.md:232`). All seven are
⛔ **independent of the PN ladder**, so they are not gated on S17/S18 and do not wait on them.
**S0.5 and the substrate-action step are prerequisites** — ⛔ nothing is banked before them.
**S22 is the deliverable.**

⚠ **Expect roughly this shape:** ~9 steps introducing genuinely new inputs, ~9 pure-consequence steps,
~7 steps whose entire result is *"this is not determined"*, and 1 (S22) that is the point of the
exercise. ⭐ **Three of those seven are Q3, Q4 and Q5** — the charge sector's sign, magnitude and
universality — ⛔ and "not determined" is their entire result, not a stalled step to come back and finish.
⛔ If S5–S8 start producing derivations instead of postulates, something is being smuggled —
that is the phase where the model is known to be weakest.

## 2. ⛔ Failure modes to watch, specific to this plan

1. **Re-deriving the PN ladder.** S17 is cite-only. It is audited and DOI'd; touching it is pure cost.
2. **Letting S22 become a census.** It is one table with loci. It is not a provenance system.
3. **Treating a clean close as validation.** See `CHARTER.md` §1.1. ⛔ A tidy v3-gravity does **not**
   mean the far field is independent of the interior — S16 does not establish that. It means the far
   field is computable **given supplied** mass, multipoles, compactness and a calibrated normalization.
4. **Back-filling `c_γ` into phase 0** to make the registry balance at S2. That is the exact error O-01
   recorded.
5. **Building registry machinery as a precondition** for banking a step. Add one quantity and one
   relation when a step needs them; nothing more.

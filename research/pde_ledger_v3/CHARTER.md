# v3 LEDGER — charter

**User decision, 2026-07-31.** Opened after step ① (the `a`-pin retirement, `407eed94`).

---

## 0. ⛔ Method history — and v3 DID take a method change

⚠ **This section said *"this is not a third method change"* until 2026-08-01. That sentence is now
FALSE and has been rewritten rather than annotated.** Two distinct events:

- **v3 OPENED as a scoping change** (2026-07-31). Same forward walkthrough, same step record, same
  checks, same registry — a **sector boundary** in place of a phase ordering. That much was true.
- ⭐⭐ **v3 then TOOK a method change (2026-08-01, user decision): REQUIREMENTS-FIRST.** Each force
  sector states what it needs to survive. **Brane and bulk are defined LAST, at the knit**, by asking
  whether one medium can satisfy every requirement at once. **A no-go between requirements IS the
  falsification.**

**Why the change:** deriving forward from the medium is **circular** here. `A13` (is `χ_B` real or
complex?) has no answer *in the medium* — it is settled by what the drain law needs — while the plan
demanded it resolved before the steps that would settle it.

⚠⚠ **Do NOT claim this is "a return to the project's original method."** The record is contested and
BOTH sides are damaged: `docs/em_gravity_requirements_inversion.md` carries a ⛔ *"SUPERSEDED … Do NOT
build on it"* banner, and `docs/em_medium_first_generative_plan.md` claims founding-method status for
**medium-first** while calling requirements-inversion *"not generative"*. ⇒ The direction rests on the
**user's decision** plus the circularity argument, ⛔ not on a historical claim.

⚠ For contrast, the change this section originally warned about: census → walkthrough (2026-07-30) cost
eleven commits and verified no physics.

⇒ The method is `docs/derivation_walkthrough_plan.md`. ⛔ This charter does not restate it and does not
supersede it. If you find yourself designing new apparatus, stop — that is the failure this project has
already had twice.

⭐⭐ **Third hard constraint, added 2026-07-31: every step is walked SIDE BY SIDE with the user.** v2 was
delegated heavily to go fast, and that is why the `a`-pin, the falsified lepton tower and the Spin
Problem all went unnoticed by the one person able to check them. ⇒ `V3_STEP_PLAN.md` §0 is the operative
description; it is not a courtesy stop-for-review but a change in who does the deriving.

**Two hard constraints, both from the decision that opened v3:**
1. ⛔ **Do not rebuild the registry machinery.** Registry, reader, dimensional gate, able-to-fail
   harness — they work, and they were built to grow one step at a time. That **is** the centralized
   variable script.
   ⚠ **AMENDED 2026-08-01 — v3 is SELF-CONTAINED (user decision).** The machinery was **copied** into
   `research/pde_ledger_v3/reduction/` (`361e8114`), and loci are repo-root-relative so v3's citations
   *into* v2 are explicit. ⛔ **Do not edit or import `research/pde_ledger_v2/reduction/`** — the
   original instruction to "reuse" it predates the self-containment rule and would split the registry
   across two trees.
2. ⛔ **`DEFECT_REGISTER.md` is a list, not machinery.** The moment it acquires a schema, a checker, or a
   manifest, it has become the census.

---

## 1. Scope — gravity · light · gravitomagnetism · charge · magnetism

**Widened 2026-07-31 (user decision): the EM/charge sector is IN scope.** ⭐ The reason is not appetite —
it is that these are not independent sectors that happen to share a deadline, they share the **same
unsolved object**. `docs/model_map.md#shared-r1-throat-solve` records that **"One nonlinear throat solve is the shared R1 for
gravity `{μ_R, ρ_br}` ... electric `bc_selection`, and magnetism `q_T` — one interior solve collapses
several knobs at once"**, and the register books charge's debts against that same solve: `Q_E` and
`bc_selection` each discharge *"under the shared sim-deferred throat solve (sibling of R10/R30)"*
(`research/pde_ledger_v2/notes/parameter_register.md:143`, `:228`; the `R63` edge row states it at `:329`
as *"a sibling of R10/R30/R33"*) — `R10`/`R30` being gravity's interior debts.

⇒ ⭐ **A debt list that names gravity's siblings of that solve while omitting charge's and magnetism's is
an INCOMPLETE listing of one object — and it would read as complete.** That is §1.1's failure mode applied
to the scope boundary itself.

⚠ **A working boundary, and its justification is CONDITIONAL — do not read it as established.** The
gravity and light sectors are **linear response on the brane**, far from the defect, and the corpus
offers a **response-side** reduction (how a compact body *moves in* a given field) that is leading-order
only and premised on compactness. ⛔ It is **not** a source-side interior-independence theorem, and
S16/S19 may show the premise fails (see §1.1). The quotations below are that conditional result, not a
licence:

> *"controlled by the worldline/worldtube multipoles rather than by arbitrary internal details of the
> defect"* — `research/4d_2_5pn/paper/4d_2_5pn.tex:613`
> *"we do not attempt to solve the fully dynamical moving-throat problem in the bulk"* —
> `research/4d_1pn_full/paper/4d_1pn_full.tex:110`

Gravitomagnetism belongs here: it is a PN-order effect inside the audited ladder, ⛔ **not** the EM
sector's "magnetism = moving throat".

⚠ **And charge is NOT insulated from the interior the way the brane-shear and linear-response sectors
are.** Its far-field **FORM** is target-blind EARNED — *"The `1/R²` falloff and the `s₁s₂` product of the
far-field shell are target-blind EARNED"*
(`research/pde_ledger_v2/notes/stages/ledger_stage031_puncture_deflection_field_identity_source.md:25`) —
but its **SIGN and MAGNITUDE route through the interior solve**. ⇒ Bringing charge in widens the debt
list; it ⛔ does **not** import a second insulated sector.

### ⛔ 1.1 What this boundary EXCLUDES — state it or it will read as completeness

Gravity looks clean **because it is insulated from the hard part** — to the extent it is insulated at
all — not because the hard part is solved.
v3 **will not claim to discharge** — and must not appear to have settled:

- the geon and the mass mechanism (**C4**)
- the drain law and `g_phys` (**C1**)
- `R*(E)`, the mass–radius relation, and the throat size (**C2**, **D1**)
- what distinguishes the leptons (**D1**)
- brane existence (**B2**, **C7**, **C8**) — ⚠ **`B2`'s scope, so this line is not over-read:** the
  *gate program* that would have grounded the brane is superseded **at the brane-existence level**; the
  no-go itself closes only *deriving* `μ_R` from a polar `P`, and light still stands on the
  **postulated** modulus (`DEFECT_REGISTER.md#B2`). ⛔ Not "light is foundationless."

⛔ **A v3-gravity ledger that closes cleanly is not evidence the model is sound.** At most it shows the
far field is computable *given supplied* mass, multipoles, compactness and a calibrated normalization — {#clean-close-not-validation}
the conditional response-side conclusion S16 actually permits. ⛔ It does **not** establish that the far
field is independent of the interior.

⛔⛔ **AMENDED 2026-07-31 — the boundary is weaker than stated above, and may not hold at all.**
Reading the worldtube theorem's actual statement rather than its slogan
(`research/4d_2_5pn/paper/4d_2_5pn.tex:605-614`), the independence is **leading-order only** and
**conditional on the defect being compact (`a ≪ r`)**, with a first correction `O(a_WT²/r²)` — controlled by the
**worldtube profile width**, which is a different object from the mouth radius (A11) and is itself
undetermined.

⚠ And the Spin Problem (**C10**) concludes you *cannot* get frame dragging from a compact defect —
*"you need the tail"*. ⇒ **If that is right, the compactness premise fails and the boundary this charter
is built on does not exist as stated.** S16 and S19 must therefore be run **before** the charter's scope
claim is treated as settled, and the charter amended to whatever they find. ⛔ Do not preserve the
boundary against the result.

### 1.2 ⭐⭐ THE DELIVERABLE HAS TWO HALVES (user decision, 2026-08-01) {#two-halves}

⛔ **This replaces the earlier framing *"a debt list, not a discovery."*** That framing was purely
negative and it understated half one.

**HALF ONE — the positive claim, and the user believes it is reachable:**

> ⭐⭐ **ONE medium supports the LINEAR part of every force — all the far-field effects of all the
> forces.**

⭐ **This is the part that can be FULLY DERIVED**, which is why it is the ledger's body. Under
requirements-first each sector states what it needs; the **knit** asks whether one medium supplies every
requirement at once. ⛔ It is **able to fail**: a **no-go between two sectors' requirements IS the
falsification** (§1.3), ⛔ not a snag to engineer around.

**HALF TWO — the hand-off, and it is BROAD (user decision: broad, not photon-specific):**

> ⭐ **What is left, that can only be done by SIMULATION** — stated precisely enough to start the work.

⛔⛔ **Nothing requiring a simulation is DONE in this ledger.** Half two **sets up future work**; it does
not attempt it. ⭐ It is **broad** — every nonlinear gap, not just light's: the throat interior, the geon,
the drain law, the nonlinear brane-shear action, the packet/soliton. ⚠ **It therefore SUBSUMES the
interior debt list** (**S22**), because S22's organising spine is already *"one nonlinear throat solve"*
and the nonlinear shear is **the same missing object as the geon** (**C6** — no closed parent action).
⇒ One closing part, ⛔ not two.

⚠⚠ **THE LINEAR/NONLINEAR SEAM IS A SEAM, ⛔ NOT A CEILING.** Do not write that the ledger "stops" at
linear response, and ⛔ do not record missing nonlinearity as a *blocker*. Linear is **half one's scope**;
nonlinear is **half two's subject**. Same fact, and the wrong framing has already caused a session to
treat the photon simulation as gated rather than as the deliverable.

**What half two collapses.** Today the debts are scattered across `R10`, `R30`, `R33`, `R36`,
`m_defect`, `J` and `INFLOW_MASS_SOURCE_MISSING` on the gravity side, and across `bc_selection`, `Q_E`,
`q_T`, `c_a` and `c_ξ` on the EM side. Collapsing that into one honest list is the win — it is a thing
that fits in one head, which is what is currently missing.

⚠ **And it is still not a discovery.** Half two names debts; it pays none of them down (§5).

### 1.3 ⭐ What would count as FAILURE {#falsification-standard}

⭐ **The point of this model is to show that circumstances CAN EXIST to support the forces we observe — ⛔
NOT to derive the universe from first principles** (user, 2026-07-31). ⇒ ⛔ **Fidelity to standard physics
is the wrong test:** superfluid hydrodynamics describes media of atoms with EM interactions, and those
equations carry that microphysics inside them; this model postulates a medium with whatever properties the
forces require. ⇒ **Departing from real-superfluid behaviour is EXPECTED, ⛔ NOT a defect.**

**What WOULD falsify:** an **internal impossibility** — ghost modes, energy unbounded below, lost
hyperbolicity, a constraint-count failure, a causality violation · a **DERIVED** value contradicting data
(⚠ a *calibrated* value cannot falsify, it was fitted; only a derived one can) · **running out of
surplus**, more knobs spent than held-out matches earned. ⛔ **NOT failure:** differing from helium-4 or
any real superfluid · failing to *derive* a legitimately postulated quantity (⚠ that leaves it postulated
— a **cost against surplus** and a §1.2 debt row, not a failure).

⭐ ⇒ **Postulating freely is not free — every postulate is a knob**, and the model is scored on surplus:
the burden shifts from *"is this like a real superfluid"* to *"how many knobs did this cost"*.

⭐⭐⭐ **STRONGER THAN "EXPECT EXTRA" — AN EXACT MATCH WOULD BE THE FAILURE** (user, 2026-08-02):

> We were never expecting exact Maxwell, and **we cannot get it, because we are IMPORTING Maxwell into
> the model.** ⭐⭐ *"If it was exact Maxwell then it wouldn't have a way to physically anchor to the
> model."*

⇒ Maxwell puts **charge in by hand** as an external source. A model reproducing it **exactly** would
**inherit that gap** — charge would have nothing in the model to attach to. ⇒ ⭐ **The extra structure IS
the anchor point**, ⛔ not a blemish on an otherwise clean match. ⚠ On record: the extra longitudinal is
what made the drum-head charge picture click (`V3_STEP_PLAN.md` § S11).
⛔⛔ **Many `FAIL_*` tokens are simply MISNAMED** — `FAIL_CAUCHY_STRAY_LONGITUDINAL` above all. ⛔ Never
read the prefix as a verdict.

⭐⭐ **EXPECT EXTRA — ⛔ do not score it as failure** (user, 2026-07-31). **The comparison targets are
themselves incomplete:** Maxwell puts **charge** in by hand, GR **matter** as a stress tensor — ⇒ neither
*derives* what this model derives. ⇒ **Deriving what the reference assumes yields structure it has no slot
for**, and a match-check reports that surplus a **defect**. ⛔ It is not one. ⭐ On record: the light sector
gave two transverse modes **and** an irremovable extra scalar; that irremovability killed exact Maxwell as
a target and forced the rethink that produced `±w` charge (**S11 · The stray-longitudinal departure**,
`research/pde_ledger_v3/V3_STEP_PLAN.md`). ⇒ ⛔ A `FAIL_*` recording *"more than the reference has"* is a
**characterized departure**, ⛔ not a falsification — ⚠ only the modes above falsify. ⚠ And the naming
misleads: `FAIL_CAUCHY_STRAY_LONGITUDINAL` names a load-bearing departure — ⛔ never read a `FAIL_` prefix
as a verdict without opening what it measured.

---

## 2. Order

**Step 0 — import the register.** ✅ `DEFECT_REGISTER.md`. Nothing else starts until the known defects
have a home; that was the stated reason for opening v3.

**Then the walkthrough**, forward, one step at a time, from the medium's defining properties into the
in-scope sectors of §1 — gravity, light, gravitomagnetism, charge and magnetism. Each step records
what the method doc specifies: *what it is · what it does · what's new · inputs by locus · the
equation(s) · class per new item · regime · departure.*

⛔ **One step at a time, with a stop for the user at each.** No fan-out across steps.

---

## 3. What carries forward from v2, and what does not

**Carried — do not re-derive:**
- the bulk GPE/Madelung system (textbook, committed)
- the PN ladder 1PN→4PN + 2.5PN, audited, GR-matched, dual-engine verified, **worldtube-reduced**
- the ℓ=2 DtN fingerprint `{1/9, 4/81, 1/27}`, `χ_Q = +1`, cross-ℓ `{1, 1/2, 1/27}`, SO(3) `λ_m = 6`
- the dimensional foundation, **minus the pin**
- ⭐ **every row of `DEFECT_REGISTER.md`** — especially section B, which does not look like progress and
  is exactly what a clean start loses
- ⭐ **the EM/charge cluster — its deferral is RETIRED (2026-07-31).** The stated reason ("that
  workstream's charge model changed recently") was checked against the record and is not supported: the
  charge mechanism changed **2026-07-21** (puncture-deflection, `520d4b1f`) and the sector's last build
  commit is **2026-07-23** (`0961b27c`). ⚠ **No later commit changed the puncture-deflection charge
  mechanism** — `087565b0` (2026-07-24) altered the shared parent-action integration and added a
  `Z_χ` draft knob, not the mechanism. What is open is labelled `R1_REQUIRED` — waiting on the
  throat solve — ⛔ **not** a model in flux.

⛔ **CONDITIONAL — neither carried nor discarded, pending S14a:** {#conditional-s14a}
- the `1/r²` law, the attractive sign, and the stage009 `RETURN_RESIDUAL_PREDICTION` falsifier.
  Correcting the drain to the dynamical `Γ_B` law severed the chain to these — they were derived under a
  Gauss **number-flux** drain. ⛔ Do not cite them as earned until the S14a bridge succeeds.

**Not carried:**
- v2's irreducible count and its headline knob claims (**E1** — understated, and v3 derives its own
  count forward rather than importing one)
- the census apparatus and manifests (already archived)

**And one distinction inside the v2 walkthrough directory:**
- ⚠ **The v2 walkthrough step records** (`research/pde_ledger_v2/walkthrough/00_*`, `01_*`) are **prior
  art, ⛔ not authority** — v3 derives S1/S2 forward. ⭐ **But `DECISIONS.md` in that directory is
  different: its numbered decisions are USER CALLS and they bind**, and the method doc
  (`docs/derivation_walkthrough_plan.md:7-10`) records that they **overturn** it where they disagree.

---

## 4. Acceptance

v3 is working if, at any point:

1. the accumulated "what's new" is traceable, each member to the step that introduced it;
2. every step passes its mandatory checks, with failures **recorded, not fixed by adjustment**;
3. a reader can follow the chain **without consulting v2** — ⚠ under requirements-first that chain runs
   **from each sector's requirement to the knit**, ⛔ not from the medium outward;
4. **HALF ONE** — the knit returns a verdict on whether **one medium supplies every sector's linear
   requirement at once**, with each requirement stated, sourced and falsifiable;
5. **HALF TWO** — everything that can only be settled by simulation is named in one place, with loci,
   across every in-scope sector.

⛔ **(4) and (5) are the deliverable, and they are different objects.** (4) is a **verdict** — and
⭐ **a no-go is a first-class result**, not a failed run. (5) is an **inventory**.

⛔ A count is not quoted without the closing certification the method doc specifies — the forward
walkthrough produces an inventory, not a certified number.

---

## 5. ⭐ The thing this charter cannot address, recorded so it is not forgotten

**B1** falsified the lepton mass tower, and with it the only family label the model had. Nothing now
explains why three leptons of wildly different mass carry *exactly* the same charge (**D1**).

That question is **out of scope here by construction** — it lives in the throat interior. ⇒ It will still
be open when v3 closes, and closing v3 will not have made progress on it.

Recorded because a debt list is only honest if it says what it is *not* paying down.

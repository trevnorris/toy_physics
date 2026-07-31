# v3 LEDGER — charter

**User decision, 2026-07-31.** Opened after step ① (the `a`-pin retirement, `407eed94`).

---

## 0. ⛔ What this is NOT

**This is not a third method change.** Read that sentence before anything else.

- census → walkthrough (2026-07-30) **was** a method change. It cost eleven commits and verified no
  physics.
- v2 → v3 is a **scoping** change. Same forward walkthrough, same step record, same checks, same
  registry. The only difference is a **sector boundary** in place of a phase ordering.

⇒ The method is `docs/derivation_walkthrough_plan.md`. ⛔ This charter does not restate it and does not
supersede it. If you find yourself designing new apparatus, stop — that is the failure this project has
already had twice.

⭐⭐ **Third hard constraint, added 2026-07-31: every step is walked SIDE BY SIDE with the user.** v2 was
delegated heavily to go fast, and that is why the `a`-pin, the falsified lepton tower and the Spin
Problem all went unnoticed by the one person able to check them. ⇒ `V3_STEP_PLAN.md` §0 is the operative
description; it is not a courtesy stop-for-review but a change in who does the deriving.

**Two hard constraints, both from the decision that opened v3:**
1. ⛔ **Reuse `research/pde_ledger_v2/reduction/`.** Registry, reader, dimensional gate, able-to-fail
   harness — they work, and they were built to grow one step at a time. That **is** the centralized
   variable script. Do not rebuild it.
2. ⛔ **`DEFECT_REGISTER.md` is a list, not machinery.** The moment it acquires a schema, a checker, or a
   manifest, it has become the census.

---

## 1. Scope — gravity, light, gravitomagnetism

⭐ **This is a principled boundary, not a convenience.** Both sectors are **linear response on the
brane**, far from the defect. That is *why* they are the solid ones:

> *"controlled by the worldline/worldtube multipoles rather than by arbitrary internal details of the
> defect"* — `research/4d_2_5pn/paper/4d_2_5pn.tex:613`
> *"we do not attempt to solve the fully dynamical moving-throat problem in the bulk"* —
> `research/4d_1pn_full/paper/4d_1pn_full.tex:110`

Gravitomagnetism belongs here: it is a PN-order effect inside the audited ladder, ⛔ **not** the EM
sector's "magnetism = moving throat".

### ⛔ 1.1 What this boundary EXCLUDES — state it or it will read as completeness

Gravity is clean **because it is insulated from the hard part**, not because the hard part is solved.
v3-gravity will not touch, and must not appear to have settled:

- the geon and the mass mechanism (**C4**)
- the drain law and `g_phys` (**C1**)
- `R*(E)`, the mass–radius relation, and the throat size (**C2**, **D1**)
- what distinguishes the leptons (**D1**)
- brane existence (**B2**, **C7**, **C8**)

⚠ **A v3-gravity ledger that closes cleanly is not evidence the model is sound.** It is evidence that
the far field does not depend on the interior — which was already known.

⛔⛔ **AMENDED 2026-07-31 — the boundary is weaker than stated above, and may not hold at all.**
Reading the worldtube theorem's actual statement rather than its slogan
(`research/4d_2_5pn/paper/4d_2_5pn.tex:605-614`), the independence is **leading-order only** and
**conditional on the defect being compact (`a ≪ r`)**, with a first correction `O(a²/r²)` — controlled
by the throat radius, which is undetermined.

⚠ And the Spin Problem (**C10**) concludes you *cannot* get frame dragging from a compact defect —
*"you need the tail"*. ⇒ **If that is right, the compactness premise fails and the boundary this charter
is built on does not exist as stated.** S16 and S19 must therefore be run **before** the charter's scope
claim is treated as settled, and the charter amended to whatever they find. ⛔ Do not preserve the
boundary against the result.

### 1.2 ⭐ Expected yield: a debt list, not a discovery

Say this up front so it does not read as failure when it arrives. The deliverable is a **precise,
followable statement of what gravity assumes from the interior**, which today is scattered across `R10`,
`R30`, `R33`, `R36`, `m_defect`, `J` and `INFLOW_MASS_SOURCE_MISSING`. Collapsing that into one honest
list is the win — it is a thing that fits in one head, which is what is currently missing.

---

## 2. Order

**Step 0 — import the register.** ✅ `DEFECT_REGISTER.md`. Nothing else starts until the known defects
have a home; that was the stated reason for opening v3.

**Then the walkthrough**, forward, one step at a time, from the medium's defining properties into the
gravity and light sectors. Each step records what the method doc specifies: *what it is · what it does ·
what's new · inputs by locus · the equation(s) · class per new item · regime · departure.*

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

⛔ **CONDITIONAL — neither carried nor discarded, pending S14a:**
- the `1/r²` law, the attractive sign, and the stage009 `RETURN_RESIDUAL_PREDICTION` falsifier.
  Correcting the drain to the dynamical `Γ_B` law severed the chain to these — they were derived under a
  Gauss **number-flux** drain. ⛔ Do not cite them as earned until the S14a bridge succeeds.

**Not carried:**
- v2's irreducible count and its headline knob claims (**E1** — understated, and v3 derives its own
  count forward rather than importing one)
- the EM/charge cluster — deferred; that workstream's charge model changed recently
- the census apparatus and manifests (already archived)

---

## 4. Acceptance

v3 is working if, at any point:

1. the accumulated "what's new" is traceable, each member to the step that introduced it;
2. every step passes its mandatory checks, with failures **recorded, not fixed by adjustment**;
3. a reader can follow the chain from the medium to any step **without consulting v2**;
4. the interior debts gravity depends on are named in one place, with loci.

⛔ **(4) is the deliverable.** ⛔ A count is not quoted without the closing certification the method doc
specifies — the forward walkthrough produces an inventory, not a certified number.

---

## 5. ⭐ The thing this charter cannot address, recorded so it is not forgotten

**B1** falsified the lepton mass tower, and with it the only family label the model had. Nothing now
explains why three leptons of wildly different mass carry *exactly* the same charge (**D1**).

That question is **out of scope here by construction** — it lives in the throat interior. ⇒ It will still
be open when v3-gravity closes, and closing v3-gravity will not have made progress on it.

Recorded because a debt list is only honest if it says what it is *not* paying down.

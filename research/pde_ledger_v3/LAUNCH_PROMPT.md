# Launch prompt — paste this after /compact or at the start of a new session

⚠ **This file is a POINTER, not a summary.** The last session's launch prompt was a lossy hand-compression
that lost a whole review round and cost the first hour. ⛔ Do not let this one become that: it names
documents and states rules; it does **not** restate findings. If it and a named document disagree, ⛔ **the
document is right.**

---

Branch `ledger-v2-rebuild`. Run `git log --oneline -5` and `git status` first —
⛔ do not trust any hash written in a doc, including this one.

▶ **READ FIRST, AND IN FULL: `research/pde_ledger_v3/SESSION_2026-08-01.md`.**
⚠ It is the **2026-08-01** record and it governs **for what it covers**. ⛔ It is NOT the latest — the **2026-08-02** work (the two halves, `K`/`n_eos`, the drum-head charge picture, S10) lives in **this file**, `CHARTER.md`, `V3_STEP_PLAN.md` and the memories. It carries the direction change, the model's foundational
postulate, the five standing decisions, what was banked, the simulation track, the open physics, and
§8 — how the previous session failed, by mechanism. ⛔ Do not start work before reading it.

▶ **THEN the orientation budget — and it is EIGHT documents, not the five `NEXT_SESSION.md` lists.**
The three it omits all govern:
- `docs/derivation_walkthrough_plan.md` — the method; **§5a** is why the registry is grown forward
- `research/pde_ledger_v2/walkthrough/DECISIONS.md` — ⛔ overturns the method doc where they disagree
- `docs/development_pipeline.md` — the operative process doc

▶ **`NEXT_SESSION.md` is still worth reading for its ELEVEN-ITEM list** of what a fresh session gets
wrong — each item is something an oriented orchestrator actually believed. ⛔ But its ordering and status
claims are stale; the session record says which.

---

## ⭐⭐ THE RULES — these are not style guidance

1. **Every step is walked SIDE BY SIDE.** Derive in the open, one move at a time, reasoning BEFORE the
   result. ⛔ Never pre-derive and present. ⭐ Flag every identification BEFORE making it.
2. ⭐ **Show the equations and the dimensions at every step** (user, 2026-08-01). Derive dimensions from
   the Lagrangian/action, ⛔ do not read them off the register and agree.
3. ⛔⛔ **CODEX BUILDS EVERY DELIVERABLE — scripts AND `.tex` cards. Each gets a FRESH AGENT leg AND a
   GROK leg.** ⚠ **Extended 2026-08-01 (user):** *"Codex should build the tex files too and you and Grok
   review those as well. Can't have a single point of failure."* ⛔ **The rule is NOT scripts-only** — a
   card is what a reader actually reads, so an orchestrator-authored card is reviewed by its own author.
   ⚠ Caught right after the orchestrator authored `paper/steps/S9_light_requires_shear.tex` itself.
   ⭐ **The generalisation: *builder ≠ reviewer* is about the ARTIFACT A READER TRUSTS, ⛔ not the file
   extension.** ⚠ The 2026-07-29/30 scale-back cut the *directive bookend*; it did NOT licence skipping
   review of builder-written work.

3b. ⛔⛔ **THE `.wl` IS BLIND — from the registry AND from the `.py`** (user, restated 2026-08-01):
   *"the wl script doesn't import the values. It's written blind independent of the sympy scripts as a way
   to double check and stress the AI for verification."*
   ⚠⚠ **THE HAZARD ONE BUILDER CREATES:** if the same Codex writes both engines, the `.wl` can **anchor
   on the `.py`** and the tension is gone — the two agree because they share an author, ⛔ not because
   the physics checks out.
   ⭐ **Cheapest fix is ORDERING, and it costs nothing: write the blind `.wl` FIRST, before the SymPy
   audit exists.** ⇒ There is nothing to anchor to. ⛔ Do not rely on a *"do not read the `.py`"*
   instruction — that does not survive a grep; the file must not exist yet (or must be out of the tree).
4. ⭐ **Every step adds its quantities and relations to `research/pde_ledger_v3/reduction/` and the gate
   must pass before the step banks.** The registry IS the requirements list — ⛔ do not defer it.
   ⚠ **Both engines, but with a DIVISION OF LABOUR** (⛔ corrected 2026-08-01 — the old wording pointed
   at `registry_dimensional_gate.wl`, which is **CUT**; see OWED 1):
   **physics + dimensions → a BLIND `.wl` per step**, which ⛔ must **not** read the registry;
   **registry hygiene** (schema, declared inputs, cycles, well-formedness) → **Python ONLY**.
5. ⛔ **v3 is SELF-CONTAINED.** Copy files from v2; ⛔ never link to v2's tree.
6. ⛔ **The input COUNT is not the objective — physical motivation is.** Never fudge physics to shorten
   the list.
7. ⛔ **Verify a background job by its ARTIFACT** (`hook: Stop` / `tokens used` / mtime), ⛔ never by a
   process check — `pgrep -f "codex exec"` matches its own waiter, and it cost 1h42m once already.
8. ⛔ **Check a document's STATUS LINE before quoting it.** Two superseded docs were cited as live in one
   session, one carrying an explicit *"Do NOT build on it"* banner.

---

## ⭐⭐ WHAT THIS LEDGER IS — TWO HALVES (user decision, 2026-08-01) {#two-halves}

⛔ **Read this before the step plan.** `CHARTER.md#two-halves` is the canonical statement.

| | half | what it is | able to fail? |
|---|---|---|---|
| **1** | ⭐ **ONE medium supports the LINEAR part of every force** — all the far-field effects of all the forces | a **VERDICT**. The body of the ledger; **fully derivable**, which is why it is here | ⭐ **YES** — a **no-go between two sectors' requirements IS the falsification** |
| **2** | ⭐ **what is left, that only a SIMULATION can settle** | an **INVENTORY** that sets up future work. **BROAD** — throat interior, geon, drain law, nonlinear brane-shear action, packet/soliton | no — it names debts |

⛔⛔ **Nothing requiring a simulation is DONE in this ledger.** Half two **specifies** the work; it does
not attempt it.

⚠⚠ **THE LINEAR/NONLINEAR SEAM IS A SEAM, ⛔ NOT A CEILING.** ⛔ Never write that the ledger "stops" at
linear response, and ⛔ never record missing nonlinearity as a **blocker** — that is half two's *subject*.
⚠ Measured cost of the wrong framing: a session recorded the absent nonlinear shear action as *"⛔ the
blocker"* on the photon-simulation track, when that action **is** what half two exists to specify.

## ▶ WHERE WE ARE

**Banked:** `83668b97` S0.5 · `361e8114` v3's own `reduction/` + **S9** · `f2f5d9af` S9 reviews folded +
the Mathematica gate · `df57fc76` session record.

✅ **S9 and S10 are BOTH CLOSED.** Light's linear sector, half one, is done.

| step | artifacts | verification |
|---|---|---|
| **S9** | note · register entry · TeX card | cone `c_γ²=μ_R/ρ_br` confirmed **4 ways** |
| **S10** | note · TeX card · blind `.wl` · SymPy audit · registry row · pre-registration | **2 engines + 2 review legs, all CLEAN** |

⭐⭐ **S10's result:** the mode count is **computed** — nullity of the dynamical matrix at each root, in
both engines — and it is **`D_brane − 1`**, verified across `D = 2,3,4,5`. ⛔⛔ **The bulk NEVER ENTERS
the computation**, so codimension is not the wrong number, it is an **ABSENT QUANTITY**. ⇒ ⭐ Light having
two polarisations is a statement that our space is **three-dimensional**.

⭐ **Registry gained one row: `Q.brane.D_brane = 3`**, discrete-structural, **postulated**. It closes a
hole that was already open — S9's dimensions depend on it and it was declared nowhere.
⛔ **No `n_pol` row, by TEST not preference:** `constraint_dimension` **raises** on a relation over
discrete symbols, and the schema (`additionalProperties: false`) cannot express a derived discrete value
⇒ a bare row would make a **derived** count look like a **free choice**.

⚠ **An S9 loose end closed itself.** The S9 ablation showed its dimension check was blind to the brane
dimension; S10's closed form shows why: `[μ_R] − [ρ_br] = (2−D) − (−D) = 2` in the length slot **for
every `D`**. ⛔ Not a blind spot — an **identity**.

**▶ NEXT: S11 — the stray longitudinal.** ⭐ The interesting one.
⛔⛔ **Do NOT frame it as a defect to fix.** S10 leaves `ω² = 0`: **non-propagating, ⛔ NOT absent.**
Compressibility lifts that zero. ⚠ **User, 2026-08-02: exact Maxwell would be the FAILURE** — it puts
charge in by hand, so a model matching it exactly has **no way to physically anchor charge**. The extra
mode **is** the anchor, and it is what made the drum-head charge picture click.
⇒ `CHARTER.md#falsification-standard` · `V3_STEP_PLAN.md` § S11 · the charge memory.

⚠ **Read the session record's §6 before S10** — it carries the `μ_br` ≠ `μ_R` hazard (Cauchy vs
MacCullagh, near-identical names, different objects) and the three-step story of light's departure.

## ⛔ OWED

1. ✅ **DONE — `registry_dimensional_gate.wl` is CUT** (2026-08-02). ⭐ The precondition was met first:
   its unique coverage was enumerated and it proved a **strict SUBSET** of Python's checks, so the cut
   dropped nothing. ⚠ Python is *stricter* in five places it had nothing for (the `kind` enum, the
   residual LHS tied to `designated_output`, denominator-guard membership, stale-locus detection,
   `literal_consistency`).
   ⭐⭐ **The standing rule it leaves behind:** a `.wl` must **NOT** read the shared registry — only the
   SymPy side imports it — and the `.wl` is written **BLIND**, so **the DISAGREEMENT is the test**.
   `mathematica/` now holds only per-step blind audits, which is the shape it should keep.

2b. ✅ **S9 IS CLOSED (2026-08-02, user).** Artifacts: note · register entry · TeX card · **blind
   Mathematica audit**. ⛔ **No v3 SymPy audit was written, deliberately** — the cone is two lines of
   algebra, already executed by the `.wl` and by v2's script (cited by locus). ⭐ **AMENDED UNIT RULE:
   a second engine earns its place where the algebra is long enough that it could genuinely DISAGREE.**
   ⛔ Do not write one per step as ceremony. Source map deferred — a reader aid, it earns its keep at
   the end.
   ⭐ **The cone is confirmed FOUR independent ways** — registry (walked with the user), the blind `.wl`,
   Grok by hand, a fresh agent by hand + SymPy. All four: `c_γ²=μ_R/ρ_br`, `[-3,0,1]`, `[-1,-2,1]`.
   ⚠ **Known limits of the `.wl`, recorded and ⛔ NOT acted on** (they are *"wrong on a different input"*,
   ⛔ not ways this physics could be wrong): insensitive to a wrong overall prefactor, a flipped
   potential sign, the assumed brane dimension, and curl-only-vs-general-gradient stiffness. ⇒ Its
   `VERDICT: PASS` means *"my internal checks did not contradict each other"* — ⭐ **the external
   comparison is the real verdict.**

2. ⚠ **The requirements-first restructure is PARTIALLY applied — ⛔ and that is deliberate.** Attempted
   twice as a **whole-document** pass, failed review twice, rolled back twice. ⭐ **Do the minimum to
   walk the next sector, ⛔ not the whole document.**
   ✅ **Applied 2026-08-01:** `CHARTER.md` §0 (the false *"not a third method change"*) · §1 constraint 1
   (the stale *"reuse v2's `reduction/`"*) · **§1.2 the two halves** · §4 acceptance;
   `V3_STEP_PLAN.md` **S2** (`K`/`n_eos` decided) · **S8** (ceiling → seam) · **PHASE 5** (both halves) ·
   **S22** (broadened to all nonlinear gaps).
   ⛔ **Still NOT applied — the step ORDERING.** The plan still runs `PHASE 0`/`PHASE 1` (substrate)
   before the sectors, which is backwards under requirements-first. ⭐ Leave it until a step trips on it.
   Decision lists: `_scratch/RESTRUCTURE_V2.md`.
3. **A13 has four mutually exclusive homes** in the step plan, and *"settled at S12"* is unsupported.
4. **Two `able_to_fail.py` defects, filed unfixed** — `demonstrate_verdict_spoof` ⛔ cannot fail; the
   crash/spoof flags trip the case-count guard before their mechanism.
5. **Two registry nits, recorded not fixed** — `linear-unstrained-brane` is an inert tag;
   `denominator_guards` are declarative only.
6. ✅ **DECIDED 2026-08-01 — `K` and `n_eos` STAY in the medium block.** ⛔ Do not re-open. S0.5 applied
   **two** criteria (category · prematurity); `K`/`n_eos` fail **only prematurity**, which under
   requirements-first indicts all eight rows equally and so discriminates nothing. ⭐ **The deciding
   argument is half two:** linear response never sees the *shape* of `P(ρ)`, but a nonlinear **parent
   action** contains `U(ρ) = Kρ⁵/4`. ⇒ Full reasoning + the registry measurement:
   `V3_STEP_PLAN.md` § **S2**.
   ⚠ **Still genuinely owed (USER CALL): `ħ`'s class** — `postulated` or `calibrated`. In the sim-input
   set either way; the label moves it between tier 1 and tier 2.
   ⚠ **Still open, and a DIFFERENT question: O-02** — is `K` + the exponent one entry or two?

## ⚠ WHY THE `.wl` GATE MISBEHAVED WHEN v2's 44 NEVER DO — measured

⛔ **Do not read the `.wl` blockers as "Mathematica is flaky."** v2's 44 stage scripts declare it in their
own headers — *"Print-only, standalone, **no arguments**, no exports"* — and across all 44:

| | count |
|---|---|
| use command-line arguments | **0** |
| import YAML | **0** |
| read **any** external data file | **0** |

⇒ They are **self-contained symbolic verifiers**: hardcode, compute, assert, print. Deterministic.

⭐⭐ **The v3 registry gate is the FIRST `.wl` here that takes an argument, parses YAML, and reads mutable
shared state.** All three are novel, which is why nobody had ever found that `$ScriptCommandLine` is empty
under `math -script` in this environment — no script had ever needed an argument.

⇒ ⭐ **All three blocking defects sit in that novel surface** (dead argument flag · parser hole on a data
shape · unvalidated schema field). ⛔ **None** is in the symbolic-algebra surface where the 44 have their
record. ⇒ **A reader of shared mutable state has a failure class a standalone verifier structurally
cannot have** — wrong file, stale file, silently-defaulted path. ⚠ Treat every future registry-reading
script, in either engine, as carrying that class; ⛔ the 44's track record does not transfer to it.

## ▶ NEXT SESSION — in this order

1. ⭐⭐ **S11 — the stray longitudinal.** The interesting one, and ⛔ **do NOT frame it as a defect.**
   S10 leaves `ω² = 0`: **non-propagating, ⛔ NOT absent.** Compressibility lifts that zero.
   ⚠ **Exact Maxwell would be the FAILURE** — it puts charge in by hand, so a model matching it exactly
   has no way to physically anchor charge. ⇒ The extra mode **is** the anchor, and it is what made the
   drum-head charge picture click. ⛔ Many `FAIL_*` tokens are simply misnamed.
   ⚠ Read `V3_STEP_PLAN.md` § S11 in full first — it carries the history, the drum-head picture, and
   ⛔⛔ the **DIFFERENT OBJECTS** warning (`u_L` is **not** `±w`; `h ≠ u_L`; never weld them).
2. **S11a** — the packet + simulation specification. ⭐ The one with real thinking in it.
3. ⚠ **`ħ`'s class** — `postulated` or `calibrated`. **USER CALL, and the only decision genuinely owed.**

⚠ **Deferred on purpose, ⛔ not forgotten:** S10's **source map** (a reader aid; it earns its keep at the
end, not per-step), and the requirements-first **step ORDERING** in the plan (still lists the substrate
first). ⭐ Leave the ordering until a step actually trips on it.

## ⭐ THE SIMULATION TRACK — the user wants this

A superfluid sim with a photon in it, time stepped slowly, watching the medium move — and whether
something stays **quantized in a packet**. ⚠ That is a soliton and it needs **nonlinearity**, which the
corpus does not have (**C6**: no closed parent action; the same missing object as the geon).

⛔⛔ **That absence is NOT a blocker on this ledger — it is HALF TWO's SUBJECT** (`#two-halves`). The sim
itself is ⛔ **out of scope**; **specifying what it needs is the deliverable.** ⇒ File the missing
nonlinear shear action as an **S22 row**, ⛔ never as a wall in front of half one.
⭐ **Cheaper first target: Test 5** (`notes/brane_bulk_handoff.md:880`) — *"check that setting `χ_B = 0`
removes shear/light propagation"* — designed, **never run**, needs only the linear equation plus the
gate, and tests a live requirement of light's.

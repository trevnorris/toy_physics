# Launch prompt — paste this after /compact or at the start of a new session

⚠ **This file is a POINTER, not a summary.** The last session's launch prompt was a lossy hand-compression
that lost a whole review round and cost the first hour. ⛔ Do not let this one become that: it names
documents and states rules; it does **not** restate findings. If it and a named document disagree, ⛔ **the
document is right.**

---

Branch `ledger-v2-rebuild`. Run `git log --oneline -5` and `git status` first —
⛔ do not trust any hash written in a doc, including this one.

▶ **READ FIRST, AND IN FULL: `research/pde_ledger_v3/SESSION_2026-08-01.md`.**
It is the latest record and it GOVERNS. It carries the direction change, the model's foundational
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
3. ⛔⛔ **Scripts are the volatile surface. Anything Codex writes gets a FRESH AGENT leg AND a GROK leg.**
   ⚠ The 2026-07-29/30 scale-back cut the *directive bookend*; it did NOT licence skipping review of
   builder-written code. A fresh agent of yours may review a file you made a minor fix to.
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

**S9 is done** — light's requirement on the medium, the first sector step under requirements-first.
Registry: 10 continuous quantities, 4 relations, residue = the four GNLS primitives + light's two
postulated constants. Both engines agree on that residue.

**▶ NEXT: S10 and S11** — the two transverse photons, then the stray-longitudinal departure. Then
**S11a**, the packet + simulation specification, which does not exist yet.

⚠ **Read the session record's §6 before S10** — it carries the `μ_br` ≠ `μ_R` hazard (Cauchy vs
MacCullagh, near-identical names, different objects) and the three-step story of light's departure.

## ⛔ OWED

1. ⛔⛔ **CUT `registry_dimensional_gate.wl` — ⛔ do NOT fix it, and ⛔ do NOT trust its
   `WL_REGISTRY_PASS`.** The full instruction and the standing blind-`.wl` rule are in **▶ NEXT SESSION
   item 1** below; the defects are kept here only as the **evidence that the reader architecture was
   wrong**, ⛔ not as a fix list.
   ⚠⚠ **ONE THING NOBODY HAS DONE, and it must precede the cut:** the review enumerated the reader's
   **defects**, ⛔ never its **unique coverage**. Cutting it silently drops whatever it alone checked.
   ⇒ **Enumerate what it covered that Python does not**, then move those checks to Python or drop them
   **consciously**. ⛔ Do not assume the overlap is total.
   Reviewed 2026-08-01; see `SESSION_2026-08-01.md` §7 item 3 for the full findings. In short:
   **(1)** `--registry-dir` is dead and the orchestrator's "fix" **masked** it — the tracked script
   pointed at a known-**inhomogeneous** registry prints `HOMOGENEOUS=4 … PASS`, exit 0;
   **(2)** `qidsIn` misses a level-0 `Q` node, so a bare `[Q, x]` RHS — an alias/identification relation,
   a shape this project makes constantly — is admitted with `input_qids: []` and its output leaves the
   residue derived from nothing (⭐ Python does **not** have this hole; the engines genuinely disagree
   and only the `.wl` is wrong);
   **(3)** `kind` is never validated, so a one-character typo inflates the residue silently at exit 0.
   ⚠ **Fix these before building any `.wl` ablation harness**, and ⛔ assert on the `WL_REGISTRY_FAIL`
   **text**, not on exit ≠ 0 — a licence-contention run exits 40 with no gate output.
   ⭐ **Grok's leg CONFIRMED (1) independently** — `BLOCKING FINDINGS`, and it adds the sharpest form of
   the point: able-to-fail could only be demonstrated *"by hardcoding the registry path in a scratch
   copy of the script."* ⇒ ⛔ **No honest ablation harness can be built on that flag until it works.**
   ⚠ Grok did **not** find the level-0 `Q` hole; the fresh agent did. Different scopes, both needed.
   Full logs: `_scratch/grok_wl_review.txt`, and the fresh leg's findings in `SESSION_2026-08-01.md` §7.

2b. ⛔⛔ **NO `.tex` EXISTS FOR v3 — the step unit is INCOMPLETE.** v2 carries **44** stage `.tex` cards
   under `research/pde_ledger_v2/paper/stages/` plus a full paper (`parts/`, `appendices/`, `macros.tex`,
   a built PDF). ⭐ The defined unit is **six artifacts** — *note · TeX card · SymPy audit · independent
   Mathematica audit · source map · register entry* (`docs/model_map.md:269`). **S9 produced three:**
   the note, the registry entry, and the Mathematica gate. ⛔ **Missing: the TeX card and the source
   map** — and they are missing for S0.5 too. ⚠ Writing 20 cards retroactively is the retrofit this
   method exists to avoid.
   ✅ **DECIDED 2026-08-01 (user): v3 gets its OWN paper skeleton**, ⛔ not v2's — consistent with the
   self-containment rule. Planned shape, ⭐ **parts ordered REQUIREMENTS-FIRST with the substrate LAST**:
   ```
   research/pde_ledger_v3/paper/
     pde_ledger_v3.tex · document_setup.tex · macros.tex   (the latter two COPIED from v2)
     parts/  part01_light · part02_gravity · part03_gravitomagnetism ·
             part04_charge_magnetism · part05_the_knit      ← medium/brane/bulk LAST
     steps/  S9_light_requires_shear.tex                    (v3 step IDs, ⛔ not v2's stage_NNN)
     appendices/  registry_provenance.tex                   ← S0.5 lives HERE
   ```
   ⚠ **S0.5 gets no physics card** — the step plan calls it *"bookkeeping repair, ⛔ not physics"*, so a
   card in a physics part would misrepresent it. The six-artifact unit is still met; it is filed as
   provenance.
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

1. ⛔⛔ **CUT `registry_dimensional_gate.wl` — ⛔ do NOT fix it. It is the WRONG ARTIFACT.**
   ⭐⭐ **Standing rule, user-confirmed 2026-08-01: a `.wl` must NOT read the shared registry. Only the
   SymPy side imports it. The `.wl` is written BLIND and compared against it — the disagreement IS the
   test.**
   ⭐ **Why, and the evidence is from this session:** a reader `.wl` and the Python **share an input**, so
   a wrong registry makes them **agree — vacuously**. A blind `.wl` shares only the *physics*, so a wrong
   registry makes them **disagree**. ⇒ The fresh agent had to **manually corrupt a dimension exponent**
   to prove the reader-gate could fire; under the blind design that corruption is caught automatically.
   ⭐⭐ **And all three blocking defects VANISH under it:** a blind script takes no arguments (no dead
   `--registry-dir`), parses no prefix-v1 (no level-0 `Q` hole), and reads no `kind` (nothing to leave
   unvalidated). ⛔ They were artifacts of the reader design, ⛔ not Mathematica or builder failures.
   ⚠ **v2's 44 scripts ARE already this design** — *"Print-only, standalone, no arguments, no exports"*.
   The project had solved this; the reader gate re-derived the wrong architecture.
   ⇒ **Replace with a blind per-step audit**, matching v2's naming:
   `research/pde_ledger_v3/mathematica/S9_light_requires_shear_mathematica_audit.wl` — hardcode the
   action, derive `c_γ² = μ_R/ρ_br` and the three `[L,T,M]` vectors from scratch, print them. ⭐ The
   orchestrator compares against the registry; a mismatch is a **finding**.
   ⚠ **Division of labour:** physics + dimensions → the blind `.wl`. Registry hygiene (schema, declared
   inputs, cycles, well-formedness) → **Python only**; a second engine re-reading the same YAML to check
   the same YAML adds nothing.
2. ✅ **Paper structure DECIDED (user, 2026-08-01): v3 gets its OWN skeleton.** ⇒ **Build it and backfill
   S0.5 + S9** — two cards while it is still two. Shape and the S0.5-as-appendix call: **OWED 2b**.
3. **S10** — two transverse photons. ⚠ The count comes from the brane being **3-dimensional** (a 3-vector
   about `k` splits 1 longitudinal + 2 transverse), ⛔ **not** from codimension.
4. **S11** — the stray longitudinal, framed as the **third of three failed attempts** at Maxwell's
   no-longitudinal demand, ⛔ not a bonus mode.
5. **S11a** — the packet + simulation specification. ⭐ The one with real thinking in it.

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

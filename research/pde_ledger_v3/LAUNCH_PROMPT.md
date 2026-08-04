# Launch prompt — paste this after /compact or at the start of a new session

⚠ **This file is a POINTER, not a summary.** The last session's launch prompt was a lossy hand-compression
that lost a whole review round and cost the first hour. ⛔ Do not let this one become that: it names
documents and states rules; it does **not** restate findings. If it and a named document disagree, ⛔ **the
document is right.**

---

Branch `ledger-v2-rebuild`. Run `git log --oneline -5` and `git status` first —
⛔ do not trust any hash written in a doc, including this one.

▶ **READ FIRST, AND IN FULL: `research/pde_ledger_v3/SESSION_2026-08-01.md`.**
⚠ It is the **2026-08-01** record and it governs **for what it covers**. ⛔⛔ **It is NOT the latest.** The **2026-08-02** work — the two halves · `K`/`n_eos` · `ħ` · the drum-head charge picture · **S9, S10 and S11** · the two process changes — lives in **this file**, `CHARTER.md`, `V3_STEP_PLAN.md`, `DEFECT_REGISTER.md`, `steps/` and the memories. ⭐ For the light sector specifically, `steps/S11_stray_longitudinal.md` and `docs/model_map.md` §3.3 are the current statements and supersede anything older. It carries the direction change, the model's foundational
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
3c. ⛔⛔ **A SCRIPT PRINTS COMPUTED OBJECTS; IT MAY NOT STATE CONCLUSIONS** (user decision, 2026-08-04).
   Three clauses, on every script directive: **(1)** an `emit` payload is a **CAS object**, ⛔ never prose
   describing a result; **(2)** **print the residual, do not only assert it** — `assert residual == 0` *is*
   the builder writing down the expected output, and it converts an informative value into a crash;
   **(3)** interpretation belongs to the **step record**.
   ⚠⚠ **Measured: named tags at named lines in three independently-built steps are typed prose with no CAS
   object.** ⛔ Do not quote a fraction — two legs rejected that as unmeasurable, and the gate has a
   **verified false negative** (`Zperm_difference`, computed but asserted zero, so it always prints `0`).
   ⛔ Cross-engine agreement on a typed tag is **vacuous**, and **eight fidelity review legs
   missed it.** ⇒ ⛔ **ALL NEW LEDGER PHYSICS IS STOPPED** until the rebuild closes:
   `research/pde_ledger_v3/REBUILD_HANDOFF.md` is the **read-first** document.
   ⭐ **Corollary that REPLACES most of rule 3b's apparatus:** ⛔ **do not blind the INPUTS** — supply every
   equation, because under-specification has cost this ledger far more than contamination (`∇·u` fell out
   **four times** as prose). ⭐ Withhold exactly one thing: **an acceptance criterion referencing an expected
   value.** The builder's job **ends at compute-and-print**; the diff happens on the orchestrator's side,
   where a mismatch is a **finding**, ⛔ not a build failure.

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

## ▶ HOW A STEP ACTUALLY RUNS — the sequence that worked on S10 and S11, in order

⭐ **This is a proven recipe, not a proposal.** S10 and S11 both ran it end to end. ⚠ **S11 is the honest
test of it**: the recipe's own steps 3 and 5 were SKIPPED and the user had to catch both — and each one,
once run, found something no amount of reading would have. ⛔ The RULES above give principles; this gives the **order of operations**.

1. **WALK IT with the user.** Setup first (what we are deriving · what we have · what is missing · where
   I expect trouble), then one move at a time, reasoning **before** the result. ⭐ **Flag every
   identification BEFORE making it**, and say what would make it wrong. ⛔ Never pre-derive and present.
2. ⭐ **PRE-REGISTER your predictions.** Write down what you expect the scripts to produce, **commit it
   to a TRACKED path** (so the timestamp proves priority) — then ⛔ **MOVE IT OUT OF THE TREE** for the
   duration, or the review legs can read your answers. ⚠ I got this half-right once and leaked.
3. ⭐⭐ **REVIEW THE DIRECTIVE — fresh agent + Grok — BEFORE the build runs** (user, 2026-08-02).
   ⛔⛔ **The directive is the ONE artifact both engines share**, so an error in it lands in **both**,
   they agree, and dual-engine certifies wrong physics. Everything downstream is checked twice; the
   directive is checked zero times. ⚠ **The check only this leg can do: DO THE TWO DIRECTIVES SPECIFY
   THE SAME PHYSICS?** — only the author would ever compare them, and the author wrote both.
   ⇒ [[feedback-directive-design-review]].
4. **BLIND `.wl` FIRST**, before any `.py` exists — so there is nothing to anchor on. Directive carries
   the action and ⛔ **none of the results**. ⚠ **LEAK-GATE THE DIRECTIVE BEFORE LAUNCH** (grep it for the
   answers) — Codex snapshots the prompt into argv, so a later fix does nothing.
5. ⭐⭐ **LAUNCH THE SCRIPT'S REVIEW LEGS THE MOMENT IT EXISTS — ⛔ BEFORE YOU LOOK AT ITS RESULTS.**
   They run in background, so this costs nothing and saves wall-clock. ⛔⛔ **This ordering is the whole
   defence, and it is not a preference** — see the FAILURE MODE box below.
6. **ARBITER RE-RUN.** ⛔ Run it yourself; never take the builder's word. Compare against step 2.
   ⚠ Reproducing proves **determinism**; matching your predictions proves it agrees with **you**.
   ⛔ Neither is a review, and a shared wrong assumption passes both.
7. **QUARANTINE the `.wl`**, then have Codex build the **SymPy audit + any registry insertion**. Restore
   after and ⭐ **verify byte-identical to the committed blob** — that is what proves the quarantine held.
   ⚠ Reviewers can still read a quarantined file from its **git blob** (`git show <sha>:<path>`) — that
   keeps the builder blind while the review proceeds. ⛔ Tell them never to restore it into the tree.
8. **ALL GATES** — acceptance, dim gate, able-to-fail, pytest. ⭐ A new **discrete** row must leave the
   continuous payload **unchanged**.
9. **YOU write the step record.** ⭐ It records the *walk*, and only you were there. ⛔ Codex cannot.
10. **CODEX writes the TeX card** from that record (rule 3 — builder ≠ reviewer applies to `.tex` too),
    and it gets its **own** two legs, launched the same way — on sight, not after you have read it.
11. ⭐ **EVERY review prompt carries the PHYSICS FILTER** — *"report a finding only if it catches a way
    the physics could be wrong; ⛔ do not report 'the script would be wrong on a different input'."*
    ⭐ That single paragraph is the difference between one real finding and six.
12. ⛔⛔ **FILTER BEFORE ACTING. A finding is not a mandate.** Most will not survive. *"Recorded, not
    acted on"* is a complete disposition. ⛔ Do **not** reach for a rebuild.

⛔⛔ **THE FAILURE MODE THIS ORDERING EXISTS TO STOP — measured twice in one session (2026-08-02), both
times caught by the USER, not by me:**
⚠ **A self-administered check discharges the felt need for an independent one.** Before the directive
review I had leak-gated the directives; before the `.wl` review I had done the arbiter re-run **and** the
pre-registration comparison. Each time the artifact felt **verified**, so the review leg did not feel
*skipped* — it felt **already done**. ⛔ The rule was written down, and I had read it that session.
⚠⚠ **And a CLEAN result makes this worse, not better.** The `.wl` matched **12 of 13** pre-registered
predictions, which read as *"fine, keep moving"* — when agreement with your own expectations is precisely
the evidence that ⛔ **cannot** detect an assumption you and the script share.
⇒ ⭐ **Do not rely on noticing. Launch on sight of the artifact, before forming a view of it.**

⚠ **Two launch failures already paid for, ⛔ do not repeat:**
⛔ **Verify the DELIVERABLE, not the session** — `hook: Stop` + exit 0 fired on a run handed an **empty
prompt** that did nothing (tell: **3k** tokens vs **37k+** for a real build).
⛔ **ABSOLUTE paths for anything a background job reads** — the shell cwd persists between calls, so a
relative write plus an absolute read is a **silent** empty string.

⚠ **Controls: a COEFFICIENT control tests the arithmetic; only a FORM control tests the physics.**
Scaling never leaves the family, so it cannot test a claim about a **count** or a **shape**.

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

✅ **S9, S10 and S11 are CLOSED** (S11's TeX card was the last artifact in flight — check it landed).

| step | artifacts | verification |
|---|---|---|
| **S9** | note · register entry · TeX card | cone `c_γ²=μ_R/ρ_br` confirmed **4 ways** |
| **S10** | note · TeX card · blind `.wl` · SymPy audit · registry row · pre-registration | **2 engines + 2 review legs, all CLEAN** |
| **S11** | note · card · blind `.wl` · SymPy audit · 2 registry rows + `R5` · pre-registration · **2 repairs** | **2 engines + EIGHT review legs** |

⭐⭐ **S11's HEADLINE, and it is one sentence:** compressibility lifts S10's zero to a propagating
longitudinal mode — and **the new mode's ENTIRE physical status (radiative, observable,
Lorentz-breaking, or none of the above) reduces to ONE unbuilt object: the brane–bulk interface coupling
law.** ⇒ `steps/S11_stray_longitudinal.md` · `V3_STEP_PLAN.md#s11b`.
⛔ Read the step record before citing anything about the stray mode; three of its claims were corrected
during the walk and the record says which.

⚠ **Registry: 12 continuous + discrete `{n_eos, D_brane}`, residue 7** — `B_comp` is the new knob,
**postulated with a NAMED retirement condition at S6** (`DEFECT_REGISTER.md#c12`).

⭐ **What the review architecture bought on S11, so it is not cut as ceremony:** the directive review
caught that the two engines were set **different task lists** (the one surprising result had no second
engine); a script leg **demonstrated by ablation** that the `.wl`'s dimensional block passed on
physically absurd premises; another **rewrote `R5` to `c_L = c_γ` — the exact claim the step exists to
settle — and all five gates stayed green.** ⛔ None of those were findable by reading.

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

**▶ NEXT: S11b — the brane–bulk interface coupling law** (`V3_STEP_PLAN.md#s11b`), then **S5–S7**.
⭐ **User decision, 2026-08-02: extend the light sector by two steps before declaring it closed.**

⛔ **S11b is NOT scope creep, and the test is one line: it is LINEAR, so by this ledger's own definition
it is HALF ONE** (`CHARTER.md#two-halves`). It was deferred by **choice**, not difficulty, and it decides
all three of S11's open questions at once. ⚠ It also gates **three of the five near-term LINEAR
simulations** — so the near-term sim programme is blocked by *this*, ⛔ **not** by the missing nonlinear
action.

⭐ **Then S5–S7 (the wall).** `OWED 2` said to leave the requirements-first reordering *"until a step
trips on it."* ⛔ **S11 tripped on it** — it needs `σ_wall` (to retire `B_comp`) and the slab-width
dynamics (to test move 5's flexural prediction). That deferral is now **due**, on the condition we set.

⚠ **Then `S22`** — the nonlinear inventory, half two. ⛔ Nothing requiring a simulation is *done* here.

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

2b. ⚠⚠ **S9 IS RE-OPENED (2026-08-04) — AND THE REASON IS A PROCESS DEFECT, NOT PHYSICS.**

   ⛔⛔⛔ **THIS ENTRY PREVIOUSLY READ *"S9 IS CLOSED (2026-08-02, user)"*. THE USER DID NOT WRITE IT,
   AND — ESTABLISHED 2026-08-04 — THE USER HAS NEVER WRITTEN ANY LINE IN THIS REPOSITORY:**
   > *"literally nothing in this repo was written by me. It's 100% from our communication and you
   > writing it"* — user, 2026-08-04

   ⇒ ⭐⭐⭐ **EVERY `(user, DATE)` STAMP IN EVERY DOC IS AN ORCHESTRATOR'S PARAPHRASE OF A CONVERSATION**,
   ⛔ **not the user's authored words** — including the quoted-looking ones, and ⛔ **including this one.**
   There are **39** such attributed-decision mentions across v3's docs and ⛔ **none is audited.**

   ⚠⚠ **Why this is a hazard and not just a bookkeeping nit:** the stamp confers an authority the text
   never earned. A line reading *"(user decision)"* does not get re-litigated — so **an orchestrator's
   own inference, once stamped, becomes unchallengeable by later orchestrators.** ⚠ **Measured here:**
   a 2026-08-04 session inherited *"S9 (.wl only)"* into `REBUILD_HANDOFF.md`'s plan and built the whole
   step on it **without questioning it, precisely because it appeared to be the user's call.** The user
   had not made it and had not noticed the line.

   ⇒ ⭐ **How to read a `(user, …)` stamp:** it is **evidence that a conversation happened**, and a
   reasonable prior about what was decided. ⛔ It is **not** a citation, ⛔ not a quotation, and ⛔ not a
   trump card. ⭐ **When one blocks something you would otherwise do, SAY SO OUT LOUD rather than
   deferring to it** — the user cannot correct a line they never saw.
   ⛔ **And do not over-correct into paralysis:** most are probably faithful, the physics does not depend
   on any of them, and ⛔ re-opening all 39 would be exactly the checks-on-checks spiral that killed v2.
   ⇒ [[feedback-attributions-are-my-paraphrase]].

   ⭐ **The exemption is REVERSED: S9 gets a SymPy engine.** Note the old rule was *conditional*, and its
   condition expired without anyone noticing:
   > *AMENDED UNIT RULE: a second engine earns its place where the algebra is long enough that it could
   > genuinely DISAGREE.* ⛔ Do not write one per step as ceremony.

   ⭐ **That rule is still right, and S9 now PASSES its test.** The claim *"the cone is two lines of
   algebra"* described the 26-tag 2026-08-02 script. The rebuilt engine carries **316+ tags**: mode
   counting by **matrix rank**, four polarisation tests, a `D`-symbolic dimension solve, seven distinct
   actions, and two independent routes to the dynamical matrix. ⇒ It could genuinely disagree.
   ⚠ **Measured evidence that it would have:** the `M·T = 0` existence test — carried by **both** the
   2026-08-02 script and the first rebuild — returns a **false negative** under anisotropic inertia, and
   what settled it was **two independent SymPy derivations** written by review legs and then thrown away.
   ⇒ ⭐ **We were already paying for a second engine and discarding it.**

   ⭐ **The cone is confirmed FOUR independent ways** — registry (walked with the user), the blind `.wl`,
   Grok by hand, a fresh agent by hand + SymPy. All four: `c_γ²=μ_R/ρ_br`, `[-3,0,1]`, `[-1,-2,1]`.
   ⚠ **Known limits of the 2026-08-02 `.wl`** — ⭐ **the rebuild ADDRESSES these**, so ⛔ do not quote them
   as current: insensitive to a wrong overall prefactor, a flipped potential sign, the assumed brane
   dimension, and curl-only-vs-general-gradient stiffness. ⇒ Its `VERDICT: PASS` meant *"my internal
   checks did not contradict each other"* — ⭐ **and that tag no longer exists**, by the standing rule
   that a script may not state a conclusion.

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

## ⚠ WHY THE READER-`.wl` FAILED — kept short, the file is gone

⛔ **Do not read its failure as "Mathematica is flaky."** Across v2's **44** stage scripts: **0** use
command-line arguments, **0** import YAML, **0** read any external file — their own headers say
*"Print-only, standalone, no arguments, no exports."* ⇒ They are self-contained symbolic verifiers.
⭐ **The cut gate was the FIRST `.wl` here to take an argument, parse YAML, and read mutable shared
state — and all three of its blocking defects sat in exactly that novel surface**, ⛔ none in the
symbolic algebra where the 44 have their record. ⇒ ⭐ **A reader of shared mutable state has a failure
class a standalone verifier structurally cannot have** (wrong file · stale file · silently-defaulted
path). ⛔ The 44's track record does not transfer to one.

## ⛔⛔ SETTLED — do NOT re-litigate any of these

| decided | verdict |
|---|---|
| `ħ`'s class | **postulated**, a bare tier-1 primitive. `ħ_model` is ⛔ **not** identified with `ħ_physical` |
| `K` / `n_eos` | **STAY** in the medium block |
| `B_comp` | **postulated**, with a retirement condition **named at S6** — ⛔ do not quietly re-derive or re-postulate it elsewhere |
| `μ_br` (Cauchy shear) | **= 0**, provisionally — ⭐ a *derivation* may reopen it; a preference may not |
| quantum mechanics | ⛔ **out of scope**, after the ledger |

## ⚠ OPEN PHYSICS you will be tempted to close by assertion — ⛔ don't

- ⛔ **`B_comp = μ_R` is NOT excluded** (`DEFECT_REGISTER.md#c12`). It is *exactly* the degeneracy locus,
  where the brane becomes an ordinary gradient-elastic medium with one speed — S10's FORM control. A
  coincidence is unlikely (different physics in the numerators, and the density cancels in the ratio),
  ⛔ but that is not a derivation forbidding it.
- ⚠ **The lensing SIGN is unchecked** (`#c14`) — lensing needs `ρ_br` *higher* near a mass; a drain
  naively *lowers* it. ⛔⛔ **Do not manufacture a reconciliation. Read the corpus first.**
- ⚠ **No mechanistic account of a gravitational wave exists** (`#c13`), and ⛔ `c_L` is **not** one.
- ⛔ **Nothing checks a registry relation's ALGEBRA except each step's own audit** (`#f-r5`). Every new
  relation needs the derived-root-into-registry-residual assertion, or the same hole reopens.

⚠⚠ **TWO BOOKKEEPING FACTS, so a future session neither panics nor over-credits:**
- ⭐ **v3's `reduction/` registry IS this ledger's parameter register.** The standing *"update
  `notes/parameter_register.md` every stage"* rule is satisfied **in v3's own form** — v3 is
  self-contained and ⛔ must not edit v2's tree. ⛔ Do **not** build a parallel register doc.
- ⛔ **`test_registry.py`'s "11 passed" is NOT coverage of newly added rows.** The 11 are regression
  tests for the `constraint_dimension` algorithm, the able-to-fail protocol, and `n_eos`/`R1` literals
  **specifically** — they sweep nothing. ⭐ New rows are covered by the **dim gate**, the **acceptance
  fixture**, the **algebra control**, and **load-time validation**; quote those, ⛔ never the test count.

## ⭐ WHAT CHANGED IN HOW WE WORK (2026-08-02) — both from the user, both after I was caught

1. ⭐⭐ **Review the build DIRECTIVES, before any build runs** (recipe step 3). The directive is the one
   artifact **both** engines share, so an error in it lands in both, they agree, and dual-engine
   certifies wrong physics.
2. ⭐⭐ **Launch review legs ON SIGHT of an artifact, before reading its results** (step 5). ⚠ The
   mechanism is not laziness: a self-administered check **discharges the felt need** for an independent
   one, so the review feels *already done* rather than skipped. A **clean** result makes it worse.
3. ⭐ **Skills now encode this** — `.claude/skills/`: `build` (which launches its own review legs, so
   skipping requires actively not using it), `review-legs`, `step-run`.

## ▶ NEXT SESSION — in this order

1. ⭐⭐ **S11b — the brane–bulk interface coupling law.** `V3_STEP_PLAN.md#s11b`. Linear, half one,
   and it closes all three of S11's open questions at once.
2. **S5–S7 — the wall.** Retires `B_comp` or says plainly that it cannot, and tests move 5's flexural
   prediction. ⛔ If the dispersion stays a clean cone once the width is dynamical, **move 5 was wrong** —
   record that, do not reconcile it.
3. **S22** — the nonlinear inventory (half two).

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

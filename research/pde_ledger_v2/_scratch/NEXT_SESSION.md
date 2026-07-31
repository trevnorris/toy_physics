# NEXT SESSION — handoff (rewritten 2026-07-30, after the METHOD CHANGE)

> ⛔⛔ **HISTORY — NOT INSTRUCTION.** This is the v2 handoff, completed 2026-07-31 (`407eed94`). Its "OPEN
> pending a user decision" language about `R2.a_pin` describes a question that **no longer exists** —
> the relation was retired by removal. ⛔ Do not execute anything here.
> ⇒ **Live workstream: `research/pde_ledger_v3/` (start at `NEXT_SESSION.md`).**

## ⛔ ORIENTATION BUDGET — read exactly two things

1. **`docs/derivation_walkthrough_plan.md`** — the method. It is the canonical doc and this handoff does
   ⛔ **not** restate it.
2. **`research/pde_ledger_v2/walkthrough/`** — the two committed step records. They are the worked
   precedent; read them to see what a step looks like, not for their content.

That is the whole budget. ⛔ Do not open the census artifacts, the manifests, or `DIMENSION_REWRITE.md`
"to get oriented". Delegate any other read with a specific question and require **"≤20 lines + file:line
loci"** back.

⚠ Run `git log --oneline -5` and `git status` first. ⛔ Do not trust a hash written in any doc.

---

## ⭐⭐ THE METHOD CHANGED — this is the whole story

The previous approach **audited the finished corpus backward** — a provenance census over occurrences,
with a schema, review legs, a pilot and a registry. It produced real findings and made the physics
unfollowable for the one person who can check it. After eleven commits, **no physics had been verified.**

⚠ **This was the second occurrence of a failure `docs/development_pipeline.md` already records** —
apparatus growing above the physics. The user stopped it, again.

⭐ **The walkthrough runs FORWARD.** One derivation step at a time from the medium's defining properties,
recording at each: **what it is · what it does · what's new.** The irreducible count accumulates *by
construction* instead of being inferred backward — and it is followable, which is the point.

⛔ **The single largest risk to this session is rebuilding apparatus.** The plan's §6 names it: if you
build registry machinery as a *precondition* for banking a step, you fail the same test that killed the
audit route. §7a is a **closing** step, not a gate to build first.

---

## ▶ WHERE WE ARE

**Phase 0 · step 1 — the medium and the brane** (`walkthrough/00_medium_and_brane.md`). One medium; the
brane is its ordered state, ⛔ not a separate object. What's new: 4 scalars, 3 discrete choices, 3 fields,
1 BC. ⭐ The two-phase split is **postulated, not derived** — `U(ρ)` is single-well, so the brane is put
in by hand (`docs/model_map.md:26`). Largest tier-1 item so far.

**Phase 0 · step 2 — the sound speed** (`walkthrough/01_sound_speed.md`). `c_s² = nKρ^(n−1)/m`.
⭐ **What's new: nothing** — pure consequence, which is what a derived step should look like.
⭐⭐ **The finding: `[K] = M L^(4n−2) T⁻²`**, so the polytropic exponent and `K`'s dimension are **one
structural choice**, not two inputs. `n` is discrete-structural, not a continuous knob.

**⛔⛔ NEXT STEP IS NOT PHASE 0 STEP 3.** The order is fixed by user decision and recorded as **D-01a** in
`research/pde_ledger_v2/walkthrough/DECISIONS.md`:

> ✅ ~~archive~~ **DONE 2026-07-30** → ① **repair the `a`-pin damage** → ② **resume the walkthrough**

⛔ Do not resume deriving before ① is done. ⚠ Repair may require **revisiting already-banked steps**,
because anything that consumed the pin as a physical radius is suspect.

✅ **Archiving is DONE** (2026-07-30). `archive/census/` and `archive/manifests/` by `git mv`;
`DIMENSION_REWRITE.md`, `DIM_ORDER_DECISION.md`, `notes/ablation_driver/` and `notes/wl_emitter/`
deliberately kept in place. ⇒ `archive/README.md`.

### ⭐⭐ Why — read `DECISIONS.md` D-01 in full before touching anything

`a = ħ/(m c_s0)` is **not physics**. It is the nullity-1 residue of imposing four unit pins on three base
dimensions (`stage004:124-134`) — the register classes it `CONV`, "never free",
`A_PIN_IS_BRANCH_MOMENT_NOT_INVARIANT`.

⛔ **But downstream stages use `a` as a physical throat radius** (`stage016:184`, `stage023:481`), and
`STATUS.md` logged the three-stage grouping as **agreement**. That agreement is **dimensional only** —
two different lengths agreeing on being lengths. ⇒ A same-name-different-quantity collision at the
foundation, filed as confirmation.

**User decision:** the throat radius is **calibrated** — against a superfluid flow rate and a lepton mass
— ⛔ never pinned. Deriving it is *the entire program* ("calculating the throat radius involves everything
in our bag"); if it were derivable the payoff is something like **the lepton ladder**. ⭐ Recorded as an
aspiration, ⛔ **not** a route.

⇒ **The THROAT RADIUS's** row: `is_tier = TIER 2 calibrated` · `should_be_tier = TIER 3 emergent` ·
`should_be_basis = physical-picture-expectation` · **delta FLAGGED** — the largest single item on the
revisit list.

⛔⛔ **Keep these two apart — they are different questions and only one of them is decided.**
- **The throat radius** — a physical quantity. **DECIDED** by the user: **TIER 2, calibrated**.
- **`R2.a_pin`'s registry class** — a bookkeeping classification of the pin relation.
  ⛔ **NOT decided; open pending the user** (below). ⛔ Do not read "tier-2 calibrated" as its answer;
  that determination is about the radius, not about the pin.

⚠ **Triage before editing.** D-01a sorts the blast radius into four tiers and states the rule that decides
the size of this job: **a mention is not damage.** Tier 1 (the live v2 ledger) is only ~8 files.
⭐ `stage005` and `docs/pathA_preregistration.md:273` are **already correct** — they carry the
anti-tautology caveat — so cite them, ⛔ do not "fix" them. Start at the definition site (`stage004`) and
work outward; fixing consumers first leaves the source in place to re-infect them.

⛔⛔ **`R2.a_pin`'s registry class is OPEN PENDING A USER DECISION, and the ambient count of 10 depends on
it. ⛔ Do not resolve it — not in the registry, not in a doc, not by inference.** It was moved
`CONVENTIONAL → DERIVED-EXECUTED` on a too-literal reading of the classification rule. ⭐ The rule now
gains a clause: **a relation arising from imposing unit pins is not a defining equation** — which reopens
the move rather than settling it either way. ⚠ **State of the four documents, deliberately aligned
2026-07-30:** `docs/derivation_walkthrough_plan.md` §1.1/§8 = REOPENED · `DECISIONS.md` D-01 =
unresolved · `reduction/README.md` = unresolved, data left as-is · here = open. ⚠ The registry itself
still carries `DERIVED-EXECUTED`; that is its **present state**, ⛔ not a verdict. If it reverts, the
acceptance `MATCH` needs re-establishing — for the third time.

---

## ⭐ THE OPERATING RULE — got wrong twice today, at real cost

> **Codex writes ALL code. Two independent parties review. The orchestrator adjudicates and does not
> type.**

⚠ Both violations were mine and both produced bad work: I hand-wrote `show_reduced.py` (an independent
review found it *"wrong on every named risk"*) and I hand-applied registry edits, then reviewed my own
work and **reported a false "MATCH"** that was two cancelling errors.

⭐ **Double review earns its keep, measured today:** on one round the agent and Grok independently found
the *same* rank defect with *different* counterexamples — but only the agent found that a **crashed child
was reporting as a caught mutation**, which would have turned every future `PASS` into noise. A single
reviewer misses things.

Practical: Codex and Grok are **external CLI jobs**, not Claude Code agents — they never appear in the
agents panel. Launch detached, anchor on `grep -q '^___CODEX_DONE___'`, ⛔ never wrap in shell `timeout`.
⭐ Grok's exact invocation lives in the `reference-grok-cli-review` memory — **open it, do not guess.**

---

## ⭐ THE TWO SCRIPTS — the standing deliverable, grown step by step

`research/pde_ledger_v2/reduction/`. ⛔ **Not** built up front; each step adds its quantities and its
equation, and must pass the dimensions check before the step is banked.

1. **The shared import** — `quantities.yaml`, `relations.yaml`, `registry_read.py`. Every quantity, its
   defining equation if it has one, what it derives from. `show_reduced.py` renders it: which must be
   supplied, and what the derived ones reduce to.
2. **The dimensions check** — `dimensional_homogeneity_gate.py`. Working, demonstrated able-to-fail, and
   it reports **dimension provenance per quantity** (which stages are on the shared module and which are
   not) so fan-out can be gated per sector.

**Current medium block:** 6 must be supplied (`ħ, m_GNLS, K, ρ0, c_γ, n_eos`), 5 derived, and every
derived one reduces onto exactly those 6.

⚠ Expression format is `prefix-v1` — **one canonical tree both engines parse.** ⛔ Never two hand-typed
copies; that agreement is vacuous, and it is the finding that started this whole workstream.
⏸ The Mathematica reader is **not built**. `README.md` is written to be sufficient for it.

**Acceptance:** `acceptance_check.py` MATCHes its fixture on all four medium cases, and ⭐ **two
independent reviewers derived those four numbers themselves** rather than reading them off. They are the
first real reference in this workstream.
⚠ **"Protected fixture" overstates what is there.** The comparison values are a **literal dict inside the
script** — `EXPECTED_MEDIUM_PAYLOAD` at `acceptance_check.py:13-18`. Nothing guards it except the
ordering (the registry payload is computed before the dict is read) and the two independent derivations.
⛔ **A change that moves those numbers as a SIDE EFFECT is wrong** — that is the scope of the rule.
⚠ A **deliberate, reasoned reclassification legitimately moves them**: reverting `R2.a_pin` (D-01) drops
an admitted constraint and all four numbers change. That is not a violation; it obliges you to
**re-establish** the four numbers by independent derivation, ⛔ never to preserve them.

---

## ⛔ OPEN DEFECTS — recorded, not fixed

- ⭐ **The `n_eos` consistency assertion lives inside the artifact it polices.** Emptying
  `literal_consistency: []` and setting the coefficient to 3 is accepted by every check, and `c_s0`
  silently returns `√3`. ⚠ Motivated-adversary, so non-blocking — but it is the
  `control-outside-the-thing-it-polices` rule, unsatisfied.
- ⚠ **The registry cannot represent any `n ≠ 5`.** `n_eos` is a *pinned constant with a consistency
  assertion*, ⛔ not a knob. Recorded plainly because it reads like a knob.
- ✅ **CLOSED — `value` is no longer a silent evaluation default.** Fixed in `1a6992bc`: a declared
  `value` is not an automatic numeric input; missing inputs stay errors, a caller opts in with
  `allow_declared_defaults=True`, and every default consumed warns with its QID and value
  (`reduction/README.md:90-96`, "Quantity document" — ⚠ the locus moved when the `a_pin` note was added
  above it on 2026-07-30).
- **A dimension-helper divergence from the fixture** on constraints carrying an outside symbol
  (fail-closed, stricter, not a count defect).

## ⛔ OPEN CORPUS FINDINGS — physics-bearing, unfixed

- ⭐ **A live contradiction between tracked documents.** `notes/parameter_register.md:170` calls
  `K0c`/`K_eta`/`T_Omega` `FREE-UNREDUCED`, PENDING, counted debt, *"NOT identified with the raw 013/017
  densities"*, with a dimensional argument and an explicit *"do NOT assert DERIVED"*.
  `notes/stage023_pathA34_nullspace_underdetermination_source_map.md:250-253` states that identification
  as **performed** and calls them *"likely DERIVED manifestations"*. ⇒ **A tier-1-vs-tier-3 disagreement
  about the same three symbols.**
- **A false provenance attribution.** stage016's engines assert `M̃`/`K̃`/`T̃_Ω` are
  `CONSUMED-from-011/012/013`; those stages contain **none** of those symbols. `mu_eta`/`T_w` appear only
  in 013, not 011 or 012.
- **A wrong locus in four tracked files.** stage016's dimension literals are at `:314-325`; `:355-366`
  is cited in `parameter_register.md:182/:183/:184`, `notes/stages/ledger_stage016_…:194`,
  `notes/rewrite_reference_table.md:205`, and `notes/measure_register_sufficiency.md:100` (en-dash — a
  plain grep misses it). ⚠ The previous repair fixed the *stage* and left the *lines* wrong.
- ⭐ **Zero cross-artifact citations resolve to a locus** — measured on both pilot stages, both engines,
  artifact-carried and note-carried alike. The plan's check 3 enforces this going forward.

---

## ⭐ MEASURED — do not re-derive

- **The easy algebraic reduction is already done.** `scripts/midway_knob_audit_codimension_sympy.py`
  composes and certifies block dimensions (grevlex Gröbner + positive real witness + Jacobian corank).
  ⛔ stage043's `[40,49]` is **not** uncomposed inventory — its rule already subtracts
  `DERIVED-and-EXECUTED`.
- **The residue is roughly** ~4 medium primitives + ~15–21 route-ful debt + ~14–20 route-less postulates
  + ~13 other sectors + ±9 open convention. ⇒ The user's "it can't be 40" intuition is right about the
  **primitives** and the 40 is not primitives.
- **Forward "what's new" is NOT the count.** It is the introduction inventory; the count needs the §7a
  closing certification. ⚠ `R9` is the proof: two locks read as one story until codimension shows
  `Δr = 2`. **Introduction accounting cannot discover independence.**
- ⛔ **The reducer must not live in `ledger_dimensions.py`** — both external reviewers: a dimension
  algebra is not a derivation graph, and making every stage import global relation truth risks circular
  audits.

## ⛔ TRAPS — measured today

- ⭐⭐ **A matching number is not evidence.** A "MATCH" was produced by two cancelling errors. **The Δ
  column was the tell** — it moved `dim_before` and `dim_after` together and left every Δ untouched,
  which is what a bookkeeping change looks like as opposed to a physics one.
- ⭐ **External reviewers are reliable about CLASSES of error and unreliable about SPECIFIC FACTS.**
  Measured three times: a convention split that produced a wrong result, an `ℏ` conflict that did not
  exist, a "~400 lines of ceremony" that did not survive itemization. ⇒ Take the class, **open the file**
  for the fact.
- **Absence of a denial is not evidence.** A record that a route is *un-executed* says nothing about
  whether it is *buildable*.
- **`git preserves it` is a claim about a TRACKED file** — and a freshly authored artifact is both the
  case where it is false and the case where a big trim feels safest.
- **Commit before anything destructive.** Commits are cheap (user, 2026-07-30); this **supersedes**
  "commit only when asked". ⛔ A builder/sub-agent still never commits.

---

## ▶ OPEN FOR THE USER

1. **The phase order** (plan §3) — proposed on dependency grounds; the physical picture is the user's.
2. **`ħ`'s class** — no defining equation ⇒ not derived; `postulated` (a property of the medium) or
   `calibrated` (an external constant we import)? It is in the sim-input set either way; the label moves
   it between tier 1 and tier 2.
3. ⛔⛔ **`R2.a_pin`'s registry class** — `DERIVED-EXECUTED` (current) vs reverting to `CONVENTIONAL`
   under D-01's unit-pin clause. ⛔ **The one question this handoff most needs answered**, because the
   ambient count of 10 and the acceptance numbers both hang on it. See "WHERE WE ARE" above.

⛔ **Archiving is NO LONGER an open user question** — it was agreed in principle and **D-01a makes it
step ①**, the next action. Its mechanics moved up into the ①→②→③ block above.

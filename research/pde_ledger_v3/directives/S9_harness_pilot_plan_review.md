# Independent review — the S9 harness pilot plan

## Artifact

`/var/projects/toy_physics/research/pde_ledger_v3/HARNESS_S9_PILOT_PLAN.md`

**Orchestrator-written.** It proposes rebuilding **how a derivation step is verified** — not the physics —
and proving the rebuild on step S9, whose values are already closed and must not move.

Context you need and should verify for yourself: this project derives toy-physics results symbolically.
Each step is built twice, independently — one Wolfram engine, one SymPy engine — from one shared
specification, and a harness (`research/pde_ledger_v3/reduction/`) compares their emitted objects. The
measured problem is that days per step go into making the two engines and the shared spec agree on **how
to report** results, while the physics was right in both engines from the first build.

## What to check

⭐⭐ **The governing question: if this plan were executed exactly as written, would the resulting harness
CATCH A WAY THE PHYSICS COULD BE WRONG that the current one catches — and would its acceptance test
actually be capable of failing?**

⛔ This is a **method** plan. A finding matters only if it changes what gets computed, what may be
claimed, or what method is adopted. ⛔ Not style, not structure, not "this could be clearer."

### ① ⭐⭐ ATTACK L3 — the base-dimension rescaling check. It is the least-vetted idea in the plan.

§3 L3 proposes: substitute `t → α·t`, `x → β·x`, and each dimensioned parameter by its declared scaling,
then require each emitted object to scale as `α^a β^b γ^c` **with the exponents it declared**. The plan
claims this **derives** the dimension from the computed expression rather than accepting a declaration,
and therefore that **a fabricated (typed-literal) value fails it** where a declared `DIM|` tag would pass.

⭐ **Attack that claim specifically:**
- Is it **circular**? It compares against *declared* exponents. Does that make it an agreement between two
  declarations rather than a derivation?
- ⭐⭐ **Can a fabricated value pass it?** ⛔ Do not answer in prose — **exhibit a counterexample.** Write a
  concrete expression (or a mutation of a real S9 emission) that is not computed from the model yet
  satisfies the rescaling test, and show the substitution.
- What does it catch that the existing `dimensional_homogeneity_gate.py` does not, and vice versa?
- Does it break on objects that are legitimately **not** homogeneous, or on dimensionless quantities?

### ② Is the acceptance test able to fail?

§4 lists A1–A10. For each, ask: **what would a wrong outcome look like, and could the criterion be
satisfied by something that is not the thing it claims to test?** In particular:

- **A1** (comparator reproduces S9's 12 AGREE / 0 DISAGREE with all 7 naming declarations deleted) is the
  stated crux. ⭐ Could A1 pass **trivially** — e.g. by a comparator that agrees too easily, by rows
  silently becoming un-compared, or by an equality test that is weaker than the one it replaces?
- **A6** lists seeded checks that must fire and must **not** fire. Is that list adequate? What defect
  would pass every one?
- ⚠ Does any criterion reference an **expected physics value**? The project rule (`CLAUDE.md` rule 5) is
  that a spec says what to compute and **never what anything equals**; an instrument claim must be stated
  on an **invariant** (a count, a stabiliser, a coverage), ⛔ never a sampled outcome. ⭐ Check every one of
  A1–A10 against that rule and say which side of the line it falls on.

### ③ ⭐ What defect class does the stack NOT catch?

§3 defines layers L0–L7. §3b promises a mutation coverage study but ⚠ **has not run it**.

⭐ Do the study's first pass yourself, on paper but concretely: for each defect class below, name **which
layer fires** — or state that none does. ⛔ Do not accept the plan's own claim about which layer covers
what; check it against the layer definitions.

wrong sign · wrong dimensionless coefficient · wrong **form** of a stiffness term · dropped dimensioned
factor · permuted root order · relabelled symbol · **a typed literal substituted for a computed value** ·
an inert premise that drives nothing · wrong root multiplicity · a wrong operator that happens to share a
kernel · a defect present identically in **both** engines because it came from the shared spec.

⭐ **The last one is the important one.** The plan asserts cross-engine comparison is structurally blind to
it. Is any other layer able to catch it, or is the stack blind to shared-spec defects end to end?

### ④ Does the sequencing prevent tuning-to-green?

§5 has Codex build the new comparator against S9's **committed** outputs with the engines untouched.
⭐ Does that actually prevent the comparator being tuned until it agrees? What in the plan stops a builder
from adding special cases until A1 passes? ⚠ The project has previously measured builders converging on
any target they can see.

### ⑤ Leak assessment — ⭐ settle this, do not hedge

⚠ The plan states a physics value in §3 L4: *"S9's gradient-elastic control has a known answer (two
transverse modes at `c² = μ_R/ρ_br`, plus a longitudinal)."*

⭐ **Is that a leak?** The plan is a method document, not an engine build directive — but it will be read
by whoever writes the spec at step 1, and that spec IS read by both engines. Judge concretely: does that
sentence make it possible for a downstream builder to converge on a result it should have computed? If
yes, say exactly what must be removed or reworded. If no, say why the distinction holds. ⭐ Also scan the
**whole plan** for any other statement of what something equals, is expected to be, or was measured to be.

### ⑥ The three-layer diagnosis and the reduction gate

§2 claims the cost sits in three negotiation layers (tag join · symbol naming · shape selection) and that
§3b's L0 "object kind" removes the third. §3b then argues that granularity is **not** grounds for
reduction — retracting an earlier false claim — and gates any reduction on a mutation coverage study.

⭐ Check both. Is the three-layer decomposition complete, or is there a fourth cost the plan misses?
⭐ Is the §3b gate real, or can it be satisfied cheaply enough to become a formality?

### ⑦ Is the pilot's premise sound at all?

⭐ Step back: **is S9 the right pilot, and is a pilot the right move?** The plan argues S9 because it is
small and has an adjudicated answer key. ⚠ Counter-argue it. Is S9 so unrepresentative — 9 actions, 12
compared pairs, no `Piecewise`-heavy objects — that a green pilot would predict nothing about S11's 3,750
tags or S11b's branch-sensitive, transcendental emissions? ⭐ If so, say what the pilot would have to
include to be predictive.

## Do not read

⭐⭐ **Reading order is the blindness mechanism.** Form your own view of what this harness should look like
from the **repository** before you read the plan's argument for its own design.

- ⛔ `/tmp/claude-1000/**` and any file named `*prior-art-review*` — **prior review outputs on the
  document this plan cites.** Reading them hands you another leg's conclusions.
- ⛔ `/var/projects/toy_physics/docs/method_prior_art_findings.md` and
  `/var/projects/toy_physics/docs/method_prior_art_brief.md` — ⛔ **not until step ④ of the method below.**
  These are the plan's own justification and its commission. ⭐ You may read them **last**, to judge
  whether the plan's use of them is honest.
- ⛔ `/var/projects/toy_physics/research/pde_ledger_v3/_scratch/` — raw build transcripts.
- ⛔ Any other `*_review*.md` under `directives/`.
- ⛔ This prompt is not evidence. Claims restated above are restated so you know what to check, ⛔ not so
  you can quote them back as established.

## Required method — DOCUMENT branch

Read the **source of truth first**, in this order:

1. ⭐ **The repository**, and form your own view of what a good harness standard for this project would be:
   - `research/pde_ledger_v3/reduction/checks_S9.yaml` (194 lines) and `checks_S10.yaml` (3,121)
   - `research/pde_ledger_v3/reduction/engine_output_checks.py` — especially `check_cross_engine`,
     `symbolic_equal`, and the naming/selection machinery
   - `reduction/harness_ablation.py`, `reduction/able_to_fail.py`,
     `reduction/dimensional_homogeneity_gate.py`, `reduction/derived_or_declared.py`,
     `reduction/measurements/declaration_load_ablation.py`
   - the engines: `scripts/S9_*.py`, `mathematica/S9_*.wl`, and for object-kind variety
     `scripts/S11bA_*.py`, `scripts/S11bB_*.py`, `mathematica/S11bB_*.wl`
   - `steps/S9_light_requires_shear.md` — ⭐ especially its VERIFICATION section and the table of what each
     mechanism caught and missed
   - `research/pde_ledger_v3/REBUILD_HANDOFF.md`, `CLAUDE.md`
2. ⭐ Write down, **before opening the plan**, what you think the check stack and the acceptance test
   should be. You will be asked to compare.
3. ⭐ **Then** read the plan and audit it against ①–⑦.
4. ⭐ **Then** read `docs/method_prior_art_findings.md` and `docs/method_prior_art_brief.md` and judge
   whether the plan's citations of them are honest. Report this as a separate, clearly-labelled section.

⛔⛔ **A PROSE ASSERTION IS WORTH NOTHING.** Every finding needs a **quotation with a locus** — a file and
line number, or a literal command and its output. Where a claim is checkable by running something, ⭐ run
it and paste the literal output. Where you claim a check can be defeated, ⭐ **exhibit the defeating
example**, concretely.

## Physics filter

> Report a finding only if it could change **what the project computes, what it may claim, or what method
> it adopts**. ⛔ No style, formatting, tone, or structure findings.

⭐ **A finding that says "this check cannot fail" or "this defect class is caught by nothing" or "this
acceptance criterion leaks the answer" is worth more than ten completeness notes.**

⚠ **A leg that returns "nothing survives the filter" is weak evidence on its own.** If that is genuinely
your conclusion, say so plainly — ⛔ but state what you checked, what you could not check, and what would
have had to be true for you to find something. ⛔ Do not manufacture findings to fill a quota.

## Operational constraints

- ⛔ **Do not modify the repository.** Read-only. Copy to `/tmp` if you need to run anything.
- ⛔ **No Mathematica kernel is required.** If you start one anyway, wrap it in `timeout 600`, and ⛔ never
  run more than one at a time — the licence has two seats.
- ⭐ Save any script you write and its literal stdout to named absolute paths, and report those paths.
- ⭐ Rank findings most-severe first. For each: the claim · the evidence (quotation + locus) · what must
  change.

## Quarantine rule

Not applicable — no parallel builder, no sibling implementation. Blindness comes from **reading order**
and the do-not-read list.

# Blocking physics review of the post-fix pilot rows

**Verdict: BLOCKING.** Reviewed `pilot/stage016_rows.md` and `pilot/stage023_rows.md` as re-run against
the eight-gap-fixed schema. The pre-fix findings were removed from the working tree before this leg ran,
so nothing here is an echo of them.

⭐ **What this leg validated** is recorded first, because it is load-bearing for the schema and should
not be lost among the defects.

---

## VALIDATED — Rule A (axis C never evicts a row from tier 1)

**CLEAN.** 30 occurrences are kept in tier 1 with an unresolved closure (`stage023` PY/WL-02, -07…-16;
`stage016` R03-R05, R20). Every one was opened. **Every one is a free symbol the artifact declares and
never assigns, or a form built over such symbols — i.e. an input the model takes.** Keeping them in is
the conservative direction §1 mandates, and none belongs outside the tier system.

⇒ The G1 decision holds under adversarial review. The failure mode on those rows is the *quality of the
axis-A evidence* (D1, D2 below), not the rule that keeps them.

**Universe boundary — spot-checks hold.** `PY:231 P0_physical` reaches only the dimension block;
`q_free`/`gain0`/`gain1`/`eta_null` are control-only, verified at their consumption sites; `PY:168 R0,R1`
are absent from `GENERATOR_DOFS`. The verdict-token adjudication is what §7.1.1 requires, not a
convenience.

---

## D1 — ⭐⭐ AXIS A IS BEING INFERRED FROM AXIS B (the structural one)

`stage023_rows.md:117-124` + E11 (`:395-396`) adopt a blanket rule — *"the reduction is performed at the
binding site itself"* — that assigns `A-REDUCED` to **every** `B-EXECUTED` row. **78 of 121 stage023 rows
are classified this way.** No computed row can then ever be `A-UNADJUDICATED`.

⛔ **This re-creates precisely the fusion §2 exists to break.** The three axes are only informative if
they are independent; a rule that derives one from another collapses the design back into the
single-label taxonomies the schema was written to replace.

**Numerically inert on stage023** (both readings land at `no-tier:unadjudicated`) — **not inert on
stage016**, where it sends **13 occurrences to `no-tier:unclassified-nonfed`, outside the tier-1 range
entirely.**

⭐ Clearest unevidenced case: `stage016_rows.md:208-216` (R11). `PY:279` is `order = list(harmonics)` — a
re-listing of dict keys — called `A-REDUCED` because it *"reduces to the basis declaration at `PY:269`"*,
while `WL:164`'s typed twin of the **identical quantity** is `A-UNADJUDICATED` and sits in the tier-1
upper span. The 5-ness is obtained in both engines from a typed 5-element basis; **the divergence is
produced by the blanket rule, not by the engines.**

⇒ **Fix: the schema must prohibit inferring axis A from axis B, and require evidence for `A-REDUCED`
that is independent of the fact that code ran.** Executing an expression is axis B's claim, not axis A's.

## D2 — Rule B's converse fails on the rows carrying stage016's whole tier-1 lower bound

`stage016_rows.md:113-117` (R03/R04/R05) and `:386-391` (R20) call route R35 (`M̃ = ∫ μ_η β₂² dV`)
**executable within this framework** on two evidences. Both were opened; neither supports it.

1. *"the integrand is already assembled in-artifact at `PY:336`/`WL:232`"* — `PY:336` is
   `mu_eta_density * beta2_dim**2 * measure`, a product of **dimension-walk free symbols** (`PY:223-227`,
   whose only assignment anywhere is an `(L,M,T)` rule at `PY:314-325`). There is no integrand:
   ⛔ **no functional form of `β₂(w)` or `μ_η(w)` exists anywhere in the ledger** —
   `ledger_stage017_…_sympy_audit.py:434-443` carries them as *names in a `calibration_inputs` list*, and
   the scalars as literals (`:41-43`, `Mtilde: 3.0`).
2. *"stage043 records the route as un-run, not un-buildable"* —
   `ledger_stage043_…_sympy_audit.py:381-395` records `moment_integral_executed=False`,
   `source_status="DERIVED-IN-FORM-UNEXECUTED"`. It says **un-executed**; it says nothing about
   buildability. ⭐ **Absence of a denial is not evidence.**

⚠ **The file contradicts itself on the same expression.** `stage016_rows.md:510-516` (R26) — the *same*
`PY:336` — records *"no code path integrates it… no symbolic derivation is present in either engine, only
the integrand product"* and sets `A-UNADJUDICATED` because *"the relation would follow from a wall
action…; the ledger posits none."* One expression cannot be both machinery-present and machinery-absent.

**Consequence: 8 of the 10 `tier1-debt` occurrences (4 of 5 QIDs) rest on this.** Under §9.0 the lower
bound goes `[10,39] → [2,39]` occurrence-level and `[5,20] → [1,20]` QID-level. Under a structural
reading the total holds but 8 occurrences move `tier1-debt → tier1-structural`. ⛔ Neither is evidenced
as written, so neither is asserted here.

⇒ **This is G5's failure mode running the other way.** G5 stopped `structural` being the cheap
assertable claim; nothing yet stops `debt` being one — and `debt` is the *optimistic* direction, because
it says the reduction is available to anyone who does the work.

## D3 — a named route pointing the wrong way

`stage023_rows.md:382-391` evidences `A-IRREDUCIBLE-STRUCTURAL` for `Omega_U, Omega_W, R_mix, g_U, g_W`
with *"Route named at `PREREG:186-195`"*. Opened: `pathA_preregistration.md:191-196` **consumes** those
five symbols to produce `Δ, Q, P, Λ, N0, P0`. **It is a route out of them, not to them.** The prereg never
defines them (`:71`, `:189`, `:192-194` only) and no `REG:` row exists for any.

⚠ The tier is probably rescued by the *other* half of its evidence — the foreclosing property
(`PREREG:180`, "requires the full closed matter/gauge solve"; `REG:300`), which §9 accepts on its own
where no route is named. ⇒ Recorded as an **evidence defect, not a tier flip**. It is ⛔ not
executable-here debt: nothing in the ledger gives these five a value or a formula.

⇒ **Fix: §9 must require a named route to be checked for DIRECTION** — a derivation that consumes X does
not reduce X.

## D4 — route-2 double-counts a single dict literal

`stage023_rows.md:204-205`, `:312-313` mint **two** constitutive propositions from one literal:
`PY:425` is `"pathA29_premise": {"Z_is_premise": True, "boundary_dof": "none"}`, and each is given
`tier1-structural` in each engine — **4 occurrences / 2 QIDs**.

Both describe a fact already counted as `Z0_ret`/`Z1_ret` (PY-10/11, same E8 evidence, same Gate-6
extension). `REG:169` states these are one 2-dof freedom — *"COORDINATE ALIASES of these SAME 2 dofs —
NOT 2 new free dofs"* — and the artifact itself prints the premise as *"PROVENANCE (documentation; not
counted)"* (`PY:814-818`).

If folded: `tier1-structural` **18 → 14**, tier-1 lower bound **26 → 22** (13 → 11 QIDs).
⇒ **Fix: route 2 needs a de-duplication rule against dofs already counted via route 1.** Widening the
universe without one lets the same freedom be counted twice.

## D5 — two smaller ones

- **SUSPECTED.** `stage016_rows.md:779-782` (S8) claims the 18 suppressed occurrences *"would land in
  tier 1 or its upper span"* — unevidenced, and at least `μ_η, T_w, K_η, T_Ω, β₂` are register-`CALIB`
  (`REG:176-178`, `:182-183`), which projects to **tier 2**.
- **DEMONSTRATED.** `stage023_rows.md:214` heads the WL block "56 occurrences"; `:318` and the grand total
  use **58** (16+37+2+3). The 58 is the one consistent with the row tables.

---

## What this means for the fan-out

⛔ Still not ready, and the reason has moved up a level. The **schema's tier-1 split is sound in
direction** (Rule A validated) but its **evidence requirements are not yet strong enough to stop a
builder asserting the optimistic branch** — `debt` over `structural` (D2), and `A-REDUCED` over
`A-UNADJUDICATED` (D1). Both defects push the same way: they make the ledger look more reducible than
the evidence supports.

⭐ That is the direction §1's asymmetry exists to forbid, and it survived two fix passes because §1 was
written about **defaults** and these are **assertions**.

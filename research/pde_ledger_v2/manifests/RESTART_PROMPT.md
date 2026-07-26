# RESTART PROMPT — paste this after a compact / at the start of the next session

> Read `research/pde_ledger_v2/manifests/HANDOFF_NEXT_SESSION.md` first, then
> `manifests/LEDGER_WIDE_PLAN.md` (plan of record) and `docs/development_pipeline.md` (process).
> Do NOT read vision docs.
>
> We are on branch `ledger-v2-rebuild`, in **Pass 1 §2.1** of the ledger-wide dimension-unification
> plan. The dimension survey has **not run yet** — the previous session built and hardened the
> apparatus for it. `research/pde_ledger_v2/schemas/` is untracked.
>
> **Do this, in order:**
>
> 1. **Round 7 is DONE and independently verified — `ACCEPT_WITH_FINDINGS`, pilot READY.**
>    Suite green (122/0). Lexical coincidence is dead: renaming every axis to alien tokens left 17/17
>    honest records accepting (under r16 that flipped to REJECT). Read
>    `_scratch/pass1_dim_survey/reviews/REPAIR7_VERIFY.md` for the two items that must be settled
>    **before the 44-way dispatch** (they do NOT block the 4-script pilot):
>    ⭐ **`DIMENSIONLESS = Dim()` has NO honest escape** (13 lines incl. `stage004:143`, a line §3.5
>    itself names) — fabrication-forcing, fix it in the DIRECTIVE not the schema; and
>    **`LENGTH = Dim(1, 0, 0)` rejects across 72 lines in 13 of 43 scripts** with an undocumented
>    narrow-span escape. Plus a minor §3.0-vs-§3.5 conflict on `basis_locus`.
>    ⛔ **Do not patch the axis grammar.** The root cause was §3.4 and §3.5 being mutually
>    unsatisfiable by construction for stage042 (registry says axis `time`, script spells it
>    `frequency`). r17 satisfies the test from the declared basis + arity instead, and that survives
>    the rename probe that killed r16.
>
> 2. **Re-run the arbiter command yourself before trusting any build report** —
>    `python3 research/pde_ledger_v2/schemas/validate_dimension_survey.py --test-examples`.
>    A prior round printed its done-marker with 45 of 110 fixtures failing.
> 3. **Then run the 4-script pilot** per `_scratch/pass1_dim_survey/directive/PILOT_DIRECTIVE.md`
>    (stage021, stage038, stage042, stage024), including the **planted-defect positive control**.
>    Bring me the three numbers in handoff §5 and a GO/NO-GO before dispatching the other 40.
>
> **Operating model: you are a THIN CONDUCTOR.** Agents return ≤12-line checkable receipts; full
> reports go to disk and you read them only on a negative verdict. Fold reviews via
> `_scratch/pass1_dim_survey/FOLD_PROTOCOL.md` — never hand-fold; every orchestrator hand-fold last
> session introduced a defect. **Pass `model: opus` explicitly on every agent dispatch.** Gate every
> Codex launch on a check that PRECEDES it (Codex snapshots its prompt into argv). Grok is the default
> tertiary reviewer, not GLM.
>
> Standing principles that decided most calls last session: **findings are the product, green is not
> the goal**; **a rule whose only honest outcome is an invented value is worse than no rule**; and
> **a check that cannot fail is worse than no check**.

---

## Quick state table (details in the handoff)

| Item | State |
|---|---|
| SymPy refactor baseline | ✅ 43/43 exit 0, reproducible across two runs |
| Mathematica `.out` determinism | ✅ 44/44 reproduce (after `$NNNNN` normalisation), 442 script-seconds |
| Survey directive | ✅ **r17** — 5 Codex rounds, Grok gate, 9 verification folds |
| Schemas + validator | ✅ green (122/0), verified — 43/43 load-bearing, 105/105 minimal violations |
| 4-script pilot | ⏭ next, then the user gate |
| Remaining 40 + module design + refactor | after the pilot |

## The three numbers the pilot must produce
1. `UNDETERMINED` leaf count per record.
2. ⭐ `CONSTRAINED_BUT_NOT_STATED` count among `REAL+ABSENT` quantities — **the Pass-1 stall signal**.
3. Wall-clock + token cost per record (survey and verification legs separately).

## Do not re-open (decided, with reasons in handoff §2)
Schema-first split · two schema artifacts · orthogonal `ownership`/`uses` · §3.11 narrowed to
dimension-valued bindings · interpreted axis names in `orders` · conditional anchor applicability ·
full case-(ii) ownership coverage · **SymPy-only module, Mathematica stays independent** · parameter
register as evidence source · **semantic dimension adjudication descoped to §4.9** (task #13).

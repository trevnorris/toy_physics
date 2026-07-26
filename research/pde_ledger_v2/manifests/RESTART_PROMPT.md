# RESTART PROMPT — paste this after a compact / at the start of the next session

> ⭐ **Read `research/pde_ledger_v2/manifests/PIVOT_TO_REWRITE.md` FIRST.** It supersedes
> `LEDGER_WIDE_PLAN.md` §2.1–§2.2. Then `HANDOFF_NEXT_SESSION.md`, then repo-root
> `docs/development_pipeline.md` (process). Do NOT read vision docs.
>
> Branch `ledger-v2-rebuild`. The approach **changed on 2026-07-26** after two measurements:
> we are **no longer surveying 43 scripts to discover idioms the module must accommodate**.
> We **rewrite** each script's dimension handling onto a designed module, gated by a
> cross-engine comparison harness plus the standing multi-AI red-team pass.
>
> **Why (do not re-litigate — evidence in `PIVOT_TO_REWRITE.md` §2):**
> - `notes/parameter_register.md` covers only **105 / 226 = 46.5 %** of the scripts'
>   dimension-bearing quantities, its identity key is **not machine-readable**, stage042's
>   stiffness-basis quantities are **absent**, stage038's four-axis quantities are normalised into
>   `{L,T,M}` **with contradictory values**, and some rows claim dual-engine verification for stages
>   that **compute nothing dimensional**. It is an ORACLE for the 105 it covers, not a seed.
> - There is **no single idiom to consolidate** — 6 basis conventions over 4 axis sets, including a
>   2-axis `(L,T)` (008) and a 4-axis (038); 8 scripts use fractional exponents.
>
> **NEXT: build the cross-engine dimension comparison harness (task #20).**
> ⭐ It does not exist — the two runners write to separate trees and **never diff**, and only 9 of 30
> `.out` files print a computed exponent triple. Every dimension-bearing stage must print its computed
> exponent triple; the harness diffs `.py` against `.wl` **on values**, not pass/fail.
> **Why it is load-bearing:** the `.wl` only reveals a `.py` error when the declared dimension and its
> assertion target *disagree*. A shared module supplies **both from one place**, so a module bug moves
> them together and both stay green. It also catches the transposition class — 22 of 30 stages need
> literals transposed between orderings, and a transposition error is symmetric, hence otherwise
> invisible.
>
> Then: fix the register's false provenance (#2 in §4) → design the module against **stage042 and
> stage038 first** (the extremes that broke the register) → rewrite stage-by-stage, **adjudicating**
> the 9 cross-stage conflicts (task #21) rather than transcribing them → red-team pass.
>
> **Operating model: THIN CONDUCTOR.** Agents return ≤12-line checkable receipts; full reports go to
> disk and you read them only on a negative verdict. **Pass `model: opus` explicitly on every agent
> dispatch.** Gate every Codex launch on a check that PRECEDES it (Codex snapshots its prompt into
> argv). ⚠ **Never wait on `pgrep` for a pattern your own waiter contains** — that cost 6h43m on
> 2026-07-26; wait on the captured PID from a pidfile, or on artifact + log quiescence.
> Grok is the default tertiary reviewer, not GLM.
>
> Standing principles that decided most calls: **findings are the product, green is not the goal**;
> **a rule whose only honest outcome is an invented value is worse than no rule**; **a check that
> cannot fail is worse than no check**; and **if each fix bans a spelling and the next probe evades
> it, the architecture is wrong.**

---

## Quick state table

| Item | State |
|---|---|
| Approach | ⭐ **PIVOTED** to rewrite-and-gate (`PIVOT_TO_REWRITE.md`) |
| Scripts surveyed | **0 / 43** — survey premise is dead, apparatus parked |
| Cross-engine harness | ⏭ **NEXT (task #20)** — does not exist |
| Shared module | not started; design against 042 + 038 first |
| SymPy baseline | ✅ 43/43 exit 0, reproducible |
| Mathematica `.out` determinism | ✅ 44/44 reproduce after `$NNNNN` normalisation |
| Survey apparatus | ⏸ parked: schemas + validator committed (`35dc6aa0`), 173 fixtures green, directive at r18 |

## Open findings carried forward (catalogue, do not silently fix)
`PIVOT_TO_REWRITE.md` §3 — **9 cross-stage dimension conflicts** (`K_eta` has three values; `r_BA`
computed in incompatible unit systems; `A_E` 038-vs-037; `ε0/ε1` vs `Z0_ret/Z1_ret`; `M0`; `μ_η`/`T_w`),
plus **tautological checks in committed work**: stage031 has 20 of 21 `PASS_UNITS_*` rows comparing an
expression to a copy of itself, and stage017 asserts a hardcoded `True` in both engines (task #22).
Register locus mis-attribution (stage017/stage015) — task in §4.2.

## Weakened gate — do not reuse "identical PASS counts" naively
16 scripts emit no `PASS tally:` (001–003, 016–028); the **ordered** marker list is not stable
(stage007 swaps two adjacent lines between baseline runs); duplicates lose 6 under set comparison
(stage013 ×5, stage042 ×1). Use a **multiset** comparison, plus the harness, plus `mathematica/out/`
byte-identity after `$NNNNN` normalisation.

## Do not re-do
`.out` determinism probe · SymPy baseline · the two measurements
(`notes/measure_register_sufficiency.md`, `notes/measure_rewrite_feasibility.md`) · repair round 8.

## Parked, only if the survey is revived
Tasks #12 (4-script pilot), #16 (pilot directive), #17 (F3 ownership on `REAL+ABSENT`), #18 (checker
basename surrogate), #19 (multi-locus `any()` acceptance), #13 (computed register adjudication).

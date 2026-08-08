# W3/W4 round 4 — three defects, then this lands

⚠ **Warning: `research/pde_ledger_v3/steps/S9_PILOT_ADJUDICATION.md` holds S9 result values. Don't read it.**

⭐ Two independent legs reviewed round 3. ⭐⭐ **The pin is sound and its independence is established by
mutation, not by reading**: 11 pairs, 0 circular. ⛔ Do not disturb it. ⭐ Three defects block the landing.

---

## ⛔⛔ R1 · A production script prints an exculpatory verdict that is not computed

`reduction/engine_dimension_pin.py:299-305` prints, **unconditionally**:

> *"the standing nonzero transcript residual is a specialisation artefact, not a physics disagreement"*

⭐ **Measured by a leg:** drive the same printer with a shifted registry and the residual becomes
`(2 - D, 0, 0)` — **−1 at the declared `D`**, a genuine disagreement — and ⛔ **the exculpatory sentence
still printed three times.**

⇒ ⭐ The printed characterisation must be **computed from the residual at the declared `D`**, ⛔ not typed.
⭐ Interpretation belongs to the step record. ⛔ A production script prints objects; it does not adjudicate.
⚠ **This one is the orchestrator's fault** — the round-3 directive said *"say so wherever it is recorded"*
and a script is not a record. ⭐ Print the residual and the value it takes at the declared `D`; ⛔ stop there.

## ⛔⛔ R2 · The pin's completeness guards are a demo, not a pin

⭐ **Measured:** of five weakenings built by a leg, **four survive the entire test population** — dropping
the covered-quantity check, dropping the population-count check, dropping the errors check, and dropping the
unmapped-symbol raise. ⭐ With both completeness guards removed:

```
EMPTY_POPULATION_GUARD pin-passed-while-comparing-nothing: True
production exit 0 · 34 passed · test-exit 0
```

⇒ ⛔ As shipped the pin **is** guarded — ⛔ but nothing would notice if it stopped being.
⭐ **Build tests that fail on each:** an empty population, a subset population that still covers the required
quantities, a non-empty error set, and an unmapped symbol.
⛔ **A test that "covers" an invariant demonstrates it on one example.** ⭐ Build each weaker implementation
yourself and show the test fails on it, with literal stdout.

## ⛔ R3 · `README.md:219-224` denies a gap this same round closed

> *"Nothing in this reduction directory witnesses the registry's `D_brane` coefficient against an
> independent symbolic engine emission… is not closed here."*

⛔ False. ⭐ The same file documents the pin at `:113-125` and lists it in the runbook at `:292`, and the
check now reports the coefficient as policed. ⚠ Both legs found this independently; ⭐ it is an **unstaged
addition from round 3**, added beside the runbook entry for the tool that refutes it.

⇒ ⭐ Rewrite it to match what is measured. ⛔ Do not replace one absolute claim with another — ⭐ the last
two attempts at this sentence were wrong in **opposite** directions.

⭐ **State the uncaught class plainly while you are there**, because it is real and both legs confirmed it:
a common-mode shift applied to the registry **and** every transcript together — one wrong shared premise —
leaves the pin, the gate and the law check all green. ⭐ `README.md:121-122` and the handoff already say
this correctly; ⛔ do not weaken it.

---

## ⚠ Two method notes — ⭐ record, ⛔ do not over-build

- `w3_duplicate_pin_ablation.py:19-22` enforces "no duplicate" by grepping **two historical spellings**.
  ⛔ That is a denylist, and a differently-named duplicate walks past it. ⭐ Its real evidence is the
  behavioural flip at `:114-115`. ⭐ Say so, ⛔ or drop the grep — ⛔ do not add more spellings.
- ⭐ Swapping the `mu_R` and `B_comp` tags still passes, because their laws are identical. ⇒ ⭐ attribution
  **between those two** carries no independent content. ⭐ Record it; ⛔ it is not a defect to fix here.
- ⚠ `Q.brane.B_comp` is corroborated by **S11 alone**. ⭐ Single-engine, not circular. ⭐ Say so wherever
  the pin's coverage is described.

## Scope

⛔ Do not modify any engine, step record, committed `.out`, or `checks_S*.yaml`. ⛔ Do not touch the pin's
comparison logic or its operand selection — ⚠ its independence is established and ⛔ re-deriving it risks
what four rounds bought.
⛔ Do not commit.

## ⭐ The regression bar

⭐ Unchanged from round 3: gate, loader, `acceptance_check.py`, `able_to_fail.py` and its cases,
`engine_output_checks.py` on both committed configs, the test population against a pristine `HEAD`, and
⭐ **every engine the import fence discovers**, re-run and diffed.
⚠ `WL_S10_RUNTIME_SECONDS` moves every run; ⚠ `wolframscript` writes config warnings to **stderr** —
⭐ separate the streams. ⛔ One Mathematica kernel at a time; ⛔ `timeout 900` on each.
⛔ **If a `D = 3` result moves, STOP and report it.**

## Rules

- ⭐ Print computed objects — both operands and the residual — then guard. ⛔ No status typed rather than
  computed. ⛔ **No expected value typed beside the value it checks.**
- ⛔ Never `git checkout`, `git stash`, or `git restore` — ⚠ the tree holds uncommitted work.
- ⭐ Absolute paths only; ⛔ never rely on the shell's working directory.

## Deliverables

1. The change.
2. ⭐ The R1 printer driven with a residual that does **not** vanish at the declared `D`, literal stdout.
3. ⭐ The four R2 weaker implementations, each shown failing its test, scripts at named absolute paths.
4. ⭐ The regression evidence.
5. ⛔ A ≤40-line report. ⛔ No conclusions about whether the physics is right.

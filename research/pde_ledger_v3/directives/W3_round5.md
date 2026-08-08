# W3/W4 round 5 — the pin's independence must be COMPUTED, not asserted

⚠ **Warning: `research/pde_ledger_v3/steps/S9_PILOT_ADJUDICATION.md` holds S9 result values. Don't read it.**

⭐ Round 4's three fixes **all hold** — two legs and the orchestrator confirmed them independently, including
with a harness authored **before** round 4 existed. ⛔ Do not disturb them.
⭐ Three defects remain. ⭐ **The first blocks the commit.**

---

## ⛔⛔ B1 · The pin's central claim is held by nothing

⭐ The pin's entire value is one sentence — **"11 pairs, 0 circular"** — and `README.md:118-119` states it as
fact. ⛔ Nothing computes it. It is true today because two legs measured it once.

⭐⭐ **Measured, in an isolated mirror.** Repoint the `S9-wl` record at the **SymPy** transcript, leaving its
label unchanged:

```
ENGINE_DIMENSION_COMPARISON S9-wl Q.brane.rho_br tag=PY_S9_MAIN_DIM_PRIMARY_INERTIA: … residual=[0,0,0]
ENGINE_DIMENSION_GUARD S9-wl Q.brane.rho_br: PASS
ENGINE_DIMENSION_PIN: PASS      exit=0
138 passed
```

⇒ ⛔⛔ **The pin compares one engine with itself and reports PASS.** It prints the contradiction — a `PY_`
tag under a `wl` label — in its own output, and nothing reads it.
⚠ A second weakening does the same: **delete the two `wl` records** ⇒ the pin becomes Python-only with 7
comparisons, still `PASS`, still `D_COEFFICIENT_POLICED_IN_REDUCTION: YES`, still `138 passed`.

⚠⚠ **This is the round-4 defect one level up.** Round 4 pinned the four **guard clauses**; ⛔ nothing pins
the **population table those guards operate on.** ⇒ we tested the check and left unchecked the thing that
decides what it checks.

### The object

⭐ **The pin's independence property, computed from the population table rather than asserted about it** —
so that a change making the population circular, or collapsing it to one engine, **fails**.

⛔ **Name the object; ⛔ do not let this become a spelling check.** ⚠ A denylist over tag names is what
round 4 was told to remove; ⛔ do not reintroduce one here.

⭐ **Both weakenings above are the acceptance surface.** ⛔ If either still passes, the change did not land.
⭐ Build them yourself and show the failure, with literal stdout. ⭐ Then build a **third** weakening the
first two do not cover, and show that too.

## ⛔ B2 · `README.md:287` is false about a command this change set added

> *"From this directory, all commands are bounded externally with `timeout 600`."*

⭐ `registry_import_fence.py --verify` was added under that heading this round. ⭐ Measured: it exits **124**
under `timeout 900`, and takes roughly **33 minutes** on an idle machine. ⭐ Its own budget contradicts the
sentence directly — `registry_import_fence.py:101` uses `3600` for `.py` and `910` for `.wl`.

⇒ ⭐ State the real bound for that entry. ⛔ Do not "fix" it by lowering the fence's budget — ⚠ the bar is
**observable progress**, ⛔ not elapsed time, and this command prints as it goes.

## ⛔ B3 · The removed denylist reappeared in the test file

⭐ Round 4 correctly removed the two-spelling source grep from `w3_duplicate_pin_ablation.py`. ⛔ It then
added the same thing over stdout:

```
test_dimension_laws.py:472:    assert "specialisation artefact" not in completed.stdout
test_dimension_laws.py:473:    assert "physics disagreement" not in completed.stdout
```

⇒ ⛔ Any re-added verdict with **different wording** walks past both. ⚠ **This is the third location this
denylist pattern has appeared in this file set.**

⭐ **The object: every non-empty line the printer emits is one the printer computes.** A leg measured that
this already holds for the shipped printer under three adverse drivers
(`lines-outside-the-computed-token-set: 0`), ⇒ ⭐ it passes today and fails on any added prose.

---

## ⭐ The regression bar

⭐ Unchanged, and **fully closed as of this round** — ⛔ do not let it reopen: gate, loader,
`acceptance_check.py`, `able_to_fail.py` and its cases, `engine_output_checks.py` on both committed configs,
the test population against a pristine `HEAD`, and **every engine** re-run byte-identical.
⭐ The orchestrator has confirmed all of these, including `S11-py` at `sha 5ed934e5…` with empty stderr.
⚠ `WL_S10_RUNTIME_SECONDS` moves every run. ⚠ `wolframscript` writes config warnings to **stderr**.
⛔ One Mathematica kernel at a time (two seats).

⛔ **If a `D = 3` result moves, STOP and report it.**

## Scope

⛔ Do not modify any engine, step record, committed `.out`, or `checks_S*.yaml`. ⛔ Do not commit.
⛔ **Do not change the pin's comparison logic or its operand selection** — ⚠ its independence is established
and ⛔ re-deriving it risks what five rounds bought. ⭐ You are adding what **maintains** that property.

## Rules

- ⭐ Print computed objects — both operands and the residual — then guard. ⛔ No status typed rather than
  computed. ⛔ No expected value typed beside the value it checks. ⛔ **No denylist.**
- ⛔ Never `git checkout`, `git stash`, or `git restore`.
- ⭐ **Absolute paths only.** ⛔ Never rely on the shell's working directory — ⚠ that has failed **four**
  times in this session, including for the orchestrator, and once it edited live files.

## Deliverables

1. The change.
2. ⭐ Both **B1** weakenings run against your build, each shown **failing**, literal stdout — plus a third
   of your own.
3. ⭐ The **B3** whitelist shown passing on the shipped printer and failing on added prose.
4. ⭐ The regression evidence.
5. ⛔ A ≤40-line report. ⛔ No conclusions about whether the physics is right.

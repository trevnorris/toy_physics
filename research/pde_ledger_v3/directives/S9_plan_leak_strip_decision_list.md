# Decision list 2 — strip the plan's leaks and create the adjudication artifact

**You are applying two decisions, ⛔ not authoring anything.** ⛔ Do not re-litigate, ⛔ do not improve
beyond this list, ⛔ do not soften a retraction.

**Files you may touch — exactly two:**
1. `/var/projects/toy_physics/research/pde_ledger_v3/HARNESS_S9_PILOT_PLAN.md` (edit)
2. `/tmp/claude-1000/-var-projects-toy-physics/36f37d88-e717-46ce-a790-6f9d1ef3d7bc/scratchpad/S9_PILOT_ADJUDICATION.md` (create)

⛔ **No other file.** ⛔ Do not `git add`, ⛔ do not commit.

**Report when done:** ≤40 lines — `D18: …`, `D19: …`, a residual-leak grep result, and the two paths.

---

## ⛔⛔ Why file 2 is outside the repository — do not "fix" this

**Blindness is enforced by ABSENCE, ⛔ never by instruction.** A do-not-read list is a denylist, and a
denylist means the architecture is wrong. The comparator builder's sandbox is the repository ⇒ a file
outside it is **unreachable**, which is the only control that actually holds.

⚠ This follows the project's existing, validated pattern: **the pre-registration lives in the session
scratchpad during the build and is folded into the step record after the reviews land.** ⛔ Do not place
file 2 in the repo. ⛔ Do not reference its *contents* from the plan.

---

## D18 · Strip every value the plan's own §3c lists as leaked

§3c enumerates predeclared outcomes a builder must not see. **The plan still contains them.** Remove each
from the plan and move it to file 2. ⛔ A value must not survive as a paraphrase in a neighbouring
sentence — that is the leak pattern this repository has already measured twice.

| # | locus | what is wrong | what to do |
|---|---|---|---|
| a | `:12`, `:492` | *"S9's values do not move"* / *"No S9 value moves"* tells the builder what a green run looks like | ⭐ Restate as an **instruction**, ⛔ not a prediction: *"If any computed value differs from the committed output, STOP and report it — ⛔ do not adjust either engine or the comparator."* |
| b | `:24`, `:32` | the **verdict** `12/12 agree, 0 disagree` (and the S10 verdict cells `545/690`, `145 no verdict`, `11 false DISAGREE`) | ⭐ Keep **inventory counts only** — config lines, declared pair counts, tags emitted, script lines. ⛔ Delete every verdict cell from the §1 table and from the surrounding prose; move them to file 2 |
| c | `:76` | `omega2 = omega**2` names S9's identity | ⭐ Make it **generic**: *"a squared-variable identity of the form `x2 = x**2`"*. The D8 point (an algebraic identity is not a spelling exception) must survive intact |
| d | `:154`, `:155`, `:161` | ⛔⛔ **the L3 counterexample block reveals the true form of S9's result** — it was introduced *while repairing a different leak* | ⭐ **Re-run the counterexamples with placeholder symbols** carrying declared dimensions (e.g. `A`, `B` with `A → γβ⁻¹α⁻²A`, `B → γβ⁻³B`), and paste the **genuine literal stdout of that placeholder run.** ⛔ Do not hand-edit the existing block into placeholders — the output must be real. ⭐ The three classes that must still be demonstrated: a wrong dimensionless **coefficient** passes · a wrong **sign** passes · a fabricated **dimensionless** quantity passes. Also keep the two stated non-run cases: `0` satisfies every declared exponent, and per-entry dimensions defeat one scalar exponent |
| e | §3c's own table | it must **describe** each leaked item by category and locus | ⛔ Verify no row restates a value. Fix any that does |

⭐ **After editing, grep the plan for residual leaks and paste the literal result in your report:**

```
grep -nE 'mu_R/rho_br|μ_R/ρ_br|12/12|12 AGREE|E2 = 2|omega2 = omega|value moves|do not move|545/690|145 no verdict|11 false' \
  research/pde_ledger_v3/HARNESS_S9_PILOT_PLAN.md
```

⚠ Some hits may be legitimate (e.g. §3c naming a category). ⭐ Report each surviving hit with a one-line
justification, ⛔ do not silently leave any.

---

## D19 · Create the adjudication artifact (file 2)

It receives everything D18 removes, plus the material §3c and §3b require to exist but never created.
⭐ **Read the values out of the repository** (`reduction/checks_S9.yaml`, the committed engine outputs
under `scripts/out/` and `mathematica/out/`, `steps/S9_light_requires_shear.md`,
`reduction/measurements/declaration_load_ablation.py`) — ⛔ do not invent any.

**Required sections:**

1. **Status header** — orchestrator-only; ⛔ never readable by a comparator builder; applied **only after
   the generic comparator is frozen**; folds into the step record after the pilot closes.
2. **Expected per-row verdicts** — S9's declared cross-engine pairs and the verdict each reaches under
   today's comparator. ⭐ This is the A1 differential target.
3. **Load-bearing declarations** — which of the six naming exceptions and the one symbol identity change a
   row when removed, and which row each moves. ⭐ Run `declaration_load_ablation.py` and paste its literal
   stdout; ⛔ do not restate it from any document.
4. **S9 result values the plan must not carry** — everything D18 removed, recorded once, here.
5. **Held-out mutants** — ⛔ **leave this section EMPTY with a stub heading.** ⚠ Held-out mutants must be
   authored **after** the comparator interface is frozen, by a party that is not the builder. ⭐ Record that
   requirement; ⛔ do not author them now — mutants written today are visible to today's plan.
6. **Probe points / seeds** — ⛔ **EMPTY stub, same reason.** ⭐ Record the rule: drawn after implementation,
   re-drawn per adjudication run.
7. **Control structural expectations** — ⭐ what each control must make **move** (a root appearing, a
   multiplicity reducing, a sign changing), ⛔ never what any quantity equals.

⛔ **Do not add an eighth section.**

---

## Standing prohibitions

1. ⛔ **Introduce no new physics value anywhere**, including in file 2 — file 2 records values **read from
   the repository**, never derived, inferred or invented.
2. ⛔ **The plan must not state, cite, quote or summarise file 2's contents** — it may say only that the
   artifact exists, where, what categories it holds, and when it applies.
3. ⛔ Add no new checks, layers or acceptance criteria of your own.
4. ⭐ Preserve existing marker conventions (⭐ ⛔ ⚠ ⇒) and every retraction block verbatim.

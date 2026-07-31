# DIRECTIVE — step ①: retire the `a`-pin from live infrastructure

**You are the builder.** Make the edits, run the acceptance commands, iterate to green, report.
Working root: `/var/projects/toy_physics`. Branch `ledger-v2-rebuild`.

⛔ **Do not commit.** The orchestrator commits.
⛔ **Do not read** `research/pde_ledger_v2/_scratch/PLAN_apin_repair_and_throat_restart.md` or anything
under `research/pde_ledger_v2/_scratch/reviews/`. They contain an author's predicted outcome for the
numbers you are about to recompute, and reading them would destroy the independence of your result.
This directive is self-contained.

---

## 1. Why

`research/pde_ledger_v2/notes/stages/ledger_stage004_gnls_action_dimensional_foundation.md:122-137`
imposes **four** natural-unit pins `{a, c_s0, ħ, m_GNLS}` on **three** base dimensions `{L,T,M}`,
reports "rank 3, nullity 1", and reads the surviving dimensionless monomial as a relation:

```
a·c_s0·m_GNLS/ħ = 1     i.e.    a = ħ/(m_GNLS·c_s0)
```

Three independent reviews confirmed the arithmetic: `{c_s0, ħ, m_GNLS}` alone is already rank 3,
nullity 0 — a complete unit system. The fourth pin is redundant, and *any* fourth pin of dimension `L`
produces exactly that monomial and no other. Its **form** is therefore forced by dimensional analysis.

But its **value** is not. Asserting `a = ħ/(m c_s0)` for a physical `a` is a substantive claim: that the
throat radius equals the medium's natural length. That claim is false in kind — `ħ/(m c_s0)` is built
only from medium primitives, so it is **one number for the whole medium**, whereas a throat radius is a
property of an individual defect, and defects form a species-indexed family.

**User decision:** the pin is retired. `a` denotes the **throat radius**, one quantity, and its value is
undetermined — there is no equation for it anywhere in the corpus.

**Scope decision:** v2 prose is becoming reference material. This directive repairs only what a future
ledger will **compute with** or **read as foundation**. Everything else is quarantined and left alone.

---

## 2. ⛔ Scope — exactly these seven files

Edit **only** these. Touching anything else is out of scope for this directive.

| # | file |
|---|---|
| 1 | `research/pde_ledger_v2/reduction/relations.yaml` |
| 2 | `research/pde_ledger_v2/reduction/quantities.yaml` |
| 3 | `research/pde_ledger_v2/reduction/acceptance_check.py` |
| 4 | `research/pde_ledger_v2/reduction/able_to_fail.py` |
| 5 | `research/pde_ledger_v2/reduction/README.md` |
| 6 | `research/pde_ledger_v2/notes/parameter_register.md` — **line 132 only** |
| 7 | `research/pde_ledger_v2/notes/stages/ledger_stage004_gnls_action_dimensional_foundation.md` |

### ⛔ Explicitly OUT of scope — do not touch, do not "fix in passing"

`parameter_register.md` anywhere except line 132 · `ledger_stage016_l2_so3_covariance.md` ·
`ledger_stage023_nullspace_underdetermination.md` · `ledger_stage005_*` · `stage005_pathA20_source_map.md` ·
`ledger_stage043_*` (`.py`, `.wl`, `.md`) · `ledger_stage028_*` · `midway_knob_audit_codimension_sympy.py` ·
`software/` (all of it) · `research/4d_*` · `docs/` (all of it) · `STATUS.md` · `archive/`.

⚠ Several of those contain statements that are wrong for the same reason. **That is known and
deliberate.** They are recorded in a defect register and are not this directive's job.

---

## 3. The edits

### 3.1 `relations.yaml` — remove the relation

Delete the `R2.a_pin` relation entry entirely. It is not a defining equation; it is the residue of an
over-specified unit system.

### 3.2 `quantities.yaml` — remove the quantity from the medium block

Remove `Q.medium.a_pin`.

**Reasoning, so you can judge whether the edit is faithful:** this registry is the medium block —
phase-0 material, the substance before any defect exists. A throat radius is a **defect-sector**
quantity. It cannot be an input to a medium-level simulation, because at that level there is no throat.
Carrying it here would repeat an error already recorded against this registry (it carries `c_γ`, a
light-sector quantity, that no walkthrough step ever introduced). The throat radius will be introduced,
with proper provenance, by whichever step first has defects in it.

⚠ Verify before deleting: confirm that nothing consumes it. Check every relation's `input_qids`, and
grep the whole `reduction/` directory. If something *does* consume it, **stop and report** rather than
forcing the deletion.

### 3.3 `acceptance_check.py` — recompute, do not preserve

The expected payload is a literal dict inside the script. Removing a relation changes what the registry
computes.

⛔ **Do not adjust the registry, the case definitions, or anything else to make the existing numbers come
back.** The old numbers described a registry that admitted a relation this directive deletes; they are
not a target.
⭐ **Recompute the payload from the edited registry and report every number you get, with your derivation
of why each is what it is.** A number that merely matches something is not evidence; the orchestrator
needs your reasoning, not just the value.

Also review `medium_cases()`: one of the mutation cases is built from the quantity being removed. If a
case can no longer be constructed, say so and propose what it should become — do not silently drop it.

### 3.4 `able_to_fail.py` — one probe's premise is destroyed

`demonstrate_entailed()` (around `:64-77`) constructs `xi_h**2 - 2*a_pin**2` and asserts that adding it
does **not** reduce the constraint dimension, because it is entailed by two relations in the registry.

Once `R2.a_pin` is deleted that entailment no longer holds, so the probe would be asserting something
false about the new registry.

⇒ Either **replace it with a genuine entailment of the relations that survive** — one you derive and can
show is identically zero on the registry's own definitions — or **retire the probe** and record in the
file why it was retired. ⛔ Do not leave it asserting a premise that no longer holds, and ⛔ do not weaken
the harness to make it pass. The harness must still demonstrate that it can catch a real defect.

### 3.5 `reduction/README.md` — the `a_pin` note

The README carries a note (around `:79-88`) describing `Q.medium.a_pin` as unresolved and pending a user
decision. That decision has been made and the entry is gone. Rewrite the note to describe the current
state and say why the quantity is absent. Keep it short.

### 3.6 `parameter_register.md` — line 132 only

That row currently presents `a` as a units-pin choice with a defining equation. It is cited as a source
locus by `quantities.yaml`, which is why it is in scope while the rest of the file is not.

Rewrite the row so it describes `a` as the **throat mouth radius** — a physical quantity with **no
defining equation** and no determined value.

⛔ **Line-count-preserving.** This file is cited by line number from at least ten places. Your edit must
leave every other line at its current line number.

⚠ Since `Q.medium.a_pin` is being removed, its citation of this line disappears with it. Check whether
any *other* live file cites `parameter_register.md:132`, and report what you find.

### 3.7 `stage004` note — the definition site

**This is the source. It is in scope precisely because a future ledger will read it as the dimensional
foundation, and leaving it intact would re-introduce the error.**

- **§4, the pin null-relation** (currently `:122-137`) — rewrite so the unit system is the three pins
  `{c_s0, ħ, m_GNLS}`, which you should verify yourself is rank 3 / nullity 0 over `{L,T,M}`. The
  relation `a = ħ/(m_GNLS c_s0)` must be **gone** as an assertion about a physical length. If you retain
  any statement about the four-pin construction, it must be explicitly historical and must not assert
  the relation.
- **`:197`** — currently lists `a` among the quantities that "follow" from the primitives, alongside
  `c_s0` and `xi_h`. `a` does not follow from them. Remove it from that list.
- Add a short banner near the top recording that the pin was retired, dated, so a reader arriving
  mid-file is not misled.

⛔ **Line-count-preserving for lines outside `122-137`.** Sixteen line-pinned ranges cite this file:
`65`, `76-83`, `85-92`, `102-120`, `122-137`. If your §4 rewrite changes the length of the `122-137`
block, you **must** update every citation of it in `reduction/quantities.yaml` in the same change, and
report that you did.

### 3.8 The engines are NOT in scope

`ledger_stage004_..._sympy_audit.py` and `..._mathematica_audit.wl` correctly print that `a` is a branch
moment and not a base invariant. They are quarantined with the rest of v2. Do not touch them.
⚠ If you find that leaving them makes some check **fail**, stop and report — do not expand scope.

---

## 4. Acceptance — run all four, from `research/pde_ledger_v2/reduction/`

These were run against the pre-edit tree and all four exited 0. Iterate until they do again, ⛔ except
where a change of result is the *correct* outcome of this directive — in which case explain it rather
than suppress it.

```
python3 test_registry.py
python3 dimensional_homogeneity_gate.py
python3 acceptance_check.py
python3 able_to_fail.py
```

⛔ **No script may exceed 10 minutes.** If one does, stop and report; do not raise a timeout.
⛔ If a check fails and the honest fix is outside this directive's scope, **stop and report**. Do not
expand scope and do not weaken the check.

---

## 5. Report — at most 60 lines

1. **Per file**: what you changed, in one or two lines each.
2. **The recomputed acceptance payload**: every number, plus your derivation of why each is what it is —
   ambient count, how many constraints are admitted, and the rank.
3. **The entailed probe**: what you did and why; if you replaced it, show the entailment is identically
   zero on the registry's definitions.
4. **Locus integrity**: whether any line-pinned citation moved, and what you did about it.
5. **Anything you found that contradicts this directive.** ⭐ This is the most valuable part of the
   report. If a premise here is wrong, say so — the directive is not authoritative over the corpus.
6. **Anything you stopped on.**

⛔ Do not report a check as passing without having run it. ⛔ Do not describe work you did not do.

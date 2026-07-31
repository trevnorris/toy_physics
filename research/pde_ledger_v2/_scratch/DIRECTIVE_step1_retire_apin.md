# DIRECTIVE — step ①: retire the `a`-pin from live infrastructure

> ⛔⛔ **HISTORY — NOT INSTRUCTION.** This is the step-① build directive, completed 2026-07-31 (`407eed94`). Its "OPEN
> pending a user decision" language about `R2.a_pin` describes a question that **no longer exists** —
> the relation was retired by removal. ⛔ Do not execute anything here.
> ⇒ **Live workstream: `research/pde_ledger_v3/` (start at `NEXT_SESSION.md`).**

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
| 6 | `research/pde_ledger_v2/notes/parameter_register.md` — **lines 132 and 269 only** |
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

⚠ Also review `medium_cases()`: one of the mutation cases resolves the quantity being removed, so it
will **crash** with a validation error rather than produce a wrong number.

⭐ **That case's job is to be an ENTAILMENT control** — it adds a residual that is already implied by the
admitted set, and the suite passes only if the dimension does **not** drop. Its replacement must do the
same job, built from relations that survive. ⛔ Do not substitute a different kind of mutation: a
vacuous or an independent one would still produce a self-consistent recomputed payload while silently
removing the suite's only entailment coverage.

⛔ Do not silently drop the case, and ⛔ do not keep the quantity in the registry merely so it still
runs. A test does not get to dictate what the registry contains.

### 3.4 `able_to_fail.py` — one probe's premise is destroyed

`demonstrate_entailed()` (around `:64-77`) constructs `xi_h**2 - 2*a_pin**2` and asserts that adding it
does **not** reduce the constraint dimension, because it is entailed by two relations in the registry.

Once `R2.a_pin` is deleted that entailment no longer holds, so the probe would be asserting something
false about the new registry.

⚠ Precisely: `resolve_qid("a")` will have no referent, so the probe **crashes** with a validation error
before any assertion runs. You will see a traceback, not a failed assertion.

⇒ **Replace it with a genuine entailment of the relations that survive.** Derive one yourself and show
it is in the ideal generated by the surviving constraints — i.e. it vanishes on the constraint variety
while being **a non-trivial polynomial**, not an expression that is identically zero as a symbol.
⛔ An expression like `x − x` satisfies "identically zero" while testing nothing; that is not acceptable.

⛔ **Do not retire the probe by removing it from `CASES`.** `able_to_fail.py:182` emits
`expected={len(CASES)}`, and `test_registry.py:62` asserts on that count with a hardcoded value. Removing
a case breaks a test in a file that is out of scope. If you conclude a valid replacement entailment
cannot be constructed, **stop and report** — do not retire, and do not edit `test_registry.py`.

⛔ Do not weaken the harness to make it pass. It must still demonstrate that it can catch a real defect.

⛔ **Three specific ways this could be faked — none is acceptable:**
- catching the resolution error and printing the expected-failure marker without ever calling
  `constraint_dimension`, so the probe reports CAUGHT without exercising anything;
- a replacement that is only *numerically* near-zero rather than an exact identity;
- a replacement that is trivially zero as a symbol.

⇒ **Your report must contain the algebraic proof** that the replacement lies in the ideal generated by
the surviving constraints — the explicit combination, not an assertion that it holds.

### 3.5 `reduction/README.md` — the `a_pin` note

Three separate places go stale, not one:

- the note around `:79-88` describing `Q.medium.a_pin` as unresolved and pending a user decision — that
  decision has been made and the entry is gone;
- `:4`, which calls this the **"ten-scalar medium block"**;
- `:108`, which enumerates `R2` as `R2.xi_h`, `R2.a_pin`, and `R2.h0`.

Fix all three so the README describes what the directory actually contains. Keep it short.
⚠ `:4` also attributes the block to a consumer script — check whether that attribution is still true and
report it; do not edit that script.

### 3.6 `parameter_register.md` — lines 132 and 269 only

**`:132`** currently presents `a` as a units-pin choice with a defining equation. Rewrite the row so it
describes `a` as the **throat mouth radius** — a physical quantity with **no defining equation** and no
determined value.

**`:269`** is the `R2` edge row. It still states `a = ħ/(m c_s0)` and describes the edge as collapsing
`ξ_h, a, h0` into primitives. ⭐ **It is in scope because `relations.yaml` cites this exact line as the
source locus for relations that SURVIVE** (`R2.xi_h` and `R2.h0`). Leaving it would point live relations
at a line asserting the relation this directive deletes — the file contradicting itself, ten rows after
you fixed `:132`. Rewrite it so it describes only the surviving collapse.

⛔ **Line-count-preserving, both edits.** This file is cited by line number from at least ten places.
One table row must become one table row. Every other line must keep its current number.

⚠ `Q.medium.a_pin`'s citation of `:132` disappears with the quantity. Other files cite `:132` in prose
(not machine-validated). **Enumerate every citation of `:132` and `:269` you find and report them** — do
not edit the citing files, they are quarantined.

### 3.7 `stage004` note — the definition site

**This is the source. It is in scope precisely because a future ledger will read it as the dimensional
foundation, and leaving it intact would re-introduce the error.**

⚠ **§4 is not the only place this file asserts the pin.** All of the following must also stop asserting
it, and all sit above line 120, so each must be reworded **in place at identical line count**:

- **`:7-8`** (Status) — lists "the pin null-relation" among the stage's exact-symbolic claims.
- **`:21-22`** (Purpose) — states `a = hbar/(m_GNLS c_s0)` outright as something the stage establishes.
- **`:43-45`** (Provenance) — says the engines "derive — do not hardcode — … the pin relation".
- the closed-claims list near `:185-190`, which carries `A_PIN_IS_BRANCH_MOMENT_NOT_INVARIANT`.

⛔ Leaving these is worse than leaving §4 intact: a reader going top-down gets the retired relation as
the stage's **Purpose** and never reaches your correction.

- **§4, the pin null-relation** (currently `:122-137`) — rewrite so the unit system is the three pins
  `{c_s0, ħ, m_GNLS}`, which you should verify yourself is rank 3 / nullity 0 over `{L,T,M}`. The
  relation `a = ħ/(m_GNLS c_s0)` must be **gone** as an assertion about a physical length. If you retain
  any statement about the four-pin construction, it must be explicitly historical and must not assert
  the relation.
- **`:197`** — currently lists `a` among the quantities that "follow" from the primitives, alongside
  `c_s0` and `xi_h`. `a` does not follow from them. Remove it from that list.
- Record the retirement — dated — **inside the rewritten §4 block**, so a reader landing there is not
  misled. ⛔ **Do not add a banner at the top of the file.** See the locus constraint below: a banner
  above line 120 silently invalidates fourteen citations.

⛔ **Nothing at or before line 120 may shift.** Sixteen line-pinned ranges cite this file — `65`,
`76-83`, `85-92`, `102-120`, `122-137` — and two of those (the `122-137` pair) disappear with
`Q.medium.a_pin`. **All fourteen survivors sit at or before line 120, i.e. above your edit.** So the
`122-137` block may change length freely, provided you insert, delete, or reflow **nothing** above it.

⛔ **A shifted locus fails silently, which is why this matters.** `registry_read.py:440-455` validates
only that the file exists, that `line_end >= line_start`, and that `line_end` is within the file. It
**does not check that the cited lines still contain the cited content.** A citation shifted by a banner
would therefore keep passing every gate while pointing at the wrong text — a silent corruption, not a
detectable failure. ⇒ Do not rely on the gates to catch a line shift; they cannot.

### 3.8 The engines are NOT in scope — and this is a deliberate narrowing

`ledger_stage004_..._sympy_audit.py` and `..._mathematica_audit.wl` do **more** than print a caveat:
they **assert** the relation as a derived result (`..._sympy_audit.py:302` and `:306`,
`expect_zero("derived pin relation a=hbar/(m_GNLS*c_s0)", …)`). So they are not harmless, and they are
not "already correct".

⛔ **Do not touch them anyway.** They are v2 stage audits; they do not feed the reduction registry, and
the future ledger will not run them. ⚠ Note this **narrows** an earlier recorded decision that listed
these engines as part of the repair — that narrowing is deliberate, made on the "fix what computes,
quarantine what only narrates" rule, and is recorded here so it is not a silent change.

⇒ **Because they stay stale, the stage note must say so.** In your `:43-45` rewrite, state plainly that
this stage's engines still assert the retired relation and are known-stale on that point. ⛔ The note
must not present them as a current certificate for a claim the note itself has withdrawn.

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

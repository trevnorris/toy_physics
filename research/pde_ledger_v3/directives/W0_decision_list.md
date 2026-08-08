# W0 — decision list for the restatement

⛔ **The orchestrator's directive `W0_emit_speed_itself.md` was reviewed by two legs and is NOT safe to
issue.** Both legs independently reached that verdict; the orchestrator verified every structural claim
below against the repository. ⇒ rule 15: **the author changes.** This file is the decision list; the
directive gets written by someone else from it.

⚠ Everything here is **measured**, with the command that measured it. ⛔ Nothing here is a conclusion about
what the physics should be.

---

## What is actually true — ⭐ verified, ⛔ replaces the old directive's diagnosis

### D1 · The old diagnosis was right for S9 and wrong for S10

| | bound tag | its committed value | what it is |
|---|---|---|---|
| S9 | `WL_S9_CANDIDATE_SPEED_SQUARED1` | `muR/rhoBr` | a **squared speed** |
| S10 | `WL_S10_MAIN_D3_Q3_DISTINCT_ROOTS` (`sequence_second`) | `((k1^2+k2^2+k3^2)*muR)/rhoBr` | ⛔ an **ω²**, carrying `|k|²` |

⇒ the S10 residual is `|k|²(μ_R/ρ_br) − √(μ_R/ρ_br)`, ⛔ **not** the form the old directive asserted.
⇒ ⭐ S9 and S10 have **different defects** and must not be described by one sentence.

### D2 · ⭐⭐ Both S10 engines ALREADY emit the squared speed — the check binds neither

⚠ **CORRECTED 2026-08-07.** ⛔ The first version of this entry showed one value against both roots. ⭐ The
measurement, in both engines and every `D` cell:

```
WL_S10_MAIN_D{2,3,4,5}_ROOT1_Q6_RATE_COEFFICIENT:               0
WL_S10_MAIN_D{2,3,4,5}_ROOT2_Q6_RATE_COEFFICIENT:               muR/rhoBr
PY_S10_MAIN_D{2,3,4,5}_ROOT1_Q6_ROOT_OVER_WAVENUMBER_NORM:      0
PY_S10_MAIN_D{2,3,4,5}_ROOT2_Q6_ROOT_OVER_WAVENUMBER_NORM:      mu_R/rho_br
```

⇒ ⭐ **W0 is primarily a BINDING defect, ⛔ not an absence of computation.** The old directive framed it as
four engines failing to emit. ⛔ That is false at S10 in both engines and in every `D` cell.

⇒ ⭐⭐ **And `ROOT1` is the ZERO root.** ⇒ selecting among candidates is **load-bearing**, ⛔ not cosmetic:
a build that selects the wrong one emits `0`, not a speed. ⭐ This is the concrete reason decision **3**
(pin the transverse certificate) cannot be left to the builder's judgement.

### D3 · R4's harness coverage is one engine, one package, one dimension

Both `registry_residual` blocks declare `engine: wl`. ⇒ ⛔ **the Python engines never run R4 at all**, at
either step. S10's binding is `MAIN` `D3` only ⇒ `D = 2,4,5` never enter R4.
⇒ ⭐ An emission change alone does **not** make R4 tested. ⛔ Do not let the directive imply it does.

### D4 · ⛔ The "additions only" regression bar is IMPOSSIBLE

```
WL_S10_RUNTIME_SECONDS: 60.62794      # AbsoluteTiming — moves every run
WL_S10_EMITTED_TAG_COUNT: 2983        # every emit increments it; adding tags MUST change this line
```
⇒ ⭐ The bar must be a **deterministic projection of physics tags**, with runtime excluded and the
inventory/count delta stated and checked exactly. ⛔ "Only additions" is also **insufficient** — every
worthless emission below is addition-only.

### D5 · ⛔⛔ A compliant build can emit a FREQUENCY instead of a speed

`sqrt(selected_root)` where the root is `ω²` yields `|ω|`, dimension `[0,-1,0]`. ⚠ It descends from the
solved dispersion relation, selects by a computed criterion, emits a premise, and **passes all three
ablations the old directive demanded.** ⇒ ⭐ the object must be pinned by its **dimension and its defining
residual**, ⛔ not by the word "speed".

### D6 · ⛔ The anti-tautology ablation does not distinguish anything

Perturbing `μ_R` upstream moves the honest implementation and the registry-recombined one **identically**.
⇒ ⭐ The distinguishing mutation must inject a **fresh sentinel after root selection and after the division
by `|k|²`**, holding `μ_R` and `ρ_br` fixed. ⭐ The honest object moves; the recombined one does not.
⚠ This is the "a test that covers an invariant demonstrates it on one example" class again.

### D7 · ⛔ The old directive leaked the target

It printed R4's residual tree, simplified the residual algebraically, named the forbidden implementation's
exact formula, and then **required the builder to implement that formula** as an exercise.
⇒ ⭐ **Split the statement in two.** The **builder** statement carries the object and its certificate; the
**integrator** statement carries R4, the residual, and the verdict. ⛔ A builder never sees the second.

### D8 · Scope outside the four builds

Only two `registry_residual` rows exist repository-wide, both R4, both `wl`. R5 has the same linear-speed
shape and **no** harness row. S11's current engines expose only an `ω²` longitudinal root.
⇒ ⭐ Either fold R5/S11 in, or ⭐ **record the deferral explicitly.** ⛔ Do not leave it unstated.

---

## Decisions the restatement must implement

| # | decision |
|---|---|
| **1** | ⭐ Split into a **builder statement** and an **integrator statement**. ⛔ The builder statement contains no relation, no residual, no expected value, no current tag name, and no selector index. |
| **2** | ⭐ Pin the object by **dimension and defining residual**, not by name — the scalar whose square equals the transverse branch's `ω²/(k·k)`, non-negative, of dimension `[1,-1,0]`. ⭐ Require the link residual to be emitted as a computed object. ⛔ This is what stops D5. |
| **3** | ⭐ Pin the **transverse certificate**: what is computed, its required fields, and that the selection cardinality is **one**. ⛔ "polarization, multiplicity, whatever this engine computes" is not a specification. |
| **4** | ⭐ Enumerate the **required cells** — which packages, which `D` — and explicitly exclude the rest. ⛔ Leaving this free is what makes four builders diverge. |
| **5** | ⭐ Replace the regression bar per **D4**: deterministic physics-tag projection, runtime excluded, count delta stated and verified. |
| **6** | ⭐ Replace ablation 1 with the **sentinel-after-division** mutation of **D6**, and state the required differing observable — ⛔ not "show what distinguishes them". |
| **7** | ⭐ State plainly that this change **does not make R4 tested**; the binding repoint and the missing `py` rows are separate reviewed changes. ⛔ Do not let W0's success be described as "R4 closes". |
| **8** | ⭐ Record the R5 / S11 deferral per **D8**. |

## ⚠ Sequencing this exposed — ⭐ the repoint may be the bigger half

⭐ Given **D2** and **D3**, a **checks repoint** is now a first-class piece of work, not a footnote:
bind the squared-speed objects that already exist, add the missing `py` rows, and decide the `D` coverage.
⚠ It is a separate reviewed change and it touches two **closed** steps' configs.
⛔ Do not fold it into the emission build — that is how one change becomes unreviewable.

## Rules for whoever writes the directive

- ⛔ Do not restate the physics. ⭐ Every factual claim above came from a command; ⛔ do not add one that did
  not.
- ⛔ Do not supply an expected value, a residual, or an acceptance number anywhere in the builder statement.
- ⭐ Name objects. ⛔ Where you find yourself writing *how*, name *what* instead.

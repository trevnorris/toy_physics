# W0 — the engines emit the wave speed, not only its square

⚠ **Warning: `steps/S9_PILOT_ADJUDICATION.md` holds S9 result values. Don't read it.**

⛔ **This is the shared object statement for FOUR engine builds** (S9-wl, S9-py, S10-wl, S10-py). Each
engine is changed in its own build, by its own builder, from this one statement. ⛔ No builder reads
another engine's change.

## The measured gap

`reduction/relations.yaml` states, at `R4`:

```
residual: Sub  [Q, Q.brane.c_gamma]  [Sqrt, [Div, [Q, Q.brane.mu_R], [Q, Q.brane.rho_br]]]
```

`Q.brane.c_gamma` is declared with dimension `[1,-1,0]` — **a speed**.

But every engine emits only the **square**, and the checks bind the squared tag to the linear quantity:

- `reduction/checks_S9.yaml:194` — `Q.brane.c_gamma: {tag: WL_S9_CANDIDATE_SPEED_SQUARED1}`
- `reduction/checks_S10.yaml` — `Q.brane.c_gamma: {tag: WL_S10_MAIN_D3_Q3_DISTINCT_ROOTS, select: sequence_second}`

⇒ the residual evaluates `(μ_R/ρ_br) − √(μ_R/ρ_br)`, which is **nonzero by construction**.
⇒ ⛔⛔ **`R4` has never been tested at S9 or at S10.** It is reported red at both, carried as a note at S9
and unmentioned at S10. `R5`/`c_L` has the same shape via `B_comp`.

## ⭐⭐ The object each engine must emit

⭐ **The speed of the transverse brane mode — as a speed — obtained as the non-negative root of the
squared speed THAT ENGINE ALREADY DERIVED**, together with the criterion by which that mode was identified
as the transverse one, and the premise under which the root was taken.

⛔ **Name the object, ⛔ not the recipe.** How you take the root, where the emission sits, and what you call
the tag are yours to decide from the existing code.

### ⛔⛔ THREE WAYS THIS EMISSION BECOMES WORTHLESS — all three have happened in this repository

1. ⛔⛔ **Constructing the speed from the registry's own formula.** If the emitted speed is built as
   `sqrt(mu_R/rho_br)`, the `R4` residual is zero **by construction** and tests nothing. ⭐ The emitted
   speed must descend from **the dispersion relation this engine solved** — the determinant, the roots, the
   `root/q` the engine already computes — ⛔ never from `μ_R` and `ρ_br` recombined by hand.
   ⭐ **Emit the squared speed and the speed as separate tags from the same object**, so a reader can see
   which one produced the other.
2. ⛔⛔ **Selecting the transverse candidate by matching a known answer.** The engines emit a *list* of
   speed candidates. ⛔ Do not pick the one that makes `R4` close. ⭐ Select by a stated **physical**
   criterion — polarization, multiplicity, whatever this engine already computes — and ⭐ **emit that
   criterion as its own tag, read out of the computation, ⛔ not typed beside it.**
   ⚠ **Measured:** engine 1 shipped **nine** defects in one class — a tag declaring what the run used,
   assembled from a literal beside the computation instead of read out of it. ⛔ Do not add a tenth.
3. ⛔ **Taking the root without a premise.** A symbolic square root needs the radicand's sign. ⭐ Emit the
   positivity premise the root was taken under, **read from the engine's own assumption set**. ⛔ If the
   engine cannot establish the sign from its declared premises, that is the honest outcome — ⭐ emit a
   pinned failure object and say so; ⛔ do not assume it.

## ⚠ What is NOT in scope

⛔ Do not change any derivation. ⛔ Do not change any existing tag's value, name, or shape. ⛔ Do not
change the dispersion relation, the action, the ansatz, or the control matrix. ⭐ This is an **additive**
emission plus whatever the addition forces.
⛔ Do not touch `reduction/registry_dimension_witness*` or `reduction/dimension_laws*`.
⛔ Do not edit `reduction/checks_S*.yaml` — the repoint is a separate reviewed change.
⛔ Do not commit.

## ⛔⛔ THE REGRESSION BAR

⭐ **Every value the engine emits today must be unchanged.** Re-run the engine and diff its full output
against the committed `.out`. ⭐ Report the diff literally; it must contain **only additions**.

⛔ **If any existing value moves, STOP and report it.** ⚠ It is either a defect you introduced or one this
change exposed — ⭐ which one it is, is the finding. ⛔ Do not adjust anything to restore the old output.

⚠ **And an unchanged baseline can be worth nothing.** Measured here: three variants — correct, unfixed,
and knowingly wrong — produced byte-identical stdout because no separator character occurred in any
identity field. ⭐ Show the new tags actually appear, with their values, in the re-run output.

## ⭐ What must become catchable

⭐ **Show it, with literal stdout, ⛔ not prose:**
1. Perturb the derived squared speed upstream — the new speed tag must move with it. ⛔ If it does not, the
   emission is disconnected from the derivation and everything above is void.
2. Perturb the selection criterion so a different candidate is chosen — show what changes.
3. ⭐ **Build the WEAKER version of your own fix** — the speed constructed from `μ_R/ρ_br` directly — and
   show what distinguishes it from yours. ⚠ Three times in this repository a half-fix passed the new test,
   the full suite, the whole battery **and** produced byte-identical output. ⛔ A test that "covers" the
   invariant on one example is a demo, not a pin.

## Rules

- ⭐ Print computed objects — both operands and the residual — then guard. ⛔ No `assert` before a value it
  guards. ⛔ No status typed rather than computed. ⛔ No conclusion as an unconditional literal.
- ⛔ **Do not tune anything to make a run green.** ⭐ A disagreement is a finding and is the deliverable.
- ⚠ Wrap every kernel run in `timeout 600`; ⛔ never raise it; ⛔ never more than one Mathematica kernel at
  a time — the licence has **two** seats.

## Deliverables

1. The engine change, and whatever it forces.
2. ⭐ The full-output diff against the committed `.out`, literal, showing additions only.
3. ⭐ The three able-to-fail runs above, with scripts at named absolute paths and their literal stdout.
4. ⛔ A ≤40-line report. ⛔ No conclusions about whether the physics is right.

# Consultation — how should the registry express a D-dependent dimension?

⛔ **This is not a build request. Do not write or modify any file.** I want your reasoning and a
recommendation with consequences. ⚠ Warning: `steps/S9_PILOT_ADJUDICATION.md` holds S9 result values —
don't read it.

## The measurement that raised this

A new instrument (`reduction/registry_dimension_witness.py`) compares each engine's **emitted** dimension
vector for a registry quantity against the **declared** exponents in `reduction/quantities.yaml`. Run
against committed S9/S10/S11 outputs it reports, among other things:

```
BRANCH_DIMENSION_MISMATCH = 4      (S10-wl and S10-py, for Q.brane.rho_br and Q.brane.mu_R)
```

with residuals walking across S10's brane-dimension branches:

```
rho_br  emitted [-D,0,1]      declared [-3,0,1]
        D=2 -> residual [ 1,0,0]
        D=3 -> residual [ 0,0,0]
        D=4 -> residual [-1,0,0]
        D=5 -> residual [-2,0,0]
```

## What I believe is true — ⭐ check me, don't accept it

- **The engines are correct.** A brane inertia density is mass per unit **D-dimensional** volume, so
  `[rho_br] = M·L^-D`. Likewise `[mu_R] = [2-D,-2,1]`.
- **The registry is correct at D=3 and silent about that.** `quantities.yaml` declares
  `exponents: [-3,0,1]` with convention `LTM-exponent-vector-v1` and no D-scope.
- ⇒ **Neither is wrong. The registry cannot express the D-dependence**, and nothing anywhere records that
  its brane dimensions are D=3 specialisations.
- ⭐ **And the D cancels in the derived quantity**: `[mu_R] − [rho_br] = [2-D,-2,1] − [-D,0,1] = [2,-2,0]`,
  so `[c_gamma] = [1,-1,0]` **independently of D**. ⚠ Verify this; if it holds it is a physics statement
  the registry currently expresses only by accident.

## The constraints any answer must live with

- `reduction/relations.yaml` holds **5** relations (R1, R2.h0, R2.xi_h, R4, R5) over **14** quantities in
  **2** sectors (`Q.medium.*`, `Q.brane.*`). ⭐ Read them.
- **~25 derivation steps remain** across light, gravity, gravitomagnetism, charge and magnetism. The
  registry must grow one step at a time, and whatever is decided here is paid on every one of them.
- `reduction/dimensional_homogeneity_gate.py` consumes these declarations and currently exits `PASS`.
  ⭐ Check what each option does to it.
- S10 varies D deliberately — the brane dimension is the thing that step derives a result *from*.
  ⚠ S9 and S11 appear to work at fixed D. ⭐ Confirm.
- This is a two-person project with no infrastructure team. ⛔ A schema change that costs a week must save
  more than a week.

## What I want from you

⭐ **Enumerate the real options and their consequences.** I can see at least these; ⭐ add any I have
missed, and ⛔ tell me if one of mine is incoherent:

1. **Exponents as expressions in a declared symbol** — `exponents: [-D, 0, 1]` with `D` bound to
   `Q.brane.D_brane`. What breaks? What has to learn to parse expressions?
2. **Scope the registry to D=3 explicitly** — add a declared specialisation field, and treat D≠3 branches
   as outside the registry's remit. ⚠ What does that cost S10-like steps that vary a parameter?
3. **Declare the functional form and let each consumer specialise** — the registry carries the law, the
   witness/gate substitutes per branch.
4. **Something else.**

⭐ **For each option, state concretely:** what changes in `quantities.yaml` · what changes in
`relations.yaml` · what each consumer (`registry_read.py`, `dimensional_homogeneity_gate.py`,
`engine_output_checks.py`, the new witness) must learn · what it costs per future step · and ⭐ **what a
wrong declaration would look like and whether anything would catch it.**

⭐ **Then recommend one, and say what would have to be true for you to be wrong.**

⚠ **Two things I especially want tested, not assumed:**
- ⛔ Is the D-cancellation in `[c_gamma]` general, or an accident of these two quantities? ⭐ Check it
  against R1, R2.h0, R2.xi_h and R5 as well.
- ⛔ Does any option make a **wrong** dimension declaration harder or easier to detect? ⚠ The one defect
  this project measured that survived two review legs and a full ablation suite was a **wrong dimension**.
  ⭐ An option that makes that class harder to see is disqualified regardless of its other merits.

⛔ Do not give me a plan, a schema draft, or code. ⭐ Give me the trade-offs and a recommendation, with the
reasoning visible.

# Phase 0 · Step 2 — the sound speed

**Upstream:** step 1 (the medium, the EOS closure `P = Kρⁿ`, the primitives).
**Sources:** `docs/model_map.md:62`, `notes/stages/ledger_stage004_gnls_action_dimensional_foundation.md:65`.

---

## What it is

The phonon / density-mode speed, obtained from the postulated EOS by the standard superfluid relation:

```
P = K ρⁿ                            (step 1, postulated closure)
c_s² = (1/m) dP/dρ = n K ρ^(n-1) / m
```

At `n = 5`: `c_s² = 5 K ρ⁴ / m` — matching `model_map.md:62` and the registry's `R1`.

## What it does

Gives the medium its **density-mode signal speed**, which is the gravity-change speed
(`model_map.md:62`). It is also the denominator of `λ_γ = c_γ/c_s`, so every later cone ratio is
measured against it.

## What's new

⭐ **Nothing.** `c_s` is pure consequence of step 1's closure plus the primitives `{n, K, ρ0, m}`.

⇒ This is what a derived step looks like: it adds an equation and **zero** inputs. Confirmed by
`show_reduced.py` — `c_s0` reduces to `sqrt(5)*rho0**2*sqrt(K/m_GNLS)`, resting on `K`, `ρ0` and `m` and
nothing outside step 1. ⚠ Note the **literal `5`**: `n_eos` is ⛔ **not** a symbol in that expression —
see the corrected consequence 3 below.

---

## ⭐⭐ The finding: `n` and `[K]` are not independent

Deriving forward exposes a coupling that a hardcoded literal hides.

```
[P]   = M L⁻² T⁻²           (4D pressure = energy / 4-volume)
[ρ]   = L⁻⁴
[K]   = [P] / [ρ]ⁿ = M L^(4n−2) T⁻²
```

At `n = 5` this is `M L¹⁸ T⁻²` — **exactly** the `[18,-2,1]` the registry carries for `K`.

⇒ ⛔ **The polytropic exponent is not a continuously-variable knob.** Changing `n` changes the
*dimension* of `K`. The two are one structural choice, not two independent inputs.

**Consequences, recorded rather than patched:**
1. `n_eos` is a **discrete structural** input, ⛔ not a continuous one. It must still be supplied for a
   simulation to run, so it stays in the sim-input set (§1.0) — but in the **discrete** partition.
   ⚠ This is why §1.0 partitions by kind: summing it into a continuous dimension would be wrong.
2. `R1`'s literal `4` **is** `n−1`, and is correct *because* `n` is fixed. ⛔ Do not "generalise" it to
   `ρ^(n−1)` — that would imply `n` is variable, which the `[K]` coupling forbids, and would make the
   dimension symbolic and undeterminable.
3. ⚠ **A latent trap for anyone who edits the registry** — ⛔ **corrected 2026-07-30; the earlier wording
   here misdescribed the registry.** It claimed `n_eos` "now appears symbolically in `R1`'s numerator".
   **It does not.** `R1`'s residual (`relations.yaml:12`) stores **integer literals** `5` and `4`, and
   `input_qids` is `[K, rho0, mass]` — `n_eos` is **absent** from it. `n_eos` appears only in `R1`'s
   `literal_consistency` block, which **asserts** that the coefficient node equals `n_eos` and the
   exponent node equals `n_eos − 1`. ⇒ The literals are pinned to `n_eos` by a **consistency assertion,
   ⛔ not by symbolic substitution**; nothing recomputes if `n_eos` changes. Setting `n_eos` to anything
   but 5 therefore trips that assertion while leaving an unchanged relation that is now both
   algebraically and dimensionally wrong. The coupling above is the reason.

---

## Checks

| # | check | result |
|---|---|---|
| 1 | dimensional homogeneity | ✅ `R1: HOMOGENEOUS dimension=[1,-1,0]` — velocity, as required |
| 2 | what's-new classification | ✅ nothing new; pure consequence |
| 3 | input provenance | ✅ every input traces to step 1, cited by locus |
| 4 | step-to-step identity | ✅ `K`, `ρ0`, `m` — `R1`'s three `input_qids` — are the same quantities step 1 introduced. ⚠ **Corrected 2026-07-30:** `n_eos` was listed here as an `R1` input and is **not** one; it enters `R1` only through the `literal_consistency` assertion on the literals `5` and `4`. Its step-1 identity is checked separately, as the discrete structural choice |
| 10 | term-by-term fidelity | ✅ `c_s² = (1/m)dP/dρ` derived from step 1's closure by hand above; the factor `n` and the power `n−1` both traced to `d(ρⁿ)/dρ = nρ^(n−1)` |
| 11–12 | dynamical health · approximation validity | n/a — no dynamical block, no truncation |
| 13 | second physics leg | ⏸ deferred to the phase-0 close; nothing derived here beyond one textbook differentiation |

⚠ **Provenance note surfaced by the gate:** `c_s0`'s dimension comes from **stage005**, which is
`shared_module=no` — a stage-local declaration, not the shared import. Recorded, not treated as wrong;
it is exactly the fan-out gate the plan describes.

---

## ⛔ Two defects this step exposed in the machinery

### A. The earlier "acceptance MATCH" was the right number for the wrong reason

`a` was marked `convention-orbit` on the **quantity** while its relation `R2.a_pin` had been set to
`DERIVED-EXECUTED`. So `a` was excluded from the ambient count, leaving 9 continuous quantities — and
adding `n_eos` as continuous brought the total back to 10. **Two errors cancelled**, and the resulting
match was reported as confirmation. ⚠ It was not.

Both entries now say derived, and `n_eos` sits on the discrete axis where it belongs. Ambient is 10 for
the right reason.

> ⛔ **REOPENED 2026-07-30, after this step was banked.** The half of that sentence that reads *"both
> entries now say derived"* is exactly what **D-01** put back in question
> (`research/pde_ledger_v2/walkthrough/DECISIONS.md`): a relation arising from imposing unit pins is not
> a defining equation, so `R2.a_pin`'s `DERIVED-EXECUTED` status — and with it "ambient is 10 for the
> right reason" — is **OPEN pending a user decision**. ⛔ Do not resolve it here. ⚠ The `n_eos`
> half stands: the discrete axis is correct and is not affected. ⭐ The **lesson** below is untouched by
> any of this and is the durable part of this record.

⭐ **Lesson: a matching number is not evidence.** The Δ column was the tell — reclassifying an axis
changed `dim_before` and `dim_after` together while leaving every Δ untouched, which is exactly what a
bookkeeping change (not a physics change) looks like.

### B. ⭐⭐ The rank computation counts an entailed relation as independent

`C-M2` adds `ξ_h² − 2a²`. Substituting the registry's own definitions:

```
ξ_h = √2·ħ/(c_s0·m),   a = ħ/(c_s0·m)   ⇒   ξ_h² − 2a² = 0   identically
```

⇒ The added relation is **entailed**; the correct `dim_after` is **5**, unchanged. The registry reports
**4** — it counts a redundant relation as an independent constraint.

⚠ **This is the off-locus Jacobian weakness already recorded** at phase 1b (*"the off-locus symbolic
Jacobian helper reports treatment-(i) C-M2 as 4; the exact satisfying-witness rank is 5"*). It is no
longer latent.

⛔ **Direction matters: it biases flattering.** A `dim_after` that is too low makes the model look more
determined — fewer free parameters — than it is. ⇒ The fix is to evaluate the Jacobian at an **exact
constraint-satisfying witness** rather than off-locus. **Not fixed here**; carried as the next
machinery task, because it is a separate piece of work from a derivation step.

**Acceptance is left at `DISAGREEMENT_UNRECONCILED`, exit 1, on C-M2 alone.** Three of four cases match.
⛔ Do not tune the registry to close it — the fixture is right and the rank helper is wrong.

> ✅ **RESOLVED 2026-07-30 in `6e79e2ec`, after this step was banked.** ⛔ **The finding above is left
> exactly as written** — it was true when banked, and the record of *how* it was found is the point.
> **What changed:** the rank is now evaluated at **exact constraint-satisfying witnesses** instead of
> off-locus. C-M2's `dim_after` is **5**, `acceptance_check.py` reports `PHASE1_ACCEPTANCE: MATCH` and
> exits **0** on all four medium cases, and ⭐ **two independent reviewers derived those four numbers
> themselves** rather than reading them off. ⭐ **The "⛔ do not tune to close it" instruction held and
> was honoured:** the rank helper was repaired; the fixture was not moved. ⇒ This is no longer an open
> defect anywhere — `docs/derivation_walkthrough_plan.md` §4 no longer carries it forward as latent.

---

## Open, carried forward

1. ✅ **CLOSED — the schema has the discrete counting axis.** This item read *"`counting_axis` is
   `[continuous-model, convention-orbit]`; `n_eos` needs `discrete-structural`. One enum value, added
   below"* — and ⛔ **nothing followed it.** The enum value **exists**:
   `registry_schema.yaml:139` is `enum: [continuous-model, convention-orbit, discrete-structural]`, and
   `Q.medium.n_eos` carries `discrete-structural` (`show_reduced.py` reports it as a *structural
   selection*, separate from the five continuous simulation inputs). Nothing outstanding.
2. ⭐ **Is `D = 4` the same kind of coupling?** Step 1 records bulk dimensionality as a discrete choice.
   `[ρ] = L⁻⁴` and `[P] = M L⁻² T⁻²` both already encode `D = 4`. ⇒ Likely the *same* situation as `n`:
   a structural choice that fixes other dimensions. **Worth checking explicitly at the next step that
   uses a volume integral.**

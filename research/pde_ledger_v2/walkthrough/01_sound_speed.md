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
(`model_map.md:34`). It is also the denominator of `λ_γ = c_γ/c_s`, so every later cone ratio is
measured against it.

## What's new

⭐ **Nothing.** `c_s` is pure consequence of step 1's closure plus the primitives `{n, K, ρ0, m}`.

⇒ This is what a derived step looks like: it adds an equation and **zero** inputs. Confirmed by
`show_reduced.py` — `c_s0` reduces to `√K·√n_eos·ρ0²/√mass`, resting on nothing outside step 1.

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
3. ⚠ **A latent trap for anyone who edits the registry:** `n_eos` now appears symbolically in `R1`'s
   numerator while `n−1` is baked into the exponent. Setting `n_eos` to anything but 5 yields a relation
   that is both algebraically and dimensionally wrong. The coupling above is the reason.

---

## Checks

| # | check | result |
|---|---|---|
| 1 | dimensional homogeneity | ✅ `R1: HOMOGENEOUS dimension=[1,-1,0]` — velocity, as required |
| 2 | what's-new classification | ✅ nothing new; pure consequence |
| 3 | input provenance | ✅ every input traces to step 1, cited by locus |
| 4 | step-to-step identity | ✅ `K`, `ρ0`, `m`, `n_eos` are the same quantities step 1 introduced |
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

---

## Open, carried forward

1. **The schema has no discrete counting axis.** `counting_axis` is `[continuous-model,
   convention-orbit]`; `n_eos` needs `discrete-structural`. One enum value, added below.
2. ⭐ **Is `D = 4` the same kind of coupling?** Step 1 records bulk dimensionality as a discrete choice.
   `[ρ] = L⁻⁴` and `[P] = M L⁻² T⁻²` both already encode `D = 4`. ⇒ Likely the *same* situation as `n`:
   a structural choice that fixes other dimensions. **Worth checking explicitly at the next step that
   uses a volume integral.**

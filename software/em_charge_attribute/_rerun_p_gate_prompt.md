# Confirm the v6 basis fix, then (only if sound) re-run the gate to a decisive verdict

The v5 build of `software/em_charge_attribute/directive_native_p_constraint_gate.md` returned **`GATE_ILL_POSED`** — both engines agreed all constraints are second-class with **no additional first-class `G`** on the regular stratum (THEORY-A 8/8, THEORY-C 12/12), but the directive's demand for a *unique exact-quartic basis modulo unrestricted nonlinear field redefinitions* is ill-defined (a nonlinear point map sends quartic `x²(x−1)²` → degree 8), so the exhaustive tuned-vs-generic classification wasn't defined. **v6 folds a fix.** Reasoning effort: **xhigh**; `--sandbox danger-full-access` is set (Mathematica may run; ONE kernel at a time).

## STEP 1 — CONFIRM the v6 fix is mathematically sound (gate; do NOT touch scripts yet)
The v6 directive (re-read it) replaces the ill-posed basis with:
- an **ORDER-BY-ORDER basis**: ≤2 derivatives, field-amplitude order-by-order from quadratic, mod IBP + **perturbative (order-by-order)** field redefinitions + holonomic constraints — which closes at each order (a redefinition `φ→φ+δφ`, `δφ` higher-order, only shifts operators upward);
- a **decisive-linearization lever**: a first-class Gauss constraint of the full nonlinear theory necessarily has a nontrivial linearization, so **ABSENCE of first-class Gauss at leading (quadratic) order ⇒ decisive `NATIVE_P_NO_EMERGENT_GAUSS`**; a leading-order PRESENCE needs nonlinear-persistence confirmation; genuine nonlinear-persistence ambiguity ⇒ `GATE_ILL_POSED` (that branch only).

**Confirm three things, concretely:** (1) is the order-by-order quotient now genuinely well-posed/finite at each order? (2) is the linearization-necessity lever correct — is quadratic-order ABSENCE of a first-class Gauss truly decisive for the full nonlinear no-go (i.e. can a nonlinear gauge symmetry exist with trivial linear part on these fields)? (3) does the fix resolve the v5 ill-posedness WITHOUT introducing a new false-pass/false-fail?

**If any is wrong → STOP, change NOTHING, and report `V6_FIX_NEEDS_WORK` + the specific problem.** Only if all three hold, proceed.

## STEP 2 — update the scripts minimally and RE-RUN (only if STEP 1 = sound)
Update `native_p_gate_sympy.py` + `native_p_gate_dual.wl` + comparator to implement the order-by-order basis + the linearization lever (keep everything else — the six controls, per-tooth ablation, `ENGINE_AGREE`, the frozen fields/theories — unchanged). Re-run `run_native_p_gate.sh` and iterate to green. **Do NOT weaken anything to force a verdict**; a genuine `GATE_ILL_POSED` or a surprise `FIRST_CLASS_*` must be reported honestly.

## Output
Update `software/em_charge_attribute/reports/native_p_constraint_gate.md` with the decisive per-theory VERDICT (expected `NATIVE_P_NO_EMERGENT_GAUSS` given the v5 all-second-class regular-stratum result + the linearization lever, but report what the computation actually gives), the order at which each verdict is decided, the six controls, and dual-engine `ENGINE_AGREE`. State clearly at the top whether STEP 1 confirmed the fix (`V6_FIX_SOUND`) or not.

# Adversarial audit — the GENUINE-REBUILD of the native-`Pᵃ` constraint-class gate

You are a FRESH, INDEPENDENT adversarial verifier. A PRIOR build of this exact gate was a **pass-by-construction RIG** — it hardcoded the decisive verdict (`additional_G_exists=False`), built **no** coupling-carrying Hamiltonian (the `{g_a}` were only display strings), FAKED the Maxwell control (`cmat = sp.zeros(2)`), and used the tautology `solve([t=−1,t=+1])=∅` as its "search." That fake produced the SAME verdict this rebuild now claims. **Your job: determine whether the REBUILD is a genuine computation, or whether it still fakes/shortcuts the decisive step.** A "can't-do" verdict that matches the builders' prior is the seductive-trap case — verify HARDER than a positive.

## Context
Toy 4D superfluid analog. Dirac–Bergmann "gate": does the native little-arrows sector (hard O(4) sigma `|P|=1` + MacCullagh curl-elastic shear `u` + brane field `χ_B`, most-general order-by-order couplings `{g_a}`) develop an ADDITIONAL first-class local `U(1)` Gauss constraint (emergent EM), or is it purely second-class (no native EM)? Decided at quadratic order (a nonlinear gauge symmetry must have a nonzero linearization → quadratic absence is decisive). Rebuild verdict claimed: **`NATIVE_P_NO_EMERGENT_GAUSS`** (both theories), with all six controls passing through a shared pipeline and both computed guards passing.

## Files (under /var/projects/toy_physics/software/em_charge_attribute/)
- directive_native_p_gate_rebuild.md   (the anti-fake requirements)
- directive_native_p_constraint_gate.md (v6 physics spec)
- native_p_gate_sympy.py                (SymPy engine — the main audit target)
- native_p_gate_dual.wl                  (independent Mathematica engine)
- native_p_gate_compare.py               (comparator)
- run_native_p_gate.sh                   (runner + ablation harness)
- reports/native_p_constraint_gate.md    (report + REBUILD NOTE claiming Q1–Q6 closed)
- reports/native_p_gate_artifacts/*.log

## Determine (concrete, cite file:line; assume nothing the report claims)
1. **Genuine coupling-carrying Hamiltonian?** Does `build_H2` actually differentiate the input Lagrangian → Hessian → Legendre map → constraints, with the `{g_a}` as FREE SYMBOLS that genuinely appear in the computed momenta/Hessian/PB matrix? Or is any of that hardcoded/short-circuited? Verify `coupling_guard` genuinely FAILS if a coupling drops (not a no-op).
2. **Do the `{g_a}` actually reach the decisive bracket/kernel?** The whole point. Confirm the couplings enter the `search_G` input, not just the constraint bookkeeping.
3. **Do ALL SIX controls run the IDENTICAL code path as native-`P`?** Is Maxwell's zero-PB / `FIRST_CLASS_GAUSS` a COMPUTED output of the same `build_H2→dirac_search→search_G` (not a literal / special case)? Are the varied control results (Maxwell `FIRST_CLASS_GAUSS`, gauged-hard-unit `MIXED`, bare-sigma `SECOND_CLASS_RADIAL`, nonconserved `INCONSISTENT_PRESERVATION`, gauge-fixed `SECOND_CLASS_NO_LOCAL_GAUGE`, global-only `GLOBAL_CHARGE_NO_LOCAL_GAUSS`) genuinely computed?
4. **Is `search_G` a real enumeration that CAN return nonzero?** Confirm it is not structurally forced to 0. `GUARD-SEARCH-CAPABLE` (Maxwell + gauged-hard-unit must yield nonzero candidates) is the safeguard — verify it is real and would fail if the search were blind.
5. **The tuned-strata finding — real or a dodge?** The rebuild reports first-class directions on tuned coupling strata (`FC=2`) but dismisses them as "transverse vector-shift primaries with zero descendants — not Gauss chains." VERIFY this classification is genuinely computed (a shift symmetry whose Hamiltonian descendant / Gauss chain is actually absent), not a hand-wave to discard a first-class direction that would otherwise be an EM candidate. **This is the most likely place a real EM signal would be swept under the rug — scrutinize it hardest.**
6. **Independent dual engine?** Does the Mathematica engine differentiate its OWN Lagrangians and run its OWN Dirac/kernel code, or read SymPy's result? Does the comparator check independently-computed outputs, not shared literals?
7. **Residual fakes?** Any remaining hardcoded verdict, `require(rank==len)` crash-guard, blanket `SECOND_CLASS` dict, tautological assert, or `X≡X`-on-literals ablation?
8. **Transliteration fidelity + artifact check:** does the coded `H₂` match the physics (hard `|P|=1` via multiplier, MacCullagh `(∇×u)²`, frozen couplings)? Any sign/convention/dropped-term bug that would spuriously zero the first-class count?

You MAY re-run the SymPy engine unchanged (`python3 native_p_gate_sympy.py ...`) and inspect output; you MAY add your OWN scratch probe that imports its functions to test whether a deliberately-gauged input yields nonzero candidates. **Do NOT run Mathematica / wolframscript / run_native_p_gate.sh** (2-seat license; a re-run is happening elsewhere). Read the committed Mathematica logs.

## Return
Verdict **VERIFICATION_CONFIRMS** (genuine computation; the no-go is trustworthy) or **VERIFICATION_ISSUES** (found a fake/shortcut/artifact), then per-item findings with file:line, any concrete defect, and a bottom line. Be adversarial and specific; do not pad. If genuine, say so plainly; if not, show the defect.

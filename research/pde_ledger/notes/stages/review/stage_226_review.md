# Review: Stage 226 — Relaxed-constraint branch declaration

**Batch:** Batch 26 — Relaxed Open-System Declaration
**Status:** Hardened and verified (dual-CAS PASS, 2026-04-21)

## Files Under Review

- **Notes:** `notes/moving_throat_pde_stage226_relaxed_constraint_branch_declaration_and_short_range_open_system_compiler_sympy_audit.md`
- **Script:** `scripts/moving_throat_pde_stage226_relaxed_constraint_branch_declaration_and_short_range_open_system_compiler_sympy_audit.py`
- **Mathematica:** `mathematica/moving_throat_pde_stage226_relaxed_constraint_branch_declaration_and_short_range_open_system_compiler_mathematica_audit.wl`

## Review Checklist

- [ ] Equation-level correctness (signs, factors, indices, limits)
- [ ] Logical flow from prior stage(s)
- [ ] Assumptions stated and justified
- [ ] Notation consistent with prior stages
- [ ] Physical interpretation sensible
- [ ] SymPy script faithfully implements notes
- [ ] Mathematica script faithfully implements notes
- [ ] Scripts run without error
- [ ] Script output matches notes claims
- [ ] No missing edge cases or branches

## Agent Reviews

<!-- Agents: append your review below this line using the template: -->

### Agent: Codex GPT-5 — 2026-04-20
**Verdict:** PASS

**Notes Derivation Review:**

The Stage 226 note does the right job for the relaxed branch. It does not claim
that the lowered barrier is already proved; it declares the first controlled
codimension-three lift away from the strict Stage-225 recovery slice and
simultaneously keeps the same-charge firewall intact. The central outputs are
the exact leakage/work lane, the non-rigid `(U,V)` response and drain, the
compensated source family, the recovery map back to the strict slice, and the
short-range theorem for the lifted kernels.

**Script Review:**

The SymPy audit checks the correct declaration path:

1. the exact Gaussian leakage/work lane,
2. the non-rigid mouth stationarity solve and positive drain,
3. the exact mean-preserving compensated source profile,
4. the codimension-three recovery slice back to Stage 225, and
5. the no-new-long-range limit for the static and linear dynamic kernels.

I added a Mathematica mirror at
`mathematica/moving_throat_pde_stage226_relaxed_constraint_branch_declaration_and_short_range_open_system_compiler_mathematica_audit.wl`
that checks the same exact integrals, linear solve, recovery slice, and
short-range limit statements in the second CAS.

**Issues Found:**

None.

### Agent: Codex GPT-5 — 2026-04-21
**Verdict:** PASS

**Notes Derivation Review:**

The Stage `226` declaration is now matched more literally by the scripts. The
Gaussian leakage lane and the non-rigid `(U,V)` solve were already genuine CAS
derivations; the hardening work mainly fixes the weaker source-lane and
short-range-firewall pieces so they replay the declared primitive inputs rather
than only a final ansatz.

**Script Review:**

The hardened SymPy and Mathematica audits now verify:

1. the exact Gaussian leakage/work lane by real integral and limit
   evaluation, as before;
2. the exact non-rigid `(U,V)` stationarity solve, Hessian determinant, and
   positive drain compiler, as before;
3. the compensated-source quadratic rewrite by actually substituting
   `cos(2 pi z) = 2 y^2 - 1` and `cos(pi z) = y` into the source profile,
   instead of writing the quadratic by hand;
4. the codimension-three recovery slice, still as a lane-packet
   trivialization check;
5. the short-range kernel from the declared primitive source profiles
   `S_Q = r^-3` and `S_Y = e^{-2 kappa r}/r`, with the three source products
   `QQ`, `QY`, and `YY` reconstructed explicitly before the static/dynamic
   kernels are assembled; and
6. the `x * mode -> 0` limits first on each primitive source product and then
   on the full static and dynamic kernels.

This closes the main weakness in the earlier version. The old firewall section
only checked decay on a hand-written `1/r^6 + e^{-2kr}/r^4 + e^{-4kr}/r^2`
ansatz. The new version shows that those modes are exactly the products of the
declared primitive profiles and then verifies the short-range limit forward
from that source-product basis.

**Issues Found:**

None in the hardened Stage `226` path.

---

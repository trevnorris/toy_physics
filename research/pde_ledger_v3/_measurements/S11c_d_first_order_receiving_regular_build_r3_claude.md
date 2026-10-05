CLEAR FOR THIS REAL-AXIS RECEIVING AND CROSS-FLUX BUILD

This is a source/build assessment only, not runtime acceptance. I read all four files (`worker.py`, `build-notes.txt`, `method-plan.txt`, `guide.txt`). I did not open any of the roughly 2,075 saved JSON evidence files. I did not see the shared `Journal` helper, so I did not check what `J.zero` and `J.nonzero` do on failure. Nothing was run.

**What holds up, from reading the code**
- **Row selection and projection:**
  - The full 5×2 `D5*Feta` is compared with the hand table before rows 2–4 are used (`worker.py:259-261`).
  - `R3` and `E3` are joined to literal matrices (`:225-226`).
  - `R3(p)*9U = 0` is checked (`:263`).
  - `R3*U = iQ/(5K2)` is checked before any cusp (`:277`).
  - I confirmed by hand that `A' = (9/2+4T)(1-T²)/L`, that the endpoint jump is 9, and that `Q·Â = -iÂ'` matches the forward-transform convention.
- **Moving projector:** the δ' product rule `fδ' = f(p)δ' − f'(p)δ` is applied with the correct sign (`:375`). The declared control entry [0,0] has a zero baseline and a mutant that drops `R3'U`.
- **Pressure and chemical normalization:** the chemical tag enters once through `10·(A_mem/ρ)` (`:345-346`). Normal signs follow the face, and `A_mem` and effective `cs` come from saved operands, not the original `c_s0`. The depth law uses `ω²/cs² − |h|² − l²`, with `q=0` at `l=±p`.
- **Domains and roots:**
  - Every original negative-power base is kept before `together` or `cancel`, transported onto `q`, and tested on both rays. The same tests cover chart K2, the entry and raw-determinant denominators, and forcing denominators.
  - The Euclid/Bezout/Sturm logic is correct: a constant gcd rules out common roots, and the count `V(0)−V(upper)` is valid because endpoints are checked nonzero separately.
  - Unresolved signs raise rather than pass.
- **Controls and cross current:** the three controls change the real operands at declared entries. Both cross contractions use the correct leg order, and the conjugacy check covers `J(k,k')^H = J(k',k)`.
- **Growth ledger:** positive exponents add across R3, the pressure factors, the eta multiplier, E3 and the multiplier. Nothing in the ledger buys decay through cancellation.
- **Basis joins:** the dual joins `D5·E5` and the transformed matrix joins `D5·physicalMatrix·E5`. The raw LEFT operand stays in the physical weak basis and is not compared with the D5/E5 basis.

None of the findings below can produce a false pass, because each one ends in a refusal. They are not blockers.

1. **Weak `l²→s` substitution (`:426`, `:484`, `:530`).** `xreplace({l*l:s})` only matches `l**2` nodes. An `l**4`, or an odd power inside a base, leaves `l` behind and trips `UNRESOLVED_EVEN_BLOCK_REDUCTION`. That would use up the single authorized run. Using `subs(l, sqrt(s))` after the evenness check, or polynomial reduction, would be sturdier. This is optional.
2. **Sturm signs at `p` (`:453`, `:455`).** `exact_sign` relies on sympy deciding `.is_positive` for expressions in `sqrt(595)`. It raises if sympy cannot decide, so it cannot give a wrong answer. This is optional.
3. **Original flux function binding (`:399-414`).** Only the six assignment ASTs are matched. The definitions of `U`, `B`, `l`, `p` and `km`/`kp` in that function are not matched literally. They are tied to the live objects only indirectly, through the `B(p)==U` join at `:148` and `:397`. A literal AST match on those definitions would complete the inherited-frame join. I do not count this as a blocker, since the guide says the completed diagonal contractions should not be replayed.
4. **Hard-coded result flags (`:547`).** `projectedUBThreeRowsZero`, `crossCurrentZero` and `realAxisC3PoleExclusion` are written as constants. They are only valid if `J.zero`, `J.nonzero` and the `require` calls raise on failure. The plan says they do, but I could not confirm that for the helper.
5. **Cross-current stop is all-or-nothing.** If either cross contraction is nonzero, `J.zero` stops the whole run before the C3 determinant is computed. That matches the plan's "save it and stop" rule, but the receiving result is then lost as well. A saved cross current plus a flagged continuation would keep both outputs.

**Still required after a pass:** nonuniform physical work balance. That means exterior radiation, memory, mass-rate/chemical work and LAB_HELD work. A pass here does not give a leakage claim.
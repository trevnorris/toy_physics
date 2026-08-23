# Independent physics review — the REPAIRED S11b SymPy engine (X-1 energy-basis fix)

## Artifact
`/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11b_interface_coupling_law_sympy_audit.py`
(a repair of committed baseline `864d6f41`; it writes `S11b_exports.py` and prints tagged objects to stdout).

## What to check
The repair was supposed to fix an energy-basis OVER-COUNT: the engine's independence test judged only
pointwise polynomial independence and omitted §5's symmetry-group bullet *"equivalence modulo total
divergences — two densities differing by a total in-plane divergence are the same term; do not count both"*
(`directives/S11b_SHARED_PHYSICS.md:286-287`). The repair had to (a) make the energy-basis independence
judgment honor that quotient, (b) eliminate any redundant invariant by REWRITING it — folding its coefficient
into the retained invariants — NOT by deleting the term, and (c) prove preservation via an emitted
`ENERGY_REEXPRESSION_RESIDUAL` whose Euler–Lagrange derivative is printed (not asserted zero) and must be able
to FAIL under a wrong fold. It must preserve the S11 carry-forward and the whole export chain. ⚠ Many
downstream symbolic COORDINATES legitimately shift when a coefficient folds — that is NOT a defect; the
invariant to preserve is the PHYSICS (EOM/energy modulo total divergence, dispersion roots, stability class,
transverse decoupling, the B2c/G13 slice, the breathing form).

⛔ You are NOT told the correct basis count, which invariant is redundant, or the fold — DERIVE them yourself
and compare to what the engine emits. A number I hand you would make your confirmation worthless.

## What you are handed
- The repaired engine (path above) and the spec `directives/S11b_SHARED_PHYSICS.md` (§5 ≈ L255-300; the
  symmetry group L282-291; note the L290 clause "judge independence as field bilinears, B1 not applied" — it
  is the B1-constraint rule, NOT an override of the divergence bullet).
- `CLAUDE.md`. You have the full repo; ⛔ there is no do-not-read list.

## Required method — SCRIPT review, and a FORM ablation is MANDATORY (rule 14), not optional
Derive independently, then ablate. ⛔ Code-reading alone has repeatedly missed real defects here.

1. **Derive the §5 basis yourself (CAS, save script + literal stdout to named absolute paths).** Enumerate
   the symmetry-allowed quadratic invariants of `(u, ∇u, θ, ∇θ, e_W, ∇e_W)` under the §5 group and compute
   the basis DIMENSION both WITH and WITHOUT the total-divergence quotient (witness: a constant-coefficient
   quadratic density is a total divergence iff its Euler–Lagrange derivative vanishes identically). Report
   BOTH numbers. Compare to the engine's emitted `ENERGY_BASIS_COUNT` / `ENERGY_BASIS` /
   `ENERGY_BASIS_INDEPENDENT_TERMS`. State whether the engine's emitted count matches the quotiented
   dimension you derive, and whether §5 mandates the quotient.

2. **FORM ablations (copy the engine to /tmp; ablate the COPY; ⛔ never touch the working tree). Report the
   LITERAL diff for each:**
   - **Quotient really implemented?** Ablate the divergence-quotient step (force the independence judgment
     back to pointwise). The emitted basis count MUST rise. If it does not move, the quotient is TYPED, not
     computed — a fatal defect. Also: which LINE computes the count? Give the line number, or report the
     count as an uncomputed literal.
   - **Divergence-equivalent invariant is not counted.** Perturb the FORM of a constructed invariant that is
     a total-divergence-equivalent of retained ones (e.g. flip a sign + an off-diagonal so it stays in the
     same divergence class). The emitted count MUST NOT rise. Then perturb a GENUINELY independent invariant;
     the count MUST rise. If the count is insensitive to both, it is not computed from the invariants.
   - **Fold, not delete (the load-bearing check).** Confirm `ENERGY_REEXPRESSION_RESIDUAL`'s EL derivative is
     identically zero as emitted, and that it becomes NONZERO under a one-sided corruption of a folded
     coefficient. Separately, construct the DELETE alternative (drop the redundant term's coefficient instead
     of folding it) and show — with your own CAS — that it changes the equations of motion (loses stiffness).
     A repair whose residual cannot fail, or that deleted instead of folded, is the defect.
   - **Two routes independent.** Corrupt ONLY the full-enumeration route of `ENERGY_REEXPRESSION_RESIDUAL`;
     the reduced re-expressed route must NOT move. If it moves, the two operands are not independent and the
     residual is decoration.

3. **Physics preserved — derive, don't read.** Independently derive the longitudinal dispersion roots /
   determinant vanishing locus, the transverse coupling (is it identically zero?), the B2c/G13 slice relation
   `Λ_p⁰ ↔ Λ_A⁰`, and the breathing-mode form, and compare to the engine's emitted objects as functions of
   the surviving free coefficients. A coordinate rename is fine; a changed root set / a nonzero transverse
   coupling / a changed slice relation is a defect (deleted or corrupted physics).

4. **Export chain (real-consumer run — ⛔ never read a cached `.out`).** `import S11b_exports` and check:
   the 1663 `S11_exports.LEDGER` keys are carried with identical values (except an `F9c` `s11b_`-prefixed
   collision); the NEW primary emissions (`ENERGY_REEXPRESSION_RESIDUAL`, its EL derivatives, the collapse
   evidence) actually appear as rows in `S11b_exports.LEDGER` (D1 — a computed-but-unexported row passes
   `import` and D3 yet is missing); every inherited object (`c_s0`, `mu_R`, `rho_br`) is bound to the
   imported `S11_exports.LEDGER` object under the same name and NOT redeclared / renamed onto a folded
   effective coefficient (§11); the bare key `v_0` never appears; digests pin `{own source, S11_exports.py,
   S11b_SHARED_PHYSICS.md}`; D3 round-trip, `_RELATIONALS`, freeze, and F6 gate are intact.

5. **No typed answer.** Grep the engine for a hardcoded count / coefficient relation / "expected". Report any
   `assert` that precedes the value it guards (an assert-before-emit hides a form defect — a perturbation
   strong enough to flip the check kills the process, so you see only PASS-or-crash). ⭐ For every emitted
   load-bearing object ask: WHICH LINE COMPUTED THIS? Give the line or report it uncomputed.

## Physics filter
Report a finding only if it catches a way the physics could be wrong, a value leak, a typed/undead check, or
a broken export — not "the script would be wrong on a different input."

## Constraints
- Copy to /tmp and ablate the COPY; ⛔ never modify the working tree.
- Save every derivation/ablation script AND its literal stdout to named absolute paths; report those paths.
  ⛔ A prose derivation is discarded.

## Output
Per item: verdict + the script path and literal stdout for every physics claim + the literal ablation diff.
End with a one-line bottom line: does the repair correctly implement §5's total-divergence quotient by
FOLDING (not deleting), with the export chain intact — or the specific defects found.

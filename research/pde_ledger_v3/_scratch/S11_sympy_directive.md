# DIRECTIVE — S11 SymPy audit + registry insertion

Two deliverables, both required:

**A.** `research/pde_ledger_v3/scripts/S11_stray_longitudinal_sympy_audit.py`
**B.** the registry insertion + gate updates described in part B below.

Iterate until every gate exits 0. Then stop and exit. ⛔ Do not commit anything.

---

## ⛔ WHAT YOU MAY NOT READ

A Mathematica script covering this same physics exists and has been **removed from the tree** for the
duration of your work. ⛔ Do not go looking for it in git history, in `git stash`, in any backup, or in
`research/pde_ledger_v2/`. The point of this exercise is that your derivation is **independent**; if you
transcribe another engine's answer, both agree vacuously and the check is worthless.

⛔ Also do not read:
- `research/pde_ledger_v3/V3_STEP_PLAN.md`
- `research/pde_ledger_v3/steps/` — any file
- `research/pde_ledger_v2/scripts/ledger_stage003_*`, `research/pde_ledger_v2/mathematica/ledger_stage003_*`

⭐ You **may and must** read `research/pde_ledger_v3/reduction/` — the SymPy side imports the registry.
That asymmetry is deliberate.

---

## THE SETUP — your physical input

A **brane** is a `D`-dimensional elastic sheet. `u` is an **in-plane** displacement `D`-vector field
`u(x,t)`, `x ∈ R^D`.

```
L  =  (1/2) rho_br (∂_t u)·(∂_t u)  −  (1/2) mu_R Curl2[u]  −  (1/2) B_comp (div u)^2

Curl2[u]  ≡  (1/2) Σ_{i,j} ( ∂_i u_j − ∂_j u_i )^2
div u     ≡  Σ_i ∂_i u_i
```

**Premises, all of them — do not supply your own:**
- `rho_br`, `mu_R`, `B_comp` are positive real constants.
- `u` is a **displacement**, so `[u] = L` (a length). ⭐ State where you use this.
- `L` is an **energy density on a `D`-dimensional spatial brane**, energy `= M L² T⁻²`.
- Wavevector `k ≠ 0`.
- Coefficients are **generic** (in particular `B_comp ≠ mu_R`) except in the one task that asks about
  coincident roots.

Separately, an unbounded medium off the brane supports a **scalar sound mode only**:
`ω² = c_s0² (k_in·k_in + k_w²)`, `k_in` along the brane, `k_w` normal to it.

---

## PART A — the audit script

Print one tagged line per result, prefix `S11_`. Compute everything; assert nothing you have not computed.

**A0 · Invariant census.** A quadratic form in the `D×D` matrix `∂_i u_j` invariant under simultaneous
rotation of both indices has some number `N` of independent coefficients. Determine `N` by explicit
construction. ⚠ **Report TWO counts and label them: `N_SO` under proper rotations `SO(D)` only, and
`N_O` under the full orthogonal group `O(D)` including reflections.** If they differ for some `D`, say
which `D` and exhibit the invariant responsible. Report both for `D = 2,3,4,5`.

**A1 · Dynamical matrix and spectrum by nullity.** For a plane wave `u = a exp(i(k·x − ωt))`, derive
the Euler–Lagrange equation, put it in the form `rho_br ω² a = M·a`, and **print `M`**. Then for
`D ∈ {2,3,4,5}` find the distinct roots `ω²` and, for **each** root, the **nullity** of the matrix at
that root plus the orientation of its kernel vectors relative to `k`.
⛔ Compute the kernel dimension. Do not read a count off an exponent in a factorised determinant.

**A2 · Cross-sector dependence test.** Determine **by computed residual** whether the
perpendicular-kernel root depends on `B_comp`, and whether the parallel-kernel root depends on `mu_R`.
Report each as a computed TRUE/FALSE with the residual shown. ⚠ Both answers are open; do not assume
either.

**A3 · Degeneracy locus.** Solve for the condition on `{mu_R, B_comp}` under which the two roots
coincide, and report the nullity structure **on** that locus.

**A4 · Dimensions, derived as closed functions of `D`.** From the premises above, derive `[rho_br]`,
`[mu_R]`, `[B_comp]` as `[L,T,M]` exponent vectors that are functions of `D`. Then, **for each root you
found in A1**, derive the dimension of that root's phase speed. Evaluate everything at `D = 3`.

**A5 · ⭐ REGISTRY CROSS-CHECK — the point of importing the registry.** Load the registry with
`registry_read` and compare your **derived** `D=3` dimensions against the **declared** ones for
`rho_br` and `mu_R`. Print `PASS`/`FAIL` per quantity with both vectors shown.
⛔ If they disagree, print the disagreement and **fail**. Do not adjust your derivation to match.

**A6 · Bulk matching — and state its SCOPE honestly.** For a brane mode of phase speed `c` coupling to
the medium above, sharing `ω` and `k_in`, solve for `k_w²`. Classify into **three** cases — `k_w` real,
`k_w` imaginary, and the threshold `k_w = 0` — and state what each means for whether the mode is
spatially bound to the brane.

⛔⛔ **Then state explicitly what this calculation does and does NOT establish.** It is **kinematic
phase matching only**: there is no brane–bulk coupling law, no interface boundary condition, and no
energy-flux calculation here. Report, as separate tagged lines, which of your conclusions are
**necessary** consequences (a channel does or does not exist) and which would additionally require a
coupling law to establish (whether an existing channel is actually used, and at what rate).
⚠ ⛔ Do not overstate the result to make it tidier. The scope limit is itself a required output.

Then attempt the same matching for the perpendicular-kernel brane mode, given that the medium supports
**only** the scalar sound mode, and report what happens.

**A7 · ⭐⭐ COMMENSURABILITY CHECK — this one exists because no gate can do it.** The bound/unbound
condition from A6 is an **inequality**, and the dimensional gate only classifies residuals equal to
zero, so nothing in the harness checks that its two sides are dimensionally comparable. Using dimensions
**imported from the registry**, verify that the two sides of that inequality have equal `[L,T,M]`
exponent vectors. Print both vectors and `PASS`/`FAIL`.

**A8 · FORM control.** Replace `(1/2) B_comp (div u)^2` with `mu_br Σ_{i,j} S̃_ij²`, `S̃` the symmetric
traceless part of `∂_i u_j`, with `mu_br > 0` and generic. Recompute A1 and A2 at `D = 3`; report the
roots, their nullities, and which roots moved and which did not.

**A9 · COEFFICIENT control.** With the original Lagrangian at `D = 3`, rescale `B_comp → 2 B_comp`.
Report the roots **and their nullities**, and which roots moved and which did not.

**A10 · Verdict.** `S11_VERDICT: PASS|FAIL`. ⛔ Define `PASS` as the conjunction of an **enumerated,
printed list** of the specific assertions you checked — print the list and each member's outcome. A
verdict not tied to enumerated assertions is worthless.
⚠ Understand its scope: it certifies your own checks did not contradict each other. It is ⛔ **not** a
verdict on the physics — the external comparison is.

---

## PART B — registry insertion and gates

⚠ **Part B runs AFTER part A and uses A's results.** ⛔ Do not let B's naming influence what A computes.

**B1 · Two new quantity rows** in `research/pde_ledger_v3/reduction/quantities.yaml`:

- `Q.brane.B_comp`, symbol `B_comp` — `kind: parameter`, `counting_axis: continuous-model`,
  `scope: [brane, light, III]`, regime matching the other brane rows.
- `Q.brane.c_L`, symbol `c_L` — `kind: intermediate`, `counting_axis: continuous-model`, same scope.
  This is the phase speed of the **parallel-kernel** root from A1.

⛔ **The declared `dimension.exponents` for each must be the values YOU DERIVED in A4**, not values you
found somewhere. The declaration is a derived result, not an assertion.

⚠ **Naming hazard, and a review leg will ask.** `B_comp` is a **brane** compression modulus. It is ⛔
not `K_br` (bulk modulus of a rejected ordinary-elastic branch), ⛔ not `B_eff` (`= ρ_B0²/χ_c`,
appearing elsewhere in the corpus), and ⛔ not the bulk medium's modulus. Do not alias it to any of
those, and do not cite their loci as its provenance.

**B2 · One new relation** `R5` in `relations.yaml`: the residual defining `c_L` from the A1 root,
prefix-v1, `designated_output: Q.brane.c_L`, with appropriate `denominator_guards` and `assumptions`.
Structurally a twin of the existing `R4`; match that row's shape and completeness.

**B3 · ⛔ The acceptance fixture — read this whole item before touching it.**
`reduction/acceptance_check.py` carries `EXPECTED_MEDIUM_PAYLOAD`, a **control fixture**. Its own comment
says: *"On any registry change, RECOMPUTE and independently re-derive; never copy forward."*

⇒ **Derive the new payload yourself first**, by reasoning about what adding one free parameter and one
eliminating relation does to the ambient count and to the constraint dimension in each of the three
cases. **Write your derivation down in your final report.** Only then run the check.

⛔⛔ **If the code's computed payload disagrees with your independent derivation, REPORT THE
DISAGREEMENT and stop. Do not set the fixture to whatever the code printed** — a fixture copied from the
code it polices is the code agreeing with itself, and it would silently destroy the only control here.

**B4 · Run every gate and report the output of each:**
```
python3 research/pde_ledger_v3/reduction/acceptance_check.py
python3 research/pde_ledger_v3/reduction/dimensional_homogeneity_gate.py
python3 research/pde_ledger_v3/reduction/able_to_fail.py
python3 -m pytest research/pde_ledger_v3/reduction/test_registry.py -q
python3 research/pde_ledger_v3/scripts/S11_stray_longitudinal_sympy_audit.py
```

**B5 · Confirm the dim gate actually COVERS the new relation.** Report `R5`'s classification line and
the `POPULATION_COUNTS`. ⛔ `R5` landing in `UNDETERMINED` is a failure, not a pass.

**B6 · Confirm the discrete payload is untouched.** `n_eos = 5` and `D_brane = 3` unchanged, and the
discrete rows must not have entered the continuous count.

---

## Report

Under 40 lines: the A0 two counts, what you derived for A4, the A5/A7 verdicts, A6's scope statement,
the B3 derivation and whether it matched, each gate's headline line, and anything that surprised you.
⛔ Do not summarise the directive back to me.

# Question-vet ONLY — the N6 build correction's TRACTABLE implementation. ⛔ Do NOT build/run; load no `.out`.

You are Codex, consulted to vet the **approach** for the N6 build correction before a directive is written. ⛔ Do
NOT write or run code; read source/spec and reason. Working dir = `/var/projects/toy_physics`; paths are under
`research/pde_ledger_v3/`.

## What is already settled (⛔ not up for re-litigation)
`S11c_c2 §5c` (committed `30d4b72d`, review-until-clear) **defines** the corrected N6 rep-invariance test: within
each FIXED anchoring `α` and density `ρ`, construct the self-energy increment by two routes and compare —
`I_E = extract(close(SLAB)) − extract(SLAB)` (Eulerian) vs `I_{M→E} = T_{M→E}[extract(close(SLAB_M)) − extract(SLAB_M)]`
(material-coordinate, mapped back by `Δρ=δρ_E+u·∇ρ⁰`), residual `R_N6 = I_E − I_{M→E}` must vanish; plus a
one-sided independence corruption. Read `directives/S11c_c2_SHARED_PHYSICS.md` §5c (≈303–369), the audit machinery
`scripts/S11c_c2_selfenergy_fold_sympy_audit.py` (`bind_inputs`, `build_case`, `extract`, `retained_shape`), and
`S11c_b_SHARED_PHYSICS.md` §5a (the material-slab route 2 template).

## The DESIGN problem (what to vet)
A sibling attempt (the F/G diagnostic) hit a **tractability wall**: a full-**symbolic** test over all retained
grades on the c2 increment operands is impractical (~11 min / ~150 MB per case;
`_measurements/S11c_c2_FG_regrounding_deferred.md`). The N6 build risks the same, since it constructs a SECOND full
slab route (`SLAB_M`) and differences two huge symbolic objects.

Proposed tractable design (VET IT):
1. **Numeric residual test (Schwartz-Zippel).** Instead of proving `R_N6 = 0` symbolically, evaluate `R_N6` at many
   independent **random numeric** assignments of the free fields/jets/parameters (per retained grade) and test it is
   numerically zero to working precision; a nonzero at any grade is the finding.
2. Construct `SLAB_M` per S11c-b §5a (material flattening + map-back) reusing `build_case`, ⛔ not a re-derivation.

## What I need from you (reason; cite source)
1. **Is the numeric Schwartz-Zippel residual test VALID for N6**, or a proxy? Does numeric-zero at randomized points
   (enough points / high enough precision, respecting the `(ε,η,σ_W)` grading and the weak/trial-test structure)
   soundly establish rep-invariance vanishing, or can it miss a structural non-vanishing? What must be randomized
   vs held (bound integration dummies, the `Z`/Fourier kernels, the trial/test potentials, the background jets) so
   the test is decisive and not accidentally zero? Name the proxy traps.
2. **Is constructing `SLAB_M` (the material route) itself tractable**, or does it hit the same wall regardless of
   the residual test being numeric? If the material *construction* is the expensive part, is there a lighter route
   (e.g. construct only the increment's slot-carrier part `S_{P,s}`, or a per-grade numeric build) that still gives
   a sound N6 residual?
3. **Where should the build live** — extend `S11c_c2_selfenergy_fold_sympy_audit.py`, or a companion script reusing
   its machinery? Does the corrected §5c residual need any object the audit does not already produce?
4. **The corrected question, if mine is a proxy**, and the decisive tractable diagnostic (numeric or otherwise) that
   soundly tests the §5c N6 residual + the one-sided independence corruption, with the proxy-success traps to forbid.

## Output
A verdict: is the numeric residual approach sound for N6 (yes / no + the corrected approach); is the material route
tractable (yes / no + the lighter route if not); where the build should live; and the decisive tractable diagnostic.
⛔ Reasoning + citations only — no code.

# DIRECTIVE — S11 blind Mathematica audit

Write **one** file: `research/pde_ledger_v3/mathematica/S11_stray_longitudinal_mathematica_audit.wl`

Iterate until it runs clean under `math -script` with exit 0. Then stop and exit.

---

## ⛔⛔ BLIND RULES — the whole point of this script

This script is an **independent** check. Its value is that it can **DISAGREE** with a SymPy script that
does not yet exist. If you seed it from anything that already contains an answer, it agrees vacuously
and the check is worthless.

⛔ **Do NOT read any of these, at any point:**

- `research/pde_ledger_v3/reduction/` — any file. The script must not import the registry.
- any `.py` anywhere under `research/pde_ledger_v3/`
- `research/pde_ledger_v3/V3_STEP_PLAN.md`
- `research/pde_ledger_v3/steps/` — any file
- `research/pde_ledger_v3/CHARTER.md`, `LAUNCH_PROMPT.md`, `SESSION_2026-08-01.md`, `DEFECT_REGISTER.md`
- `research/pde_ledger_v2/scripts/ledger_stage003_*` and
  `research/pde_ledger_v2/mathematica/ledger_stage003_*` — these compute a related problem
- `docs/model_map.md`, `docs/conceptual_foundation.md`,
  `software/stage1_solver/decisions/15_em_medium_native_physical_picture.md`

⭐ Everything you need is in **THE SETUP** below. If you think you need more, you do not — derive it.

**Script conventions** (matching the 44 `.wl` audits already in this repo): standalone, **no
arguments**, imports **no files**, **print only**, no exports. `$ScriptCommandLine` is empty under
`math -script` here — do not use it.

---

## THE SETUP — this is your only input

A **brane** is a `D`-dimensional elastic sheet. `u` is an **in-plane** displacement, a `D`-vector field
`u(x,t)` with `x ∈ R^D`. Its Lagrangian density is

```
L  =  (1/2) rhoBr (D[u,t] . D[u,t])  −  (1/2) muR Curl2[u]  −  (1/2) Bcomp Div[u]^2
```

with

```
Curl2[u]  ≡  (1/2) Sum_{i,j} ( D[u_j, x_i] − D[u_i, x_j] )^2
Div[u]    ≡  Sum_i D[u_i, x_i]
```

`rhoBr`, `muR`, `Bcomp` are positive real constants. Work in general `D` where asked, and specialise
where asked.

---

## TASKS — compute each, print the tagged line

**T1 · Invariant census.** A quadratic form in the `D×D` matrix `D[u_j, x_i]` that is invariant under
simultaneous rotation of both indices has some number `N` of independent coefficients. Determine `N` by
explicit construction (decompose the matrix into pieces that do not mix under rotation and count them).
Report `N` for `D = 2,3,4,5`.
→ `WL_S11_INVARIANT_COUNT: {D -> N, ...}`

**T2 · Dynamical matrix.** For a plane wave `u = a Exp[I(k.x − w t)]`, derive the Euler–Lagrange
equation and put it in the form `rhoBr w^2 a == M . a`. Print `M` with its `k`-dependence explicit.
→ `WL_S11_DYNAMICAL_MATRIX: <expression>`

**T3 · Spectrum by NULLITY.** For each `D ∈ {2,3,4,5}`: find the **distinct** roots `w^2` of
`Det[M − rhoBr w^2 IdentityMatrix[D]] == 0`. For **each** root, compute the **NULLITY** of
`M − rhoBr w^2 I` — the dimension of its kernel — and the orientation of the kernel vectors relative to
`k` (parallel / perpendicular).
⛔ **Do not read a count off an exponent in a factorised determinant. Compute the kernel dimension.**
→ `WL_S11_SPECTRUM_D<n>: {root -> {nullity, orientation}, ...}`
→ `WL_S11_NULLITY_SUM_D<n>: <integer>`

**T4 · Cross-sector independence.** Determine whether the root whose kernel vectors are
**perpendicular** to `k` depends on `Bcomp`, and whether the root whose kernel vectors are **parallel**
to `k` depends on `muR`. Decide each by a computed residual, not by inspection.
→ `WL_S11_PERP_ROOT_DEPENDS_ON_BCOMP: TRUE|FALSE`
→ `WL_S11_PARA_ROOT_DEPENDS_ON_MUR: TRUE|FALSE`

**T5 · Degeneracy locus.** Solve for the condition on `{muR, Bcomp}` under which the two roots coincide.
→ `WL_S11_DEGENERACY_CONDITION: <expression>`

**T6 · Dimensions — DERIVED, not assumed.** Using `[L,T,M]` exponent vectors and the single premise
that `L` above is an **energy density on a `D`-dimensional spatial brane** (energy `= M L^2 T^-2`),
derive each coefficient's dimension as a **closed function of `D`**, then evaluate at `D = 3`.
→ `WL_S11_DIM_RHOBR`, `WL_S11_DIM_MUR`, `WL_S11_DIM_BCOMP`, `WL_S11_DIM_SPEED_LONGITUDINAL`
   (each printed both as a function of `D` and at `D = 3`)

**T7 · Bulk matching.** A separate unbounded medium fills the space off the brane and supports a
**scalar sound** mode only, with dispersion `w^2 == cs0^2 (kIn.kIn + kw^2)`, where `kIn` is the
wavevector component along the brane directions and `kw` the component normal to the brane. A brane
mode with in-plane wavevector `k` and frequency `w` couples to this medium sharing **both** `w` and
`kIn = k`.

Solve for `kw^2` in terms of the brane mode's phase speed and `cs0`, and determine the condition under
which `kw` is real versus imaginary. State what each case means for whether the brane mode is spatially
**bound** to the brane or carries energy away from it.
→ `WL_S11_KW_SQUARED: <expression>`
→ `WL_S11_BOUND_CONDITION: <inequality in Bcomp, rhoBr, cs0>`

Then attempt the **same** matching for the brane mode whose kernel vectors are perpendicular to `k`,
given that the medium above supports **only** the scalar sound mode. Report what happens.
→ `WL_S11_TRANSVERSE_MATCHING: <what happened>`

**T8 · FORM control.** Replace the `(1/2) Bcomp Div[u]^2` term with `muBr Sum_{i,j} Stilde_{ij}^2`,
where `Stilde` is the **symmetric traceless** part of `D[u_j, x_i]`. Recompute T3 and T4 at `D = 3` and
report which roots moved and which did not.
→ `WL_S11_FORM_CONTROL_D3: <roots, nullities, what moved>`

**T9 · COEFFICIENT control.** With the **original** Lagrangian, rescale `Bcomp -> 2 Bcomp` at `D = 3`.
Recompute T3 and report which roots moved and which did not.
→ `WL_S11_COEFFICIENT_CONTROL_D3: <roots, nullities, what moved>`

**T10 · Verdict.** `WL_S11_VERDICT: PASS` if every internal consistency check you built agrees, `FAIL`
with a reason otherwise.
⚠ Be clear what this means: it certifies **your own checks did not contradict each other**. It is ⛔
**not** a verdict on the physics — the external comparison is.

---

## ⛔ Output discipline

⛔ Do **not** print any expected value, target, reference result, or comparison to a known answer.
Print **what you compute**. If a result surprises you, print it anyway — a surprising computed value is
the most valuable thing this script can produce.

## Environment

- Mathematica: **2-seat licence**. Never run more than 2 `math -script` processes at once.
- Do not wrap anything in a shell `timeout`. Keep every run under 10 minutes; if something is slower
  than that, reformulate it rather than waiting.

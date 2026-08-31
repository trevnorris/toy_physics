# S11c-b — §1d variable-coefficient energy-basis REFERENCE (build directive, #86, rev 2)

> ⚠ **OUTCOME (see `_measurements/S11c_b_86_reference_result.md`).** This directive built a reference, but its
> "object to compute" left the uniform coefficient carrier ambiguous, and the builder read it as **per-source
> (2 carriers)** → uniform 22, combined 50, frozen 36 ≠ engine 26. The Claude build-leg caught it (the Grok
> leg missed it). The §1d fix itself (spurion Hessian, per-source 8→15, **+14**) is carrier-independent and
> correct. The settled corrected basis is **40** (1-carrier, matching the engine's frozen 26), **established by
> four consistent computations** (engine anchor 26, Codex decision-leg, a Claude build-leg's from-scratch
> 1-carrier reconstruction, and the orchestrator crux CAS check) — not merely a user call. The defective
> 2-carrier deliverable named in "Task and authority" below was **removed** (archived at
> `~/.s11_build/S11c_b_86_probes/DEFECTIVE_2carrier_reference.{py,out}`). Kept as the process record; ⛔ do not
> rebuild from this directive without pinning the uniform to the engine's 1-carrier convention (frozen must
> reproduce the engine's emitted basis, count AND selections).

## Task and authority
Build a **standalone reference** that computes the correct §1d variable-coefficient invariant quotient for the
S11c-b stored-energy density — the **dimension** and a **spanning representative set** (with divergence
certificates) of the DOF-bilinear energy invariants admissible on the non-uniform background, judged
independent **modulo the variable-coefficient total in-plane divergence**. Deliverable:
`research/pde_ledger_v3/scripts/S11c_b_energy_basis_reference.py`, stdout to
`research/pde_ledger_v3/scripts/out/S11c_b_energy_basis_reference.out`. Those two are the only writes; every
other file is read-only.

`CLAUDE.md` binds. `research/pde_ledger_v3/directives/S11c_b_SHARED_PHYSICS.md` is the sole physics authority
(§1a, §1c, §1d, §3a govern) and wins every conflict. This is a **reference / third route**: **import nothing**
from `scripts/S11c_b_brane_operator_sympy_audit.py`; re-derive the enumeration and the quotient from the spec
(the `S9` blindness control). It is a computation, not a repair — touch no engine file.

⛔ **Withhold the answer (rule 5).** The corrected quotient dimension, the per-family sub-dimensions, the
per-source selected/omitted sets, and the number of newly independent invariants are the computed result.
Write **no** count, expected value, "should be N", or acceptance criterion into the script or its tags. Every
count, rank, and selection must be **reached by computation from the constructed candidates**, never typed.
The script computes and prints; the diff happens on our side. (A build leg will form-ablate to confirm no
count is hard-coded.)

## The object to compute (a finite quotient; named, not a recipe)
Construct one **finite** candidate seed of parity-even `O(3)` scalar bilinears in the DOF data, then quotient
by the variable-coefficient total in-plane divergence. The seed has two families:

- **Uniform family:** zero-spurion bilinears of the DOF first-jet data
  `{u_i, ∂_j u_i, θ, ∂_iθ, e_W, ∂_i e_W}`, each carried with a **single** background-scalar coefficient
  (`W_bg` or `μ_R,bg`) at fixed background-amplitude order — **not** an arbitrary function library of the
  background (do **not** admit `{W_bg^n · …}` as independent classes; the coefficient is one background scalar
  at one amplitude order, per §3a "its own free symbolic coefficient"). On the uniform background `u` enters
  only through its gradients (§1a translation-invariance), so undifferentiated `u` does **not** appear here.
- **Spurion family:** bilinears **linear in exactly one** background first jet `∂W_bg` or `∂μ_R,bg`,
  contracted with the DOF data. Here §1a permits `u` (either part) to enter **undifferentiated** contracted
  with the background gradient. Whether a given undifferentiated-`u` contraction survives the quotient is
  computed, **not** assumed.

Impose the exact thickness map `e_W,bg=(W_0/W_bg) e_W` (S11c-a §2a) **before** enumeration/rank so the
thickness sector is not double-counted. Independence is judged as field bilinears (B1 not applied); the
unbroken `w→−w` reflection is enforced; no positivity/T-reversal.

The **equivalence relation** is §1d, verbatim: "the uniform quotient does **not** lift trivially … integrating
by parts a variable coefficient generates first-background-jet terms (`c∇·F ≡ −(∇c)·F` modulo a boundary
term) … representatives equivalent uniformly differ by a first-jet invariant that is **physics in the
operator/kernel**." §3a, verbatim: "a second spatial derivative of `W_bg` is still first order in background
amplitude … and is **not dropped**." ⇒ Quotient in the jet algebra where the background scalars carry their
**full jet tower** — `D_i W_bg = g_i`, `D_i g_j = H_ij = H_ji`, and the analogous chain for `μ_R,bg` — so the
product rule `∂_i(c F_i) = (∂_i c) F_i + c ∂_i F_i` holds with **both** terms present. Backgrounds are
differentiated but are **not** dynamical DOFs (the quotient's total-divergence currents may depend on them).

**Grade projection (§3a).** The background bookkeeper order bounds **powers** of `η`/`σ_W`, not spatial
derivative order. Retain higher spatial background jets (`∂²W_bg`, …) at the amplitude order of their
originating factor, but **project out** any candidate or generated term whose background-amplitude grade
exceeds `(η,σ_W) ≤ 1` — e.g. the exact `e_W,bg` map expands to `σ_W²` terms that must be removed. Emit the
result both before and after this projection.

## Diagnosed defect this reference measures against (context — evidence, not recipe)
In the frozen engine the §3a independence test (`basis_euler_signatures` / `quotient_independent_indices`,
`scripts/S11c_b_brane_operator_sympy_audit.py`) builds its Euler–Lagrange total-derivative map **only from the
DOF fields**; the background jet is absent, so `∂W_bg` is treated as constant and the `∂²W_bg` term §3a
retains is never generated. The reference must judge the object above correctly; **how** (extend the
derivative map, an explicit divergence-current basis, or any other construction) is yours — name the object,
not the recipe. ⛔ Do **not** read
`directives/_measurements/S11c_b_coupling_84_basis_verification.md` or any file that estimates the corrected
count; the frozen selection you need for the regression is the engine's own emitted `.out` (public, and only
the *frozen* value).

## Controls and diagnostics (compute and PRINT each; ⛔ assert no expected outcome)
Every item below emits a computed object. ⛔ Do **not** write "must move", "must change", "must appear", or any
target — that is a rule-5 leak. We read the deltas on our side.
1. **Frozen-quotient regression.** Provide a switch that freezes the background jet (drops `∂²c`, treats the
   spurion as constant). Emit its per-source selection. It is the regression anchor against the engine's
   *emitted frozen* selection (the only public value).
2. **Background-Hessian ablations.** Emit the quotient dimension under: (a) `∂²W_bg→0` only; (b) `∂²μ_R,bg→0`
   only; (c) both `→0`; (d) full retention. Also emit a `W_bg↔μ_R,bg` source-exchange consistency object
   (the two single-source sub-dimensions).
3. **Map-ordering diagnostic.** Emit the dimension under (a) the `e_W,bg` map imposed **before** rank, and
   (b) the map imposed **after** rank / two independent thickness coordinates.
4. **Undifferentiated-`u` diagnostic.** Emit whether the spanning set contains undifferentiated-`u` spurion
   representatives, and an example representative if so — as a computed object, not a requirement.
5. **Certificates (the real independence evidence).** For **every omitted** candidate, emit an explicit
   total-divergence / linear-combination certificate proving its dependence. For the selected span, emit an
   independent rank certificate. Emit that the enumeration is `O(3)`-even and that the divergence-free
   bilinear currents (the would-be null-Lagrangians) are parity-**odd** and therefore excluded.
6. **Combined-vs-separated diagnostic.** Emit the combined-quotient rank and the per-family ranks
   (uniform, `W_bg`-spurion, `μ_R,bg`-spurion) so their relationship is a **measured** object, not asserted.

## Three script clauses — non-negotiable (build skill)
1. PRINT computed objects (expressions, ranks, certificates, booleans). ⛔ Never a typed conclusion.
2. PRINT the dimension/residual; ⛔ do **not** `assert dim == N` — that is the withheld answer written down.
3. Interpretation belongs to the step record, not the script.

⭐ The only place physical symbols may be combined by hand is in **constructing the candidate bilinears and
the divergence currents**. Every rank, selection, certificate, and control dimension must be **reached by
computation**; every control re-enters at the construction, never at a result.

## Emit (RESULTS tags — objects, not verdicts)
```text
S11CB_REF_BASIS_DIMENSION              corrected combined-quotient dimension (post grade projection)
S11CB_REF_BASIS_DIMENSION_PREPROJECT   same, before grade projection
S11CB_REF_BASIS_SPAN                   spanning representative invariants (expressions; note: non-unique)
S11CB_REF_FAMILY_RANKS                 per-family ranks {uniform, W_spurion, muR_spurion} and their sum
S11CB_REF_OMITTED_CERTIFICATES         per omitted candidate: its divergence/linear-combination certificate
S11CB_REF_FROZEN_SELECTED_PER_SOURCE   control-1 frozen selection (regression vs the engine's emitted set)
S11CB_REF_HESSIAN_ABLATION_DIMS        control-2 dims {W_off, muR_off, both_off, full} + source-exchange
S11CB_REF_MAP_ORDER_DIMS               control-3 dims {map_before_rank, map_after_rank}
S11CB_REF_UNDIFF_U                     control-4 boolean + example representative
S11CB_REF_PARITY_EXCLUSION             control-5 evenness of the seed; parity of the divergence-free currents
```

## Report (≤15 lines)
Report only: the deliverable + `.out` paths; that it imports nothing from the engine and does not read the
count-estimating record; whether control-1 reproduced the engine's emitted frozen selection; and any
construction ambiguity resolved against the spec. ⛔ Do not state whether the corrected dimension is "right" —
that adjudication is ours.

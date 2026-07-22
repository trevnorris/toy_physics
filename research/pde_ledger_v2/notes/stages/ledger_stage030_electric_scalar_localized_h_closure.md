# ledger_stage030_electric_scalar_localized_h_closure

## Status

**Part IV — Charge. IV-1 (build-order 030; the FIRST Part-IV stage).** LIGHT-with-scoping reshape of
the electric-scalar subset of the monolithic G0 closure checker
`software/em_charge_attribute/g0_closure_card_v0_checks.{py,wl}` (source card
`g0_closure_card_v0.md` §§2.3–2.4). Surviving headline verdict token:

- **`ELECTRIC_SCALAR_CLOSURE_STATIC`** — the static algebraic closure of the charge electric-scalar
  sector: the localized parent scalar `H`, its unique gapless zero mode `f₀`, the normalized reduced
  mode `h = P₀H`, the reduction relations `{M_h = N₀M₄, K_h = N₀K₄ = M_h c_E²}`, and the coupled
  `(u_L, h)` kernel with independent Sylvester-positivity and positive generalized wave speeds.

**EARNED** (what this stage genuinely proves, *given* the postulated action): the parent operator
`O_⊥ = A†A ≥ 0`; the unique normalizable **zero** mode `f₀` and the positive continuum gap `4/ℓ²`;
the zero-mode norm/projection and its idempotence; the three reduction identities; static kernel
positivity (`D > 0`) and two positive squared wave speeds; reduced-`h` masslessness (`Q_s(0,0)=0`);
conservative-Hessian symmetry; and the units-restored `[L,T,M]` dimensional homogeneity of the
electric-scalar-block terms. **POSTULATED** (the honest landing, first-class): the entire G0-v0
closure — the choice of Pöschl–Teller operator, the parity representation, and every numeric witness
value — is a labeled **`DRAFT v0` postulated closure**, not a derivation.

Ledger-local earned-label (NOT a source verdict token): `ELECTRIC_SCALAR_STATIC_PREREQUISITES_VERIFIED`.

This stage lays the **substrate** the charge sector is built on; it does **not** contain a force, a
sign, or a mouth source. Charge = the static ±w throat: a particle punctures the brane into `±w`, the
bend is the geometric field `ξ_w = ℓ h`, and `h` is the localized zero mode assembled here. The
puncture-deflection **mechanism** (the mouth `χ↔h` source, the exterior `1/R²` far-field FORM) is
stage 031, and the four BC ensembles + the honest `R1_REQUIRED(bc_selection)` landing are stage 032
(4-stage split, 2026-07-22) — both *dress* this substrate.

## Purpose

The charge sector needs a distinct electric scalar that (i) is genuinely **localized** on the brane
(so charge is a defect, not a delocalized bulk field), (ii) reduces to a single **gapless** propagating
mode `h` (so the far field is massless `1/R`, not Yukawa/pinned), and (iii) couples to the longitudinal
brane displacement `u_L` through a **stable** two-field kernel. G0-v0 postulates such a sector; stage
030 is its formal ledger home and the place the static algebraic prerequisites are paid for in public,
dual-engine, before any force language. Every predicate here is a *static prerequisite* for the
assembled Green mode, **not** a solved assertion of it (the assembled one-body / two-body / far-field
claims are Part VII; the mouth source is stage 031).

## 1. The localized parent `H` and its Pöschl–Teller operator (EARNED spectral structure)

The distinct electric scalar is a parent field `H(x, w, t)` of dimension `L⁻¹` living in the full
`y = (x, w) ∈ R³ × R`. Its transverse operator is the reflection-even Pöschl–Teller operator

```text
O_⊥ = −∂_w² + [4 − 6 sech²(w/ℓ)] / ℓ²,
S_H = ½ ∫ dt d⁴y { 2 M₄ (∂_t H)² − 2 K₄ [ |∇_x H|² + H O_⊥ H ] }.
```

**Factorization / nonnegativity (`PASS_TRANSVERSE_FACTORIZATION`).** With the superpotential
`A = ∂_w + 2 tanh(w/ℓ)/ℓ`, the engines verify term-by-term on a symbolic probe that

```text
O_⊥ = A†A,     A† = −∂_w + 2 tanh(w/ℓ)/ℓ,
```

so `O_⊥` is formally nonnegative (`⟨Aψ, Aψ⟩ ≥ 0`).

**The zero mode (`PASS_PARENT_TRANSVERSE_ZERO_EIGENVALUE`).** Direct substitution gives

```text
O_⊥ f₀ = 0,     f₀(w) = 1 / [ ℓ cosh²(w/ℓ) ].
```

**Uniqueness of the normalizable kernel (`PASS_UNIQUE_NORMALIZABLE_KERNEL`).** Because `O_⊥ = A†A`,
the zero mode is exactly `ker A`. `A f = 0` is a **first-order** ODE, so its solution space is
one-dimensional: the integrating-factor solution is `f ∝ cosh⁻²(w/ℓ)`, and the engines confirm its
ratio to `f₀` is `w`-independent (`= ℓ`) with finite positive `L²` norm. Hence `ker A = span{f₀}` —
`f₀` is the **unique normalizable zero mode**.

**Continuum gap (`PASS_POSITIVE_CONTINUUM_GAP`).** As `w → ±∞`, `V_H = [4 − 6 sech²(w/ℓ)]/ℓ² → 4/ℓ²`
(verified by symbolic `limit` in both directions), so the continuum threshold is `4/ℓ² > 0` and the
zero mode sits a **positive gap** below it.

> **Scope of the spectral claim (no over-claim).** The engines claim exactly: a **unique normalizable
> zero mode** `f₀ = ker A`, below a continuum starting at `4/ℓ²`. They do **not** enumerate the full
> bound-state spectrum, and `f₀` is **not** the only bound state: as an *analytic PT-theory aside* (NOT
> part of this stage's dual-engine verification), this well has `s(s+1) = 6` (`s = 2`), so it additionally
> supports one **excited** bound state at `E = 3/ℓ²` (still below the `4/ℓ²` continuum). That excited mode
> is real but is neither computed nor relied on here; the reduction below keeps only the `E = 0` zero mode.

## 2. Zero-mode norm, projection, and the reduction relations (EARNED)

**Norm (`PASS_ZERO_MODE_NORM`).** With the parent inner weight `2` fixed by `S_H`,

```text
N₀ = ∫ dw 2 f₀² = 8 / (3ℓ)  > 0     (Integrate / integrate, closed form).
```

**Normalized projection (`PASS_NORMALIZED_ZERO_MODE_PROJECTION`).** Writing `H = f₀ h + H_⊥` with
`∫ dw 2 f₀ H_⊥ = 0`, the reduced mode is the normalized zero-mode projection

```text
h = P₀H := N₀⁻¹ ∫ dw 2 f₀ H.
```

The engines verify `P₀ f₀ = 1`, `P₀(f₀ h) = h`, and idempotence `P₀² = P₀` on the zero mode. `[h] = 1`
(dimensionless): `[N₀⁻¹] · [f₀] · [H] · [dw] = L · L⁻¹ · L⁻¹ · L = 1`.

**Reduction relations (separate teeth, each genuinely able-to-fail).** The zero-mode reduction of `S_H`
gives, each as its own dual-engine predicate (not folded into `M_h` alone), with dimensions checked:

| Relation | Predicate | Dependency |
|---|---|---|
| `M_h = N₀ M₄` (witness lands `M_h*=1`) | `PASS_REDUCED_H_INERTIA` | `M₄` is an **input**; `M_h` DERIVED (`M L⁻¹`) |
| `K_h = N₀ K₄` (witness lands `K_h*=1`) | `PASS_REDUCED_H_STIFFNESS` | `K_h` DERIVED (`M L T⁻²`) |
| `[K₄] = [M₄ c_E²]` (dimensional) | `PASS_PARENT_STIFFNESS_DIMENSIONAL_CONSISTENCY` | `{M₄, c_E}` **inputs** → `K₄` DERIVED (`E`) |
| speed preservation `c_E² = K₄/M₄ = K_h/M_h` | `PASS_REDUCED_SPEED_PRESERVATION` | tests the reduction is speed-preserving |

The physical positivity `M₄ N₀ > 0` (with `[M₄ N₀] = [M_h]`) is its own predicate (`PASS_PHYSICAL_H_NORM`).

> **On the two `c_E`-relation teeth (genuinely able-to-fail, not tautological).** The *value* relations
> `K₄ = M₄ c_E²` and `K_h = M_h c_E²` are **definitional** ({M₄, c_E} are the inputs, so `K₄` is *defined*
> as `M₄ c_E²`); a value self-comparison would be a vacuous `X≡X`. So the two teeth verify their genuine,
> production-live content instead: `PASS_PARENT_STIFFNESS_DIMENSIONAL_CONSISTENCY` checks `[M₄ c_E²]=[K₄]`
> (fails under a `[K₄]` dimension mutation), and `PASS_REDUCED_SPEED_PRESERVATION` checks that the reduced
> speed equals the parent speed, `c_E² = K_h/M_h = K₄/M₄`, using the **separately-reduced** `K_h=N₀K₄`,
> `M_h=N₀M₄` — so an inconsistent reduction (e.g. `K_h=N₀²K₄`) makes `K_h/M_h=160/3 ≠ 1=K₄/M₄` and the
> tooth fires on a *real* error, not a decoupled copy. The genuine derived witness landings `M_h*=K_h*=1`
> are carried by `PASS_REDUCED_H_INERTIA` / `PASS_REDUCED_H_STIFFNESS` (the Integrate-derived `N₀*=160/3`
> times the calibrated `M₄*`/`K₄*`).

> **Dependency direction (load-bearing).** `M₄` (a postulated `ACTION` primitive) and `c_E` (input)
> are the free data; `K₄ = M₄ c_E²`, `M_h = N₀ M₄`, and `K_h = N₀ K₄ = M_h c_E²` are all **derived**
> from them. The reduced normalization `M_h* = 1` is a **`CONV`** choice imposed **on `M₄`**, not a
> reduction of `M_h`.

## 3. The coupled `(u_L, h)` scalar block (EARNED — two INDEPENDENT positivity facts)

After eliminating the brane phase `θ_B` (its Schur shift absorbed into the inertia
`A_eff = ρ_br + C_J²/κ_phase`), the reduced scalar action is

```text
S_Lh = ∫ dt d³x [ ½ A_eff (∂_t u_L)² + ½ M_h (∂_t h)²
                  − ½ B_eff |∇u_L|² − ½ K_h |∇h|² − C_hu ∇u_L · ∇h ],
```

with inertia matrix `M = diag(A_eff, M_h)` and stiffness matrix `K = [[B_eff, C_hu], [C_hu, K_h]]`.
The stage splits kernel health into **two teeth that are genuinely independent** (this split is used
in BOTH engines):

**(a) Sylvester positivity (`PASS_STABILITY`).** Static stability is positive-definiteness of the
stiffness matrix `K`, i.e. `B_eff > 0` and

```text
D = B_eff K_h − C_hu²  > 0     (physical [D] = M² T⁻⁴).
```

The engines confirm `K` is positive-definite (SymPy `is_positive_definite`; Wolfram
`PositiveDefiniteMatrixQ` + positive `Eigenvalues`) and that `[B_eff K_h] = [C_hu²] = [D]`.

**(b) Positive generalized wave speeds (`PASS_POSITIVE_GENERALIZED_WAVE_SPEEDS`).** The two coupled
squared characteristic speeds `z = c²/c_E²` solve the generalized eigenproblem `det(K − z M) = 0`. In
the explicitly-normalized star units (`A_eff* = M_h* = 1`, so `M* = I`; `B_eff* = 2`, `K_h* = 1`,
`C_hu* = 1/2`):

```text
det(K* − z* M*) = z*² − 3 z* + 7/4,     z*_± = c_±*² = (3 ± √2)/2  >  0.
```

Both roots are positive, so both coupled modes propagate. (The *physical* `det(K − z M)` is
dimensionful and is **never** equated to this numeric polynomial — the determinant is stated only
through the explicitly-normalized star form.)

> **Why (a) and (b) are independent facts.** `D` depends on the **stiffness** matrix only; the wave
> speeds depend on **both** `K` and the inertia `M`. The able-to-fail tooth for (b) mutates the inertia
> (`A_eff*: 1 → 2`): this **keeps** the Sylvester margin `D* = 7/4` intact yet **changes**
> `det(K − zM)` and its roots — driving (b) to exit 1 at its *own* assert while (a) stays green. So the
> roots check is not a corollary of the margin check.

## 4. Reduced-`h` masslessness and conservative Hessian symmetry (EARNED)

**Masslessness (`PASS_REDUCED_H_MASSLESSNESS`).** The time-dependent quadratic operator is

```text
Q_s(ω, k) = [[ A_eff ω² − B_eff k²,   − C_hu k² ],
             [ − C_hu k²,             M_h ω² − K_h k² ]].
```

At star values `Q_s(0,0) = 0`: there is **no `k⁰ h²` term** (no bare `h` mass), the static prerequisite
against `FAIL_PINNED_BRANON` / `FAIL_YUKAWA`. The tooth reinstates a `k⁰ h²` coefficient (`0 → 1`) and
`Q_s(0,0) ≠ 0` fires.

**Hessian symmetry (`PASS_CONSERVATIVE_HESSIAN_SYMMETRY`).** The gradient-energy Hessian of the
`(u_L, h)` block equals the displayed Euler operator and is symmetric:

```text
∂²/∂(∇u_L, ∇h)² [ |∇u_L|² + ½|∇h|² + ½ ∇u_L·∇h ]  =  [[2, 1/2], [1/2, 1]]  =  its transpose.
```

This is the action-derived integrability prerequisite; the tooth desymmetrizes the displayed cross
coefficient (`1/2 → 3/4`) so `Hessian ≠ displayed operator`.

## 5. Dimensional firewall (EARNED — units-restored `[L,T,M]`)

The electric-scalar-block terms are checked exponent-triple homogeneous, units restored, able-to-fail.
The four-dimensional parent-`H` terms target the **bulk** density `E/L⁴`; the three-dimensional reduced
`(u_L, h)` terms target the **brane** density `E/L³`:

```text
{ H_kin, H_grad, H_pot }                         → E/L⁴  (bulk)
{ u_kin, u_stiff, h_kin, h_stiff, mix }          → E/L³  (brane)
```

with `[H] = L⁻¹`, `[u_L] = L`, `[h] = 1`, `[M₄] = M`, `[K₄] = E`, `[A_eff] = M L⁻³`,
`[B_eff] = M L⁻¹ T⁻²`, `[M_h] = M L⁻¹`, `[K_h] = M L T⁻²`, `[C_hu] = M T⁻²`, `[c_E] = L T⁻¹`. The
`mix`-term tooth corrupts `[C_hu]` (`M T⁻² → M T⁻¹`) and the homogeneity residual fires. *(The mouth
Robin/source term dimensions travel with the mouth SOURCE in stage 031, not here.)*

## Numeric witness slice (star `*` = the §7 nondimensional witness units)

All numeric witnesses carry a star; decimals are quarantined to this slice.

```text
ℓ* = ℓ/a = 1/20,     a* = 1
N₀* = 160/3
M₄* = K₄* = 3/160          (M₄ input; K₄ = M₄ c_E² derived)
M_h* = K_h* = c_E* = 1      (M_h* = 1 is a CONV normalization ON M₄)
A_eff* = 1  (⇒ M* = I),  B_eff* = 2,  C_hu* = 1/2
D* = B_eff* K_h* − C_hu*² = 7/4        (physical [D] = M² T⁻⁴, kept DISTINCT from the number 7/4)
c_±*² = (3 ± √2)/2   [decimals, labeled]:  0.7928932188…  and  2.2071067812…
```

## Scope, deferrals, and first-class departures (recorded honestly)

1. **The whole G0-v0 closure is a POSTULATED `DRAFT v0` closure**, not a derivation and not a final
   model. Every witness value (`M₄*`, `c_E*`, `B_eff*`, `C_hu*`, …), the Pöschl–Teller operator choice,
   and the `R_w`-odd parity representation of `{H, h, u_L}` are labeled postulates. What is EARNED is
   that, *given* this action, the static prerequisites of §§1–5 hold.

2. **Mouth SOURCE mechanism DEFERRED to stage 031.** `PASS_BARE_MOUTH_SOURCE` (the orientation-odd
   `χ↔h` mouth functional `η_i(k_m h − g_χh s_i)`, `Q_χ[r_Σ, s_i] = s_i`, `f₀(0) = 1/ℓ`, the nonzero
   projected source) is **not** folded here. It is the puncture-deflection mechanism (IV-2). ⚠ Binding
   requirement carried to 031: stage 031 must own a **full able-to-fail** `PASS_BARE_MOUTH_SOURCE`
   reconstruction (annulus `∫η = 1`, `f₀(0) = 1/ℓ` from the actual profile, `J_m → g_χh = J_m/ℓ`, the
   orientation projection, and the nonzero projected source) — not merely a re-typed assigned map.

3. **Out-of-scope, UNDISCHARGED cross-sector G0 predicates (Part VII de-dup obligation).** The source
   checker is monolithic across all seven G0 gates. Stage 030 folds ONLY the charge-new electric-scalar
   subset. The remaining G0 predicates are **NOT re-folded and NOT already covered** by the built
   Part-I/II stages (those verify physically-adjacent but *different* facts) — they are
   **`OUT-OF-SCOPE / UNDISCHARGED CROSS-SECTOR G0 PREDICATE`**, whose reconciliation is a Part VII
   G0-vs-Part-I de-dup obligation:
   - `PASS_PLANAR_WALL_FACTORIZATION` (G0 wall Hessian; stage006 verifies a kink EL + tension, not this);
   - `PASS_BULK_U_SECOND_DERIVATIVE` (only partially overlapped by stage005),
     `PASS_BULK_FOURIER_DENSITY_HESSIAN`, `PASS_PHASE_CONSTANT_MODE_QUOTIENT` (G0 bulk witness
     `{K*=4/5, ħ*=√2/10, μ*=1, U''(ρ₀)*=4}` + quantum Fourier Hessian);
   - `PASS_DRAIN_KERNEL_NORMALIZATIONS`, `PASS_DRAIN_MASS/MOMENTUM/ENERGY_CONTROLLER` (G0 `D_i/R_0`
     normalization + controllers; stages 008–010 are a return-target constraint spec, not these).

4. **Excluded assembled claims → Part VII.** The class-(2) assembled one-body and class-(3)
   two-body/far-field claims (ghost-free constrained Hessians, the dressed `1/R` monopole, pair-force
   integrability) are **not** discharged here; they are solved-response objects for Part VII.

5. **`c_E` has no committed cone lock.** `c_E` is the first *dimensional* entry of the electric sector
   and is `FREE-UNREDUCED`; `c_s/c_E = 2` was chosen expressly to avoid imposing a cone lock. Its cone
   relation to `c_γ`/`c_s` is re-adjudicated in Part VI (`pathA_40`, `PENDING`).

## Parameter-ledger labels (honest earned/derived/imposed split; per-stage register update carries these)

- `H` — `ACTION` (postulated PT operator), `L⁻¹`. `f₀`, `N₀ = 8/(3ℓ)` — `DERIVED` (`L⁻¹`).
  `h = P₀H` — `DERIVED`, dimensionless (the normalized zero-mode projection).
- `M₄` (`ACTION`) and `c_E` (input) are the inputs; **`K₄ = M₄ c_E²` is DERIVED**; `M_h = N₀ M₄`
  (`DERIVED`, `M L⁻¹`) and `K_h = N₀ K₄ = M_h c_E²` (`DERIVED`, `M L T⁻²`). `M_h* = 1` is a **`CONV`**
  normalization on `M₄`.
- `ℓ` — the length-**ratio** `ℓ/a = 1/20` is `IMPOSED`/`CALIB` (frozen handoff scale); only the unit
  pin `a` is `CONV`. (Do not tag `ℓ` ambiguously.)
- `A_eff = ρ_br + C_J²/κ_phase` — the *relation* is `DERIVED` (Schur elimination of `θ_B`); the numeric
  witness `A_eff* = 1` and the partition `ρ_br/A_eff = 3/4`, `(C_J²/κ_phase)/A_eff = 1/4` are
  `[POSTULATE]` (IMPOSED edge — not a reduction).
- `D = B_eff K_h − C_hu²` — `DERIVED`, `M² T⁻⁴`, dimensionless witness `D* = 7/4`.
- `B_eff = ρ_B0²/χ_c` — **CONSUMED** from built stage003 (R16), `DERIVED`; do NOT re-count. Its numeric
  witness `B_eff* = 2` is an IMPOSED edge (consumption does not make the *value* derived).
- `C_hu` — **stage 030 is its FIRST built dual-engine dimensional verification** (stage 041 is NOT
  built; any earlier "consumed from 041" attribution was wrong). Keep `C_hu` `FREE-UNREDUCED`/`PENDING`
  (reduction route Part VI/`pathA_41`); its now-verified dim is `M T⁻²`; witness `C_hu* = 1/2` is an
  IMPOSED edge.

No IMPOSED/CALIB value is dressed as DERIVED; consumption ≠ derivation of the numeric witness.

## Verification

- **Dual-engine, both exit 0, genuinely independent routes.**
  `scripts/ledger_stage030_electric_scalar_localized_h_closure_sympy_audit.py` — **SymPy 16 PASS**.
  `mathematica/ledger_stage030_electric_scalar_localized_h_closure_mathematica_audit.wl` —
  **Mathematica 16 PASS**. Standalone, print-only, assert-zero (`raise SystemExit(1)` / `Exit[1]` on
  failure), no argparse harness, no JSON/YAML payload I/O, zero file-I/O between engines.
- **The `.wl` is a genuinely independent route, not a transliteration.** It reorders the derivation
  (native matrix route first) and uses native Wolfram machinery where the `.py` uses SymPy: `Integrate`
  for `N₀`, `DSolve` (counting one integration constant) for the `ker A` uniqueness, `Limit` for the
  `4/ℓ²` continuum threshold, and `PositiveDefiniteMatrixQ` / `CharacteristicPolynomial` / `Eigenvalues`
  for the coupled kernel — versus SymPy's integrating-factor kernel, `sp.limit`, `is_positive_definite`,
  and `sp.solve`.
- **16 per-tooth ablations, each `FIRED_AT_OWN_ASSERT`** (env switch `LEDGER_STAGE030_MUTATION`): a
  card-primitive mutation for every folded predicate, including the two **distinct** kernel teeth
  (finding 5: `C_hu*` → margin fires `PASS_STABILITY`; inertia `A_eff*` → roots fires
  `PASS_POSITIVE_GENERALIZED_WAVE_SPEEDS`), the genuine reduction teeth (the `[K₄]` dimension mutation on
  `PASS_PARENT_STIFFNESS_DIMENSIONAL_CONSISTENCY`; the `K_h=N₀K₄→N₀²K₄` reduction-inconsistency mutation on
  `PASS_REDUCED_SPEED_PRESERVATION`; the norm/projection weight on `N₀`/`h=P₀H`; `M₄*` on the inertia/norm),
  the spectral teeth (factorization coefficient, `f0_power`, `V_H`-asymptote → gap), and the mundane dim
  tooth (`[C_hu]`). Out-of-scope predicates carry no teeth (dropped, not silently green). Float-free
  symbolic payload throughout (machine-`Float` atoms raise). *(A pre-commit adversarial re-verify replaced
  two originally-tautological `c_E`-relation value self-checks with the genuine dimensional/speed-preservation
  teeth above — see §2.)*

## Downstream consumers

- **`ledger_stage031` (IV-2, puncture-deflection MECHANISM):** consumes `f₀` (with `f₀(0) = 1/ℓ`), `N₀`,
  `h = P₀H`, the reduced kernel `S_Lh`, and the positive block `D` (its witness `D* = 7/4` and the response
  objects `z_g`, `m_gg = B_eff z_g²/D`, `det m = z_g²/D` are IV-2 objects). Stage 031 adds the mouth `χ↔h`
  source and the exterior `h(r) = h_A a/r` ⇒ the target-blind EARNED far-field FORM (verdict
  `THROAT_H_SOURCE_1_OVER_R2`). *(4-stage split, 2026-07-22: the `R1_REQUIRED(bc_selection)` sign landing is
  the SEPARATE **stage 032**, which consumes 031's `m`/`m_gg`/`S_gg` + the far-field shell.)*
- **Parameter register:** rows added/re-homed as in the label section above; `c_E` seeded here (no cone
  lock, Part VI re-adjudicates); `B_eff` consumed (stage003 R16), not re-counted; `C_hu` gets its first
  built dim verification here.
- **Part VII (integration):** owns the G0-vs-Part-I de-dup of the out-of-scope G0 predicates (§Scope 3),
  the assembled/solved claims (§Scope 4), and the shared throat-solve that would convert the electric
  sign/magnitude `R1` to resolved.

## Provenance

- **Physics source:** `software/em_charge_attribute/g0_closure_card_v0.md` §§2.3–2.4 (the localized
  parent `H`, the zero-mode reduction, and the `(u_L, h)` scalar block), with dims from §§1, 7 and the
  committed static-check record in §10.4. The card is `DRAFT v0` — a postulated closure for the
  Codex→Grok→Codex gauntlet, not a derivation.
- **Source checker (monolithic):** `g0_closure_card_v0_checks.{py,wl}` spans all seven G0 gates; stage
  030 folds **only** the electric-scalar subset (finding 1). The bulk/wall/drain predicates are the
  undischarged cross-sector obligations of §Scope 3, not covered here.
- **Ratified split (context):** `research/pde_ledger_v2/notes/part4_charge_atomic_split.md` (IV-1 row +
  register preview). Reshape directive + review trail:
  `research/pde_ledger_v2/_scratch/stage030_reshape_directive.md` (v2, folds 7 Codex findings).
- **Governing:** `notes/ledger_v2_blueprint.md` §5 (reshape) + §6 (verification protocol);
  `docs/model_map.md` §3.4 / §5.1.

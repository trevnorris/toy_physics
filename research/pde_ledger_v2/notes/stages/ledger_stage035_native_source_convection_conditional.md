# ledger_stage035_native_source_convection_conditional

## Status

**Part V — Magnetism. V-2 (build-order 035; the SECOND stage of the 6-stage Part-V split,
user decision 2026-07-22).** The self-contained native source-law + parity-census cluster of the
moving-throat (magnetism) sector. Where stage 034 (V-1) built the *action row* — importing the
pathA_36 transverse-vector kinetic/gradient row into G0 and writing the one finite-profile moving
coupling `q_T Σ_i s_i η_a V_i·u_T` — stage 035 asks the complementary *source* question: what current
does the moving ±w throat actually source, when you translate its **actual signed dent** under local
continuity? — and records the honest answer:

- **`CONVECTION_LIKE_CONDITIONAL`** — the native moving-throat source is the continuity-derived
  convective current, its tensor form is a profile-smeared `sV`, and its magnitude is conditional on
  the unfixed R1 throat charge `q_T`.

  ```text
  σ_i = s_i η_a(x − X_i(t)),     ∂_t σ_i + ∇·(σ_i V_i) = 0
        ⇒ continuity residual = (α − 1)·Σ_k V_{i,k} ∂_{x_k} σ_i = 0  ⇔  α = 1  (UNIQUE),
  I_i = s_i η_a V_i,     J_{T,i} = q_T I_i = q_T s_i η_a V_i,     q_T = λ_T τ_d.
  ```

  This native source is **DERIVED from signed-dent continuity, not imported** from the barred
  `j ∝ sV` / pathA_39 current — enforced at runtime by the `SOURCE_NOT_IMPORTED` free-symbol-disjointness
  guard (the computed source vector's free-symbol set is disjoint from the barred amplitudes
  `{N_u, a_T, a′_T, a_L, q_A^T, q_L}`). The stage also derives the full **24-cell parity census** of
  `{s, V, τ_d/q_T, J_T, u_T, b_T}` under `R_w` / `P_w` / rotation / time-reversal. Verdict token
  **`CONVECTION_LIKE_CONDITIONAL`** (both engines, both exit 0).

**⭐ Scope class: EARNED (source law + parity census, target-blind).** This is an EARNED result — the
continuity translation fixes the unique isotropic flux coefficient `α = 1` and yields the native source,
and the parity census is a target-blind structural fact. It is a **CurrentIdentity CLASSIFICATION**
(the source-identity token), **NOT a sealed landing** — the sealed 1152-cell first-match landing
(`R1_REQUIRED(electric_bc_selection)` + co-blockers) is **stage 038 (V-5)**. It is **NOT** a knob,
**NOT** a reduction/codimension edge, and it does **not** shrink the irreducible count. It carries no
electric branch/sign; the source's MAGNITUDE rides the unfixed `q_T` (R67), which is exactly why the
compact-limit classification is *conditional*.

**⚠ CRITICAL TYPING (import vs earn).** `q_T = λ_T τ_d` (the moving-throat charge), `τ_d` (the
active-drain time-arrow), and the transverse shear field `u_T` with its finite throat profile `η_a`
(`[η_a] = L⁻³`) are **CITED from stage 034 (V-1)** (register edges R67/R68) and, upstream of that, from
**Part III / stage003** (`u_T`, `∇·u_T = 0`, `c_γ² = μ_R/ρ_br`). Stage 035 **re-instantiates** these
imported symbols for its own dual-engine check but does **NOT** re-derive them, does **NOT** re-count
`{q_T, τ_d, u_T, η_a}` in the irreducible set, and does **NOT** re-open the action row. The
magnetism-NEW content of 035 is **exactly** `{the continuity-derived native source law, the parity
census}`.

**⚠ CRITICAL BOUNDARY (`b_T` parity is RAW here, not the departure).** The census records `b_T = ∇×u_T`
as an **axial** vector and **time-reversal EVEN** as **RAW computed facts**. Their interpretation as the
"magnetism is NOT exact Maxwell" DEPARTURE (`B_TIME_REVERSAL_EVEN` — a physical Maxwell `B` is T-odd) is
**stage 039 (V-6)**, which CITES these census facts. Stage 035 does **NOT** frame `b_T`'s parity as the
departure and does **NOT** run the `b_T`-vs-Maxwell comparison.

## Purpose

Record, as an EARNED source-law + parity result, that magnetism's native current is the
continuity-derived convective source and that its magnitude is conditional on the R1 throat charge. The
decisive object is the **continuity residual**: translating the actual signed dent
`σ_i = s_i η_a(x − X_i(t))` — the ±w throat's finite mouth profile moving rigidly with the throat
velocity `V_i` — and imposing the local conservation law `∂_t σ_i + ∇·(σ_i V_i) = 0` leaves a residual
`(α − 1)·Σ_k V_{i,k} ∂_{x_k} σ_i` for an isotropic local flux with coefficient `α`; it vanishes for
arbitrary translating profiles **only** at `α = 1`. The stage certifies at runtime that (i) `α = 1` is
the UNIQUE root of the continuity residual, (ii) the resulting native source is `I_i = s_i η_a V_i`,
`J_{T,i} = q_T I_i = q_T s_i η_a V_i` and reduces exactly to that continuity-flux basis (no
doubled/rescaled/ansatz basis), (iii) the source is DERIVED — its free-symbol set is disjoint from the
barred pathA_39 amplitudes (the `SOURCE_NOT_IMPORTED` surviving-solution guard, a computed free-symbol
object, NOT a grep), (iv) the six-object parity census reproduces (`s`, `V`, `τ_d`/`q_T`, `J_T`, `u_T`,
`b_T` under `R_w`/`P_w`/rotation/T), and (v) the whole-stage dimensional firewall closes on the real
expressions. Both engines run the same seven Q-CURRENT teeth plus the build-global
target-blindness / dual-engine / units-restored guards, the computed-verdict re-derivation, and the
source-to-stage manifest — and each independently reaches `CONVECTION_LIKE_CONDITIONAL`.

Consumes **nothing numeric** from 030/031/032/033 — it is an independent source-law + census construction
that CITES `q_T`/`τ_d`/`u_T`/`η_a` from stage 034 / stage003 as provenance (see §10). It is the source
twin of 034's action row.

## 1. The signed-dent continuity translation → the unique flux coefficient (`SOURCE_TRANSLATION_CONTINUITY`)

The ±w throat is a finite signed dent of the brane order parameter,

```text
σ_i(x, t) = s_i η_a(x − X_i(t)),        X_i(t) = X_i(0) + V_i t,        ∫ d³x η_a = 1,
```

where `s_i = ±1` is the throat's `±w` orientation (the charge-sector sign) and `η_a` is the normalized
finite mouth profile (`[η_a] = L⁻³`), rigidly translating with the throat velocity `V_i`. The
translation-derivative identity for a rigidly moving profile is

```text
∂_t σ_i = s_i ∂_t η_a(x − X_i(t)) = −s_i Σ_k V_{i,k} ∂_{x_k} η_a(x − X_i(t)) = −V_i·∇σ_i.
```

Imposing local continuity with a candidate **isotropic local flux** `α σ_i V_i` (coefficient `α` to be
fixed) gives

```text
∂_t σ_i + ∇·(α σ_i V_i) = −V_i·∇σ_i + α V_i·∇σ_i = (α − 1)·Σ_k V_{i,k} ∂_{x_k} σ_i.
```

Since `V_i` is arbitrary and `σ_i` is a genuine translating profile (its gradient is not identically
zero), the residual `(α − 1)·Σ_k V_{i,k} ∂_{x_k} σ_i` vanishes for all admissible configurations
**only** at

```text
α = 1        (the UNIQUE root; solved as the root of the continuity-residual object).
```

The unique local flux is therefore `σ_i V_i` — the convective flux — and there is no free flux
coefficient. (Selecting a wrong coefficient — `α = 2`, the build's `source_scale = 2` mutation — or
dropping a continuity term makes the residual recompute NONZERO and the tooth fires at its own assert.
⚠ This is the primary computed object the verdict re-derives from.)

## 2. The native moving-throat source law (`SOURCE_BASIS`)

With `α = 1` fixed by continuity, the native per-defect flux current and transverse field source are

```text
I_i = s_i η_a V_i,        J_{T,i} = q_T I_i = q_T s_i η_a V_i,        q_T = λ_T τ_d,
```

where `q_T = λ_T τ_d` is the moving-throat charge (coupling strength `λ_T` times the active drain's
time-arrow `τ_d`), CITED from stage 034. This is exactly the current that appears in 034's coupling
`q_T Σ_i s_i η_a V_i·u_T` (its `J_T·u_T` structure) — 035 derives it from continuity rather than
positing it. The source reduces **component-wise** to precisely this continuity-flux basis:

```text
J_{T,i,k} = q_T s_i η_a V_{i,k}     for each component k,
```

not a doubled, rescaled, or ansatz basis. (Doubling/rescaling the source — `J → 2J`, the `SOURCE_BASIS`
mutation — makes the component-wise zero-difference against the derived `q_T s η V` recompute nonzero and
the tooth fires.)

## 3. The surviving-solution guard: DERIVED, not imported (`SOURCE_NOT_IMPORTED`)

The native source **REPLACES** the barred `j ∝ sV` / pathA_39 current — it is DERIVED from signed-dent
continuity, **NOT** imported from the pathA_39 amplitudes. This is the surviving-solution discipline: the
ledger shows the surviving native source clean; the barred pathA_39 current is NOT imported and its
post-mortem scope goes to the failures-paper backlog. The discipline is enforced at runtime as a computed
**free-symbol-disjointness** object:

```text
freesymbols(J_{T,i})  ∩  BARRED_SOURCE_MARKERS  =  ∅,
BARRED_SOURCE_MARKERS = { N_u, a_T, a′_T, a_L, q_A^T, q_L }   (code: Nu, aT, aTp, aL, q_A_T, q_L),
```

i.e. the free-symbol set of the computed source vector is disjoint from the barred pathA_39 amplitude
markers. **⭐ This is a real runtime free-symbol object, NOT a source grep** — a grep for a forbidden
token is dodgeable (string-concat / `chr()`), whereas a free-symbol-set intersection is enforced on the
computed source expression itself. (Injecting a code-native barred marker — e.g. `Nu`, from
`BARRED_SOURCE_MARKERS` — into the computed source's free-symbol set makes the intersection nonempty and
the tooth fires; the property is enforced on the COMPUTED source, so no string dodge defeats it.)

## 4. The derived parity census — 24 cells (`PARITY_RW`, `PARITY_PW`, `PARITY_ROTATION`, `PARITY_TIME_REVERSAL`)

The census is a **target-blind structural fact**: six objects × four discrete transformation columns
(`R_w`, `P_w`, rotation type, time reversal) = **24 cells**, computed as invariants of the derived source
and field (not a re-typed lookup):

| object | `R_w` | `P_w` | rotation | time reversal |
|---|---:|---:|---|---:|
| `s` (±w orientation) | odd | odd | scalar | even |
| `V` (throat velocity) | even | even | polar vector | odd |
| `τ_d`, hence `q_T` | even | even | scalar | odd |
| `J_T = q_T s η_a V` | odd | odd | polar vector | even |
| `u_T` (transverse shear displacement) | odd | odd | polar vector | even |
| `b_T = ∇×u_T` | odd | odd | axial vector | even |

The four census teeth read off this table as computed invariants:

- **`PARITY_RW`** — `J_T` and `u_T` are both `R_w`-**odd** (odd under `w`-reflection).
- **`PARITY_PW`** — `J_T` and `u_T` are both `P_w`-**odd** under the `w`-type reflection.
- **`PARITY_ROTATION`** — `u_T` is a **polar** vector and `b_T = ∇×u_T` is an **axial** vector (the curl
  of a polar vector is axial).
- **`PARITY_TIME_REVERSAL`** — `τ_d` and `V` are **T-odd**, while `J_T` and `u_T` are **T-even**.

**⚠ The `P_w` w-reflection caveat (fold verbatim).** Here **`P_w` denotes a `w`-type reflection of the
transverse coordinate**, NOT ordinary three-dimensional spatial parity `x → −x`. Under that reading the
`s`, `∇`, and `b_T` assignments are self-consistent. (This caveat is load-bearing: `s` is the ±w
orientation, which flips under a `w`-reflection but not under ordinary `x → −x`.)

**⚠ `b_T` parity is RAW here.** The census records `b_T`'s **axial + T-even** parity as RAW derived
facts. A passive time-reversal-EVEN throat would not supply this `O(V)` action row, so magnetism requires
G0's active-drain time-arrow `τ_d` (T-odd). The interpretation of `b_T`'s T-even parity as the "not exact
Maxwell" DEPARTURE (`B_TIME_REVERSAL_EVEN`) is **stage 039 (V-6)**, NOT here — see §8. (Flipping any
census entry — `u_T`'s `R_w` or `P_w` to even, demoting `u_T` to scalar or making `b_T` polar, or
flipping `J_T`/`τ_d`/`V` in the T-column — makes the corresponding parity object fail and the tooth
fires; the verdict then re-derives to `PARITY_INCONSISTENT`, see §11.)

## 5. Why `CONVECTION_LIKE_CONDITIONAL` — a classification, not a landing

The source law is **fully DERIVED in tensor form**: `α = 1` is fixed, and `J_{T,i} = q_T s_i η_a V_i`
reduces in the compact limit to a profile-smeared `s V` — a convective current. But the **magnitude** is
conditional: `q_T = λ_T τ_d` rides the R1 nonlinear moving-throat solve and is **not fixed**. Hence the
compact-limit token is `CONVECTION_LIKE_CONDITIONAL` — "convection-like" (the tensor form is the
convective `sV`), "conditional" (its magnitude is conditional on the unfixed `q_T`).

This is a **CurrentIdentity CLASSIFICATION** — the source-identity token — **NOT a sealed landing**. The
build-native `CurrentIdentity` enum is
`{CONVECTION_LIKE_CONDITIONAL, CHARACTERIZED_SOURCE_DEPARTURE, NULL_SOURCE, R1_SOURCE_BASIS}`; stage 035
classifies the native source into `CONVECTION_LIKE_CONDITIONAL`. The sealed §4 first-match landing (the
1152-cell truth table → `R1_REQUIRED(electric_bc_selection)` + the three co-blockers) is a **downstream**
object, **stage 038 (V-5)** — 035 emits no landing. An **unresolved magnitude is NOT an unresolved source
basis**: production's R1 magnitude still yields `CONVECTION_LIKE_CONDITIONAL`, NOT `R1_SOURCE_BASIS` (see
§11).

## 6. The dimensional identity — the whole-stage firewall (`UNITS_RESTORED`)

Units RESTORED (natural-units-independent, able-to-fail, free-carrier-independent). The preregistered
field is the transverse shear **displacement** `u_T` (`[u_T] = L`); the curl candidate `b_T = ∇×u_T` is
**dimensionless**. The dimensions of every real expression in scope:

```text
[σ]     = L⁻³           (signed dent density)
[η_a]   = L⁻³           (normalized mouth profile, ∫η_a = 1)
[V]     = L T⁻¹
[q_T]   = M T⁻¹
[I]     = L⁻² T⁻¹       (= [η_a][V]; the flux current s η_a V)
[J_T]   = M L⁻² T⁻²     (= [q_T][I])
[u_T]   = L
[b_T]   = 1             (= [∇×u_T])
```

The source identity `J_T = q_T s η_a V` closes:

```text
[q_T][s][η_a][V] = (M T⁻¹)(1)(L⁻³)(L T⁻¹) = M L⁻² T⁻² = [J_T].     ✓
```

Every parity assignment is dimensionless/structural (a discrete ±1 eigenvalue, no `[L,T,M]`). (Corrupting
one dimension — e.g. `[J_T] = ML⁻²T⁻² → ML⁻¹T⁻²`, or `[q_T] = MT⁻¹ → MLT⁻¹`, or `[b_T] = 1 → L⁻¹` —
makes the dimensional-homogeneity object inhomogeneous and the units tooth fires; free-carrier-independent.)

## 7. Cite-not-count discipline: what 035 consumes vs earns

Stage 035 **imports** and re-instantiates, but does NOT re-derive or re-count:

```text
CITED (provenance, NOT re-counted; de-dup deferred to Part VII):
  q_T = λ_T τ_d   (stage 034, register R67 — the moving-throat charge)
  τ_d             (stage 034, register R68 — the active-drain time-arrow, T-odd)
  u_T, η_a        (stage 034 / stage003 — the transverse shear field + finite throat profile)

EARNED-NEW (the magnetism-NEW content of 035):
  { the continuity-derived native source law  J_{T,i} = q_T s_i η_a V_i,
    the 24-cell parity census }
```

The source law introduces **NO new irreducible knob** — its only free magnitude is the already-registered
R67 `q_T`. `q_T`, `τ_d`, `u_T`, `η_a` are re-instantiated for the dual-engine check but are cited, not
re-counted (they are NOT re-added to the irreducible set; de-dup against stage 034 / Part III is deferred
to Part VII). The parity census is a structural/dimensionless fact — it does NOT shrink the irreducible
count.

## 8. Scope — what 035 does NOT do (deferred to 036–039)

035 is the **source-law + parity-census gate only**. Named explicitly as OUT OF SCOPE:

- **The transverse-vector action row** `S_{T+move}` (the `½ρ_br|u̇_T|² − ½μ_R|∇×u_T|² + q_T Σ s_i η_a V_i·u_T`
  amendment, `internal_inconsistency = none`, transverse stability, `G0_DAMAGE`) is **stage 034 (V-1) —
  DONE**. 035 CITES the coupling term's `J_T·u_T` structure but does NOT re-derive the action, the
  Hessian, or the G0-damage accounting.
- **Route A (the Maxwell–Darwin reference kernel `(δ_ij + n_i n_j)/8πR`, boost of the electric
  interaction)** is **stage 036 (V-3)** — 035 derives the source but computes NO far field.
- **Route B (direct moving-throat shear, blind to A), the structural COMPARISON, and the ratios**
  (`r_BA = q_T²/(ρ_br A_E)`, `δ_BA = r_BA − 1`, `r_cone = c_E²/c_γ²`, `ΔU`) are **stage 037 (V-4)** — NOT
  here.
- **The SEALED §4 first-match landing** — the 1152-cell truth table → `R1_REQUIRED(electric_bc_selection)`
  + the three co-blockers (`R1_REQUIRED(direct_moving_throat)` / `(magnitude)` / `(consistency)`) and the
  active-`F_flux` caveat — is **stage 038 (V-5)**. 035 emits NO landing; its verdict is the source-identity
  classification only.
- **⚠ CRITICAL boundary — the `b_T = ∇×u_T` time-reversal-EVEN DEPARTURE interpretation**
  (`B_TIME_REVERSAL_EVEN`, the "not exact Maxwell" characterization — a physical Maxwell `B` is T-odd) is
  **stage 039 (V-6)**. Stage 035 **COMPUTES the parity census, including `b_T`'s axial + T-even parity, as
  RAW derived facts**, but does **NOT** interpret them as the Maxwell departure and does **NOT** run the
  `b_T`-vs-Maxwell comparison. The boundary is explicit: **035 = the census facts; 039 = the departure
  characterization that CITES those facts.**
- **The imported `q_T = λ_T τ_d`, `τ_d`, `u_T`, `η_a`** are **CITED provenance** from stage 034 / stage003
  — NOT re-derived, NOT re-counted in the irreducible parameter set (de-dup deferred to Part VII).
- **`q_T`'s value / the throat normalization** is the sim-deferred R1 throat solve — the source law's
  MAGNITUDE is conditional on `q_T` (this is exactly why the verdict is `CONVECTION_LIKE_CONDITIONAL`,
  "conditional"). 035 does NOT resolve `q_T`.

## 9. Source-to-stage predicate manifest

Completeness certificate: **no silently-dropped source claim.** Every source tooth in scope from the
source build (`magnetism_moving_throat_check.{py,wl}`, the 35-tooth `TOOTH_ORDER`, commit `53cf049f`)
lands as **PRESERVED** (folded as-is / a stronger reconstruction owned here), **REPLACED_BY_STRONGER** (a
stronger whole-stage reconstruction), or **SCOPED_OUT** (deferred to its own Part-V stage, each named
LITERALLY — NO wildcards). The partition is disjoint + exhaustive, computed at runtime in **both** engines
from the same canonical `(id, disposition, owner)` triples:

```text
partition = { PRESERVED: 8, REPLACED_BY_STRONGER: 2, SCOPED_OUT: 25 }   (35 total),
manifest digest (SHA-256) = d85de3d8f6b7d615900c8ead24f1589eed3c681e74cf546adb999f2cc7fa5b36.
```

- **PRESERVED (8):** the seven Q-CURRENT teeth owned by 035 — `SOURCE_TRANSLATION_CONTINUITY`,
  `SOURCE_NOT_IMPORTED`, `SOURCE_BASIS`, `PARITY_RW`, `PARITY_PW`, `PARITY_ROTATION`,
  `PARITY_TIME_REVERSAL` — plus the build-global `TARGET_BLINDNESS` (folded as the source-side-dependencies
  blindness guard).
- **REPLACED_BY_STRONGER (2):** `UNITS_RESTORED` (→ whole-stage dimensional firewall on the real
  expressions) and `DUAL_ENGINE_TERMS` (→ canonical term inventory across both engines).
- **SCOPED_OUT (25)** — the source teeth of V-1 (done upstream) and V-3…V-5, each named verbatim (same
  partition both engines):
  - **V-1 (034, owned UPSTREAM — DONE):** `FIELD_IDENTITY_UNITS`, `ACTION_KINETIC`, `ACTION_COUPLING`,
    `ACTION_STABILITY`, `G0_DAMAGE`, `LEDGER_READY_ROW` (6)
  - **V-3 (036):** `BOOST_PROJECTOR`, `BOOST_GENERAL_VELOCITIES`, `BOOST_NEXT_ORDER` (3)
  - **V-4 (037):** `ROUTE_INDEPENDENCE`, `BOOST_COMMON_VELOCITY`, `DIRECT_SOURCE`, `DIRECT_PROJECTOR`,
    `DIRECT_EXCHANGE_SIGN`, `DIRECT_FALLOFF`, `DIRECT_VELOCITY_ORDER`, `COMPARE_COMPUTED`, `DELTA_RATIO`,
    `CONE_RATIO`, `QMAG_R1` (11)
  - **V-5 (038):** `ACTIVE_FLUX_CAVEAT`, `HOOK_LORENTZ`, `TRUTH_TOTALITY`, `TRUTH_PRECEDENCE`,
    `LANDING_OWNERSHIP` (5)

  **V-6 (039) owns NO additional source-build tooth** — it REUSES 035's parity results
  (`PARITY_TIME_REVERSAL` + `PARITY_ROTATION`) and authors its own NEW `B_TIME_REVERSAL_EVEN` departure
  assertions, cited-not-owned.

Both engines assert the identifier set (35 unique), the three-way disposition set, the exact counts, and
the committed digest at runtime and agree (`SOURCE_TO_STAGE_MANIFEST`; the mutation drops the final
partition row — the `DUAL_ENGINE_TERMS` `REPLACED_BY_STRONGER` row — which changes the identifier set /
count / digest and fires the tooth; dropping any one of the 35 rows fires it).

## 10. Consumes / cites

- **Consumes NOTHING numeric from 030/031/032/033.** 035 does not touch the puncture-deflection response
  matrix `m`, `m_gg`, the four `A_X`, the BC-ensemble data, or the native-`P` constraint signatures — it
  is an independent source-law + census construction.
- **Cites (provenance, NOT consumed numerically, NOT re-counted):** the moving-throat charge
  `q_T = λ_T τ_d`, the active-drain time-arrow `τ_d`, and the transverse shear field `u_T` with its finite
  mouth profile `η_a` (`∫η_a = 1`, `[η_a] = L⁻³`) — from **stage 034 (V-1)** (register R67/R68); and
  upstream of that, `u_T` (`∇·u_T = 0`, `[u_T] = L`) and `c_γ² = μ_R/ρ_br` from **Part III / stage003 /
  pathA_36**. De-dup of `{q_T, τ_d, u_T, η_a}` against stages 034/003 is deferred to Part VII; they are
  NOT re-added to the irreducible count here. Also cites the charge-sector `±w` orientation `s_i`
  (Part IV) as the source's sign carrier.

## Verification

- **Dual-engine, both exit 0, genuinely independent routes.**
  `scripts/ledger_stage035_native_source_convection_conditional_sympy_audit.py` — **SymPy 12 teeth**.
  `mathematica/ledger_stage035_native_source_convection_conditional_mathematica_audit.wl` —
  **Mathematica 12 teeth**. Standalone, print-only, assert-zero (`raise SystemExit(1)` / `Exit[1]`), no
  argparse harness, no JSON/YAML payload, **zero file-I/O between engines**. Each engine reaches
  `CONVECTION_LIKE_CONDITIONAL` on its own and prints its own tokens — cross-engine agreement is that they
  independently produce the same tokens, not a compare pass.
- **The `.wl` is a genuinely INDEPENDENT route** — materially distinct construction, not a line-port of
  the `.py`. The **SymPy** route builds the continuity residual via the **abstract translation identity**
  `∂_t σ = −V·∇σ` and `sp.solve`s the residual for the flux coefficient `α`; the **Mathematica** route is
  independent — it constructs `σ` with an **explicit normalized Gaussian mouth profile with time-dependent
  centers** `X_i(t)`, forms the continuity residual by native `D[·,t]` / `Div` differentiation, `Reduce`s
  the vanishing condition for `α`, generates the parity transforms via **operator-generated reflection /
  rotation / time operators** applied to the tabulated objects, classifies polar vs axial via
  `Det` / `Cross`, and derives the dimension ratios by **scale-substitution**. Each engine
  differentiates/derives its OWN source. **No dual-engine disagreement.**
- **Per-tooth ablation** (env switch `LEDGER_STAGE035_MUTATION`): all **12 mutations FIRED_AT_OWN_ASSERT
  per engine** — the seven Q-CURRENT teeth (`SOURCE_TRANSLATION_CONTINUITY`, `SOURCE_NOT_IMPORTED`,
  `SOURCE_BASIS`, `PARITY_RW`, `PARITY_PW`, `PARITY_ROTATION`, `PARITY_TIME_REVERSAL`), the build-global
  `TARGET_BLINDNESS` / `DUAL_ENGINE_TERMS` / `UNITS_RESTORED`, the `VERDICT_REDERIVATION` tooth, and the
  `SOURCE_TO_STAGE_MANIFEST` tooth. `SOURCE_NOT_IMPORTED` fires when a **code-native barred symbol** (`Nu`,
  from `BARRED_SOURCE_MARKERS`) is injected into the computed source's free-symbol set — a **computed
  free-symbol intersection, not a grep** (no string-concat/`chr()` dodge defeats it).
- **The verdict tooth is non-tautological** — it mutates a **COMPUTED** object, not the final token. The
  precedence it re-derives (from the computed continuity residual, the basis residuals, the barred-symbol
  intersection, and the census invariants — never a duplicated literal):

  | Computed state | Re-derived token |
  |---|---|
  | valid continuity (`α = 1`) + basis + provenance (empty barred intersection) + consistent census | **`CONVECTION_LIKE_CONDITIONAL`** (production) |
  | nonzero continuity residual (`α ≠ 1`), wrong unique coefficient, failed `SOURCE_BASIS`, **or** nonempty barred-symbol intersection (`SOURCE_NOT_IMPORTED` fails) | **`R1_SOURCE_BASIS`** (the native basis is NOT earned) |
  | valid native source but a failed census invariant (a `PARITY_*` object fails) | **`PARITY_INCONSISTENT`** (the sole authored stage-local token) |

  Concretely: the verdict mutation sets `flux_coefficient = 2` → the continuity residual recomputes
  NONZERO → the pipeline RE-DERIVES **`R1_SOURCE_BASIS`** (a source resting on a broken/imported basis has
  NOT earned its native basis); a census flip re-derives **`PARITY_INCONSISTENT`**; and — crucially — an
  **unresolved MAGNITUDE (`q_T`, `MagnitudeFact.R1`) does NOT change the source-identity token**:
  production's R1 magnitude still yields **`CONVECTION_LIKE_CONDITIONAL`**, NOT `R1_SOURCE_BASIS` (an
  unresolved magnitude is not an unresolved source basis). Flipping the `CONVECTION_LIKE_CONDITIONAL`
  literal tests nothing (the 030 `X ≡ X` lesson). `PARITY_INCONSISTENT` is the SOLE authored stage-local
  addition to the build-native `CurrentIdentity` enum (same token/partition in both engines).
- **Tri-review outcome (falsification-first — recorded transparently, not hidden).**
  - **FIDELITY:** found **2 NITs — both can't-fail code** (neither could ever fail, so neither weakened a
    real check): (1) a guaranteed-true `structural_cells_ok` / `ParityStructural` conjunct ANDed into
    `UNITS_RESTORED` (a tautological sub-clause riding the real dimension check); (2) a vestigial
    `deriveVerdict` dead branch in the `.wl`. **Both REMEDIATED** — the tautological conjunct was removed
    so `UNITS_RESTORED` now asserts **only** the able-to-fail dimension-equality (and its mutation still
    fires), and the dead `deriveVerdict` branch was removed; the **manifest digest is unchanged**
    (`d85de3d8…`) → **CLEAN**.
  - **ADVERSARIAL: CLEAN with 3 DOCUMENTED non-blocking notes** (verified-safe, NOT remediated):
    1. The census **INPUT** rows `s` / `V` / `τ_d` are **fixed generators** while the **DERIVED** rows
       `J_T` / `u_T` / `b_T` are propagated from them and mutation-tested; the `u_T`-polar sub-conjunct is
       **true-by-construction** (same class as the stage-034 transversality note) — it cannot fail
       individually but rides a genuinely able-to-fail census invariant.
    2. The compound teeth carry **one mutation covering their sub-checks** — within the directive's
       permissive per-tooth-ablation "e.g." scope.
    3. A **cosmetic ablation-string label** in the manifest tooth (a descriptive string, not a computed
       object).

  Arbiter re-runs of both engines reproduce the tokens and the manifest digest `d85de3d8…`; the tri-review
  leaves the source law and census EARNED and the verdict `CONVECTION_LIKE_CONDITIONAL`.

## Downstream consumers

- **Part V continues** after 035: 036 (Route A — Maxwell–Darwin kernel) → 037 (Route B + boost-consistency
  comparison) → 038 (the sealed `R1_REQUIRED(electric_bc_selection)` landing) → 039 (the `b_T` T-even
  departure, which CITES 035's parity census). 035 is the source-law + census foundation those stages
  consume by citation.
- **Parameter register:** edge **R69** (`J_{T,i} = q_T s_i η_a V_i`, native moving-throat source,
  `[ML⁻²T⁻²]`, **DERIVED** from defect-continuity with the unique `α = 1`; value conditional on R67 `q_T`;
  REPLACES the barred `j ∝ sV`; introduces NO new knob — its magnitude rides R67) + the parity census as a
  **structural/dimensionless fact** (no `[L,T,M]`, no knob, not a reduction/codim edge; `b_T`'s axial +
  T-even parity recorded RAW, its departure interpretation booked at 039). `q_T`/`τ_d`/`u_T`/`η_a` are NOT
  re-counted.
- **`docs/model_map.md` §3.5:** the Q-CURRENT source-law bullet — the moving throat sources the
  continuity-derived convective current `J_T = q_T s η V`; source-identity `CONVECTION_LIKE_CONDITIONAL`,
  magnitude R1 via `q_T`, `b_T` T-even departure downstream.
- **Part VII:** the cited `{q_T, τ_d, u_T, η_a}` edges enter the de-dup (counted ONCE, with stage 034 /
  Part III); the `q_T` obligation enters the honest R1 debt alongside the electric `A_E`.

## Provenance

- **Physics source:** `software/em_charge_attribute/magnetism_moving_throat_result.md` (the "Q-CURRENT and
  field identity" section — the signed-dent continuity, the derived source `I_i = s_i η_a V_i` /
  `J_{T,i} = q_T I_i`, the `CONVECTION_LIKE_CONDITIONAL` classification, and the derived parity census
  table) + `software/em_charge_attribute/magnetism_moving_throat_check.{py,wl}` (the Q-CURRENT block, the
  35-tooth dual-engine source build, commit `53cf049f`). Stage 035 extracts the **source-law +
  parity-census cluster** into a focused standalone dual-engine pair; the source `argparse`/`--compare`
  harness, the `--ablate-tooth`/`--out-dir` payload plumbing, the Route A/B/comparison/truth-table
  machinery (out of 035's scope), and all JSON/log file writing are STRIPPED (print-only / zero-file-I/O /
  independent-tokens contract). The stage `.wl` is a materially distinct Wolfram route (explicit Gaussian
  mouth profile + native `D`/`Div`/`Reduce` + operator-generated parity transforms), not a line-port of
  the `.py`.
- **Consumes:** nothing numeric from 030/031/032/033 (independent construction).
- **Cites (provenance):** `q_T = λ_T τ_d`, `τ_d`, `u_T`, `η_a` from stage 034 (register R67/R68) and
  `u_T` / `c_γ² = μ_R/ρ_br` from Part III / stage003 / pathA_36 — NOT re-derived, NOT re-counted; the
  Part IV `±w` orientation `s_i`.
- **Governing:** `notes/ledger_v2_blueprint.md` §5 (reshape spec) + §6 (per-tooth ablation);
  `notes/part5_magnetism_atomic_split.md` (V-2 = the native source law + parity census; the
  tooth-allocation table); `docs/model_map.md` §3.5 (Q-CURRENT source-law bullet). Reshape directive +
  review trail: `research/pde_ledger_v2/_scratch/stage035_reshape_directive.md`. ⛔ **Not retained** — it lived in gitignored `_scratch/` and no copy survives; this line records that a directive existed, it is not an auditable citation.

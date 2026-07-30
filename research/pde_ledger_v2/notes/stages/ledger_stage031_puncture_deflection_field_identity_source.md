# ledger_stage031_puncture_deflection_field_identity_source

## Status

**Part IV — Charge. IV-2 (build-order 031; the EARNED-mechanism stage of the 4-stage Part-IV split,
user decision 2026-07-22).** Reshape of the puncture-deflection source build
`software/em_charge_attribute/puncture_deflection_electric_sign_result.md` (+ its `check.{py,wl}` /
`independent_recompute.py`), with SEVERAL genuine new builds (not a wrap): the full bare-mouth
reconstruction, the independent generic `I₊>0` certification, the exterior stationary derivation, and the
full completed-square response matrix. Surviving headline verdict token:

- **`THROAT_H_SOURCE_1_OVER_R2`** — the target-blind EARNED puncture-deflection mechanism: the field
  identity `ξ_w=ℓh`; the orientation-odd mouth source `η_i(k_m h − g_χh s_i)` with `Q_χ[r_Σ,s]=s` PROVED
  (odd kernel + the independently nonzero sleeve integral); the exterior stationary field `h(r)=h_A a/r`;
  the full completed-square response matrix `m`; and the target-blind neutral far-field shell
  `A=m_gg·C` with `U=s₁s₂A/(4πR)`, `F=s₁s₂A/(4πR²)`.

**EARNED** (within the postulated G0 closure): the field identity `ξ_w=ℓh` and the reduced coupling
`g_χh=J_m/ℓ`; the annulus normalization `∫η=1`; `f₀(0)=1/ℓ` evaluated (not assigned) from the actual
profile; the orientation projection `Q_χ=s` through the odd kernel AND the **generic `I₊>0`
reflection-dominance certification** (not `N_χ/N_χ`); the nonzero bare forcing `−g_χh s`; the exterior
`h(r)=h_A a/r` with positive holding curvature; the FULL response matrix `m` (`κ`, `Z`, `m_uu`, `m_ug`,
`m_gg`, `det m`, star witness) with `m_gg≥0`/`det m≥0` earned via `z_g²`; the self-response `S_gg` (defined
+ positive for positive `L_h`); and the units-restored `[L,T,M]` firewall over EVERY mechanism primitive.
The `1/R²` falloff and the `s₁s₂` product of the far-field shell are target-blind EARNED.

**POSTULATED** (honest, first-class): the bounded frozen `r_{Σ,+}∈(0,1]` representative with the `+w`
sleeve of frozen length `L_s` is a `[POSTULATE]` **completing G0's already-postulated `r_Σ` profile class**
(NO new knob class); strict `z_g>0` (and hence strict `m_gg>0`/`det m>0`) is a POSTULATED
Robin-admissibility witness (`z_g=1` at G0) unless independently derived.

**DEFERRED to stage 032** (`R1_REQUIRED(bc_selection)`): the class-dependent numerator `C` (which fixes the
force SIGN), the four BC-ensemble coefficients `A_V/A_J/A_M/A_MIXED`, the `internal_inconsistency=none`
(Q-AMEND) verification, the sealed 23040-cell landing table, and the `Q_E`/magnitude re-home. 031 emits
`m`, `m_gg`, `S_gg`, the named fact `NONZERO_HA_REQUIRES_CORE_HOLDER`, and the neutral shell `A=m_gg C`;
it does NOT land the R1 and picks NO sign.

Ledger-local earned label (NOT a source verdict token): `PUNCTURE_DEFLECTION_MECHANISM_TARGET_BLIND`.

Charge = the static ±w throat. This stage builds the MECHANISM that *dresses* the static substrate of
stage 030: a particle punctures the brane into `±w`, the bend is `ξ_w=ℓh`, `h` is the localized zero mode
assembled in 030, and the mouth functional sources it orientation-oddly so the charge SIGN is the ±w Z₂
orientation. 032 turns the neutral shell into the honest sign landing.

## Purpose

Stage 030 paid for the static prerequisites of the electric scalar `h` (the localized parent `H`, its
gapless zero mode `f₀`, the stable `(u_L,h)` kernel, the positive block `D`). It contained no force, no
sign, and no mouth source. Stage 031 is the puncture-deflection MECHANISM: it (i) identifies the geometric
field with the reduced scalar (`ξ_w=ℓh`), (ii) reconstructs — not assigns — the bare mouth source that
makes the charge sign the ±w orientation, (iii) solves the exterior stationary field, (iv) specifies and
protects the FULL completed-square response matrix `m`, and (v) assembles the target-blind neutral
far-field shell. It stops at the EARNED, target-blind structure: WHICH BC class holds, and hence the force
sign, is stage 032's `R1_REQUIRED(bc_selection)`.

Consumes (by citation, NOT re-derived): from stage 030 (`ELECTRIC_SCALAR_CLOSURE_STATIC`) the bundle
`{f₀, f₀(0)=1/ℓ, N₀=8/(3ℓ), h=P₀H, S_Lh, D, D*=7/4}` (engines print `CONSUMES_STAGE030=…`); from stages
003/030 the consumed blocks `{B_eff, C_hu}` (engines print `CONSUMES_STAGE003_030=…`).

## 1. The field identity `ξ_w=ℓh` and the reduced coupling (EARNED)

The normal embedding displacement of the punctured brane IS the reduced scalar:

```text
ξ_w(x) = ℓ h(x),     [ξ_w] = L,   [h] = 1.
```

The engines verify `ξ_w/ℓ − h = 0` for `ξ_w=ℓ·h` (`PASS_FIELD_IDENTITY`). Because `H = f₀ h + H_⊥` and
`f₀(0)=1/ℓ` (consumed from 030), the live parent coupling reduces on the mouth,
`−J_m Q_χ H → −g_χh Q_χ h` with

```text
g_χh = J_m / ℓ,     [g_χh] = E = M L² T⁻².
```

`g_χh` is a **decided** committed coupling (nonzero), NOT the sign R1.

## 2. The full bare-mouth-source reconstruction (EARNED — the centerpiece)

Every piece is reconstructed and able-to-fail; nothing is assigned.

**Annulus normalization (`PASS_ETA_ANNULUS_NORMALIZATION`).** The actual annular kernel is

```text
η_i = 3 · 1_A / [4π((a+ℓ)³ − a³)],     [η] = L⁻³,
∫_A η d³x = ∫_a^{a+ℓ} 4π ρ² η dρ = 1     (native radial Integrate).
```

**Mouth value evaluated (`PASS_F0_MOUTH_VALUE_EVALUATED`).** From the actual profile
`f₀ = 1/[ℓ cosh²(w/ℓ)]`, direct evaluation gives `f₀(0)=1/ℓ` (not assigned); the tooth also protects the
consumed `sech²` shape through the profile curvature at `w=0` (a wrong exponent shares the same value at
zero but not the same curvature).

**Reduced coupling nonzero (`PASS_REDUCED_COUPLING_NONZERO`).** `g_χh = J_m/ℓ ≠ 0`.

**Orientation projection `Q_χ=s` (two independent parts).**
*(a) Oddness ⇒ antisymmetry (`PASS_REFLECTION_ANTISYMMETRY`, `PASS_Q_CHI_ORIENTATION`).* With the odd
kernel

```text
o_ℓ(w) = w e^{−w²/(2ℓ²)} / (√(2π) ℓ³),     [o_ℓ] = L⁻²,   o_ℓ(−w) = −o_ℓ(w),
```

the even background `r_bg(w)` (`r_bg(−w)=r_bg(w)`), and the explicit body reflection
`r_{Σ,−}(w)=r_{Σ,+}(−w)`, one has `I_− = −I_+` for `I_± = ∫η ∫dw o_ℓ(w)[r_{Σ,±}² − r_bg²]`, so
`Q_χ = I_s/N_χ` gives `Q_+ = +1`, `Q_− = −1` as ratios.

*(b) The generic `I₊>0` reflection-dominance certification (`PASS_BOUNDED_PLUS_SLEEVE`,
`PASS_I_PLUS_REFLECTION_DOMINANCE`, + controls; option A, completing the G0 postulate).* Oddness alone
does NOT give `I₊≠0`. Supply a concrete **bounded, peak-normalized** signed-distance representative faithful
to G0's cylinder/cap geometry — the tanh slab with its `+w` edge extended by the FROZEN sleeve length
`L_s`, reaching `W_B/2+L_s` on `+w` and `W_B/2` on `−w`, mapped into `(0,1]`. In centered form

```text
r_{Σ,+}(w) = (1 + cosh(2·halfspan/ℓ)) / (cosh(2(w−center)/ℓ) + cosh(2·halfspan/ℓ)),
center = L_s/2,   halfspan = (W_B+L_s)/2,   peak r_{Σ,+}(center) = 1.
```

(Range certificate `cosh(2u/ℓ) − 1 = 2 sinh²(u/ℓ) ≥ 0` gives the `(0,1]` bound.) Define the
`η`-averaged sleeve difference `D(w) = ⟨r_{Σ,+}² − r_bg²⟩(w)`. The one-sided sleeve makes `D`
reflection-dominant. The reflected denominator difference is structural and sign-fixed by `w>0` and `L_s>0`:

```text
den(−w) − den(+w) = 2 sinh(2w/ℓ) sinh(L_s/ℓ) > 0     (for w>0, L_s>0),
D(w) − D(−w) = (r₊ − r₋)(r₊ + r₋) > 0                (both factors proved > 0),
```

so, pairing `w ↔ −w` with `o_ℓ(−w)=−o_ℓ(w)`,

```text
I₊ = ∫_ℝ o_ℓ(w) D(w) dw = ∫_{w>0} o_ℓ(w)[D(w) − D(−w)] dw > 0
```

parameter-generically (`o_ℓ(w)>0` and positive kernel mass `∫_0^∞ o_ℓ > 0`). This is STRUCTURAL to the
one-sidedness, holding over the whole admissible class — not a tuned profile. Hence `N_χ = I₊ ≠ 0` is
certified BEFORE any division (the anti-`N_χ/N_χ` ordering, `PASS_N_CHI_NONZERO_GUARD`), and `Q_χ=s`
follows from (a)+(b). The positive orientation SIGN of `N_χ` is a labelled CONVENTION; the certified
positive MAGNITUDE is the earned content.

*One-sidedness teeth (load-bearing, not the amplitude).* A reflection-symmetric deformation `D(w)=D(−w)`
→ integrand odd → `I₊=0` (`PASS_I_PLUS_EVEN_DEFORMATION_CONTROL`); `L_s=0` (bare symmetric slab ⇒ `D`
EVEN, `D(w)=D(−w)`, not `≡0`) → `I₊=0` (`PASS_I_PLUS_ZERO_SLEEVE_CONTROL`); the inadmissible `−w`-sleeve
control flips `I₊<0` (`PASS_I_PLUS_NEGATIVE_SLEEVE_CONTROL`). *STOP fallback:* if the generic certificate
fails symbolically over the family, the engines print `STAGE031_STOP: I_plus_not_generic` and exit 1 (the
retained anti-rescue guard).

**The live nonzero bare forcing (`PASS_BARE_FORCING_LIVE_VARIATION`).** Varying the reduced mouth
functional gives the `h`-Euler density `δΩ/δh = η_i(k_m h − g_χh s_i)` with `k_m = K_m/ℓ²`; the integrated
Euler amplitude is `k_m h − g_χh s`, whose bare (`h=0`) component is

```text
bare source = − g_χh s ≠ 0,
```

obtained from the live functional-to-Euler wiring (including `k_m`), not a hand-formed `−g_χh s`.

## 3. Exterior stationary solution + the named core-gap fact (EARNED)

Define the held mouth datum `h_A = ξ_w|_A/ℓ = P₀H|_A`, `[h_A]=1` (distinct from `u_L`). The exterior
stationary field (3D radial Euler `d/dr(r² dh/dr)=0`, decaying branch, `h(a)=h_A`) is

```text
h(r) = h_A a / r     ⇒ the h-mediated 1/R potential,
E_ext = 2π κ a h_A²,   κ > 0 the exterior h-gradient stiffness (a GENERIC positive constant
                        in the 031 exterior block — NOT identified with D/B_eff or the response κ=D/b),
holding curvature ∂²E_ext/∂h_A² = 4π κ a > 0
        (PASS_EXTERIOR_ONE_OVER_R, PASS_POSITIVE_HOLDING_CURVATURE).
```

**Named neutral fact (`PASS_NONZERO_HA_REQUIRES_CORE_HOLDER`).** A nonzero signed `h_A` is NOT a
stationary point of the core-less exterior energy — the only stationary point of `2π κ a (h_A − shift)²`
at `shift=0` is `h_A=0`. Export this as `NONZERO_HA_REQUIRES_CORE_HOLDER`: a guaranteed-nonzero amplitude
requires the deferred core holder (the sim-deferred nonlinear throat solve). 031 does NOT resolve it and
does NOT claim a guaranteed-nonzero interaction; stage 032 CONSUMES this named fact.

## 4. The full completed-square response matrix `m` (EARNED; consumes 030's `D`)

With `b=B_eff`, `c=C_hu`, `K_h` the reduced stiffness, and `D = bK_h − c²` (consumed from 030), the
response Robin scale is `κ = D/b` (distinct from the exterior stiffness of §3):

```text
κ = D/b,     Z = [[1, 0], [−(c/b) z_b, z_g]]     (typed; lower-left entry [·]=L),
m = Zᵀ · diag(b⁻¹, κ⁻¹) · Z,
m_uu = (D + c² z_b²)/(bD),   m_ug = − c z_b z_g / D,   m_gg = b z_g²/D,   det m = z_g²/D.
```

Predicates: `PASS_KAPPA_REDUCTION`, `PASS_Z_B_WIRING`, `PASS_RESPONSE_M_UU`, `PASS_RESPONSE_M_UG`,
`PASS_RESPONSE_M_GG`, `PASS_RESPONSE_SYMMETRY`, `PASS_RESPONSE_DETERMINANT`.

**Star witness (`PASS_RESPONSE_STAR_WITNESS`).** At `b=2, K_h=1, c=1/2, z_b=z_g=1`, `D*=7/4` (consumed
from 030) and

```text
m* = (1/7) [[4, −2], [−2, 8]],
m_uu = K_h/D = 4/7  (NOT 1/b = 1/2; m_uu=1/b holds ONLY when c z_b = 0),
m_ug = −c/D = −2/7,   m_gg = b/D = 8/7.
```

**Positivity honesty.** `m_gg ≥ 0` and `det m ≥ 0` are EARNED via `z_g²` given `B_eff>0`, `D>0`
(`PASS_M_GG_NONNEGATIVE`). Strict `z_g>0` (and hence strict `m_gg>0`/`det m>0`) is a **POSTULATED**
Robin-admissibility witness (`0<z_g≤1`, `z_g=1` at G0; `PASS_Z_G_POSTULATED_WITNESS`) — the source only
sets a positive symbol, and the 031 engine carries `z_g` as an opaque escape factor (it does NOT verify a
`z_g=1−k_m S_gg` reduction).

**Self-response (`PASS_S_GG_SELF_RESPONSE`).** `S_gg = ⟨η, L_h⁻¹ η⟩` is given an explicit Galerkin
quadratic-form definition; for positive `L_h` the form is positive (DERIVED for positive `L_h`), `[S_gg]=E⁻¹`.
It enters 032's Dirichlet coefficient.

## 5. The target-blind far-field FORM (EARNED neutral shell, no sign)

From items 2–4 the two-throat interaction has the FORM (`PASS_NEUTRAL_FAR_FIELD_FORM`)

```text
U = s₁ s₂ A / (4πR),     F = s₁ s₂ A / (4πR²),     A = m_gg · C,
[A] = E L,   [C] = E²   (C an unspecified real).
```

The `1/R²` falloff, the `s₁s₂` product (from the bilinear response + the orientation-odd source), and
`[A]=EL` are **target-blind EARNED**. The class-dependent `[E²]` numerator `C` is NOT specified and NO sign
is picked — the earned result is the **allowed/conditional leading form, possibly with zero coefficient**
(`PASS_ALLOWED_ZERO_COEFFICIENT`: `U|_{C=0}=0` is admissible). A nonzero stationary amplitude requires the
deferred core holder (`NONZERO_HA_REQUIRES_CORE_HOLDER`). The four BC-class numerators, any `A_V/A_J/A_M`
formula, and the `internal_inconsistency=none` claim are NOT written here — they are stage 032.

## 6. Dimensional firewall (EARNED — units-restored `[L,T,M]`)

Every mechanism primitive + reduction + the far-field shell is exponent-triple homogeneous, units restored,
able-to-fail (21 `PASS_UNITS_*` predicates + live reduction chains):

```text
[ξ_w]=L, [h]=[h_A]=1, [K_m]=M L⁴T⁻², [J_m]=M L³T⁻²,
[k_m]=[g_χh]=E=M L²T⁻², [η]=L⁻³, [o_ℓ]=L⁻², [N_χ]=L⁻¹,
[m_uu]=L³/E=M⁻¹LT², [m_ug]=L²/E=M⁻¹T², [m_gg]=L/E=M⁻¹L⁻¹T²,
[κ]=E/L=M LT⁻², [Z_21]=L, [det m]=M⁻²T⁴, [S_gg]=E⁻¹=M⁻¹L⁻²T²,
[C]=E², [A]=EL=M L³T⁻², [U]=E, [F]=E/L.
```

Live reduction chains checked: `k_m=K_m/ℓ²`; `g_χh=J_m/ℓ`; `N_χ=∫η dx · ∫o_ℓ dw`; `A=m_gg C`; `U=A/R`;
`F=A/R²`.

## Scope, deferrals, and first-class departures (recorded honestly)

1. **EARNED within the postulated G0 closure.** What is earned (items 1–6) is target-blind and holds
   *given* the postulated G0 action + `r_Σ` profile class. The whole G0-v0 closure remains a labeled
   postulate (per stage 030).
2. **The bounded `r_{Σ,+}` representative is a `[POSTULATE]`** completing G0's already-postulated `r_Σ`
   (the `+w`-extended peak-normalized signed-distance form; frozen sleeve `L_s`) — NO new knob class beyond
   G0. The generic `I₊>0` certification is structural to the one-sidedness, not to the profile magnitude.
3. **Strict `z_g>0` is POSTULATED** (Robin-admissibility witness); `m_gg≥0`/`det m≥0` are earned via `z_g²`.
4. **The SIGN R1 + magnitude are stage 032 / deferred.** The class numerator `C`, the BC-ensemble
   coefficients, `internal_inconsistency`, and the R1 landing are stage 032; the `Q_E`/magnitude re-home
   (with `c_a`,`c_ξ` unbounded at tier-A) is 032's; a guaranteed-nonzero amplitude needs the sim-deferred
   core holder (`NONZERO_HA_REQUIRES_CORE_HOLDER`, R61).

## Joint 031↔032 ownership (summary)

| Block | 031 (this stage) | 032 |
|---|---|---|
| Field identity `ξ_w=ℓh`, `g_χh=J_m/ℓ` | ✅ owns | — |
| Bare-mouth source (`∫η=1`, `f₀(0)=1/ℓ`, `Q_χ=s` via odd kernel + generic `I₊>0`, `−g_χh s`) | ✅ owns | — |
| Exterior `h(r)=h_A a/r`; named fact `NONZERO_HA_REQUIRES_CORE_HOLDER` | ✅ owns / exports | consumes the named fact |
| Response `m` (`κ`, `Z`, `m_uu/m_ug/m_gg`, `det m`, `z_g/z_b`, star); `S_gg` (defined) | ✅ owns | consumes `m`, `m_gg`, `S_gg` |
| Neutral shell `A=m_gg·C` (`1/R²`, `s₁s₂`, `[A]=EL`) | ✅ owns (neutral) | consumes the shell |
| Class numerator `C`; `A_V/A_J/A_M/A_MIXED`; force sign; MIXED admissibility | — | ✅ owns (`R1_REQUIRED(bc_selection)`) |
| Sealed 23040-cell landing table; `internal_inconsistency=none` (Q-AMEND); `Q_E`/magnitude re-home | — | ✅ owns |
| Units/dimension checks | mechanism primitives + shell (`[L,T,M]`) | ensemble/landing dims |

## Source-to-stage predicate manifest (all 36 atomic teeth + Q-block/App-A–E claims)

Completeness certificate: **no silently-dropped source claim** — every atomic tooth of the source engine
`software/em_charge_attribute/puncture_deflection_electric_sign_check.py` (the 36-entry `TOOTH_ORDER` at
`:1082`, predicates `:1203–1248`) + the Q-block/App-A–E claims lands in 031, 032, or is explicitly
scoped-out with a defensible reason. The two partitions are disjoint and exhaustive; the 031 engine prints
`DEFERRED_TO_STAGE032: BC-class numerator C, force sign, ensembles, internal_inconsistency, R1 landing`.

### The 36 atomic teeth (`TOOTH_ORDER`)

| Source tooth | Description | Disposition |
|---|---|---|
| `FIELD_IDENTITY_UNITS` | ξ_w=ℓh, [ξ_w]=L | PRESERVED-031 (`PASS_FIELD_IDENTITY`, `PASS_UNITS_XI_W`) |
| `FIELD_PARENT_MAP` | ℓ·f₀(0)=1; N₀=8/(3ℓ), h=P₀H | REPLACED-031 (`PASS_F0_MOUTH_VALUE_EVALUATED` adds sech² curvature; N₀/P₀ consumed from 030) |
| `FIELD_LIVE_QH` | −J_m Q_χ H → −g_χh Q_χ h | REPLACED-031 (full bare-mouth reconstruction: `PASS_BARE_FORCING_LIVE_VARIATION`, `PASS_REDUCED_COUPLING_NONZERO`, `Q_χ=s` chain) |
| `ACTION_TRANSCRIPTION` | S_Lh action + Hessian [[b,c],[c,k]] | physics PRESERVED-031 (D=bk−c², `PASS_RESPONSE_STAR_WITNESS` D*=7/4); the source-**document** string-faithfulness harness SCOPED-OUT (build-time source verification; action is a cited stage-030 input) |
| `AMEND_REPLACE` | REPLACE re-types h-source BC, new_rows=0 | OWNED-032 (Q-AMEND) |
| `AMEND_ADD` | ADD keeps source + one core h-holding row | OWNED-032 (Q-AMEND) |
| `SHOLD_SCOPE` | S_hold constrains r_B−½ only | OWNED-032 (Q-AMEND consistency) |
| `ZERO_LEDGER` | 13 preserved [POSTULATE] zero rows | OWNED-032 |
| `MATRIX_DETERMINANT` | det m = z_g²/D | PRESERVED-031 (`PASS_RESPONSE_DETERMINANT`) |
| `BC_ACTUAL_GAP` | classifier → UNDETERMINED_ANALYTICALLY | OWNED-032 (consumes 031's `PASS_NONZERO_HA_REQUIRES_CORE_HOLDER` + exterior evidence) |
| `BC_FREE_CONTROL` | free mouth → FIXED_SOURCE | OWNED-032 (Q-BC 1/4) |
| `BC_VALUE_CONTROL` | imposed value → DIRICHLET_VALUE | OWNED-032 (Q-BC 2/4) |
| `BC_MONOPOLE_CONTROL` | imposed conormal → FIXED_MONOPOLE | OWNED-032 (Q-BC 3/4) |
| `BC_MIXED_CONTROL` | imposed relation → MIXED | OWNED-032 (Q-BC 4/4) |
| `FORCE_V_FUNCTIONAL` | A_V=m_gg φ²/S_gg² | OWNED-032 (consumes 031's m_gg + `PASS_S_GG_SELF_RESPONSE`) |
| `FORCE_M_HWORK` | A_M=m_gg(q²−g²) | OWNED-032 |
| `FORCE_J_FUNCTIONAL` | A_J=−m_gg(j+g)² | OWNED-032 |
| `FORCE_MIXED_FUNCTIONAL` | A_MIXED=m_gg[(1−2λ)q²−2λqg−g²] | OWNED-032 |
| `MIXED_FULL_RANGE` | endpoints A_M, −m_gg(q+g)² | OWNED-032 (MIXED three-regime) |
| `FALLOFF` | force power = 1/R² | PRESERVED-031 (`PASS_NEUTRAL_FAR_FIELD_FORM`, `PASS_EXTERIOR_ONE_OVER_R`) |
| `UNITS_RESTORED` | [A]=EL, [U]=E, [F]=EL⁻¹ on the A_X | OWNED-032 (ensemble units); neutral-shell [A]/[U]/[F] PRESERVED-031 (`PASS_UNITS_SHELL_*`) |
| `COMBINE_REPLACE` | REPLACE totals per class | OWNED-032 |
| `COMBINE_ADD` | ADD totals per class | OWNED-032 |
| `NO_DOUBLE_COUNT` | no double-counted held-h work | OWNED-032 |
| `RANGE_SIGN_FLIP` | sign-flip → outcome_not_invariant | OWNED-032 (range control) |
| `RANGE_TOUCH_ZERO` | zero-touch-without-flip control | OWNED-032 (range control) |
| `RANGE_SUBDOMINANT` | derived subdominant constant control | OWNED-032 (range control) |
| `MAG_FREE_FACTOR` | magnitude_free_factor (c_a, c_ξ unbounded) | OWNED-032 (Q-MAG / Q_E magnitude; co-blocker R1_REQUIRED(magnitude)) |
| `DENSITY_HOOK` | density NO / no_local_prediction | OWNED-032 (Q-MAG 1/4) |
| `MONOPOLE_HOOK` | radial_monopole UNDETERMINED | OWNED-032 (Q-MAG 2/4) |
| `MODULUS_HOOK` | universal_quantization NO | OWNED-032 (Q-MAG 3/4) |
| `VERDICT_TOTALITY` | truth-table exact over 23040 cells + digest | OWNED-032 (sealed landing) |
| `VERDICT_PRECEDENCE` | class gap fires before variant/tier/mag | OWNED-032 |
| `LANDING_OWNERSHIP` | landing recomputed from typed upstream, injection rejected | OWNED-032 (directive item 3) |
| `TARGET_BLINDNESS` | neutral labels; target only in §4 adjudicator | OWNED-032 |
| `DUAL_ENGINE_TERMS` | SymPy vs Wolfram payload agree | split: landing/ensemble terms OWNED-032; the mechanism terms (`D`, `κ`, `m_uu/ug/gg`, `det m`) PRESERVED-031 **conditional on OPAQUE escape factors** (031's own dual-engine .wl); the source escape-factor reductions (`z_g=ζ`, `z_b=1−k_m·reb`) are SCOPED-OUT/deferred — 031 carries `z_g`,`z_b` opaque, verifying no formula for them |

### Q-block / App-A–E claims not 1:1 with a `TOOTH_ORDER` entry

| Source claim | Description | Disposition |
|---|---|---|
| `variant_unresolved` | realization fact feeding §4 landing | OWNED-032 (landing input) |
| Q-MAG close_range hook | UNDETERMINED / out-of-scope R≈r_e | OWNED-032 (Q-MAG 4/4) |
| S_gg self-response (App-B) | self coefficient in A_V | PRESERVED-031 (`PASS_S_GG_SELF_RESPONSE`, DERIVED S_gg>0) |
| G0 witness m=⅐[[4,−2],[−2,8]] (App-A) | numeric response witness | PRESERVED-031 (`PASS_RESPONSE_STAR_WITNESS`) |
| m_uu, m_ug, m_gg, κ, z_g, z_b (App-A) | completed-square entries | **PARTIAL** — the response-matrix WIRING + star values PRESERVED-031 (`PASS_RESPONSE_M_*`, `PASS_KAPPA_REDUCTION`, `PASS_Z_B_WIRING`, `PASS_Z_G_POSTULATED_WITNESS`, `PASS_M_GG_NONNEGATIVE`); the source `z_g`/`z_b` **formulas** (`ζ` / `1−k_m·reb`) are NOT verified by 031 — carried as OPAQUE escape factors, outside the verified 031 result |
| Exterior E_ext=2πκa h_A², curvature 4πκa>0 (App-C) | exterior energy + holding curvature | **CORRECTED-031** — the energy form + positive holding curvature PRESERVED with a GENERIC exterior stiffness `κ>0`; 031 does NOT adopt the source/directive identification `curvature=4π(D/B_eff)a` (REJECTED — the 031 exterior block uses a fresh generic positive `κ`, not `D/B_eff`) |
| 23040-cell digest `7627417a…ac49` + counts (App-D) | exhaustive landing digest/totals | OWNED-032 (committed literal, directive item 3) |
| App-E production checks (both exit 0; 36 teeth; units; range) | build-level acceptance envelope | split: landing acceptance OWNED-032; mechanism half 031 |

## Verification

- **Dual-engine, both exit 0, genuinely independent routes.**
  `scripts/ledger_stage031_puncture_deflection_field_identity_source_sympy_audit.py` — **SymPy 50 PASS**.
  `mathematica/ledger_stage031_puncture_deflection_field_identity_source_mathematica_audit.wl` —
  **Mathematica 50 PASS**. Standalone, print-only, assert-zero (`raise SystemExit(1)` / `Exit[1]`), no
  argparse harness, no JSON/YAML payload, zero file-I/O between engines; float-free / machine-real-free
  payload throughout. Verdict `THROAT_H_SOURCE_1_OVER_R2`, earned label
  `PUNCTURE_DEFLECTION_MECHANISM_TARGET_BLIND`.
- **The `.wl` is a genuinely independent route.** It reorders the derivation (native response
  matrix/`Inverse` first, then `DSolve` for the exterior branch, then the native `Integrate` mouth route
  last) and uses native Wolfram machinery where the `.py` uses SymPy: `Integrate` for the annulus/η
  normalization, the certified `I₊` sleeve integral and `N_χ`; `DSolve` for `h(r)=h_A a/r`; native matrix
  construction/`Inverse` for `m`.
- **50 per-tooth ablations, each `FIRED_AT_OWN_ASSERT`** (env switch `LEDGER_STAGE031_MUTATION`): a
  primitive mutation for every folded predicate, including the SEPARATED norm/normalization teeth
  (`∫η` vs the `N_χ=0` pre-division guard vs the generic `I₊>0` certification), the one-sidedness
  load-bearing teeth (even-deformation, `L_s=0`, `−w`-sleeve), the per-entry response teeth (`m_uu`,
  `m_ug`, `m_gg`, `z_b` wiring, symmetry, `det m`, the `z_g` witness, the `m_gg≥0`-via-`z_g²` fact), the
  exterior teeth (`1/r`, holding curvature, the exported named fact), the neutral-shell teeth (`1/R²`+`s₁s₂`,
  `[A]=EL`, the allowed-zero-coefficient control), and the complete units firewall.
- **Tri-review CLEAN:** fidelity + adversarial fresh-agent review both CLEAN (50/50 mutations fire at own
  assert in both engines; the `I₊>0` certification a parameter-generic symbolic proof, not single-profile;
  the `N_χ/N_χ` trap defeated; no vacuous/subsumed tooth; zero remediation).

## Downstream consumers

- **`ledger_stage032` (IV-3, ensembles + R1 landing):** consumes 031's `m`, `m_gg`, `S_gg`, the neutral
  shell `A=m_gg·C`, and the named fact `NONZERO_HA_REQUIRES_CORE_HOLDER`; selects the class numerator `C`
  per BC class (`A_V/A_J/A_M/A_MIXED`), runs the sealed 23040-cell landing, computes
  `internal_inconsistency=none`, re-homes `Q_E`/magnitude, and lands `R1_REQUIRED(bc_selection)`.
- **Parameter register:** rows + edges R54–R61 added; the `A_X`/class-`C`/ensemble/`internal_inconsistency`/
  `Q_E` rows are 032's, and the plan-doc register-preview (`part4_charge_atomic_split.md`) is re-tagged
  accordingly.
- **Stage 033 (IV-4):** `NATIVE_P_NO_EMERGENT_GAUSS` (independent).

## Provenance

- **Physics source:** `software/em_charge_attribute/puncture_deflection_electric_sign_result.md`
  (field identity, mouth source, exterior solution, response Q-blocks, App A–E) + its
  `puncture_deflection_electric_sign_check.{py,wl}` + `..._independent_recompute.py`; the mouth functional
  normative provenance is G0 card `g0_closure_card_v0.md` §§1.1, 3.1–3.2, the `f₀` profile §2.3, the
  `r_Σ`/`r_bg` profiles §2.3 (`:47,49,55,57`). The source engine *assigns* `f₀(0)=1/ℓ` and sets a positive
  `zeta`; 031 reconstructs these as new able-to-fail builds.
- **Consumes:** stage 030 `{f₀, f₀(0)=1/ℓ, N₀=8/(3ℓ), h=P₀H, S_Lh, D, D*=7/4}` + stages 003/030
  `{B_eff, C_hu}` (cite, no re-count).
- **Governing:** `notes/ledger_v2_blueprint.md` §5/§6; `notes/part4_charge_atomic_split.md` (IV-2 row +
  register preview); `docs/model_map.md` §3.4. Reshape directive + review trail:
  `research/pde_ledger_v2/_scratch/stage031_reshape_directive.md`. ⛔ **Not retained** — it lived in gitignored `_scratch/` and no copy survives; this line records that a directive existed, it is not an auditable citation.

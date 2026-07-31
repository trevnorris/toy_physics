# ledger_stage005_sound_speed_light_ratio

> ⛔⛔ **WHOLLY HISTORICAL — DO NOT CONSUME (2026-07-31).** This stage note carries two claims that HEAD has
> retired, and they are woven through the body, not confined to a header:
> 1. the pin relation `a = ħ/(m_GNLS c_s0)` — **deleted** (`407eed94`);
> 2. `λγ` as a single underived calibration input — HEAD distinguishes **three** things: `c_γ` free/debt
>    until v3 **S9**, the **ratio** `λγ = c_γ/c_s` **derived** at v3 **S20a**, and only **`λγ = 1`**
>    calibrated/uncommitted.
>
> ⇒ **v3's S2 derives the sound speed FORWARD from S1/S1.5 — it does NOT re-bank this artifact.** That is
> also what v3 acceptance criterion (3) requires: the chain must be followable without consulting v2.
> ⛔ Cite this only for provenance of the historical result, never as a source of classification.

## Status

EARNED (sound speed + wave-sector ceiling, within the imposed EOS closure) with a
CALIBRATED honest landing -- `C_GAMMA_RATIO_UNDERDETERMINED`,
`PASS_WITH_NAMED_RESIDUALS`.

The sound-speed law `c_s^2 = 5 K rho^4 / m_GNLS` is DERIVED from the imposed EOS
`P = K rho^5`; the terminal ceiling `c = c_gamma` is DERIVED from the wave sector
(not from `E = m c^2`). The **headline** is the honest CALIBRATED landing: the
light/sound ratio `lambda_gamma = c_gamma/c_s` is **UNPINNED by the parent
action** -- a **free calibration input** -- and the tail `(c/c_s)^3 =
lambda_gamma^3` is carried symbolic (NOT set to 1). This is not a FAIL; it is a
genuinely-free knob, proven free by a negative control that is reversible to
`C_GAMMA_EQUALS_C_S` only when a source equation is inserted (which the parent
action does not supply).

This is the second stage of **Part I (The medium)**. It **consumes**
`ledger_stage004` (I-1) -- the `{L,T,M}` dictionary and the `EOS_FROM_GNLS_FACTOR`
handoff (`h0`, `xi_h`) -- as cited inputs, and OWNS the sound-speed derivation
proper that I-1 deferred here.

## Purpose

Fix the medium's velocity structure: derive the phonon sound speed from the EOS,
pin the three velocity scales `{v_b, c_s, c_gamma}`, derive the `c = c_gamma`
wave-sector ceiling (with the bound-mode clock giving time dilation), and land the
light/sound ratio `lambda_gamma = c_gamma/c_s` as a free calibration input
(`C_GAMMA_RATIO_UNDERDETERMINED`). The `pathA_20b` coupled linearization is folded
in as reference support: it sharpens the ratio into a two-layer bulk/brane verdict
and supplies the able-to-fail negative control.

## Provenance

Earned source (provenance, not content):
`software/stage1_solver/reports/pathA_20_velocity_constants.md` (verdicts
`C_GAMMA_RATIO_UNDERDETERMINED`, `STATIONARY_PROFILE_UNDERDETERMINED_BY_BRANCH_DATA`,
`HBAR_PROVENANCE_UNDETERMINED`; `PASS_WITH_NAMED_RESIDUALS`) and
`software/stage1_solver/reports/pathA_20b_cgamma_cs_linearization.md`
(`C_GAMMA_BULK_UNDERDETERMINED` / `C_GAMMA_RATIO_STILL_UNDERDETERMINED`; the coupled
principal symbol and the `FORCED_EQUALITY_REJECTED_WITHOUT_SOURCE_EQUATION`
negative control). The SymPy logic reshaped here was extracted from the shared
harness `software/stage1_solver/src/stage1_solver/dimensional_check.py`
(`run_patha20_velocity_constants`, `--patha20-velocity`, 21 dim + 5 alg checks;
`run_patha20b_cgamma_cs`, `--patha20b-cgamma-cs`, 11 dim + 7 alg checks; already
passing). The derivation below is inlined so a reader never needs to open those
reports or the harness.

Script-backed status: the exact-symbolic claims are checked by
`scripts/ledger_stage005_sound_speed_light_ratio_sympy_audit.py` and independently
by `mathematica/ledger_stage005_sound_speed_light_ratio_mathematica_audit.wl`. The
`.wl` is a genuinely independent route (native `D`/`Det`/`Factor`/`Solve`/`Reduce`
and `UnitDimensions` construction), not a transliteration of the `.py`. Both
engines derive -- do not hardcode -- the sound-speed law, the state slope, the
dispersion determinants, the ratio, and the negative-control verdict, and both
carry genuine able-to-fail teeth (a reversibility flip, dimensional/derivation
firewalls, and a consuming-handoff citation-integrity break). Coverage counts of
the harness surface (21/5/11/7) are asserted with genuine live counters, not
`N-N` placeholders.

## 0. Why needed

The medium supports three distinct speeds, and the whole program's falsifiability
turns on their ratios. `c_s` is the phonon/gravity-change speed; `c_gamma` is the
light-cone speed; `v_b` is the condensate flow. Downstream sectors need (i) the
sound speed as an explicit function of state, and (ii) an honest accounting of
whether the light/sound ratio `lambda_gamma` is fixed by the action or is a free
input. `lambda_gamma` enters the gravity-sector normalization to the fifth power
(`(c/c_s)^3 = lambda_gamma^3`, and further factors downstream), so mis-stating it
as `1` would silently smuggle a calibration. This stage settles both questions.

## 1. The sound speed from the EOS (owned derivation)

The medium carries the imposed stiff-polytropic EOS from the GNLS potential
`U(rho) = (K/4) rho^5`, so `P = K rho^5` (`EOS_CLOSURE_IMPOSED` -- postulated, not
derived from a deeper substrate here). The barotropic sound speed is

```text
c_s^2(rho) = (1/m_GNLS) dP/drho = (1/m_GNLS) d/drho (K rho^5) = 5 K rho^4 / m_GNLS,
[c_s] = L T^-1.
```

The factor `5` FALLS OUT of the derivative (it is not written by hand). The state
dependence is derived by the log-slope: since `c_s(rho) proportional to rho^2`,

```text
d ln c_s / d ln rho = rho * (d c_s/drho)/c_s = 2.
```

This law is EARNED **relative to the imposed EOS** (an `hbar`-free microscopic EOS
derivation is outside this stage). The stationary quantum-Bernoulli balance
(`(1/2) m_GNLS v_b^2 + h(rho) + V_conf + Q`) and the bulk continuity equation
`d_t rho + div(rho v_b) = 0` are dimensionally homogeneous; no stationary profile
is solved by dimensions alone.

## 2. Three velocity scales and the `c = c_gamma` ceiling

The three velocity scales, dimensionally pinned (`[.] = L T^-1` each):

```text
v_b = (hbar/m_GNLS) grad(theta)   condensate flow (theta dimensionless; pin consumed from I-1)
c_s                               phonon/density (gravity-change) speed, c_s(rho) with asymptote c_s0
c_gamma                           photon/gauge-wave (light-cone) speed; value relative to c_s is dynamical
```

**The `c = c_gamma` ceiling is derived from the wave sector, non-circularly.** The
massless dispersion `omega^2 = c_gamma^2 k^2` gives group speed `c_gamma`. A
trapped transverse mode has `omega^2 = c_gamma^2 (k_par^2 + k_perp^2)`, so

```text
d omega/dk = c_gamma * k / sqrt(k^2 + k_perp^2),
c_gamma^2 - (d omega/dk)^2 = c_gamma^2 k_perp^2 / (k^2 + k_perp^2)  >= 0,
```

strictly positive for the trapped mode (`k_perp != 0`): the group velocity is
bounded by `c_gamma` and approaches it at high drive. The degenerate `k_perp = 0`
case is the free wave at exactly `c_gamma` (the ceiling itself), consistent.

**Bound-mode clock (time dilation).** The trapped-mode rest oscillation is
`omega0 = c_gamma k_perp`. A boosted wave-operator solution has phase
`exp[-i omega0 gamma (t - v x/c_gamma^2)]`; along the packet center `x = v t` the
internal clock advances at `omega0/gamma` (with `gamma = 1/sqrt(1 - v^2/c_gamma^2)`).
**No `E = m_defect c_gamma^2` or Compton premise is used** -- time dilation is a
wave-kinematic result.

## 3. The flux law (carried undetermined)

The candidate sonic number fluxes are dimensionally consistent
(`[rho_* c_{s,*} A_{3,*}] = [rho_{3,*} c_{s,*} A_{2,*}] = T^-1`, and
`[P0] = [K rho0^5] = [P]`), but no `J_crit` law is accepted here. A conditional
ideal Euler-nozzle algebra (upstream `v0 = 0`, Bernoulli
`(1/2) v_*^2 + c_{s,*}^2/4 = c_s0^2/4`, `c_s proportional to rho^2`) gives

```text
c_{s,*}/c_s0 = 3^(-1/2),   rho_*/rho0 = 3^(-1/4),   Jcrit/(rho0 c_s0 A_*) = 3^(-3/4),
```

but these are **recorded, not accepted** (`CONDITIONAL_NOT_ACCEPTED_AS_BRANCH_LAW`).
Verdict `STATIONARY_PROFILE_UNDERDETERMINED_BY_BRANCH_DATA`; the actual stationary
profile and DtN data remain branch-dependent (pathA_21 must consume this verdict,
not an unconditional choked flux). `NO_NET_ACCRETION_BC_UNDERIVED`.

## 4. Role-catalog dimensions and the consumed pin/healing

Owned dimensional checks: circulation `kappa = closed-integral v_b.dl = h n/m_GNLS`
(`[kappa] = L^2 T^-1`), phase momentum `p = hbar grad(theta)` (`[p] = M L T^-1`),
quantum pressure `Q = -(hbar^2/2 m_GNLS) laplacian(sqrt rho)/sqrt rho`
(`[Q] = ENERGY`), and both mass-bridge candidates (Section 6).

Consumed from `ledger_stage004` (I-1), **cited, not re-derived**: the `{L,T,M}`
dictionary and the pin/healing handoff `a = hbar/(m_GNLS c_s0)`,
`xi_h = sqrt(2) hbar/(m_GNLS c_s0)`, `h0 = (m_GNLS c_s0^2)/4`
(`EOS_FROM_GNLS_FACTOR`). The stage asserts the dimensional consistency of these
(`[m_GNLS c_s0 a] = [hbar]`, `[m_GNLS c_s0 xi_h] = [hbar]`) AND an **exact-value
citation-integrity** cross-check `h0 - (m_GNLS c_s0^2)/4 = 0`,
`xi_h - sqrt(2) hbar/(m_GNLS c_s0) = 0` -- expressed via the sound speed **derived
in Section 1** -- tying this stage's `c_s0` to I-1's exported handoff. It does NOT
re-derive them from the core balance (that is I-1's owned content). The
anti-tautology caveat is carried: `hbar = m_GNLS c_s0 a` is a pin rearrangement,
not an `hbar`-emergence proof (`HBAR_PROVENANCE_UNDETERMINED`).

## 5. The coupled principal symbol and both dispersions (pathA_20b support)

On the neutralized homogeneous background (`psi0 = sqrt(rho0) e^{-i mu t/hbar}`,
`v_b0 = 0`, `A_M0 = 0`, with an explicit neutralizing external source so the
background Maxwell equation is `0 = 0`,
`LEGAL_WITH_EXPLICIT_NEUTRALIZING_EXTERNAL_SOURCE`), the coupled principal symbol
block-diagonalizes:

```text
Madelung phonon block  M = [[omega, -(rho0 hbar/m_GNLS) k^2], [-h', hbar omega]],  h' = 5 K rho0^3,
det M = hbar (omega^2 - c_s0^2 k^2),   c_s0^2 = rho0 h'/m_GNLS = 5 K rho0^4/m_GNLS,
gauge transverse       P_T = C_E omega^2 - C_B k^2 = C_E (omega^2 - c_bulk^2 k^2),  c_bulk^2 = C_B/C_E,
coupled                det P = det M * P_T^2   (two physical transverse polarizations).
```

`c_s0^2` FALLS OUT of the Madelung determinant (via `h'`), reproducing Section 1
on the homogeneous background. The off-diagonal GNLS-gauge principal terms
`VANISH` on the neutralized background (the linearized current / London /
source-coupling terms `J_psi0^0 = q_* rho0`, `(q_*/m_GNLS) rho0 dA_i`,
`A_M dJ^M` are all LOWER-ORDER than the second-derivative principal operator, so
they do **not** set the cone). Two dispersions read off: phonon
`omega^2 = c_s0^2 k^2` (the Bogoliubov `k^4` quantum-pressure correction is
dispersive and does not set the cone), gauge `omega^2 = c_bulk^2 k^2`, with the
transverse branch `BULK_PRINCIPAL_TRANSVERSE_BRANCH_ESTABLISHED`. That
`[c_s] = [c_gamma] = L T^-1` is explicitly **NON-EVIDENTIARY for equality**.

## 6. The CALIBRATED landing: `lambda_gamma` is a free calibration input

`lambda_gamma = c_gamma/c_s` is introduced as a **labeled free calibration input**
(dimensionless), with the tail `(c/c_s)^3 = lambda_gamma^3` carried symbolic and
**NOT reduced to 1**. The parent action does not pin it; `c_gamma = c_s` from
shared dimensions or from legacy weak-field prose is **REJECTED**. The `pathA_20b`
sharpening is two-layer:

- **bulk** `C_GAMMA_BULK_UNDERDETERMINED` (`BULK_METRIC_SPEED_NORMALIZATION_UNSPECIFIED`);
  conditional `c_bulk/c_s0 = sqrt((C_B/C_E) m_GNLS/(5 K rho0^4))`, carried symbolic;
- **brane** `C_GAMMA_RATIO_STILL_UNDERDETERMINED` (`BRANE_ZERO_MODE_REDUCTION_UNDERIVED`,
  sub-residual `BRANE_PHOTON_CONE_REQUIRES_PROFILE`).

If `c_bulk` is `rho0`-independent then `lambda_gamma proportional to rho0^-2`
(derived: `d ln lambda_gamma / d ln rho0 = -2`, since `c_s0 proportional to rho0^2`).

The mass-bridge candidate `m_defect = alpha_J hbar J/c_gamma^2` (or
`alpha_J h J_nu/c_gamma^2` for a cycle rate) is a **dimensional candidate only**
(`[hbar J/c_gamma^2] = M`): it does not collapse `M`, `alpha_J` is not derived,
and it is deferred to pathA_21. `h` and `hbar` share dimensions; the `2 pi`
placement is not decided by dimensions.

## 7. The negative control and reversibility (able-to-fail)

The `C_GAMMA_RATIO_UNDERDETERMINED` landing is **computed**, not asserted. With
`C_E, C_B, K, rho0, m_GNLS` independent symbols,

```text
equality_residual = C_B/C_E - 5 K rho0^4/m_GNLS   (not identically zero: the gauge
                    metric and the acoustic metric are independent),
forced_equals_valid = source_metric_equation_present AND (equality_residual == 0).
```

With `source_metric_equation_present = False` (no such equation is found in the
parent action), `forced_equals_valid = False` -> `FORCED_EQUALITY_REJECTED_WITHOUT_SOURCE_EQUATION`
and the headline `C_GAMMA_RATIO_UNDERDETERMINED`. **Reversibility (able-to-PASS):**
registering a genuine source equation (`source_metric_equation_present = True` AND
substituting `C_B -> 5 K rho0^4 C_E/m_GNLS`) drives `equality_residual -> 0` and
flips the verdict to `C_GAMMA_EQUALS_C_S` (`lambda_gamma = 1`). So the landing is
genuinely reversible -- not rigged toward `UNDERDETERMINED` -- and `EQUALS` is
reachable **only** via an actually-inserted source equation, never from a bare
dimensional match. (The reversibility branch is a labeled counterfactual; it does
NOT close the carried `C_GAMMA_RATIO_UNDERDETERMINED`.)

The dimensional/derivation firewalls all fire: corrupting the EOS exponent
(`P = K rho^4`) breaks both `[c_s]` and the log-slope; dropping `m_GNLS` breaks
`[c_s]`; corrupting the Madelung off-diagonal breaks the determinant factorization.
The consuming-handoff teeth fire on a dimensionless corruption: dropping the
`sqrt(2)` in `xi_h`, or the `1/4` in `h0`, breaks the exact-value
citation-integrity check.

## 8. Carried residuals (first-class, not repaired)

- `EOS_CLOSURE_IMPOSED`: `c_s` is earned only relative to the imposed
  stiff-polytropic EOS.
- `C_GAMMA_RATIO_UNDERDETERMINED` (`BLOCKS_NUMERIC_C_GAMMA_OVER_C_S`): carry
  `lambda_gamma` and the `lambda_gamma^3` tail to pathA_21/22.
- `STATIONARY_PROFILE_UNDERDETERMINED_BY_BRANCH_DATA`, `NO_NET_ACCRETION_BC_UNDERIVED`:
  the flux law and no-net-accretion BC are branch/boundary data, not derived here.
- `HBAR_PROVENANCE_UNDETERMINED`, `HBAR_FREE_SUBSTRATE_RELATION_MISSING`
  (`BLOCKS_HBAR_EMERGENT`), `H_2PI_RATE_CLASSIFICATION_UNDERDETERMINED`: `hbar`
  remains an explicit action coefficient; the `h`-vs-`hbar` `2 pi` split is
  winding-bookkeeping only.
- `BULK_METRIC_SPEED_NORMALIZATION_UNSPECIFIED`,
  `PARENT_METRIC_ACOUSTIC_IDENTIFICATION_MISSING` (`BLOCKS_BULK_EQUALS_C_S` -- no
  `EQUALS` verdict absent a source equation), `BRANE_ZERO_MODE_REDUCTION_UNDERIVED`,
  `BRANE_PHOTON_CONE_REQUIRES_PROFILE`: the bulk and brane cone normalizations are
  not fixed here.

## What this achieves physically

The medium's velocity structure is settled: a phonon sound speed `c_s^2 = 5 K
rho^4/m_GNLS` that scales as `rho^2`, a light-cone ceiling `c = c_gamma` earned
from wave kinematics (with time dilation), and -- the load-bearing honesty -- a
light/sound ratio `lambda_gamma` that the parent action leaves **free**. Fixing
`lambda_gamma` (and hence `c_gamma = c_s`) would require a source equation the
action does not contain; this is the first appearance of the calibration knob that
the knit (Part VI, cone-lock) later pins by calibration, not derivation.

## What is still missing

The numeric value of `lambda_gamma = c_gamma/c_s` (bulk and brane), the flux/`J_crit`
law, and the provenance of `hbar` are all deferred, carried as named residuals.
This stage derives the sound speed and pins the ceiling; it does not derive the
light/sound ratio (by design -- it is a free input) or close the flux/`hbar` gaps.

## Next step

Stage006 (I-3) specifies the two-phase material-state ontology (the order field
`chi_B`) -- the postulated microstructure that distinguishes the ordered
(brane) phase from the bulk, with its recovery reduction and the carried
`theta`-as-Maxwell-`phi` no-go -- consuming this stage's velocity structure and
the dimensional dictionary.

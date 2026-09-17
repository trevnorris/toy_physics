**CLEAR**

Packet identifier (as supplied, not recomputed): `MANIFEST.json` revision **S11 nonlinear-pencil finite-core NP1–NP4 fidelity contract v1**, `aggregate_sha256` `9c7e8e41dd1820647c3dc787ad171043944038c70d37b873699d414f49396e57`, author Codex, `created_utc` 2026-09-17T03:44:41.853371+00:00.

This is a statement-fidelity clearance of the bounded finite core. It is not completion of NP1–NP4 (a second independent review is still required) and not a physical S11c pole, scattering, or Fredholm/Keldysh existence result.

---

## Checked assumptions and conventions

- Field/source/modal spaces `X`, `Y`, `K` are independently typed modules over a field. Maps are `A : X → Y`, `V : K → X`, `W : Y → K`, with a supplied linear equivalence `D : K ≃ K` satisfying `D = W A V`.
- The structure field `PairingData.derivative` is a linear map. The name does not prove `A = L'(z*)`; that identification is an application premise.
- `W L₀ = 0` is **not** part of the basic projector algebra. It appears only in `residue_of_inverse_coefficients`.
- Circles are centered at `0`, positively oriented, normalized by `1/(2πi)`. Radius `R > 0` keeps `0` off the path. The `z²−1` example uses radius `2`, so `±1` lie inside the disk and not on the contour.
- Lean field inversion is totalized at `0`. Inverse identities carry an explicit `z ≠ 0` (or `z ≠ ±1`) hypothesis; contour identities rewrite only on the circle.
- Matrix-valued integrals use Mathlib’s finite elementwise matrix norm. That is a finite-dimensional integral convention, not an outgoing operator-norm or domain claim.
- Coordinates and synthetic operators are dimensionless. Observation/forcing order is `O L⁻¹ B`.
- `nonlinearPoleV2` (`directives/S11c_d_NONLINEAR_POLE_CONTRACT.md`) supersedes the unrestricted residue/projector wording in pinned v10 `S11c_d_SHARED_PHYSICS.md` §3b. The v10 file is unchanged.
- Literature hypotheses in `POLE_ASSESSMENT.md` (analytic Fredholm/Keldysh existence, general argument principle, Riesz calculus, physical-sheet tests) are **not** Lean theorems of this increment.

---

## NP1 — typed pairing algebra (`S11NonlinearPole/Modal.lean`)

`PairingData` supplies separately typed `X,Y,K` and an actual `LinearEquiv` `pairing` with `pairing_eq : ∀ k, left (derivative (right k)) = pairing k`.

| Object | Formal identity | Theorem |
|---|---|---|
| `R` | `V ∘ D⁻¹ ∘ W` | `residueCandidate` |
| `RAR = R` | `residueCandidate.comp (derivative.comp residueCandidate) = residueCandidate` | `residue_sandwich` |
| `P_X = RA` | `fieldProjection` | definition |
| `P_Y = AR` | `sourceProjection` | definition |
| Idempotency | `P_X² = P_X`, `P_Y² = P_Y` | `fieldProjection_idempotent`, `sourceProjection_idempotent` |
| Full ranges | `range P_X = range V`, `range P_Y = range(A ∘ V)` | `range_fieldProjection`, `range_sourceProjection` |
| Ranks | `finrank(range P_X) = finrank K = finrank(range P_Y)` when `K` is finite-dimensional | `finrank_fieldProjection`, `finrank_sourceProjection` |

`V` and `A ∘ V` are injective from invertibility of `D` (`right_injective`, `derivative_right_injective`). A singular pairing cannot inhabit this structure: `pairing` is a `LinearEquiv`. The zero map on `ℂ` is separately shown not injective (`zero_pairing_not_invertible` in `Controls.lean`).

Zero-dimensional `K` is allowed: there is no `Nontrivial K` hypothesis. The algebra does not assert that a pole exists.

**Uniqueness** (`residue_unique`): if `range R ≤ range V` and `W A R = W` (stated as `∀ y, left (derivative (R y)) = left y`), then `R = V D⁻¹ W`. Proof extracts `R y = V k` and uses `D k = W y`.

**Coefficient identification** (`residue_of_inverse_coefficients`): from

- `range V = ker L₀`
- `W L₀ = 0`
- `L₀ R = 0`
- `L₀ H + A R = I`

the theorem concludes `R = residueCandidate`. These are supplied algebraic leading/constant identities for a simple inverse expansion. The statement does not construct an analytic expansion and does not prove a Keldysh semisimplicity/existence criterion. The docstring matches the conclusion.

No stronger analytic interpretation is in the formal conclusions. Suggestive names (`derivative`, `residueCandidate`) remain naming, not theorems.

---

## NP2 — finite Laurent moments (`Moments.lean`, `Laurent.lean`)

`moment R f` is the actual Mathlib circle integral `∮_{C(0,R)} f`, scaled by `(2πi)⁻¹`. `moment_zpow` extracts exponent `-1` for every `n : ℤ`. On a complete complex normed space, `moment_finite_laurent` does the same for finite sums `∑ z^{n i} • v i`. `moment_polynomial` is the nonnegative-exponent case and vanishes by that integration, not by definition.

Positive orientation is in the formal identity itself: `moment R (z ↦ z⁻¹) = 1` (`moment_zpow`, `simple_residue`), and the sign mutant `= -1` is rejected.

`doublePrincipal C₁ C₂ z := z⁻¹ • C₁ + (z²)⁻¹ • C₂` is a **supplied** principal part, not an inverse-existence theorem.

The ordered product with affine `A + z B` expands to eight monomials with exponents `[-1, 0, -2, -1]` on `[C₁A, C₁B, C₂A, C₂B]`. The moment is `C₁ A + C₂ B` (`double_log_expansion`, `double_log_moment`). Interpreting `A+zB` as a pencil derivative (comment: `B` as second derivative) is an application premise, not a proved identification.

Rectangular response: `O₀,O₁ : Matrix o n ℂ`, `B₀,B₁ : Matrix n u ℂ` (maps `n→o` and `u→n`). The eight-term expansion has exponents `[-1,0,-2,-1,0,1,-1,0]` on the ordered products `Oᵢ Cⱼ Bₖ`. The residue is exactly

`O₀ C₁ B₀ + O₀ C₂ B₁ + O₁ C₂ B₀`

(`double_response_expansion`, `double_response_moment`). Factor order is retained. The `z⁰` term `O₁ C₁ B₀` is present in the expansion and correctly absent from the residue. No holomorphic remainder is introduced or discarded; the theorems quantify over exactly these finite operands.

---

## NP3 — actual small pencils (`Scalar.lean`, `Jordan.lean`)

**`L = z²`.** `deriv = 2z`. Inverse identity requires `z ≠ 0`. Inverse moment `0` (`square_residue`). Logarithmic moment `2` (`square_log_moment`), and `2*2 ≠ 2` (`square_log_not_idempotent`). Weighted higher-coefficient moment `moment(z · z⁻²) = 1` (`square_higher_coefficient`). Inverse is nonzero at `z=1`.

**`L = z²−1`.** Zeros are `±1` (`twoRoot_zeros`). Radius-`2` path excludes those points (`mem_sphere`) and encloses them (`integral_sub_inv_of_mem_ball`). Logarithmic formula `(z-1)⁻¹+(z+1)⁻¹` requires `z ≠ ±1`. The actual contour moment is `2` and is not idempotent (`twoRoot_log_moment`, `twoRoot_log_not_idempotent`). This is a contour calculation, not a named algebraic count.

**`N = [[0,1],[0,0]]`.** Pencil is the genuine affine `z I − N` (`jordanPencil`; native AST `[[z,-1],[0,z]]` matches). Actual derivative `I` (`jordan_derivative`). Two-sided inverse `I/z + N/z²` for `z ≠ 0` (`jordan_inverse_left`, `jordan_inverse_right`). Determinant `z²`. Kernel at `0` is the first coordinate axis (`v 1 = 0`). Chain `N e₂ = e₁`, `N e₁ = 0`. Inverse/logarithmic moment is `I`, idempotent, trace `2`. Higher coefficient is the nonzero `N`. Defectiveness does not invalidate this affine state projection.

Faithful scalar transfer is inverse entry `(0,1)`, proved equal to `z⁻²` (`jordan_transfer_exact`). Residue `0`, value `1` at `z=1`. Observation `2+5z` and forcing `1+3z` give `2/z² + 11/z + 15` with residue `11` (`affine_response_expansion`, `affine_response_residue`). Frozen maps `2 · (inverse residue) · 1 = 0` (`frozen_response_residue`).

These examples do not claim general Riesz calculus, Smith/chain classification, or an abstract operator argument principle.

---

## NP4 — controls, compact identification, native instrument

**17 paired rejections and 21 positives**, as recorded in `_measurements/S11_lean_pole_contract_checks.json` and generated from `_measurements/S11_lean_pole_contract_check.py`. Every inspected mutant has `expected/outcome REJECTED`, exit status `1`, exactly one diagnostic in `contract_control`, and exactly one `⊢ False`. Positives have empty diagnostics and exit `0`. These are isolated statement files, not replacements in canonical sources.

| Stem | True statement | Rejected mutant | What it tests |
|---|---|---|---|
| `pairing_normalization` | `residueCandidate 1 = 1/2` | `= 1` | `D⁻¹` scaling |
| `projection_normalization` | `fieldProjection 1 = 1` | `= 2` | `RA` vs raw `A` |
| `full_modal_rank` | rank `2` | rank `1` | full 2D pairing |
| `singular_pairing` | `¬ injective (0·)` | injective | singular exclusion |
| `ordered_log` | ordered product entry `1` | `0` | noncommuting order (`norm_num only [ordered_log_entry]`) |
| `contour_sign` | `moment z⁻¹ = 1` | `= -1` | orientation |
| `double_zero_residue` | inverse moment `0` | `1` | residue vs pole |
| `double_log_count` | log moment `2` | `1` | algebraic count |
| `unjustified_idempotency` | `2² ≠ 2` | equality | scalar log not a projector |
| `higher_coefficient` | weighted moment `1` | `0` | `z⁻²` coefficient |
| `nonlinear_cluster` | two-root log not idempotent | equality | cluster ≠ single projection |
| `Jordan_residue_entry` | `I₀₀ = 1` | `0` | defective affine residue |
| `Jordan_higher_entry` | `N₀₁ = 1` | `0` | nilpotent coefficient |
| `zero_residue_response` | transfer at `1` is `1` | `0` | zero residue ≠ zero response |
| `frozen_maps` | residue `11` | `0` | frozen O,B |
| `omit_observation_derivative` | `11` | `6` | drop `O₁ C₂ B₀` |
| `omit_forcing_derivative` | `11` | `5` | drop `O₀ C₂ B₁` |

The truncated targets `0`, `6`, `5` are the three contributions in `2·0·1 + 2·1·3 + 5·1·1 = 11`. Extra positives: Jordan idempotency, actual inverse at `z=1`, physical transfer residue `0`, constant polynomial moment `0`.

Admissible scaled/full pairings in `Controls.lean` are inhabited and consistent (`scaled_residue`, `scaled_projection`, `full_projection_rank`).

**Native instrument** (`_measurements/S11_lean_pole_source_check.py`, report `S11_lean_pole_source_checks.json`): 19 checks and four translation controls, all recorded `passed/true`. It parses selected original AST assignments and reads recorded synthetic operators; it does not import or execute S11c scripts. Identified operands:

- `scalarDouble`: `z²`
- `scalarBothRoots`: `z²−1`
- `affineJordan`: `[[z,-1],[0,z]] = zI−N`
- realization polynomial `z²`, forcing `1+3z`, observation `2+5z`
- state operator `N`, injection `e₂`, recovery `e₁ᵀ`

Translation controls reject wrong Jordan sign `zI+N`, same-kernel rescaling `2z²`, forcing `e₁`, and frozen observation `2`. Six historical native records are read as already-satisfied evidence; the 89-check repair is explicitly historical and not kernel-certified by this increment.

Pinned v10 SHA256 `fd76447db9c6f2b31d7a72219f4705052023c45788392ec70d761f03aea626c5` matches the addendum’s baseline pin. v10 §3b still writes unrestricted `𝓟_* ≡ (1/2πi) ∮ L⁻¹ ∂_ω L` as a Riesz projector with scalar pairing `⟨l, ∂_ω L r⟩ = 1`. The Lean core supports the addendum’s correction; it does not re-certify that unrestricted wording.

---

## Finite proved core versus literature / application

Proved here: typed full-pairing projector algebra; uniqueness and coefficient identification **from supplied identities**; actual finite Laurent circle moments; explicit scalar and Jordan contour counterexamples; compact synthetic identification and mutation controls.

Not proved, and not required for these bounded claims: analytic inverse existence; Fredholm/Keldysh theory; general argument principle or Riesz calculus; Smith/root-chain classification; physical S11c pole search, scattering, homotopy/displacement, D5, or a systematic CAS bridge. Arbitrary pencils still need actual expansions, derivative identification, domains, and remainder hypotheses before the algebraic identities apply. A zero-dimensional modal space still supplies no pole.

---

## Verification limits

I did **not** run Lean, native Python, hash recomputation, or package rebuilds. Builds, 55 axiom lists, 45 check records, object hashes, fifteen package pins, four direct Mathlib source/object pairs, preservation of forty VC plus five T1 objects, and run1–run3 lineage are **author evidence** (`POLE_VERIFICATION.txt`, `S11_lean_pole_validation.json`, `S11_lean_pole_contract_checks.json`). Mathlib source and compiled objects are not in the packet.

I independently read the six modules, audit root, coverage/fidelity/assessment/policy documents, both instruments, native addendum/v10 §3b, and the recorded control sources/diagnostics. Orientation `+1` is in the Lean statements (`moment_zpow`, `simple_residue`), not only in an unread Mathlib file. The axiom parser in the current instrument accepts wrapped/empty lists via a multi-line `[^\]]*` match and an allowlist `{propext, Classical.choice, Quot.sound}`; I did not re-execute its scratch regression.

No required fidelity correction. Optional stronger results, not needed for this clearance: uniqueness/identification mutants; a separate theorem that truncated observation/forcing formulas equal `6` and `5`; general holomorphic-remainder vanishing; any analytic existence theorem.

**CLEAR** for the bounded NP1–NP4 finite core in this packet. Bounded completion still waits on the second independent review and is not claimed here.

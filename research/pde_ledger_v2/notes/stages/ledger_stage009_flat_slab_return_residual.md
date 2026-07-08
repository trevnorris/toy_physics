# ledger_stage009_flat_slab_return_residual

## Status

**Part II — Gravity. II-B2 (build-order 009).** Reshape of **Check A** of gate **`pathA_29` v3** (the gate splits
into stages 009 + 010 per the finalized Part-II split). Source top-line, verbatim:
**`RETURN_RESIDUAL_PREDICTION`** — a JOINT Check-A + Check-B verdict; this stage computes the **Check-A
component** (the bounded residual-radiation prediction) and prints the qualifier plainly: *Check A component
computed here; Check B localization = `ledger_stage010`.* No standalone headline is faked.

- **POSTULATED (labeled):** the executable family — a flat finite slab, brane at `w=0`, return/absorber at `w=d`.
- **PREMISE vs ACCOUNTING (the v3 relabel, carried exactly):** `Z<0` (drain admissibility) is the **PREMISE**
  (`Z_is_premise = true`); the signed local accounting `Z = −M0·ε0/(1+ε0)` is **ACCOUNTING** derived from the
  channel split. The formula never derives the premise.
- **EARNED (Check A):** the solved round-trip transport, the DC continuity fractions + steady balances, the
  bounded residual orders `p_res(ℓ0)=1`, `p_res(ℓ1)=3`, the Z accounting + sign certificate, and the per-channel
  perfect-return reversibility controls.

Ledger-local earned-label (NOT a source verdict token): `RETURN_TRANSFER_DERIVED_RESIDUAL_BOUNDED`.

## Purpose

Stage 008 specified WHAT any brane↔bulk return must cancel; this stage derives — on the postulated flat-slab
family — what a finite DC-sink return actually DELIVERS: a transfer `T_ℓ(ω)` whose deviation from perfect
cancellation is O(ω⁰), leaving a **bounded monopole/dipole `c_s`-radiation residual tied to the drain strength** —
the falsifiable departure from GR (which forbids these channels). Open-item #9 is **sharpened, NOT closed**.

## The derivation (both engines, own routes)

- **Round-trip transport (EARNED):** from the bidirectional Helmholtz basis
  `Φ_ℓ = A_ℓ e^{iωw/c_s} + B_ℓ e^{−iωw/c_s}`, the forward × return phase ratios give the round-trip phase; the
  transit time is SOLVED from it (`τ = d(log phase)/dω / i`) and asserted: `e^{iωτ} ≡ phase`, `τ = 2d/c_S`.
- **DC continuity fractions:** returned fraction `α_ℓ = 1/(1+ε_ℓ)`, neighbor (absorber) fraction `ε_ℓ/(1+ε_ℓ)`,
  with independent per-channel transmissions `ε0, ε1 > 0`. Return transfer `T_ℓ(ω) = α_ℓ·e^{iω·2d/c_s}`; the DC
  values follow by LIMIT and are asserted against the continuity fractions. Steady circulation balances
  (`J_leak = J_return + J_neighbor`, both moments) assert-zero.
- **The bounded residual (the headline content):** a genuine ω-order scan computes
  `ν_ℓ = ord_ω(1 − T_ℓ) = 0` — a finite DC sink leaves a nonzero deviation at ω⁰ (exactly
  `1 − T_ℓ(0) = ε_ℓ/(1+ε_ℓ)`). Return moments `R_ℓ(ω) = −moment·T_ℓ(ω)`; residual amplitudes
  `A_res = kernel_ℓ · moment · (1−T_ℓ)`; **residual orders `p_res(ℓ) = p_raw(ℓ) + ν_ℓ = 1 and 3` COMPUTED** —
  bounded but lower order than the quadrupole. The Check-A classification (`A_strict_pass` False,
  `A_residual_pass` True) is computed from the orders, never a stamp.
- **Z accounting + sign certificate:** channel sum `Z_throat = −M0`, `Z_return = +M0·T0(0)` (+ declared-zero
  replenishment/boundary-DOF channels) reduces to `Z ≡ −M0(1−T0(0)) ≡ −M0·ε0/(1+ε0)`; the certificate
  `−Z(1+ε0)/(M0ε0) ≡ 1` (Z<0 GIVEN ε0>0, M0>0 — computed).
- **Per-channel perfect-return reversibility (EARNED; the contingency proof):** at `ε0→0⁺`: `T0 → e^{iωτ}`,
  `strict_ν0 = 1` (computed) → strict `p_res0 = 2`, and `Z → 0`; at `ε1→0⁺`: `strict_ν1 = 1` **computed directly
  from the ℓ=1 limit** (the source reused ν0 here — the reshape computes it independently) → strict `p_res1 = 4`.
  The residual prediction is contingent on the per-channel transmissions, not baked in.

## Consumed inputs

**From `ledger_stage008` (II-B1), cited — no file reads:** the three radiation kernels (`i·aω/c_S`,
`i·a³ω³/(2c_S³)`, `i·a⁵ω⁵/(27c_S⁵)`), each tied by a **genuine dual-site citation-integrity check** (the consumed
block types the kernels independently; integrity = consumed − pipeline ≡ 0 — remediated from a
set-then-compare-to-self defect the adversarial leg caught, see Verification). The DC relation to stage008's
targets printed: `R_ℓ(0) = −moment·α_ℓ → −moment` as `ε_ℓ→0`. `kernel₂`/`Q2` carried INERT
(`T2_applied = false` — nothing here is derived ℓ=2 physics).

## Exports

- The **falsifiable residual-radiation prediction** (drain-strength-tied, `ε_ℓ`-parameterized) → stage 023
  (pathA_34 residual matching must agree in form/sign/order) + the conceptual falsifier (GR forbids these
  channels; the drain predicts them, bounded).
- `Z` premise-vs-accounting labels; the `ε_ℓ` transmissions (register rows; `FREE-UNREDUCED`, Gate-6-routed).
- Register edges: R24 (`Z` accounting DERIVED, collapsing `Z` into `{M0, ε0}`); the slab-family rows `{d, ε0, ε1}`.
- Stage-010 pointer: the localization leg (`p=2` dsolves, counterfactual guard, DC-sink classifier, NOGO warp) —
  the radiation/Sommerfeld boundary is recorded `ac_check_a_only` in the source, not a Check-B branch.

## Verification

- **Reshape (blueprint §5) — the heaviest bridge-severing of the cluster:** stripped the `.py`→`.wl` JSON/digest
  agreement bridge, the **runtime YAML reads of pathA_28's results** (replaced by the cited-inputs-from-stage008
  pattern), the own-results YAML/report writers, and ALL SHA-256 trace bookkeeping. Both engines standalone,
  print-only, **zero file I/O**, float-free, ledger idioms; the `.wl` derives everything natively (own basis/phase
  solve, own `SeriesCoefficient` order-scan loop, own limits/balances/teeth) with the arity self-check; its
  `Off[Limit::alimv]` suppression was adversarially unmasked and shown benign (warning chatter from `$Assumptions`
  on limit variables; every silenced limit's RESULT independently asserted).
- **Dual-engine:** SymPy **48 PASS / 0 FAIL** · Mathematica **52 PASS / 0 FAIL**, both exit 0, CWD-independent;
  runner transcripts under `scripts/output/` + `mathematica/output/`.
- **Tri-review (fresh agents):** arbiter re-run via the runners (9/9 both engines); **`FIDELITY_CLEAN`** (full
  Check-A coverage diff, every value hand-re-derived incl. all 8 tooth residuals; 009/010 boundary confirmed
  clean; the strict_ν1 independent computation confirmed); **adversarial `ADVERSARIAL_CONCERNS` → 1 BLOCKING
  CAUGHT:** the consumed-kernel citation-integrity block was **set-then-compare-to-self** (a single-source
  `27→11` corruption exited 0) — a genuine rig-class defect. **Remediated** (dual-site kernels; acceptance test:
  all four one-site corruptions now exit 1 in both engines) along with 3 nits (T2 stamp de-counted; the two
  tautological dim checks composed from the formula structure; dead-code hygiene) → fresh-agent re-verify.
- **Teeth (8, all fire):** returning-basis sign flip; continuity-fraction swap; the DC-term-subtraction order
  tooth; Z-composition drop; sign-certificate flip; strict-limit corruption; consumed-kernel corruption
  (dual-site); the `p_res0→0` classification flip. 16-mutant adversarial matrix on top (incl. the dual-corruption
  class failing both engines independently and a planted arity mismatch caught at a real call site).

## Provenance

- Source gate: `software/stage1_solver/tools/pathA_29_brane_bulk_return_{sympy.py,.wl}` (Check-A slice; sources
  unchanged); `software/stage1_solver/reports/pathA_29_brane_bulk_return.md` + `pathA_29_results.yaml`.
- Reshape directive + tri-review artifacts: `research/pde_ledger_v2/_scratch/ledger_stage009_*` +
  `_scratch/adv_stage009/` + `_scratch/reverify_stage009/`.
- Split row: `research/pde_ledger_v2/notes/part2_gravity_atomic_split.md` (id 009; its fan-out slip
  `ε0=1−T0(0)` corrected to `1−T0(0)=ε0/(1+ε0)` — a fidelity-leg catch).

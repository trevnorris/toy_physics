# ledger_stage008_monopole_dipole_return_spec

## Status

**Part II — Gravity. II-B1 (build-order 008; first new Part-II stage).** Reshape of gate **`pathA_28`**. Top-line
verdict, verbatim: **`MONOPOLE_DIPOLE_RETURN_CONDITIONAL`**.

**⚠ Honest scope (load-bearing, carried verbatim-class):** this is a **VERIFIED CONSTRAINT-SPECIFICATION, not a
falsifiable test of no-monopole-radiation.** Verified-solid: the raw ℓ=0/1/2 outgoing amplitudes and their
radiation orders, and the exact moments/orders at which any brane↔bulk return must cancel. NOT established: that
monopole/dipole radiation is suppressed or unavoidable. The cancellation condition (`R0=−M0`, `R1=−D1`) is
algebraic **bookkeeping** — its symbolic "derivation" is an `x−x` identity. `cancellation_possible` is a **literal
parameter flag**, printed as such (`SCOPE: parameter, not computed — track-3 decides`). The real falsification
lives downstream: whether an admissible return can actually deliver the targets (stages 009/010; the Gate-6
`Z0_ret/Z1_ret` selector).

Ledger-local earned-label (NOT a source verdict token): `DTN_LADDER_DERIVED_RETURN_TARGETS_SPECIFIED`.

## Purpose

GR forbids monopole/dipole gravitational radiation (Birkhoff / mass conservation); the drain breaks brane mass
conservation, so the medium generically radiates in ℓ=0/1 unless the brane↔bulk return cancels the source moments.
This stage makes the target precise — the old ledger's open-item #9 rendered as exact amplitudes, orders, and
cancellation conditions — and is the formal home of the targets that pathA_29 (stages 009/010) and pathA_34
(stages 022/023) consume. Register edge R23 records the obligation.

## The DtN ladder (EARNED — scanned, not typed)

From the spherical Hankel closed forms `h_ℓ = j_ℓ + i·y_ℓ`, form `Λ_ℓ = k·h_ℓ′/h_ℓ |_{z=ka}`, series-expand,
normalize the admittance by its static value, and **scan** for the first k-power whose coefficient is purely
imaginary. Both engines find (each by its own scan loop):

| ℓ | first radiating power p_raw | radiation-phase coefficient | kernel |
|---|---|---|---|
| 0 | **1** | `a` | `i·(ω/c_S)·a` |
| 1 | **3** | `a³/2` | `i·(ω/c_S)³·a³/2` |
| 2 | **5** | `a⁵/27` | `i·(ω/c_S)⁵·a⁵/27` |

Raw amplitudes (nonzero, exact closed forms): `A_raw(ℓ0)=i·M0·aω/c_S`, `A_raw(ℓ1)=i·D1·a³ω³/(2c_S³)`,
`A_raw(ℓ2)=i·Q2·a⁵ω⁵/(27c_S⁵)`. Steady control: `lim_{ω→0} raw₀ = 0` (a static drain does not radiate).
Dominance: `p(ℓ0) < p(ℓ2)` from the scanned powers — the monopole channel is ω-dominant where it exists.

## The cancellation bookkeeping (labeled, not earned)

- ℓ=0: the return must deliver `R0(ω) = −M0(ω)`, `M0 = ∫_brane S_leak d³x` — cancels the raw O(ω¹) coefficient.
- ℓ=1: `R1_i(ω) = −D1_i(ω)`, `D1_i = ∫ x_i S_leak d³x + ∫ j_i d³x` (including any carried odd wake) — cancels O(ω³).
- The with/without residual pair: `without` (R→0) symbolically nonzero; `with` (R0→−M0, R1→−D1) ≡ 0, **printed with
  the label: an `x−x` bookkeeping identity, NOT an earned cancellation result.**
- Derivative outlet vertex `g_W0(ω)=ηω`: recorded as `branch_assumption` (adds two powers of ω in the two-vertex
  bookkeeping); NOT verdict-bearing.

## The S_leak anchors (live-source consistency; stage006 boundary)

Both engines derive, by direct exact integration, the two frozen closed forms the moments ride on: the stage-243
Gaussian-derivative identity `∫W′ j^w dw = −√2·ℓ_w·j₀/4` and the stage-244 one-mode law
`S_leak = √2·μ_w·ρ0·E0/(2√π·λ³)` — evidence the leakage lane is genuinely nonzero (no `S_leak=0` shortcut).
Two computed observations back the strict-recovery guard: `sleak243|_{ℓ_w→0} = 0` and `sleak244|_{E0→0} = 0`
(observed_but_not_used — the strict-recovery limit exists but is NOT taken as a basis). **Boundary:**
`ledger_stage006` (I-3) owns the recovery-reduction derivation of the projected law; here the closed forms are
cited frozen anchors only.

## The computed verdict + guards

`compute_verdict(raw_present, condition_works, cancellation_possible)` with the first two booleans COMPUTED:
baseline → `MONOPOLE_DIPOLE_RETURN_CONDITIONAL`; the synthetic control `(True, False, False)` →
`MONOPOLE_RADIATION_UNAVOIDABLE` (the rung is reachable — machinery, not a stamp); `INCONCLUSIVE` otherwise.
The 9 source controls carried honestly: computed where they compute (steady-no-radiation, M0→0 kill,
derivative-vertex-not-basis, quadrupole-survives p=5), printed declarations where they declare (no `S_leak=0`,
no strict-recovery basis, no projection-locking, no track-3 bulk kill) — with the report's own framing: *beyond
keeping the source live, these pass by construction; guards against obvious tautologies, NOT evidence of
suppression-vs-unavoidable.* Declarations are printed provenance, **not counted checks** (tally hygiene).

## Consumed inputs + frozen slice

- **Consumed from `ledger_stage005` (I-2):** the sound-speed law `c_s² = 5Kρ⁴/m_GNLS` — cited, with the exact-value
  citation-integrity check `5·(1/500)·(√10)⁴/1 − 1 = 0` (the frozen slice's `K_eos=1/500`, `m=1`, `ρ=√10` →
  `c_s²=1`). Not re-derived.
- **Frozen slice (CALIBRATED, cited as printed provenance):** `G=c=c_s=1`, `K_eos=1/500`,
  `(a*,L*)=(4731/2500, 18121/10000)` — exact rationals; benchmark calibration, not register knobs.

## Exports (downstream consumers)

- Targets `R0=−M0`, `R1=−D1` + raw amplitudes/kernels → stages 009/010 (pathA_29) and 022/023 (pathA_34).
- **`Q2` = FREE ANCHOR** (NOT derived ℓ=2 physics) → stage 026 (pathA_43 ancestry consumes it as an anchor only).
- Provenance: the DtN expansions reuse `research/4d_2_5pn` spherical-Hankel machinery (cite; stage 029 is the
  corpus's formal DOI-cite home); `M0`/`D1` are the gate's own constructions consistent with the old-ledger
  Part-VIII projected continuity — NOT verbatim Part-VIII objects.

## Verification

- **Reshape (blueprint §5):** the source pair's engine bridge severed in both directions (the `.py`'s JSON/YAML
  writes + report regeneration; the `.wl`'s sympy-JSON Import + agreement payload). Both engines standalone,
  print-only, **zero file I/O** (this stage reads nothing at runtime), no argparse, float-free (frozen slice exact
  rationals), ledger idioms. The `.wl` keeps + extends its own native route (own Hankel closed forms, own
  `Coefficient`/`ComplexExpand` radiating-power scan, own `Integrate` anchors, own `Which` verdict) — no shared
  payloads; def/call arity verified (an in-script arity self-check + the review-side scan; the stage007 lesson).
- **Dual-engine:** SymPy **54 PASS / 0 FAIL** · Mathematica **57 PASS / 0 FAIL** (3 extra = arity self-checks),
  both exit 0, CWD-independent; runner transcripts under `scripts/output/` + `mathematica/output/`.
- **Tri-review (fresh agents):** arbiter re-run via the runners (8/8 both engines); **`FIDELITY_CLEAN`**
  (coverage-diff vs the source — no dropped load-bearing check; the DtN ladder, amplitudes, both anchors, and the
  frozen-slice identity independently re-derived by hand); **`ADVERSARIAL_CLEAN`** (14/14 mutation matrix: scan
  genuineness proven per-engine via Hankel-form corruption, verdict machinery computed + synthetic control real,
  dual-corruption class fails both engines, teeth no-op probes all raise, arity scan clean incl. a
  planted-mismatch detection proof, zero file I/O confirmed).
- **Remediation (nits folded):** 8 vacuous tally entries de-counted (6 frozen-slice self-comparison stamps + 2
  declaration stamps → printed provenance); the 2 recovery-slice observations restored as computed checks; tooth 7
  re-routed through the pipeline's scanned kernel; `.py` assumptions aligned to `real=True`; the `.wl` dead
  parameter wired. Post-remediation fresh-agent re-verify: `REVERIFY_CLEAN`.

## Provenance

- Source gate: `software/stage1_solver/tools/pathA_28_monopole_{sympy.py,.wl}` (reshaped; sources unchanged);
  `software/stage1_solver/reports/pathA_28_monopole.md` + `_results.yaml` + `pathA_28_cancellation_condition.yaml`.
- Reshape directive + tri-review artifacts: `research/pde_ledger_v2/_scratch/ledger_stage008_*` +
  `_scratch/adv_stage008/` + `_scratch/reverify_stage008/`.
- Split row: `research/pde_ledger_v2/notes/part2_gravity_atomic_split.md` (id 008).

---
unit_id: 129
batch: IV.4
verifier_model: claude-opus-4-7[1m]
verify_date: 2026-05-27T00:00:00Z
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 0
findings_total: 0
material_change: false
---

# Verification — unit 129

## Per-finding outcomes

The original auditor report (`redteam/reports/stage_129.md`) returned
`verdict: clean` with `findings_count: 0`. No directive was generated for
stage 129 (no `redteam/directives/stage_129.md` exists). The batch-IV.4
orchestrator-direct mass-fix (`redteam/resolutions/batch_IV4_paper_alignment.md`,
Cluster A / mechanical groups M1+M2) was applied to two cosmetic strings
that the auditor explicitly flagged as non-finding-grade:

- `.wl` banner at `mathematica/moving_throat_pde_stage129_mouth_boundary_layer_mathematica_audit.wl:26` now reads
  `banner["STAGE 129 — EXPLICIT GNLS + LOCALIZED-MAXWELL MOUTH BOUNDARY LAYER"]`
  (previously `STAGE 112`). Confirmed.
- Notes H1 at `notes/stages/moving_throat_pde_stage129_mouth_boundary_layer.md:2` now reads
  `# Moving-Throat PDE — Stage 129: Explicit GNLS + Localized-Maxwell Mouth Boundary Layer`
  (previously `Stage 231`). Confirmed.

Neither edit touched any assertion, substitution, integral, derivative,
or symbolic identity. The three load-bearing checks (D2 normalization,
D3 zero-flux current, D4 boundary-layer ODE residual under
`V_1 → Π Θ/L`) are byte-identical to the audited state.

## Exec log assessment

**SymPy:** exit=0. Notable lines from `scripts/output/moving_throat_pde_stage129_mouth_boundary_layer_sympy_audit.txt`:

- `sigma_Pi(z) = Pi*exp(Pi*(L - z)/L)/(L*(exp(Pi) - 1))` (L9)
- `Normalization = 1` (L10)
- `Zero-flux current J_sigma = 0` (L11)
- `Stationary zero-flux ODE residual = 0` (L12)
- `EXIT_CODE: 0` (L17)

**Mathematica:** exit=0. Notable lines from `mathematica/output/moving_throat_pde_stage129_mouth_boundary_layer_mathematica_audit.txt`:

- `STAGE 129 — EXPLICIT GNLS + LOCALIZED-MAXWELL MOUTH BOUNDARY LAYER` (L3, banner now correct)
- `PASS: profile normalization` (L8)
- `PASS: zero-flux current` (L11)
- `PASS: boundary-layer ODE residual` (L14)
- `Stage 129 Mathematica audit passed.` (L19)

**Output freshness:** The mathematica output transcript was refreshed
post-banner-fix per the batch-IV.4 resolution memo (M6 stale-output
refresh) and now displays `STAGE 129` as the banner header at L3,
consistent with the corrected `.wl` source. The sympy output (dated
2026-05-11) was not regenerated because no `.py` content changed; the
audited state is preserved. Both engines exit 0; both still PASS all
three assertions.

## Material-change assessment

`material_change`: false.

Only two cosmetic strings changed (a Mathematica `Print` banner label
and a notes H1 line). No derived result, no assertion, no closed form,
no constant, no integration, no substitution, no symbol assumption,
no numerical value was affected. Downstream units have no dependency
on the banner string or the notes H1 numbering.

## Side observations (non-blocking)

None. The mass-fix was strictly cosmetic, the auditor's "clean" verdict
was already in place pre-fix, and both engines continue to verify the
three operational deliverables (D2, D3, D4) corresponding to the paper
card's exponential zero-flux profile.

## Verdict justification

The auditor returned `findings_count: 0` with `verdict: clean` and
flagged only a cosmetic banner mislabel as a non-finding-grade
observation. The batch-IV.4 orchestrator-direct mass-fix corrected
both the `.wl` banner (`STAGE 112` → `STAGE 129`) and the notes H1
(`Stage 231` → `Stage 129`) without touching any math. The refreshed
Mathematica transcript exits 0 with the corrected banner; SymPy exits 0
unchanged. All three load-bearing assertions (normalization, zero-flux
current, ODE residual under `V_1 = Π Θ/L`) still pass on both engines.
Verdict: verified.

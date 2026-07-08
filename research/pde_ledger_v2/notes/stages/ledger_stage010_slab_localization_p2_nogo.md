# ledger_stage010_slab_localization_p2_nogo

## Status

**Part II — Gravity. II-B3 (build-order 010).** Reshape of **Check B** of gate **`pathA_29` v3** (the gate splits
into stages 009 + 010 per the finalized Part-II split). **This is the COMPLETING stage of the joint verdict.**
Source top-line, verbatim: **`RETURN_RESIDUAL_PREDICTION`** — a JOINT Check-A + Check-B verdict. Stage 009 computed
the Check-A component (the bounded residual-radiation prediction); **this stage computes Check B — the localization
leg — and emits the COMPLETED joint verdict**, composing it honestly from *(Check A: `A_residual_pass=True`, cited
from `ledger_stage009`)* ∧ *(Check B: `p=2` both completions, computed here)*. The `RETURN_NOGO` branch stays
genuinely reachable (the anti-localizing warp control proves it).

- **POSTULATED (labeled):** the executable family — a flat finite slab, brane at `w=0`, return/absorber at `w=d`
  (stage 009's geometry, carried). Localization is DERIVED on it.
- **BOUNDARY SELECTIONS (labeled — NOT derived):** the outgoing-spherical-branch choice (`C1/C2` substitution in
  the dynamic radial solve) and the constant-branch removal + `1/(4πd)` overlap normalization in the static solve.
- **EARNED (Check B):** two normalizable `m=0` transverse zero modes → a genuine 3D-radial `dsolve` giving flow
  ∝ `1/r²` (**`p=2`**) on BOTH admissible DC-sink completions; the static–dynamic consistency; the gapped Yukawa
  contrast; the counterfactual guard; the anti-localizing warp NOGO control (**`p=3`**); the computed classifier.

Ledger-local earned-label (NOT a source verdict token): `RETURN_LOCALIZATION_P2_DERIVED_NOGO_REACHABLE`.

## Purpose

Stage 009 derived the finite-slab return *response* (a bounded residual prediction, given the drain premise). This
stage answers the complementary question: does Newtonian `1/r²` gravity SURVIVE the finite slab? **Yes, on the
localizing flat-slab family** — both admissible DC-sink completions carry a normalizable transverse `m=0` zero mode
whose 3D-radial Green function gives flow ∝ `1/r²`. The gate is genuinely able-to-fail: a *delocalizing* return (the
anti-localizing warp) gives `p=3` and the same classifier returns `RETURN_NOGO`. Open-item #9 is **sharpened, not
closed** — localization holding for the FAMILY is not the family being SELECTED by dynamics (R19/`W_slab`).

## The derivation (both engines, own routes)

- **Two DC-sink transverse spectra (EARNED).** On `0 ≤ w ≤ d`, each admissible completion is a transverse
  eigenproblem `f'' + m²f = 0` with a normalizable, constant `m=0` zero mode:
  - *destructuring/absorbing* (compact cell, Neumann): `f0 = 1/√d` (`m0=0`), `f1 = √(2/d)·cos(πw/d)` (`m1=π/d`);
    `f'(0)=f'(d)=0`.
  - *bloch_stack* (q=0 periodic cell): `f0 = 1/√d` (`m0=0`), `f1 = √(2/d)·cos(2πw/d)` (`m1=2π/d`); value and
    derivative match at `w=0,d`.
  Normalizations `∫₀^d f² dw = 1`, the ODE residuals, and the BCs are all asserted. The load-bearing fact: **each
  completion has a normalizable `m=0` zero mode** (the localization principle: a normalizable transverse zero mode
  → `1/r²`).
- **The 3D-radial zero-mode solve → `p=2` (EARNED; genuine dsolve, two routes).**
  - *Dynamic route:* seed with the computed `m=0` eigenvalue (a guard raises if seeded with `m≠0`); `dsolve`
    `g'' + (2/r)g' + ((ω/c_s)²−m²)g = 0` FIRST; select the outgoing spherical branch via the `C1/C2` substitution
    (a boundary SELECTION; normalization = the compact zero-mode overlap `1/d`); assert operator residual ≡ 0; then
    `ω→0` → the limit Green `1/(4πdr)`; flow `= −d(Green)/dr = 1/(4πdr²)`; the large-`r` exponent (extracted, not
    typed, via `−lim_{r→∞} r·d(log flow)/dr`) is **`p_dynamic = 2`**.
  - *Static route:* set `ω=0` first → `g'' + (2/r)g' = 0` → general `C1 + C2/r`; select `C1=0, C2=1/(4πd)` (SELECTION)
    → `g = 1/(4πdr)`; residual ≡ 0; flow `= 1/(4πdr²)`; **`p_static = 2`**.
  - *Static–dynamic consistency (computed):* `p_dynamic − p_static = 0`, **and** the stronger equality
    `dynamic_limit_green − static_green ≡ 0` (both `= 1/(4πdr)`; the `ω→0` limit of `e^{iωr/c_s}/(4πdr)` is
    `1/(4πdr)`).
- **The gapped/Yukawa contrast (EARNED).** The massive radial solve `g'' + (2/r)g' − μ²g = 0` with `μ=m1>0` selects
  the decaying `e^{−μr}/r` branch (residual ≡ 0): the higher transverse modes are Yukawa-suppressed, so ONLY the
  normalizable `m=0` zero mode sets the far-field `1/r²`. (An illustrative superposition `Green = Z·(zero) + C·(massive)`
  with `Z` a *cited* Check-A symbol and `C` free shows the zero mode dominates the far field — not a load-bearing
  assert.)
- **Both completions agree + the counterfactual guard (EARNED; the existing tooth).** Running the static solve for
  both completions gives `p_abs = p_bloch = 2` (asserted `p_abs − p_bloch = 0`). The counterfactual guard multiplies
  the solved Green `1/(4πdr)` by `r⁻⁴` (→ `1/r⁵`) and computes the radial-operator residual (with `κ²=0`):
  **`5/(π·d·r⁷) ≠ 0`** → REJECTED. (Hand-check: `∂²_r(r⁻⁵)=30r⁻⁷`, `(2/r)(−5r⁻⁶)=−10r⁻⁷`, sum `20r⁻⁷`, ×`1/(4πd)`
  = `5/(πd)r⁻⁷`.) The correct candidate's residual ≡ 0 is asserted alongside — the solve genuinely has teeth.
- **The NOGO warp control (EARNED; the mandatory able-to-fail companion).** The anti-localizing half-line warp
  `μ(w) = exp(2·k_warp·w)`: the zero-mode norm `∫₀^W exp(2 k_warp w) dw → ∞` as `W→∞` (COMPUTED non-normalizability);
  the continuum Green `∫₀^∞ e^{−mr} dm /(4πr) = 1/(4πr²)`; flow → exponent **`p_delocalizing = 3`**; the SAME
  classifier returns **`RETURN_NOGO`**. The falloff-tension witness (`p_abs=2` AND `p_delocalizing≠2`) fires
  (`tension_status = "witnessed"`).
- **The computed classifier + the completed joint verdict.** `classify_dc_sink_gate(branch_ps, quadrupole_survives)`
  is a function of the computed exponents (target `=2`; `quadrupole_survives=False → NOGO`; all `==2 →
  RETURN_RESIDUAL_PREDICTION`; none `==2 → NOGO`; else `BC_DEPENDENT`). With `p_abs=p_bloch=2` and
  `quadrupole_survives=True` → the Check-B headline `RETURN_RESIDUAL_PREDICTION`. The COMPLETED joint verdict is
  printed as the composition `(Check A cited) ∧ (Check B computed)`, with `RETURN_NOGO` reachable if the return
  delocalizes (warp → `p=3`) OR the quadrupole channel does not survive.

## Consumed inputs

**Cited — no file reads; genuine DUAL-SITE citation-integrity (two independently-typed sites, `consumed − pipeline
≡ 0`, plus an exact-value anchor):**
- **From `ledger_stage008` (II-B1):** `p_raw(ℓ2) = 5` (the ℓ=2 quadrupole raw outgoing order). Used to DERIVE
  `quadrupole_survives := (p_raw(ℓ2) is a finite non-negative integer order) → True` (the GR-allowed quadrupole
  channel radiates). An exact-value anchor (`p_raw(ℓ2) == 5`, the frozen stage008 export) pins it against a
  coordinated both-sites corruption. NOT recomputed ℓ=2 physics (`T2_applied = false`).
- **From `ledger_stage009` (II-B2):** the Check-A component `A_residual_pass = True`, used ONLY to compose the joint
  headline print. Anchored `cited True`. Not recomputed.
Every one-site corruption of either value — and a coordinated both-sites corruption of the stage008 value — fails
BOTH engines.

## Exports

- The localizing **bulk Helmholtz `m=0` zero-mode** + the `p=2` survival on the slab family → stage 024/026
  (pathA_43 consumes the `Φ_ℓ(w,r)` bulk mode + the projected-continuity operator lineage; this stage is where the
  bulk mode's localization is formally earned in the v2 ledger).
- The completed joint verdict + the reachable-NOGO able-to-fail → the conceptual falsification record.
- Register: **zero new counted knobs**; `k_warp` (`[k_warp]=L⁻¹`) is a control-construction symbol (tracked, not
  counted); edge R25 (localization `p=2` EARNED-within-family, able-to-fail via the warp NOGO — does NOT discharge
  R19 or R23). The R19/`W_slab` caveat is printed, not upgraded.

## Verification

- **Reshape (blueprint §5) — same bridge-severing as stage 009:** stripped the `.py`→`.wl` JSON/digest bridge,
  ALL SHA-256 static/dynamic/per-branch trace ids + `structure_id` + `expr_digest`, and the YAML/report writers.
  Both engines standalone, print-only, **zero file I/O**, float-free, ledger idioms; the `.wl` derives everything
  natively (own eigenfunctions + `Integrate` norms, own `DSolve` + branch selection + `Limit[…,ω→0]` + native
  large-`r` exponent, own counterfactual, own warp `Integrate` divergence, own classifier) with the arity
  self-check. No `Limit`/`Integrate` messages are silenced.
- **Dual-engine:** SymPy **71 PASS / 0 FAIL** · Mathematica **81 PASS / 0 FAIL** (81 = 71 shared + 10 `.wl` arity
  self-checks), both exit 0, CWD-independent; runner transcripts under `scripts/output/` + `mathematica/output/`.
- **Tri-review (fresh agents):** arbiter re-run via the runners (both engines, repo root + foreign CWD);
  **`FIDELITY_CLEAN`** (full Check-B coverage diff, every value hand-re-derived incl. the counterfactual
  `5/(πdr⁷)`, `p=2`/`p=3`, transverse norms/BCs, the Green-equality, and the independent MMA basis; the 009/010
  boundary confirmed clean — Check-A content cited not recomputed); **`ADVERSARIAL_CLEAN`** (22-mutant matrix — a
  hardcoded-wrong-Green mutant proved the counterfactual guard is genuine not decorative; the classifier genuinely
  reaches `RETURN_NOGO`; the warp divergence is computed; the `.wl` is a genuinely independent DSolve basis; teeth,
  arity, and tallies reconcile). The adversarial leg's one substantive find — a *coordinated both-sites* corruption
  of the stage008 citation to another finite integer escaped (magnitude is not load-bearing for
  `quadrupole_survives`) — was folded to parity with stage 009 by an exact-value anchor.
- **Remediation (4 nits) → fresh-agent `REVERIFY_CLEAN`:** added the stage008 exact-value anchor (proven the sole
  gate closing the both-sites escape); de-counted 4 vacuous/set-then-compare dim stamps + the m=0 seed restatement
  (genuinely vacuous — physics coverage untouched; retained derived dim checks keep their teeth, proven by a `dim_r`
  ablation); removed a cosmetic `rewrite(sp.exp)` no-op. Tallies dropped 76/86 → **71/81** accordingly.
- **Teeth (8, all fire):** counterfactual-guard; warp-normalizable flip; radial-operator coefficient; `m=0` seed
  guard; both-completions-agree; classifier able-to-fail (two prongs — `p=3` and `quadrupole_survives=False` both →
  `RETURN_NOGO`); transverse zero-mode normalization; the dual-site consumed-input corruption (all four one-site +
  the coordinated both-sites, now anchored).

## Provenance

- Source gate: `software/stage1_solver/tools/pathA_29_brane_bulk_return_{sympy.py,.wl}` (Check-B slice; sources
  unchanged); `software/stage1_solver/reports/pathA_29_brane_bulk_return.md` (Check-B content :11–17,:24–25) +
  `pathA_29_results.yaml`.
- Reshape directive + tri-review artifacts: `research/pde_ledger_v2/_scratch/ledger_stage010_*` +
  `_scratch/adv_stage010/` + `_scratch/reverify_stage010/`. Running-start source map:
  `research/pde_ledger_v2/notes/stage010_pathA29_checkB_source_map.md`.
- Split row: `research/pde_ledger_v2/notes/part2_gravity_atomic_split.md` (id 010). Completes the pathA_29 fold
  begun at stage 009.

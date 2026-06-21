# Directive pathA_C0e — Is the τ≈0.029 wall gauge? (gauge-invariant curl gate + mechanism ID)

**Status:** READY (Claude-authored 2026-06-20; design-review = SOUND-WITH-FIXES → all 13 fixes applied [gate the 4
Maxwell modes not the phase mode; corrected `curl∘grad`→`div∘grad` incompatibility; mode 2 enters the C0f basis when
gauge; pinned saved-state/matrix artifact policy; added AMBIGUOUS/band-overlap outcomes; dropped the coupled-capture
monotonicity claim + added the real-lane δψ formula; trial-residual-eval scope note; explicit curl-norm spec;
non-circular multi-sample stream-function negative control; residual-budget blockwise Jv; C5 as an order estimate needing
row/block evidence + λ-preimage classification; primary/secondary mechanism + unexplained budget; full plan path]). Then
the user EXECUTION gate. Three-reviewer-converged plan:
`software/stage1_solver/_scratch/pathA_C0e_agreed_plan_for_GLM.md` (§6 = Claude+Codex consult; §8 = GLM 5.2 review folded
in + Claude's C5 odd-even refinement). Follows `pathA_C0d` (verdict MIXED_GAUGE_PLUS_RESIDUAL, fidelity-verified) which
adversarial + Codex + GLM review judged OVERSTATED: modes 1/3/4 are ~99.9% gradient (gauge), only mode 2 is a partial
exception, because C0d's gate used a hand-picked weighted-divergence threshold instead of the gauge-INVARIANT field
strength. Step of option C, task #78.

**Date:** 2026-06-20
**Owner:** Codex (codes + iterates until exit 0; extends the C0c/C0d diagnostics in
`patha_c0_conditioning_spike.py`). Claude reviews after (fidelity + adversarial clean agents).
**Trigger:** C0d classified gauge-vs-stiffness on `P_G` + a weighted-divergence-RESIDUAL threshold (0.1), which is a
gauge-fixing-operator response, NOT a physical-content measure; it never measured the gauge-INVARIANT field strength
`F_rw = curl A`. Three reviewers converged: the decisive test is the per-mode curl fraction, with **mode 2 as the swing
vote**. This directive runs that gate, IDs the conditioning mechanism, and recommends the fix — it does NOT implement the
fix (that is the gated follow-up C0f).

## Scope (READ FIRST — this is the first step that touches the SOLVER, so the boundary is strict)
**DIAGNOSIS + framing-check ONLY.** C0e MUST NOT modify the physical residual, the faithful PDE operators, frozen
physics, or the `physical_export_permitted` guard. C0e MUST NOT implement any fix: no deflation, no gauge-fix, no stencil
change, no re-crawl, no changed-ξ Jacobian re-assembly. The only solver interaction permitted is a SINGLE read-only
linear solve `J·δx=−F` (C0e-0) plus a few trial RESIDUAL EVALUATIONS at trial points `x+δx`; these run only on temporary
tensors — NO accepted state, crawl, line search, or Newton update is ever written, and the state is not advanced.
Evaluate everything else on the EXISTING C0b/C0c saved states + saved Jacobian matrices
(`runs/pathA_C0b_wall_diagnosis/matrices/*.npz`; recompute the SVD modes by fresh dense SVD as C0c/C0d do).

**Artifact policy (PIN THIS — verdict integrity).** Use the SAME saved state + assembled Jacobian for C0e-0..C0e-3.
Prefer the DEEPEST saved stall attempt for which an assembled `J` exists (the τ closest to the 0.029 wall). If only the
converged τ=0.03 matrix C0d used (`attempt_tau_0p03_bt_0.npz`) is available, use it AND say so; if both a stalled-τ≈0.029
and a converged-τ=0.03 matrix exist, run the gate on BOTH and LABEL each verdict with its τ source (do not silently mix).
State in the report exactly which artifact every number came from.

The fix itself (Phase B) is the SEPARATE, gated directive **C0f**, whose FORM (minimal adaptive deflation vs a
consistent/staggered stencil) is DETERMINED by C0e's mechanism finding (C0e-3). CPU; `timeout 600` per script (split if
needed; timeout → NOT_MEASURED); standalone `python3`; no commentary `python3 -c`; YAML/markdown human output, JSON only
for machine artifacts; chunk-1a/1b/1c gates must still pass.

## Key physics (confirmed in the consult + GLM review + Codex design-review against the actual code)
- The Maxwell A-sector residual block = a PHYSICAL curl term + a SOFT gauge-fixing penalty (`operators.py` ~429–438):
  `F_rw = d_r(aw) − d_w(ar)` (the gauge-INVARIANT field strength); penalty `(1/ξ)·grad(Z·D_A·A)` with
  `D_A = axisymmetric_vector_divergence(ar,aw)`, `Z = localization_weight`, **ξ = 1.0** (verified: BranchSmokeConfig
  default, not overridden through `default_closed_branch_config`→`frozen_patha_b2a_branch`), and
  `sponge_gauge_strength = 0.0` (no zeroth-order boundary gauge fixing).
- This is a CHARGED-MATTER system: the local U(1) gauge generator acts on BOTH fields,
  `(δψ = i·(q/ħ)·λ·ψ, δa0 = 0, δar = ∂_r λ, δaw = ∂_w λ)` for scalar λ — in real lanes (ψ=psi_real+i·psi_imag):
  `δψ_re = −(q/ħ)·λ·psi_imag`, `δψ_im = +(q/ħ)·λ·psi_real` (δa0=0 because the ansatz is stationary so −∂_t λ=0; a0 is
  dynamical/elliptic via Gauss and contributes no near-null modes). A pure A-lane `∇λ` (C0d's basis) is an INCOMPLETE
  gauge transformation — the coupled generator is the correct one. Source `q/ħ` (the gauge charge/coupling) from
  `coupled_branch.py` (~386), don't guess.
- **Discrete-operator structure (Codex design-review correction):** `F_rw` and `∇λ` use the SAME separable centered
  gradient operators, so `curl(grad λ) ≈ roundoff` — a pure gradient is essentially curl-free, so the curl fraction IS a
  clean gauge discriminator. The INCOMPATIBLE operator is the gauge penalty's `div∘grad` (face-average divergence ∘
  center gradient ≠ a compatible discrete Laplacian). Still measure `F_rw` DIRECTLY (do NOT infer curl from `P_G`): the
  small non-gradient remnant (`1−P_G≈0.1%`) and high-k content can carry curl; the boundary-vs-interior split is a CHECK,
  not a presumed failure.
- **Scale clue (Claude C5 — an ORDER ESTIMATE, not a verdict).** A rough Maxwell-block expectation: a smooth-gauge-mode
  k² penalty response would land at σ~O(1–10) on this grid — ~9 orders ABOVE the measured σ≈1e-11. That gap suggests the
  near-null gauge modes are NOT smooth but live in the `div∘grad` penalty Laplacian's near-null space — most likely
  ODD-EVEN/CHECKERBOARD decoupling from the mismatched stencils (and/or a small assembled Maxwell-row scaling). This is a
  HYPOTHESIS to be settled by C0e-3's row/block-norm + spectral evidence (NOT asserted up front), because if it is
  odd-even decoupling a consistent/staggered `div∘grad` stencil is a STRUCTURAL fix that generalizes across τ (where a
  one-shot deflation may not).

## Work items (revised ordering from plan §8; lead with the decisive test)

### C0e-0 — framing insurance: does the full Newton step reduce ‖F‖?
At the pinned saved stall state, solve `J·δx = −F` once (read-only; existing assembled J). Report `‖F(x+δx)‖` vs `‖F(x)‖`
for the FULL step (no line search), and `‖F(x+δx_phys)‖` with the near-null/gauge component of `δx` removed; report each
component's norm. (These are trial residual evaluations on temporary tensors — no state advance, no line search, no
Newton update written.) Interpretation: if the full step's failure to reduce ‖F‖ is dominated by the amplified
(~1/σ≈1e11) gauge-direction component, the CONDITIONING framing is confirmed; if ‖F‖ is irreducible even with the gauge
component removed, GLOBALIZATION/geometry is an independent contributor (flag it).

### C0e-1 — THE GATE: gauge-invariant curl fraction (the 4 Maxwell-lane modes; mode 2 = hard swing vote)
For the **4 Maxwell-lane modes** (recomputed by fresh dense SVD of the pinned saved matrix), measure the gauge-INVARIANT
curl fraction `‖Z·F_rw[v]‖ / ‖A[v]‖` using the REAL `F_rw` operator (`operators.py` ~431) and weight Z. **Norm spec:**
report BOTH the unweighted code-space Euclidean norm AND the cell-volume-weighted physical norm, and classify using the
SAME norm used for the controls. Report the boundary-vs-interior split of each mode's curl (a check for boundary-closure
artifacts). **Mode 0 (the C0c phase mode) is handled SEPARATELY** — it has ~no A-lane energy so `‖Z·F_rw‖/‖A‖` is
denominator-noise; report its `A_norm` and `Z·F_rw_norm` for the record and classify it by its (already-confirmed) phase
capture (C0c) + the C0e-2 coupled capture, NOT by the curl ratio.

Calibrate the gauge-vs-transverse decision against CONTROLS, not a hand-picked absolute (the C0d sin): (a) several
held-out smooth + high-k discrete gradients → curl fraction ~0 (the "gauge" band); (b) several independently-constructed
transverse fields (see acceptance 6) → curl fraction O(1) (the "physical" band). Classify each Maxwell mode with an
EXPLICIT outcome and the numeric band separation/margin: `GAUGE` (curl in the gauge band) / `TRANSVERSE` (curl in the
physical band) / `AMBIGUOUS` (between bands) / `CONTROL_BANDS_OVERLAP` (controls don't separate — gate inconclusive) /
`LOW_CAPTURE_LOW_CURL_OTHER` (neither gradient-captured nor curl-carrying — investigate). **Mode 2 is the decisive
gate:** curl in the gauge band ⇒ all Maxwell modes gauge; curl in the physical band ⇒ genuine transverse stiffness in
exactly ONE mode.

### C0e-2 — coupled-gauge projection (confirmatory; answers the A-only-incompleteness objection)
Build the COUPLED local-gauge generator basis `G_cpl = span{ (δψ_re=−(q/ħ)λ_k·psi_imag, δψ_im=+(q/ħ)λ_k·psi_real,
δa0=0, δar=∂_r λ_k, δaw=∂_w λ_k) }` over a multi-λ scalar basis (reuse the C0d gradient machinery for the A-lanes; add the
matter-phase partner from the saved state's ψ; source `q/ħ` from `coupled_branch.py` ~386, shown not guessed).
Orthonormalize. Re-project the 5 modes onto `G_cpl`; report per-mode coupled capture BESIDE C0d's A-only `P_G` — do NOT
require coupled ≥ A-only (the coupled subspace is a graph over λ, not a superset of pure A-gradient space; capture may go
either way). High coupled capture confirms a mode is a full gauge direction, not just an A-lane gradient.

### C0e-3 — mechanism ID (NOT merely confirmatory): residual budget + λ-preimage + spectral structure
Identify WHY the gauge modes are near-null at σ≈1e-11 despite the O(1) penalty coefficient, with EVIDENCE before any label:
- (a) **Blockwise `J·v` residual budget** for each mode — split `J·v` by residual block: Maxwell curl rows, Maxwell
  gauge-penalty rows, Maxwell current-source rows, matter covariant/A-coupling rows, Gauss/charge rows, wall/mass-constraint
  rows, plus an explicit UNEXPLAINED remainder. Report each block's norm to show which terms cancel and at what scale
  (cancellations may live OUTSIDE the ar/aw rows).
- (b) **Assembled Maxwell-row scaling:** read the assembled Jacobian's ar/aw row scaling relative to the other blocks (is
  there a small row weight making the penalty's J-contribution tiny?).
- (c) **λ-preimage spectral classification:** reconstruct the scalar λ-preimage of each near-null gauge mode (the λ whose
  coupled/gradient generator best matches the mode), and classify ITS grid-frequency content SMOOTH (low-k) vs
  CHECKERBOARD (odd-even, high-k) — classify the λ-preimage, not only the raw A lanes.
Conclude a `primary_mechanism` (ROW_SCALING / ODD_EVEN_DECOUPLING / SMOOTH_K2 / OTHER) + `secondary_flags` (mixed causes
are realistic: odd-even + boundary closure + row scaling/current cancellation) + the unexplained budget. Do NOT conclude
ODD_EVEN_DECOUPLING without the row/block-norm + λ-preimage evidence. If ODD_EVEN is the primary mechanism, state whether
a consistent/staggered `div∘grad` stencil is a viable STRUCTURAL fix (DESIGN NOTE ONLY — do NOT implement a stencil change
in C0e).

### C0e-4 — VERDICT + Phase-B (C0f) recommendation
Exactly one cluster verdict from the gate (C0e-1) + mechanism (C0e-3):
  - **WALL_IS_ALL_GAUGE** — all 4 Maxwell modes in the gauge band (incl. mode 2) AND the phase mode confirmed ⇒ recommend
    C0f = minimal ADAPTIVE deflation of the FULL near-null gauge set {phase + coupled-gauge modes 1,2,3,4 + the 2-D
    G_harm} (mode 2 IS included here since it is gauge in this branch) — and, if mechanism=ODD_EVEN, also/instead the
    stencil fix — + re-crawl, with re-SVD at each new stall.
  - **GAUGE_PLUS_ONE_TRANSVERSE** — modes 1,3,4 gauge but mode 2 in the physical band ⇒ recommend C0f deflate the 4 gauge
    directions (phase + 1,3,4 + G_harm) + only mode 2's gradient PROJECTION; the surviving transverse remnant is the
    narrow A-sector question (NOT a production solver).
  - **GAUGE_FRAMING_REFUTED** — modes 1/3/4 unexpectedly show real curl (gauge band test fails) ⇒ STOP; the gauge reading
    is wrong; report the per-mode curl evidence and reconsider before any fix.
  - **DIAGNOSTIC_INCOMPLETE** — required evidence NOT_MEASURED, controls overlap (gate inconclusive), or modes land
    AMBIGUOUS so the cluster cannot be classified.
Include the falsifiable numeric `verdict_support` (per-mode curl fraction + control bands + margins; coupled capture vs
A-only P_G; mechanism primary/secondary + its row/block + λ-preimage evidence) and the concrete recommended C0f design
(deflation basis dims incl. the explicit mode-2 treatment, or the stencil change), plus the Single-Arbiter compliance
plan for C0f (deflation lives ONLY in the Newton linear-solve/preconditioner path; `patha_closed_branch_residual` stays
the sole convergence arbiter; merit/line-search on the original ‖F‖) and the acceptance caveats (the no-deflation
reference exists ONLY at τ≈0.029; deflation must be adaptive across τ).

## Acceptance criteria (PASS/FAIL; exit-0 NECESSARY not SUFFICIENT)
1. C0e-0 reports the full-Newton-step ‖F‖ ratio AND the gauge-vs-complement decomposition of δx, from a single read-only
   `J·δx=−F` solve at the pinned saved stall (trial residual evals on temp tensors only; no state advance).
2. C0e-1 reports the gauge-INVARIANT curl fraction `‖Z·F_rw[v]‖/‖A[v]‖` for the 4 Maxwell-lane modes in BOTH norms
   (unweighted code-space + cell-volume-weighted), calibrated against BOTH control families (held-out gradient band ~0;
   independent transverse band O(1)) with a boundary-vs-interior split, and an EXPLICIT outcome per mode
   (GAUGE/TRANSVERSE/AMBIGUOUS/CONTROL_BANDS_OVERLAP/LOW_CAPTURE_LOW_CURL_OTHER) with numeric band separation. Mode 0 is
   reported separately (A_norm, Z·F_rw_norm) and classified by phase/coupled capture, not the curl ratio. Mode 2's
   outcome is called out as the gate.
3. C0e-2 builds a genuine COUPLED multi-λ generator basis `G_cpl` (matter-phase partner via the real-lane δψ formula + the
   A-lane gradient; q/ħ sourced from code, shown not guessed) and reports per-mode coupled capture BESIDE A-only P_G (no
   monotonicity assumed).
4. C0e-3 reports the blockwise `J·v` RESIDUAL BUDGET (all listed blocks + unexplained remainder) + the assembled
   Maxwell-row scaling + the λ-preimage smooth-vs-checkerboard classification, and concludes a primary_mechanism +
   secondary_flags + unexplained budget WITH that evidence (no ODD_EVEN label without it).
5. Exactly one C0e-4 verdict with falsifiable numeric `verdict_support` + the concrete C0f recommendation (incl. the
   explicit mode-2 treatment, the Single-Arbiter plan, and the two acceptance caveats).
6. Anti-hardcode controls (multiple samples each, as unit tests): POSITIVE — several held-out smooth + high-k discrete
   `∇λ` ⇒ curl fraction ~0, high coupled capture, non-A remainder ~0. NEGATIVE (independent, NON-circular) — construct
   transverse A fields from a stream function (pick smooth and high-k scalars χ, build an axisymmetric divergence-free /
   `F_rw`-carrying A from χ) ⇒ curl fraction O(1), LOW gauge capture; do NOT construct the negative control by subtracting
   a vector's own gradient projection. Report the control band separation.
7. Faithful operators + frozen physics + export guard UNTOUCHED (diff); NO fix implemented (no deflation/gauge-fix/
   stencil change/re-crawl/ξ re-assembly); no Newton iteration / state advance / accepted-state write beyond the C0e-0
   read-only solve + trial residual evals on temp tensors; chunk gates pass; report + machine JSON.

**Fail conditions:** classifying gauge-vs-stiffness on anything but the gauge-INVARIANT curl fraction calibrated to
controls; a hand-picked absolute curl threshold not tied to the control bands; applying the `‖ZF‖/‖A‖` ratio to the
A-energy-free phase mode; a circular or single-sample negative control; an A-only "coupled" generator (missing the
matter-phase partner); asserting a mechanism (esp. ODD_EVEN) without the residual-budget + λ-preimage evidence; a
WALL_IS_ALL_GAUGE recommendation that omits mode 2 from the C0f basis; implementing ANY fix or advancing the Newton
state; altering operators/frozen/export; masking NOT_MEASURED; raising the timeout cap.

## Out of scope (→ C0f, gated on this verdict)
The deflation/gauge-fix/stencil fix itself; the re-crawl; adaptive re-SVD during a crawl; the full-budget crawl; the
production-solver decision; `pathA_22`.

## Review (orchestrator, after Codex)
Fidelity agent: is the curl fraction the REAL `Z·F_rw` operator in the stated norms; is `G_cpl` a genuine coupled
generator (matter-phase partner present via the real-lane δψ formula, q/ħ from code); are the controls non-circular
(independent stream-function transverse construction, multiple samples); modes recomputed not cold-loaded; the
residual-budget Jv / row-scaling / λ-preimage classification genuinely computed? Adversarial agent: is the
GAUGE/TRANSVERSE call SOUND (control-band separation real with margin, not a re-introduced arbitrary threshold; AMBIGUOUS
handled honestly); is the mechanism conclusion (esp. ODD_EVEN vs ROW_SCALING vs SMOOTH_K2) supported by the evidence or
over-read; is mode 2 classified honestly and treated correctly in the C0f basis; is the C0f recommendation correct +
minimal + Single-Arbiter-safe? Diff-check operators untouched + no fix implemented. Then gate C0f (the chosen fix +
re-crawl).

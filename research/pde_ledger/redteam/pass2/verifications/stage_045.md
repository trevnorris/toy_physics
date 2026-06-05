---
unit_id: 045
batch: III.1
role: verifier
verifier_model: claude-opus-4-8[1m]
verification_date: 2026-06-05T00:00:00Z
verification_status: verified
material_change: false
fix_applied: false
engines_agree: true
outputs_fresh: true
---

# Audit unit 045 — pass-2 verification

## Scope

Independent clean-context confirmation of the pass-2 audit's clean verdict for
unit 045 (coherent local D/N kernel tracking). No fix was applied. The audit's
sole finding is a low-severity `stale_output` about a CROSS-ref label
(`Stage-27` → upstream stage 044) deliberately deferred to the dedicated
SCRIPT/OUTPUT-band numbering pass. Read-and-reason only; nothing executed.

## Source self-labels (canonical vs. deferred)

Confirmed the canonical self-labels are already in place:

- SymPy docstring (`scripts/...sympy_audit.py:3`) `Stage 045 SymPy audit.` — canonical
- SymPy banner (`:31`) `STAGE 045 — COHERENT LOCAL D/N KERNEL TRACKING AUDIT` — canonical
- SymPy closing (`:202`) `All Stage-045 symbolic checks passed.` — canonical
- Mathematica banner (`.wl:26`) `STAGE 045 — COHERENT LOCAL TRACKING` — canonical
- Mathematica closing (`:153`) `Stage 045 Mathematica audit passed.` — canonical
- Mathematica section-4 import comment (`:119,121`) uses `Stage-044` correctly

The ONLY stale labels are the SymPy `Stage-27` strings at lines 9, 125, 128, 134
(echoed in the committed transcript at output line 35). These describe the
imported continuum-selected quadratic branch equation, which is **Stage 044's**
— confirmed by the same script's correct `Stage-044` import comment (`:173`),
the `D_cont_stage044`/`F_cont_stage044` symbol names (`:174,178`), and the
`expect_zero` labels at `:189,199`. `Stage-27` is the +17 pre-renumber drift
(27→044). The Mathematica engine carries no such label. This is the deferred
SCRIPT/OUTPUT-band class — intentional, not a verifier-blocking defect.

## Outputs internally consistent and passing

- SymPy output: `STAGE 045` banner (`:3`), every residual `= 0` (lines 9–16, 24–25,
  39, 41–44), `All Stage-045 symbolic checks passed.` (`:47`). Line 35 echoes the
  stale `Stage-27` section banner — the single label artifact, deferred.
- Mathematica output: `STAGE 045` banner (`:3`), all `PASS:` lines, no FAIL,
  `Stage 045 Mathematica audit passed.` (`:42`). No stale labels.

Both transcripts internally consistent (canonical headers) and all checks pass.

## Independent identity spot-checks (non-tautology + paper alignment)

1. **Coherence identity `g_B g_R = g_W g_S`.** Re-derived from
   `coupling_density = (λ_W W + λ_φ φ)(η − γ U)`: cross-coefficients give
   `g_R = γλ_W/√(μ_U μ_W)`, `g_S = γλ_φ/√(μ_U μ_φ)`, etc. Both products equal
   `γ λ_W λ_φ/√(μ_η μ_φ μ_U μ_W)`. Equality holds **only** by the shared `γ` and
   same-source `(η − γU)` structure — genuinely non-tautological; the
   extracted-vs-reference cross-check (A1) guards sign/coefficient transcription.
   Matches eq -coherent-condition. ✓

2. **Quadratic→tracking collapse.** With `R_φ → R_U` the tracking numerator
   `(Mmix+Msupp)(δ+ξ+λR_U²ξ) − ξ(δ+ξ)` cancels term-by-term against
   `collapsed_num` at `M_tr = Mmix+Msupp` (`−δξ−ξ²` pairs with `+δξ+ξ²`; the
   `M_tr` terms pair off). Sum = 0 — a real identity, not a vacuous pass.
   Mathematica reaches `numTrack` via the independent `Series[…,{rPhi,rU,0}]`
   route. ✓

3. **D/N closed forms.** `M_tr_req = ξ(δ+ξ)/(δ+(1+λ₀R_U²)ξ)`; at `λ₀=2/9`,
   clearing 9 gives `9ξ(δ+ξ)/(9δ+(9+2R_U²)ξ)` = notes `G_tr` D/N form (out line
   42 residual 0). `F_tr` D/N form matches likewise (out line 44 residual 0). ✓

4. **Range identities.** `1−R_tr = χ₀δ_U/[(1+χ₀)(1+δ_U)] > 0` and
   `R_tr − 1/(1+δ_U) = δ_U/[(1+χ₀)(1+δ_U)] > 0` ⇒ `1/(1+δ_U) < R_tr < 1`,
   matching eq -Rtr range. ✓

## Value reconciliation

The audit's 12-value reconciliation table (0 misaligned) is consistent with the
script lines and output lines on my independent read. `M_mix`/`M_supp`/`M_tr`
are printed carry-forwards (Stages 022/026), legitimately not asserted. The
`coherent normalization residual = R_target − F_tr` is correctly printed (not
asserted to zero) — `R_target` is the free physical target, so this is a test
quantity, not an identity. No genuine `paper_misalignment`.

## Engine cross-check

Both engines independently agree at the level claimed: every `expect_zero`/
`expectZero` residual is 0, and printed intermediates match modulo cosmetic
factoring (e.g. `(eps−1)` vs `−(1−eps)` sign grouping; the same `R_tr`,
`M_tr_req`, tracking num/den). The `.wl` uses genuinely different routes
(`D[D[...]]` cross-derivatives; `Series` expansion) — not a transliteration.
`engines_agree: true`.

## Verdict

Confirmed clean. Math is sound and non-tautological on independent re-derivation;
source self-labels are canonical; both committed outputs are fresh, consistent,
and all-passing. The only finding is the low-severity `Stage-27` → 044 CROSS-ref
label drift, which is the user-deferred SCRIPT/OUTPUT-band numbering class — no
directive issued, intentionally left for the dedicated content-keyed numbering
pass (`redteam/NUMBERING_SCRIPT_OUTPUT_BAND_PLAN.md`). No fix was applied.

`verification_status: verified`
`material_change: false`

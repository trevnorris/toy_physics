# Directive pathA_17 — B2c comparison-validity audit: is the ~9-order MISS real, or an artifact of the simplification / matching dictionary / units?

**Date:** 2026-06-19
**Owner:** Codex (verification/audit; may write read-only check scripts). Claude reviews the findings afterward.
**Trigger:** User concern (correct, and a known repeat-risk): the B2c verdict is a **~7-to-9 order-of-magnitude**
deficit of the toy model's `P0` vs the GR target `54/5=10.8`. A gap that large is exactly the size produced by a
**units / normalization / dimensional-reduction artifact** (dropped volume factor, wrong `c_s`/`a` power, a
mode-normalization convention, a zeroed channel) rather than genuine missing physics (which tends to be O(1)–O(10)).
**Before we accept "missing physics," rule out "we are comparing two quantities that are not the same thing."**

## This is a VERIFICATION-ONLY audit
Do **NOT** modify the model, the frozen family, the extraction formulas, or "fix"/rescue the miss. The deliverable
is a **finding**: is the deficit REAL or ARTIFACT, and if artifact, the specific suspect factor and its order of
magnitude. You MAY write read-only dimensional/sanity-check scripts under `_scratch/` or `runs/` (gitignored). Do not
touch `research/pde_audit/simulation/` or `physical_export_permitted`.

## What is being compared (the object under audit)
- `R_norm = mhat0²·S_port·P0 − 54·G·c_s⁵/(5·a⁵·c⁵)` (see `patha_extraction.py::observable_residuals` and
  `extract_branch`). With `G=c_s=a=c=mhat0=S_port=1`, the target is `54/5=10.8`.
- `P0 = N0/D0`, `D0 = K − B0 − Z0`. `N0` is the Maxwell self-energy transfer normalization (B2b;
  `Γ_port = 4a⁵/(27c_s⁵)`).
- The matching dictionary: `gamma_eff = mhat0²·S_port·P0·a⁵/(27c_s⁵)` vs `gamma_GR = 2G/(5c⁵)`; `R_norm=0` ⇔
  `gamma_eff = gamma_GR` ⇔ `P0 = 54G c_s⁵/(5 a⁵ c⁵ mhat0² S_port)`.
- The simplifications in play: a **rest** defect (B2a background has `A_r0 ≡ A_w0 ≡ 0`), **l=2-only** wall mode,
  the **frozen `homogeneous_isotropic_hooke_v1`** wall (strict-harmonic, g=0), **pinned `m̂0=S_port=1`**, the
  DtN/port normalization, gauge charge `q=0.35`.

## Audit items (each: state the source you checked + the conclusion)
### A1 — Target provenance (is 54/5 the right number?)
Derive `P0_target` from first principles. Confirm `gamma_GR = 2G/(5c⁵)` is the correct GR leading-order
radiation-reaction / quadrupole anchor — find and cite where it is derived (paper/`notes/`, the `graph/` atlas,
`decisions/08`/`decisions/11`). Confirm the `27` in `gamma_eff` traces correctly to `Γ_port=4a⁵/(27c_s⁵)` and that
re-deriving `gamma_eff=gamma_GR` reproduces `54/5`. Flag any inconsistency or unjustified constant.

### A2 — Dimensional consistency (apples-to-apples?)
Assign physical dimensions to `P0, N0, D0, K, gamma_eff, gamma_GR, a, c_s, c, G, m̂0, S_port`. Verify the two terms
of `R_norm` share the same dimension, and that `P0` carries exactly the dimension the target assumes (likely
dimensionless). Verify the natural-unit settings `G=c_s=a=c=1` are applied consistently on BOTH sides. **Critically:
check whether pinning `m̂0=S_port=1` (GATE-A) silently absorbs a large dimensionful or dimensionless factor** that, if
restored, would move `P0_target` or `gamma_eff` by orders of magnitude.

### A3 — Reduction fidelity (does the simplified P0 equal the full-model quantity?)
Audit whether the `§5 field→{K,M,B,Z,N}→P0` reduction faithfully represents the full-model radiative normalization
the GR target is meant to match. Is `P0=N0/D0` the right object? Does `N0` (Maxwell self-energy transfer, with
`Γ_port`, the `ω⁵` scaling, gauge charge `q=0.35`, the χ-overlap/`χᵀWχ=1` normalization) carry the normalization
`gamma_eff` expects? **Hunt specifically for a dropped or convention-buried factor** — a 4D→reduced volume/Jacobian,
a `(c/c_s)` power, a `4π`, a per-mode vs total normalization, an `a⁵`/`c_s⁵` mismatch — anything that could make the
reduced `P0` smaller than the full-model value by many orders. Quantify any suspect factor's order of magnitude.

### A4 — Simplification robustness (does the miss survive relaxing the simplifications?)
For EACH simplification (rest defect `A_r0=A_w0=0`; l=2-only; frozen strict-harmonic Hooke wall; pinned `m̂0=S_port=1`;
DtN/port normalization), assess whether the ~9-order miss is **robust to it** or **contingent on it**. In particular:
does the rest-defect `A_r0=A_w0=0` zero out a vector-current / radiative channel that would dominate `N0`/`P0` for a
moving or excited defect? Is the l=2-only truncation dropping a channel that carries the GR-scale coupling? Could any
single relaxation move `P0` by orders of magnitude? (This tells us whether the miss is a property of the toy model or
of the *rest/minimal* slice of it.)

### A5 — Verdict
State plainly: is the 6.7–9.6-order `P0` deficit a **REAL physical under-production** of the (rest, minimal-Hooke)
toy model, or **(partly/wholly) an ARTIFACT** of the matching dictionary, units, or the simplification? If artifact:
name the factor, its order, and which item (A1–A4) it lives in. If real: state exactly what is robustly established
and which simplifications it does/doesn't depend on (so κ_PV targets the right thing).

## Acceptance
A1–A5 each answered with the source checked and a concrete conclusion; any suspect factor quantified by order of
magnitude; an explicit REAL-vs-ARTIFACT verdict with its dependency on the simplifications. No model/formula changes.

## Review (orchestrator, after Codex)
Claude reviews the findings (clean adversarial agent + own read), cross-checking the dimensional analysis and the
target provenance independently before we either (a) commit B2c as a validated real miss, or (b) treat the deficit as
an artifact and revise the verdict / matching dictionary.

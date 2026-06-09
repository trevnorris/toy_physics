---
unit_id: 185
batch: V.2
created_at: 2026-06-08T00:00:00Z
findings_count: 1
stop_cold: null
applied: true
applied_at: 2026-06-08T22:44:51-06:00
findings_applied: 1
findings_blocked: 0
verification_status: pending
needs_user_resolution: false
is_checkpoint: true
---

# Codex directive — unit 185 (CHECKPOINT — higher bar)

Apply the finding below. After applying, append an `## Applied: F1` block with:
`files_changed`, `summary` (one sentence), and `deviation` (or "none"). If the change is
genuinely ambiguous or unsafe, append `## Blocked: F1` with a specific question.

Do NOT touch paper.tex, notes/, or any prose document — scripts only. Do NOT change the
SymPy `.py`. After editing, RUN `math -script <wl>` (under `timeout 600`) and iterate until it
exits 0 with every in-file check PASS. A timeout (exit 124) is a FAILURE — reformulate, never
raise the cap. The orchestrator independently re-runs afterward.

This is a CHECKPOINT stage: the higher bar applies — both engines required, assertions must be
substantive and non-tautological, and the load-bearing quantities must be genuinely re-derived
(not trusted as literals).

## F1 — mathematica_transliteration (USER-AUTHORIZED full re-author)

**Target:** `mathematica/moving_throat_pde_stage185_microscopic_monomials_mathematica_audit.wl`

**Status:** The user has AUTHORIZED a full re-author of this `.wl`. NOT a paper_misalignment —
value reconciliation is clean (17/17). The pass-2 audit confirmed the checkpoint deliverables
are correct; the defect is that the second engine is a near-line-by-line port.

**Issue:** The `.wl` mirrors the `.py`:
- the slope operator `firstRatioDrift[r] = (D[r,epsVar]/.epsVar->0)/lamScale` is the exact
  transcription of py:36-38 `sp.diff(ratio,eps).subs(eps,0)/lam`, applied to the same
  exponential-ratio objects;
- the load-bearing monomial-exponent vectors `ctrRatioPrimitive` (wl:148-154) and
  `cntRatioPrimitive` (wl:163-171) are hand-coded character-identically to py:157-162 /
  py:172-179, then `firstRatioDrift`-differentiated the same way. A transcription error in
  those exponent vectors would pass both engines.
The det-`M_*` minor check (wl:217-222) and the zero-defect `Solve` (wl:223-226) are
substantive and may be kept.

**Requirement (you design the route):** Re-author the `.wl` so the load-bearing
**monomial-exponent compilation** is derived INDEPENDENTLY rather than hand-coded and compared
to a SymPy mirror — i.e. the `C_tr,*` / `C_nt,*` / `ε_η` monomial drifts must be obtained from
the primitive-ratio structure by a route the `.py` does NOT use (the `.py` hand-writes the
primitive exponent vector and differentiates; the `.wl` must instead establish those exponents
from the monomial definitions so a wrong exponent fails). You choose the route. Keep the
det-`M_*` minor and the `Solve` zero-defect block (they are already independent-ish), but make
the monomial compilation that feeds them genuinely derived. Do NOT change the `.py`. Preserve
EXACTLY (re-derive, don't trust the literal):
- D1: the 13 primitive/microscopic ratio drifts (`Σ_chi`, `Σ_delta`, `Σ_eps`, `Σ_Z`, `Σ_eta`
  and the 8 primitive `δln`'s);
- D2: tracking monomial `C_tr,* = χ₀^(1+δ_U*) δ_U^(1+χ₀*)` ⟹ `Σ_tr`;
- D3: nontracking monomial `C_nt,* = (Z_W/Ω_W²) ε_W^(E_*) δ_U^(-F_*)` ⟹ `Σ_nt`;
- D4: dressing monomial `ε_η = c_ηU²/(K_U K_η)` ⟹ `Σ_eta`;
- D5: observable triangular law (`Θ₁`, `Ξ₁`, `R₁+Ξ₁`);
- D6 (load-bearing checkpoint quantity): `det ∂(Σ_tr,Σ_nt,Σ_eta)/∂(τ₁,κ_η,μ₁) = 1+χ₀*`;
- D7: the exact zero-defect compatibility solve for `(τ₁, κ_η, μ₁)` and the three
  back-substitution checks.

**Acceptance:** the verifier confirms (a) the monomial-exponent compilation is no longer a
hand-coded mirror of the `.py` differentiated by the same `firstRatioDrift` — an independent
derivation of the monomial exponents is visibly present and could fail on a wrong exponent;
(b) all seven deliverable groups still verify, with `det = 1+χ₀*` re-derived; (c) `math -script`
exits 0; (d) `.py` unchanged. Higher bar: the re-author must not leave any load-bearing
monomial as an unexercised literal.

## Applied: F1

- files_changed:
  - `mathematica/moving_throat_pde_stage185_microscopic_monomials_mathematica_audit.wl`
- summary: Re-authored the Mathematica monomial compiler to derive primitive exponent vectors from formal log-coordinate monomial definitions and routed the determinant/zero-defect solve through those compiled drifts.
- deviation: none

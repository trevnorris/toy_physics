---
unit_id: 186
batch: V.2
created_at: 2026-06-08T00:00:00Z
findings_count: 2
stop_cold: null
applied: true
applied_at: 2026-06-09T00:15:03-06:00
findings_applied: 2
findings_blocked: 0
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 186

Apply each finding below in order. After applying, append an `## Applied: F<n>` block under
that finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").
If a change is genuinely ambiguous or unsafe, append `## Blocked: F<n>` with a specific
question and continue with the rest.

Do NOT touch paper.tex, notes/, or any prose document — scripts only. After editing, RUN the
affected scripts (`python3 <py>`, `math -script <wl>`, under `timeout 600`) and iterate until
they exit 0 with every in-file check PASS. A timeout (exit 124) is a FAILURE — reformulate,
never raise the cap. The orchestrator independently re-runs afterward.

## F1 — mathematica_transliteration (USER-AUTHORIZED full re-author)

**Target:** `mathematica/moving_throat_pde_stage186_similarity_orbit_closure_mathematica_audit.wl`

**Status:** The user has AUTHORIZED a full re-author of this `.wl`. NOT a paper_misalignment.

**Issue:** The `.wl` shares the same hand-coded `M_*` matrix and monomial-drift exponent lists
as the `.py` (wl:35-62 ↔ py:44-56). Its "M_* row N matches paper" checks (wl:69-71) compare the
hand-coded `mMat` against `mMatDerived`, where `mMatDerived` is `Coefficient`-extracted from
`monomialLogDrifts` — itself built from the SAME hand-coded exponent lists (wl:35-48). The check
compares two hand-coded copies of the same numbers and cannot detect a transcription error in
the exponent vectors. The `.py` builds `dlog_Ctr`/`dlog_Cnt`/`dlog_eta` as explicit algebraic
expressions; the `.wl` should not just re-encode them as exponent dot-products.

**Requirement (you design the route):** Re-author the `.wl` so the central object `M_*` is
DERIVED from the physical monomial definitions (in the notes / Stage-185 carry-forward) rather
than hand-coded and compared to a second hand-coded copy. The derived rows must then be checked
against the hand-coded `mMat` reference, so a wrong exponent fails. You choose the route. Do NOT
change the `.py` (its construction is the algebraic-expression reference). Preserve EXACTLY:
- D1: the 3×8 `M_*` with rows = `δln C_tr`, `δln C_nt`, `δln ε_η`;
- D2: the convenient 3×3 minor det = `1+χ` (already substantive — keep);
- D3: the linear compatibility solve for `(τ₁, κ_η, μ₁)` and the three closed-form matches;
- D4: the finite five-parameter orbit preserving `C_tr`, `C_nt`, `ε_η`, and the linearization
  reproducing the Stage-185 ledger.

**Acceptance:** the verifier confirms (a) the row-match residuals now come from
monomial-DERIVED rows (not a hand-coded-vs-hand-coded `Coefficient` round-trip); (b) all
deliverables still verify; (c) `math -script` exits 0; (d) `.py` unchanged.

## Applied: F1

- files_changed:
  - `mathematica/moving_throat_pde_stage186_similarity_orbit_closure_mathematica_audit.wl`
- summary: Re-derived the Mathematica `M_*` rows from physical direct monomial definitions under exponential primitive-variable rescaling before checking them against the reference matrix.
- deviation: none

## F2 — insufficient_verification (both engines)

**Target:**
- `scripts/moving_throat_pde_stage186_similarity_orbit_closure_sympy_audit.py:117-124`
- `mathematica/moving_throat_pde_stage186_similarity_orbit_closure_mathematica_audit.wl:139-148`

**Issue:** The block commented "Non-tautological ground check" is a round-trip: it defines
`eps_eta_logdrift = 2C - U - eta_scaling`, solves `==0` for `eta_scaling` (trivially `2C-U`),
then asserts `solved_eta - (2C-U) == 0` — cannot fail. The companion
`eps_eta_logdrift.subs(eta_scaling, Eta_exp)` is guaranteed because `Eta_exp := 2C-U`. The
comment's non-tautology claim is the inverse of the truth.

**Required change (both engines):** Either (a) make the block a genuine ground check — derive
the preserving `K_η` scaling from the physical `ε_η = c_ηU²/(K_U K_η)` log-drift, reading the
`2`/`-1` coefficients from the monomial exponent vector rather than hand-typing `2C-U`, then
solve and confirm `= 2C-U`; OR (b) delete the round-trip block entirely (the genuine content is
already exercised by `finite orbit preserves eps_eta`). Pick (a) if it is cleanly expressible
from machinery already in the script, else (b). In either case remove/correct the
"Non-tautological ground check" comment so it does not assert non-tautology for a tautology.

**Acceptance:** either the round-trip block is gone, or the printed scaling is derived from the
monomial exponent vector (not a `2C-U`-shaped literal); both scripts exit 0.

## Applied: F2

- files_changed:
  - `scripts/moving_throat_pde_stage186_similarity_orbit_closure_sympy_audit.py`
  - `mathematica/moving_throat_pde_stage186_similarity_orbit_closure_mathematica_audit.wl`
- summary: Replaced the eta-scaling round-trip with a derivation from the physical `eps_eta = c_etaU^2/(K_U K_eta^eff)` monomial in both engines.
- deviation: none

---
unit_id: 009
batch: I.1
created_at: 2026-06-04T00:00:00Z
findings_count: 1
stop_cold: null
applied: true
applied_at: 2026-06-04T19:16:20-06:00
findings_applied: 1
findings_blocked: 0
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 009 (RESOLVED; user-approved 2026-06-04)

The single finding is a `paper_misalignment` the user has resolved via direction (a): ADD a script-side check (both engines) that exercises the card's GENERIC-kernel first-moment claim, so the audited verification matches the card's stated generality. This is a SCRIPT coverage addition (no paper edit). After editing, append `## Applied: F1` (`files_changed`, `summary`, `deviation`), RUN both scripts under `timeout 600`, and iterate until each exits 0 with all checks passing.

## F1 — paper_misalignment (script verifies only the special kernel) → RESOLVED direction (a)

**Paper claim (generic, currently under-exercised):**
- `paper/stages/stage_009.tex:16-21`: `W_\ell(s)=\frac1\ell w(s/\ell)`, `\int_0^\infty w(u)\,du=1` — kernel left GENERIC.
- `paper/stages/stage_009.tex:23-29`: `\langle X\rangle_\ell = X(0)+\ell\,\mu_1 X'(0)+O(\ell^2)`, `\mu_1=\int_0^\infty u\,w(u)\,du` — the leading correction carries the FREE first moment `μ₁`.

**Script side (currently special-cased):**
- `scripts/moving_throat_pde_stage009_projected_maxwell_near_throat_sympy_audit.py:108`: `Wexp = exp(-w/ell)/ell` (so `μ₁=∫₀^∞ u e^{-u} du = 1`; the free `μ₁` never appears).
- `mathematica/moving_throat_pde_stage009_projected_maxwell_near_throat_mathematica_audit.wl:21`: `Wel = Exp[-w/ell]/ell`.

**Requirement (Codex designs the route — do NOT just hardcode):** ADD, in BOTH engines, a check that reproduces the card's generic mouth expansion leading term with `μ₁` carried as a FREE symbol — i.e. demonstrate `⟨X⟩_ℓ = X(0) + ℓ·μ₁·X'(0) + O(ℓ²)` with `μ₁ = ∫₀^∞ u w(u) du` for a kernel that is NOT the exponential special case. Acceptable approaches (Codex's choice): (i) a generic normalized one-sided kernel with symbolic moments (carry `μ₁` as a free parameter and assert the leading coefficient equals `μ₁ X'(0)`), or (ii) a SECOND representative normalized kernel (e.g. one-sided box / one-sided Gaussian / Gamma-shape) whose `μ₁ ≠ 1`, asserting the leading coefficient matches its own `∫₀^∞ u w(u) du`. KEEP the existing exponential-kernel checks intact (as the `μ₁=1` representative). The new check must FAIL if the leading coefficient were the kernel-independent `X'(0)` (i.e. it must genuinely depend on `μ₁`).

**Acceptance criteria:**
- A new assertion in each engine ties the O(ℓ) coefficient to a FREE/representative `μ₁ = ∫₀^∞ u w(u) du`, not the collapsed value 1.
- Both scripts still exit 0; all pre-existing checks still pass; no derived value of the existing exponential-kernel results changes.
- This is additive coverage. Flag `material_change` only if a previously-asserted value changes (it should not); the verifier will assess downstream impact (expected: none).

**Verification:** orchestrator re-runs `redteam exec-sympy 009` and `redteam exec-mathematica 009` (refresh outputs); a clean verify agent confirms the new check genuinely exercises the free `μ₁` first-moment (not a tautology / not collapsed to 1) and that all prior checks still hold.

## Applied: F1

- files_changed:
  - `scripts/moving_throat_pde_stage009_projected_maxwell_near_throat_sympy_audit.py`
  - `mathematica/moving_throat_pde_stage009_projected_maxwell_near_throat_mathematica_audit.wl`
- summary: Added a normalized Gamma-shape one-sided mouth-kernel check in both engines whose computed first moment is μ₁=2 and whose O(ell) Taylor coefficient is asserted as ell·μ₁·X'(0).
- deviation: none

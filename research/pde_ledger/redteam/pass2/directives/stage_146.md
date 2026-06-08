---
unit_id: 146
batch: IV.5
created_at: 2026-06-07T00:00:00Z
findings_count: 2
stop_cold: null
applied: true
applied_at: 2026-06-08T10:43:11-06:00
findings_applied: 1
findings_blocked: 0
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 146

Apply F1 below (it now subsumes the former F2). After applying, append an `## Applied: F1` block with:
`files_changed`, `summary` (one sentence), and `deviation` (or "none").

Do NOT introduce new features, refactors, or stylistic changes beyond what F1 requires. Do NOT touch
paper.tex, notes/, or any prose document (the value reconciliation is clean — no paper edit). Edit only
the affine-law check blocks in the two scripts named.

After editing, RUN both scripts (`python3 <py>`, `math -script <wl>`) under `timeout 600` and iterate
until both exit 0 with all in-file checks passing. (The orchestrator independently re-runs afterward.)

## F1 — tautological_check / insufficient_verification (BOTH engines; user-authorized de-taut 2026-06-08)

**Targets:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage146_positive_deformation_expansion_sympy_audit.py:96-127` (the "Exact affine law" block)
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage146_positive_deformation_expansion_mathematica_audit.wl:88-121` (the "Affine laws tested via integral form" block)

**Issue:**
The convex-family affine moment laws are deliverables (notes:76-79): `g_eps = g_* + eps*(gbar_v - g_*)`
and `S_eps = S_* + eps*(Sbar_v - S_*)`, with `g_* = g(Pi_*)`, `S_* = S(Pi_*)`. Both current checks
fail to exercise them:
- **g-side** (py:114, wl:103): the residual subtracts `gminus` as the intercept, so it collapses by
  linearity to `(1-eps)*(g(Pi_*) - gminus)` — the Pi_* ROOT-closure residual, mislabeled "affine law."
  It tests that Pi_* is a root, not the affine moment structure.
- **S-side** (py:115, wl:113): the residual subtracts `Sformula(Pi_*)` (= `sStar`), the SAME quantity
  that forms the Sigma-part of `Sbar_phys`/`sDirect(Pi_*)`. The residual is `(1-eps)*(integral - Sformula(Pi_*))`,
  which is the kernel identity already checked at py:53 / wl:51 — and the eps-slope `(Sbar_v - S_*)`
  cancels entirely, so the affine (eps-dependent) content is never exercised. As written it is
  effectively `x - x` and would PASS with a wrong slope.

**Required change (requirement + acceptance — apply symmetrically to g and S, in BOTH engines):**
Rewrite both affine-law checks so each genuinely exercises the affine moment law with the closed-form
intercept and an independently-computed nonzero slope:

1. **Intercept = the closed-form moment at Pi_*, not gminus / not a round-trip.** Use
   `g_* := gPi(Pi_*)` (the closed form `gFormula` evaluated at Pi_*) and `S_* := Sformula(Pi_*)` as the
   `eps=0` intercepts. (For g this replaces `gminus`; for S keep `Sformula(Pi_*)` but see step 2 so it
   is not a self-cancellation.)
2. **Compute the deformed moment by ONE direct quadrature of the assembled profile**, distinct from the
   symbolic-linearity decomposition: build `Sigma_eps = (1-eps)*Sigma(Pi_*) + eps*varsigma_test` and
   evaluate `gbar_eps := Integrate[/ integrate Sigma_eps*cos(pi x/2)]` and
   `Sbar_eps := Integrate[/ integrate Sigma_eps*K_q]` as single integrals (high precision). Compare each
   to the affine RHS `g_* + eps*(gbar_v - g_*)` resp. `S_* + eps*(Sbar_v - S_*)`, where `gbar_v`,
   `Sbar_v` are the test-profile moments `Integrate[varsigma_test*cos]`, `Integrate[varsigma_test*K_q]`
   computed by their OWN quadrature. Require the residual `< 1e-25` at `eps = 1/10` and `eps = 1/2`
   (keep the existing tolerance/precision). With this construction the residual reduces to
   `(1-eps)*(integral(Sigma(Pi_*)) - closed_form(Pi_*))`, i.e. it now FAILS if the moment closed form
   `gFormula`/`sFormula` is transcribed wrong — the kernel identity becomes the falsifiable content,
   evaluated AT Pi_* (the old check could not catch this).
3. **Add a non-vacuity guard** asserting the eps-slopes are nonzero: `Abs[gbar_v - g_*] > 1e-3` and
   `Abs[Sbar_v - S_*] > 1e-3` (the test profile `varsigma_test = 6 x (1-x)` differs from the canonical
   source, so both slopes are genuinely nonzero). This proves the `eps`-term is actually present, so the
   affine check is not passing trivially via a vanishing slope.
4. **Relabel honestly** — the printed labels/PASS strings should describe the check as the convex-family
   affine moment law verified via direct quadrature against the closed-form intercept `g_*`/`S_*` and a
   nonzero slope; drop nothing-burger phrasing. No label may claim more than the assertion proves.

**Acceptance criteria (the verifier will check ALL):**
- Neither residual subtracts `gminus` as the affine intercept (the g-side intercept is `g(Pi_*)`).
- The deformed moments `gbar_eps`/`Sbar_eps` are computed by a direct integral of the assembled
  `Sigma_eps`, and `gbar_v`/`Sbar_v` by their own quadrature (so a wrong `K_q` or wrong closed-form
  intercept surfaces).
- Non-vacuity asserts on both slopes are present and pass (`Abs[gbar_v - g_*] > 1e-3`,
  `Abs[Sbar_v - S_*] > 1e-3`).
- Both residuals are `< 1e-25` at `eps = 1/10` and `eps = 1/2`.
- The printed labels no longer over-claim and accurately name the affine-moment-law-via-quadrature check.
- The independent moment-formula verifications at py:33-53 / wl:44-51 (each engine's own `integrate`/
  `Integrate` of `Sigma*cos` and `Sigma*K_q` vs the closed forms) are LEFT INTACT — they remain the
  cross-engine-independent core.
- Both scripts exit 0; no deliverable VALUE changes (Pi_*, g_*, S_*, the retuning-law forms unchanged).

**Verification command:**
`redteam exec-sympy 146` and `redteam exec-mathematica 146`; confirm the new intercept/quadrature/
non-vacuity structure appears in both engines, both exit 0, and no printed deliverable value moved.

## Applied: F1

- files_changed:
  - `scripts/moving_throat_pde_stage146_positive_deformation_expansion_sympy_audit.py`
  - `mathematica/moving_throat_pde_stage146_positive_deformation_expansion_mathematica_audit.wl`
- summary: Reworked both affine moment checks to use closed-form `g_*`/`S_*` intercepts, direct assembled-profile quadrature at `eps = 1/10` and `eps = 1/2`, and nonzero-slope guards.
- deviation: none

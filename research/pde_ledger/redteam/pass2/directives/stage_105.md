---
unit_id: 105
batch: IV.2
created_at: 2026-06-06T00:00:00Z
findings_count: 2
stop_cold: null
applied: true
applied_at: 2026-06-06T07:47:09Z
findings_applied: 2
findings_blocked: 0
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 105 (CHECKPOINT)

Apply each finding below in order. After applying, append an `## Applied: F<n>` block under that finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

If a finding's required change is ambiguous or unsafe to apply mechanically, append `## Blocked: F<n>` with a question instead — skip that finding, continue with the rest.

You DESIGN the route. This directive states the requirement and the acceptance criteria only; choose the actual constructions yourself. Do NOT introduce unrelated features, refactors, or stylistic changes.

After editing, RUN the affected scripts (`timeout 600 python3 <path>` for SymPy, `timeout 600 math -script <path>` for Mathematica) and iterate until they exit 0 with all in-file checks passing. A timeout (exit 124) is a FAILURE — reformulate the math, never raise the cap. Getting the scripts to run cleanly is your job; the orchestrator independently re-runs afterward.

Do NOT touch paper.tex, notes/, or any prose documents. Do NOT change any numbering label/cross-reference.

## Context

This is the checkpoint that owns the canonical pin `chi_Q = 1` — the outgoing DtN l=2 fingerprint carried forward by 097/100/106. The value `chi_Q = 1` is CORRECT (independently re-derived) and must NOT change. The fix changes only HOW the pin is forced, not the result. All emitted deliverable values and the committed Mathematica/SymPy outputs must be preserved (the change is method-only).

## F1 — insufficient_verification (checkpoint tautology) + F2 (transliteration, subsumed)

**Target:**
- `scripts/moving_throat_pde_stage105_chiQ_fix_from_outgoing_dtn_sympy_audit.py:49-55`
- `mathematica/moving_throat_pde_stage105_chiQ_fix_from_outgoing_dtn_mathematica_audit.wl:52-59`

**Issue:** The `chi_Q = 1` pin solves `[retarded i*omega^5 coeff] = a^5/(27 c_s^5)` where the LHS is by construction `chi_Q * a^5/(27 c_s^5)` (because `sigma_can` is defined as `4 a^5/(27 c_s^5)`) and the RHS is a TYPED COPY of the same literal. The equation collapses to `chi_Q * K = K → chi_Q = 1` for any `K`; it never evaluates the actual outgoing l=2 DtN fingerprint. For a checkpoint that owns this pin, that is an insufficient-verification defect. (F2: the retarded half of the `.wl` also mirrors the `.py`'s choreography on the same typed RHS — re-grounding the match in genuinely-derived, per-engine-distinct fingerprint constructions resolves this too.)

**Requirement:** In BOTH engines, the target outgoing z^5 (→ omega^5) coefficient that `chi_Q` is matched against must be DERIVED from the exact outgoing l=2 spherical-Hankel DtN operator, not typed. Concretely, before the `chi_Q` solve, construct the outgoing l=2 DtN operator `Lambda_2^out(z) = z d/dz ln h_2^(1)(z)` from the spherical Hankel function, series-expand it, and assert it equals `-3 + z^2/3 + z^4/9 + i z^5/9` (this is the card's Check #3 fingerprint exercise — a genuine, can-fail check); normalize to `Y_2^out = -3/Lambda_2^out`; read its imaginary z^5 coefficient as a DERIVED quantity (it should come out `1/27`); carry that through the `z = a*omega/c_s` map to obtain the derived `a^5/(27 c_s^5)`; then solve `[retarded i*omega^5 coeff] = [DERIVED coefficient]` for `chi_Q` and assert `chi_Q - 1 == 0`. The literal `a^5/(27 c_s^5)` must no longer appear as a TYPED RHS of the canonical `chi_Q` match — it must arrive from the Hankel derivation.

**Acceptance criteria:**
1. A new in-script check derives `Lambda_2^out` from the spherical Hankel function and asserts its series equals `-3 + z^2/3 + z^4/9 + i z^5/9` (passes, and is can-fail: a wrong Hankel form would break it).
2. The `chi_Q` solve references the DERIVED outgoing coefficient; no typed `a^5/(27 c_s^5)` literal remains on the RHS of the canonical match.
3. `chi_Q - 1 == 0` still passes in both engines; `chi_Q` stays a free symbol up to the solve (not pre-substituted to 1).
4. The two engines reach `Lambda_2^out` / its z^5 coefficient by VISIBLY DIFFERENT constructions (e.g. one via a native special-function primitive, the other via an explicit closed form or the spherical Bessel `j_2`/`y_2` pair) — so the retarded/canonical half is no longer a line-by-line port. The verifier will check this independence.
5. The deformed-branch section (py:57-66 / wl:61-78, the polynomial-inversion `xi_Q` route) is left UNCHANGED — it is already independent and correct.
6. The pre-existing retarded coefficient asserts (`omega^2`, `omega^4`, imag `omega^5` carrying `chi_Q`) and the `sigma_Q^can` check still pass; all emitted values unchanged; both scripts exit 0.

## Applied: F1

- files_changed:
  - `scripts/moving_throat_pde_stage105_chiQ_fix_from_outgoing_dtn_sympy_audit.py`
  - `mathematica/moving_throat_pde_stage105_chiQ_fix_from_outgoing_dtn_mathematica_audit.wl`
- summary: Replaced the typed canonical `chi_Q` match target with a coefficient derived from the exact outgoing `l=2` spherical-Hankel DtN fingerprint in both engines.
- deviation: none

## Applied: F2

- files_changed:
  - `scripts/moving_throat_pde_stage105_chiQ_fix_from_outgoing_dtn_sympy_audit.py`
  - `mathematica/moving_throat_pde_stage105_chiQ_fix_from_outgoing_dtn_mathematica_audit.wl`
- summary: Made the two retarded-match witnesses visibly independent by using an explicit spherical `j_2 + i y_2` construction in SymPy and native `SphericalHankelH1` in Mathematica.
- deviation: none

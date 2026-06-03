---
unit_id: 219
batch: VII.1
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-06-02T00:00:00-06:00
verdict: findings
stop_cold: null
findings_count: 1
paper_alignment: aligned
scripts_checked:
  sympy: present
  mathematica: missing
  engines_agree: n/a
  outputs_fresh: true
docs_read:
  paper_stage_tex: present
  notes_stage_files: [moving_throat_pde_stage219_one_port_mixed_bundle_static_kernel_and_square_law_suppression_test_sympy_audit.md]
  paper_appendix: present
---

# Audit unit 219 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_219.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage219_one_port_mixed_bundle_static_kernel_and_square_law_suppression_test_sympy_audit.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part07.tex` (section "One-port same-charge static kernel", lines 230-364, plus index row line 50)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage219_one_port_mixed_bundle_static_kernel_and_square_law_suppression_test_sympy_audit.py`
- mathematica: `(missing)`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage219_one_port_mixed_bundle_static_kernel_and_square_law_suppression_test_sympy_audit.txt`
- mathematica output: `(missing)`

## What the paper claims

`\stagefield{Output}` (stage_219.tex:15) states verbatim: "Exact square-law suppression theorem: the one-port static mixed bundle produces only $x^{-6}$, $e^{-2\kappa x}x^{-4}$, and $e^{-4\kappa x}x^{-2}$ product families, and no new $-1/x$ or Yukawa-strength long-range same-charge attraction." The `\stagefield{Derivation ledger}` (line 13) and notes Section 9 enumerate the supporting deliverables: (1) the reduced static `3x3` kernel `K_red` and the determinant identity `det K_red = Delta*D0`; (2) the Schur-complement formula for `D0`; (3) the six exact inverse-susceptibility entries `chi_qq..chi_WW`; (4) the on-shell shift `delta_V_mix = -1/2 J^T K_red^{-1} J`; (5) the collinear-source factorization `delta_V_mix = -1/2 chi_s S^2` with the stated `N_s`; (6) the outgoing-prefactor bridge `chi_qW = Lambda/D0` and `chi_qW^2 = P0/D0` with `Lambda = P/Delta`, `N0 = Lambda^2`, `P0 = Lambda^2/D0`; (7) the primitive-source product-kernel decomposition into exactly the three families `x^{-6}`, `e^{-2\kappa x}/x^4`, `e^{-4\kappa x}/x^2`. The appendix Theorem (line 359-362) restates the family-membership / suppression verdict. Admissible branch: `Omega_U^2>0, Delta>0, D0>0`.

## What the script claims to verify

The SymPy script builds `K_red` from the paper matrix (lines 37-41), and asserts: det identity `det K_red = Delta*D0` (line 46, det computed natively via `K_red.det()`); Schur complement equals `D0` (line 59, computed via the internal-block inverse `K_int.inv()`); the six inverse entries match their closed forms (line 85, against `K_red.inv()`); the collinear factorization `delta_V_col = -1/2 chi_s S^2` (line 110); the two outgoing-load identities (lines 130-131); and the product-kernel family decomposition `delta_V_primitive = -1/2[C6/x^6 + 2 C4 e^{-2kx}/x^4 + C2 e^{-4kx}/x^2]` (line 164). It also verifies a concrete admissible numeric slice is positive-definite via leading principal minors and eigenvalues (lines 199-202). All assertions are genuine equality-to-zero checks on independently-computed objects.

## Paper ↔ script cross-check

| Paper deliverable | Script check | Status |
|---|---|---|
| (1) `det K_red = Delta*D0` | line 46 | match |
| (2) Schur complement = `D0` | line 59 | match |
| (3) six inverse entries `chi_*` | line 85 (loop over 6 entries) | match |
| (4) `delta_V_mix = -1/2 J^T K_red^{-1} J` | constructed line 92; exercised via (5)/(7) | match (printed + used downstream) |
| (5) collinear `delta_V = -1/2 chi_s S^2`, `N_s` | line 110 | match |
| (6) `chi_qW = Lambda/D0`, `chi_qW^2 = P0/D0`, `Lambda=P/Delta` | lines 130-131 (+119-121) | match |
| (7) product-kernel = exactly 3 families (suppression verdict) | line 164 | match |
| admissibility `Delta>0, D0>0, PD` | lines 199-202 (numeric sample) | match (constructive witness) |
| Independent second-engine (Mathematica) verification | none | missing |

`paper_alignment: aligned` — every paper-side deliverable has a faithful, non-tautological SymPy check; no value/target/sign disagreements found.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 46 | `simplify(det_identity - Delta*D0) == 0` | deliverable 1 | yes |
| A2 | sympy | 59 | `simplify(schur - D0) == 0` | deliverable 2 | yes |
| A3 | sympy | 85 | `simplify(together(Kinv_ij - chi_ij_expected)) == 0` x6 | deliverable 3 | yes |
| A4 | sympy | 110 | `simplify(together(delta_V_col - (-1/2 chi_s S^2))) == 0` | deliverables 4,5 | yes |
| A5 | sympy | 130 | `simplify(together(chi_qW - Lambda/D0)) == 0` | deliverable 6 | yes |
| A6 | sympy | 131 | `simplify(together(chi_qW**2 - P0/D0)) == 0` | deliverable 6 | yes |
| A7 | sympy | 164 | `simplify(together(delta_V_primitive - 3-family decomp)) == 0` | deliverable 7 (suppression) | yes |
| A8 | sympy | 199-202 | numeric `Delta,D0,minors,eigenvalues > 0` on sample | admissibility witness | yes (partial: one point) |

All SymPy rows are "yes": each compares an independently-computed object (det, Schur via inverse, full matrix inverse, the quadratic form built from `Kinv`) against the paper's closed form, so each can fail if the paper claim is wrong. None is tautological (the `expected` closed forms are NOT the objects being computed; e.g. `Kinv = K_red.inv()` is computed before being compared to `P_U/(Delta*D0)` etc.). A8 is a single-point positive-definiteness witness, which is appropriate scope (the paper states admissibility as an assumed branch condition, not a theorem the script must prove for all parameters).

## Findings

### F1 — missing_verification_script (subtype: missing_mathematica)

**Severity:** medium
**Files:**
- `(missing)` — no `mathematica/moving_throat_pde_stage219_*.wl` exists; stage_219.tex:11 states "Mathematica audit: none yet."

**What's wrong:**
This stage is `is_checkpoint: False` but `is_status_only_candidate: False`, so it is a substantive computational stage that requires both engines wherever a second engine CAN independently verify. Every deliverable of stage 219 is a finite, closed-form symbolic computation that Mathematica verifies trivially with native primitives: a 3x3 symbolic determinant (`Det`), a Schur complement via an internal-block inverse (`Inverse`), the full `3x3` symbolic inverse (`Inverse`), a scalar quadratic form, two rational-function identities (`Together`/`FullSimplify`/`Simplify`), the product-kernel family separation (extractable via `Coefficient`/`Exponent` on the `1/x` and `Exp[-2 kappa x]` structure), and a numeric positive-definiteness witness (`Eigenvalues` / `PositiveDefiniteMatrixQ`). Independent Mathematica verification is clearly POSSIBLE here, therefore (per the dual-engine rule — the test is "is it possible," not "is it necessary") a `.wl` is REQUIRED. There is no genuine impossibility argument: nothing in this stage is numerically-only, asymptotic-only, or otherwise out of Mathematica's reach.

**Why this matters:**
Single-engine coverage on a non-status, non-checkpoint computational stage leaves the load-bearing square-law suppression theorem (deliverable 7) without an independent cross-check. A second engine that derives the result through a *different* decomposition is the only protection against a SymPy-specific simplification artifact silently passing.

**Required change:**
Author a NEW Mathematica script (NOT a transliteration of the `.py`) at the exact target path in the directive. It must derive the claims independently — see the claim manifest M1..M7 and the anti-transliteration guard in the directive.

**Verification:**
After Codex applies, the verifier runs `redteam exec-mathematica 219` and confirms the new `.wl` exists, exits 0, and emits the manifest's per-claim PASS lines.

## Independent-derivation check (Mathematica)

No `.wl` exists, so transliteration is not yet assessable. The directive's anti-transliteration guard (different decomposition + native primitives) is mandatory so that when the `.wl` is authored it is a genuine independent re-derivation, not a port of the SymPy choreography.

## Engine cross-check

Not applicable — only the SymPy engine is present. This is itself finding F1.

## Verdict justification

The SymPy script is faithful and adversarially sound: every paper deliverable maps to a non-tautological assertion that compares an independently-computed object (native `det`, Schur via matrix inverse, full `K_red.inv()`, the quadratic form) against the paper's stated closed form. I attacked the determinant identity (det is computed natively, not constructed from `Delta*D0` — passes), the inverse entries (signs and the `Q - K_star*Delta = -(Delta*D0)` denominator flip check out against the saved output), the symbol assumptions (`nonzero` only, no over-strong `positive` hiding branch issues; positivity is checked numerically on a sample, not assumed symbolically), and the suppression theorem (line 164 forces the full computed expression to equal *exactly* the three families, so any spurious `1/x` term would break it — genuine encoding of the family-membership claim). No `paper_misalignment`, `tautological_check`, `hardcoded_result`, `symbol_assumption_error`, `missing_branch`, or `insufficient_verification` defect was found. The sole finding is the missing second engine (F1, `missing_mathematica`), which the dual-engine rule requires because independent Mathematica verification is plainly possible for this entirely closed-form, finite stage.

## Self-test notes

(1) Variable independence: the new `.wl` performs no `D[..,var]` derivatives — all claims are algebraic identities and a quadratic form — so the zero-derivative trap does not apply. (2) Symmetry/parity: no unbounded-domain integrals are involved (the `x`-families are evaluated as algebraic functions, not integrated), so the parity trap does not apply. (3) Trivial-case pre-check: I mentally substituted the script's own admissible sample (`K_star=11, Omega_U=3, Omega_W=4, R=2, G_U=1, G_W=2`) — `Delta=140>0`, `D0=74/7>0`, `det=1480=Delta*D0` (140*74/7=1480 ✓), confirming the det identity numerically and that the manifest checks reduce correctly. (4) Path: the directive Target uses the exact `mathematica/...stage219...mathematica_audit.wl` string with the `_mathematica_audit` suffix. (5) Paper round-trip: F1 adds an engine, changes no constant or identity, so it introduces no new paper_misalignment.

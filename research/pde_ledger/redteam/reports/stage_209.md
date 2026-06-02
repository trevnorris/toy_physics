---
unit_id: 209
batch: VI.1
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-06-01T00:00:00-06:00
verdict: findings
stop_cold: null
findings_count: 2
paper_alignment: partial
scripts_checked:
  sympy: present
  mathematica: missing
  engines_agree: n/a
  outputs_fresh: true
docs_read:
  paper_stage_tex: present
  notes_stage_files: [moving_throat_pde_stage209_pairwise_ratio_optimizer_and_mixed_ray_winner_theorem.md]
  paper_appendix: present
---

# Audit unit 209 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_209.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage209_pairwise_ratio_optimizer_and_mixed_ray_winner_theorem.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part06.tex` (rows: line 49 manifest row; lines 878-918 the "Pairwise ratio optimizer" subsection; line 938 boundary-edge reference; line 1306 anchor MTDC-T10.3)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage209_pairwise_ratio_optimizer_and_mixed_ray_winner_theorem_sympy_audit.py`
- mathematica: (missing)
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage209_pairwise_ratio_optimizer_and_mixed_ray_winner_theorem_sympy_audit.txt`
- mathematica output: (missing)

## What the paper claims

The paper card (`stage_209.tex:15`) states the `\stagefield{Output}` verbatim: "Finite pairwise candidate set, optimized pairwise bracket, special diagonal-neutral and pair-symmetric reductions, pairwise promotion theorem, and mixed-pair winner theorem." The derivation ledger (line 13) says the stage "writes each certified objective as an algebraic function of the ratio, differentiates it, clears radicals and denominators, reduces stationary candidates to quartic equations, and splices those candidates with endpoints and admissibility boundaries." The notes enumerate seven deliverables: (1) the explicit square-root-rational algebraic form of `tau_{ij,*}(r) = 2H0 sqrt(1+r^2)/(k_i + r k_j + sqrt(A + B r + C r^2))` with `A = k_i^2 - 2H0 u`, `B = 2 k_i k_j - 4 H0 v`, `C = k_j^2 - 2 H0 w`; (2) the stationary numerator theorem `N_{ij,*}(r) = 2(k_j - k_i r) sqrt(...) + B + 2(C-A) r - B r^2 = 0` via the denominator functional `Phi = (k_i + r k_j + sqrt(...))/sqrt(1+r^2)` with derivative law `dPhi/dr = N/(2(1+r^2)^{3/2} sqrt(...))`; (3) the quartic elimination theorem `Q_{ij,*}(r) = [B + 2(C-A) r - B r^2]^2 - 4(k_j - k_i r)^2 (A + B r + C r^2) = 0`; (4) the finite candidate-set theorem (Q quartic => at most four real roots, at most six evaluations per envelope, twelve per pair); (5) the optimized pairwise bracket `tau_min^lo <= tau*^best <= tau_min^hi`; (6) two special reductions — diagonal-neutral (`u=w=kappa_*, v=0` => optimizer collapses to `r_grad = k_j/k_i`) and pair-symmetry (`k_i=k_j, u=w` => `tau(r)=tau(1/r)` so `r=1` is critical); (7) the pairwise promotion theorem and mixed-pair winner theorem (bracket-comparison inequalities `tau_{ij,min}^hi < min_m tau_{m,lo}^prim`, etc.). Claim status is "Mixed: ExactClosure, Numerical." (Note: the notes file repeatedly self-labels as "Stage 243"/"Stage 192" in places; this is an internal derivation-order artifact — the paper-restored stage number is 209, and the appendix subsection at line 881 unambiguously attributes this content to "Stage 209".)

## What the script claims to verify

The SymPy script (banner mislabels itself "STAGE 192", line 35, but content is unambiguously stage 209) verifies, by `expect_zero` identities: §I the explicit algebraic form of `tau` against `tau_expected` and the discriminant-numerator reduction `(1+r^2)(k_ij^2 - 2H0 kappa) = A + B r + C r^2` (deliverable 1); §II the identity `tau = 2H0/Phi` and the exact `Phi` derivative law `dPhi/dr - N/(2(1+r^2)^{3/2} S) = 0` (deliverable 2); §III that `Q` has degree 4, the factorization `Q = (J - 2(k_j - k_i r)S)(J + 2(k_j - k_i r)S)`, and `N = J + 2(k_j - k_i r)S` (deliverable 3, and the algebraic kernel of deliverable 4); §IV that on the diagonal-neutral branch `kappa(r) = kappa_*` is constant and `d(tau_diag)/dr` vanishes at `r = k_j/k_i` (deliverable 6a); §V that `tau_sym(r) = tau_sym(1/r)` and `d(tau_sym)/dr` vanishes at `r=1` (deliverable 6b); §VI prints the explicit quartic coefficients in original variables (no assertion). All checks pass with exit 0; output (`...sympy_audit.txt`) mtime 2026-05-11 is newer than the script (2026-04-14), so it is fresh.

## Paper ↔ script cross-check

| Paper/notes deliverable | Script-side check | Status |
|---|---|---|
| (1) Algebraic form of `tau`, A/B/C, discriminant reduction | §I `explicit algebraic tau form`, `discriminant numerator reduction` | match |
| (2) Stationary numerator `N` + `Phi` derivative law | §II `tau = 2H0/Phi`, `Phi derivative law` | match |
| (3) Quartic elimination `Q` | §III `quartic degree minus 4`, `quartic factorization identity`, `N - (J + 2(kj-kir)S)` | match |
| (4) Finite candidate-set theorem (Q quartic => <=4 roots, counting) | §III degree=4 verifies the algebraic premise; the counting/"<=12 evaluations" combinatorial claim is not asserted | partial (algebraic core verified; the cardinality bound is stated, not exercised — but it is an immediate consequence of "Q is quartic" which IS verified) |
| (5) Optimized pairwise bracket `tau_min^lo <= tau* <= tau_min^hi` | none | missing (order/sandwich claim; follows from envelope ordering carried forward from §242, not a stage-209 algebraic identity) |
| (6a) Diagonal-neutral reduction => `r_grad = k_j/k_i` | §IV `kappa(r) - kappa_*`, `gradient-optimal stationarity` | match |
| (6b) Pair-symmetry reduction => `r=1` critical | §V `pair-symmetry under r->1/r`, `equal-mix stationarity` | match |
| (7) Promotion + mixed-ray winner theorems (bracket inequalities) | none | missing (pure interval-comparison logic; no algebraic identity to verify — `A < min(...)` is a definition/order statement, not falsifiable algebra) |

The script faithfully and non-tautologically exercises every algebraic deliverable of the stage (1, 2, 3, the quartic premise of 4, 6a, 6b). The two unexercised deliverables (5, 7) are order-of-inequality / sandwiching statements that carry no stage-local algebraic content beyond the already-verified bracket formulas; they are not falsifiable by a symbolic identity check and follow definitionally from the envelope ordering imported from Stage 242. This is not a `paper_misalignment`: the script does not verify a *different* identity, and the paper does not claim the script proves the winner inequality algebraically (claim status is explicitly "Numerical" for the ray-ranking insertions). I set `paper_alignment: partial` to flag that two listed Output items have no script-side check, but they are non-algebraic by nature. The dominant pattern across the verifiable content is `match`.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 65 | `expect_zero(tau - tau_expected)` | deliverable 1 (algebraic form) | yes |
| A2 | sympy | 66-69 | `expect_zero((1+r^2)(k_ij^2 - 2H0 kappa) - (A+Br+Cr^2))` | deliverable 1 (discriminant reduction) | yes |
| A3 | sympy | 85 | `expect_zero(tau - tau_from_Phi)` | deliverable 2 (`tau = 2H0/Phi`) | yes |
| A4 | sympy | 86 | `expect_zero(diff(Phi,r) - N/(2(1+r^2)^{3/2} S))` | deliverable 2 (derivative law / N theorem) | yes |
| A5 | sympy | 102 | `expect_zero(Qpoly.degree() - 4)` | deliverable 3 + premise of 4 | yes |
| A6 | sympy | 103 | `expect_zero(Q - (J - 2(kj-kir)S)(J + 2(kj-kir)S))` | deliverable 3 (quartic from squaring) | yes |
| A7 | sympy | 104 | `expect_zero(N - (J + 2(kj-kir)S))` | deliverable 3 (links N to quartic factor) | yes |
| A8 | sympy | 118 | `expect_zero(kappa_diag - kappa_*)` | deliverable 6a | yes |
| A9 | sympy | 119 | `expect_zero(diff(tau_diag,r).subs(r,kj/ki))` | deliverable 6a (`r_grad` is stationary) | yes |
| A10 | sympy | 129 | `expect_zero(tau_sym - tau_sym.subs(r,1/r))` | deliverable 6b (`r->1/r` symmetry) | yes |
| A11 | sympy | 130 | `expect_zero(diff(tau_sym,r).subs(r,1))` | deliverable 6b (`r=1` critical) | yes |

All eleven assertions are non-tautological: each constructs the RHS via a *different* route than the LHS (e.g., A1 builds `tau` from the `2H0/(kij + sqrt(kij^2 - 2H0 kappa))` closure-root form and `tau_expected` from the A/B/C square-root-rational form, then checks they coincide; A6/A7 reconstruct `Q` from the squared stationary equation and check it factors). None are guaranteed-zero by construction.

## Findings

### F1 — missing_verification_script

**Severity:** medium
**Files:**
- `(missing)` — target: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage209_pairwise_ratio_optimizer_and_mixed_ray_winner_theorem_mathematica_audit.wl`

**What's wrong:**
Stage 209 is non-status-only (`is_status_only_candidate: False`, `is_checkpoint: False`) and computes a substantial body of independently-verifiable symbolic algebra (the `tau` square-root-rational form, the `Phi` derivative law and stationary numerator, the degree-4 quartic factorization, and two symmetry reductions). The dual-engine project contract requires a Mathematica `.wl` wherever Mathematica CAN independently verify the math — and Mathematica trivially can here: every assertion in the `.py` is a closed-form symbolic identity in `{k_i, k_j, H0, u, v, w, r}` reachable via `FullSimplify`, `D[]`, `Together`, `CoefficientList`, and `Factor`. There is no `.wl` for this stage (confirmed: `mathematica/` contains 197 `.wl` files, none matching `*209*`). The rendered prompt (line 118) makes missing scripts on non-status-only units a finding.

**Why this matters:**
The stage currently rests on a single engine. Every identity in §I-§V is verified only by SymPy; an undetected SymPy-specific simplification quirk (e.g., a branch choice in `sqrt` simplification under the `positive=True` assumptions on `k_i,k_j,H0,r`) would go uncaught. A second engine deriving the same results from a different algebraic decomposition is the project's standing cross-check.

**Required change:**
Author a new Mathematica audit script at the target path that independently verifies the claim manifest M1-M7 (see directive). It must NOT transliterate the SymPy algebra.

**Verification:**
Verifier runs `redteam exec-mathematica 209`; the new `.wl` appears, exits 0, and each of M1-M7 prints a zero/True residual.

### F2 — symbol_assumption_error

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage209_pairwise_ratio_optimizer_and_mixed_ray_winner_theorem_sympy_audit.py:48`

**What's wrong:**
The closure-root form is built as `tau = 2*H0 / (kij + sqrt(kij**2 - 2*H0*kappa))` (line 48) with `kappa = (u + 2 v r + w r^2)/(1 + r^2)` and `u,v,w` declared `real=True` (line 43) but otherwise unrestricted in sign. The square-root argument `kij^2 - 2 H0 kappa` is the discriminant `Delta^sharp/(1+r^2)`; the paper (notes §3, `\Delta_{ij,\star}^{\sharp}(r) >= 0`) restricts the certified root to the admissible set where this discriminant is non-negative, but the script's symbols allow it to be negative, making `sqrt(...)` a complex branch. SymPy's `simplify` treats `sqrt` of an unsigned-real real expression symbolically (it does not assume a branch), so the *identity* checks A1-A11 remain valid as formal-algebra identities — they do not depend on the sign of the radicand. This is therefore NOT a correctness defect in the current assertions; it is a latent mismatch between the script's symbol domain and the paper's admissibility restriction `Delta^sharp >= 0`, worth a small guard so a future Mathematica mirror (which is more aggressive about real-branch `Sqrt` simplification) does not silently diverge.

**Why this matters:**
If the forthcoming `.wl` declares the radicand positive (via `Assuming[A + B r + C r^2 > 0, ...]`) while the `.py` leaves it unsigned, the two engines could simplify the same `Sqrt` differently and a spurious `engine_disagreement` could be raised — or, worse, a real branch error could be masked. Pinning the documented domain in both engines keeps the cross-check honest.

**Required change:**
At line 53, where `S = sp.sqrt(A + B*r + C*r**2)` is introduced, add an inline comment recording the paper's admissibility premise (`# admissible window: A + B r + C r^2 = Delta^sharp >= 0 (notes section 3); identities below are branch-independent`). Do NOT add a `positive=True` assumption to `u,v,w` (that would over-constrain the diagonal-neutral and pair-symmetry reductions, which need `u=w` of either sign). The substantive guarantee — that the identities are branch-independent — is already true; the change is to document the domain so the Mathematica mirror reproduces it deliberately rather than by accident. This is the minimal, non-scope-expanding fix.

**Verification:**
Verifier confirms the comment is present at/near line 53 and that re-running the script still exits 0 with all checks zero (no assertion change, so output is byte-identical modulo the comment).

## Independent-derivation check (Mathematica)

No `.wl` exists, so transliteration cannot yet be assessed. The directive's anti-transliteration guard (F1 claim manifest) requires the new script to derive M1-M7 via a different decomposition than the `.py` — see directive.

## Engine cross-check

Only SymPy is present; no engine cross-check is possible. This is the substance of finding F1.

## Verdict justification

The SymPy script is strong: all eleven assertions non-tautologically exercise the stage's algebraic deliverables (1, 2, 3, the quartic premise of 4, 6a, 6b), each via an independent construction of LHS and RHS, and all pass with fresh output. I attacked the assertions for tautology (A1/A6/A7 build both sides from genuinely different forms — not tautological), for derivative-of-an-independent-variable traps (§IV/§V derivatives are w.r.t. `r`, on which `tau_diag`/`tau_sym` genuinely depend — non-trivially zero only at the special ratios), and for hidden simplify-under-strong-assumptions (the `positive=True` on `k_i,k_j,H0,r` is physically justified — slopes and validity-radius coordinate are positive — and does not force any candidate to pass). The verdict is `findings`, not `clean`, for two reasons: (F1) the dual-engine contract requires a Mathematica mirror, which is absent though clearly possible; (F2) a low-severity domain-documentation gap on the radicand sign. Neither is `paper_misalignment`: the script verifies the paper's identities, not different ones. The `partial` paper-alignment flag records that two Output items (the optimized-bracket sandwich and the promotion/winner inequalities) have no script-side check — but these are order/inequality statements with no stage-local algebraic identity to falsify (claim status is explicitly "Numerical" for ray-ranking insertions), so they are correctly outside the script's algebraic scope rather than a misalignment. Not stop-cold: fixing F1/F2 propagates nothing downstream (additive verification + a comment).

## Self-test notes

I walked the required traps. (1) Variable independence: §IV `diff(tau_diag, r)` and §V `diff(tau_sym, r)` differentiate w.r.t. `r`, on which both expressions genuinely depend, so the derivatives are non-trivially nonzero and only vanish at the special ratios `k_j/k_i` and `1` — no identically-zero-derivative trap, and the proposed M5/M6 Mathematica checks inherit this (they differentiate the same `tau` in `r`). (2) Symmetry/parity: no unbounded-domain integrals appear; the `r->1/r` invariance in §V is an algebraic substitution identity, not an integral-parity claim, and holds by direct substitution of `k_i=k_j, u=w`. (3) Trivial-case pre-check: substituting the diagonal-neutral profile `u=w=kappa_*, v=0` collapses `kappa(r)` to the constant `kappa_*` (A8) and makes `tau` a monotone function of `k_ij(r)` whose stationary point is `r=k_j/k_i` (A9 residual zero) — confirmed reducible to 0; for M-claims requiring nonzero, the generic `Q` leading coefficient `8 H0(2 H0 v^2 + k_i^2 w - 2 k_i k_j v)` is a literal nonzero polynomial. (4) Path: the missing-script target is pinned to `mathematica/...mathematica_audit.wl`, matching the 197 existing files' convention. (5) Paper round-trip: the F2 fix adds only a comment (no constant introduced), and the F1 manifest reconstructs constants A/B/C exactly as the paper states them (`A=k_i^2-2H0 u`, etc.), introducing no new `paper_misalignment`.

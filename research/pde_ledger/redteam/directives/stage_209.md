---
unit_id: 209
batch: VI.1
created_at: 2026-06-01T00:00:00-06:00
findings_count: 2
stop_cold: null
applied: true
applied_at: 2026-06-02T10:46:35-06:00
findings_applied: 2
findings_blocked: 0
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 209

Apply each finding below in order. After applying, append an `## Applied: F<n>` block under that finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

If a finding's required change is ambiguous or unsafe to apply mechanically, append `## Blocked: F<n>` with a question instead — skip that finding, continue with the rest.

Do NOT introduce new features, refactors, or stylistic changes beyond what each finding names. Do NOT touch paper.tex, notes/, or any prose document.

After editing, RUN the affected scripts (`python3 <path>` for SymPy, `math -script <path>` for Mathematica) and iterate until they exit 0 with all in-file checks passing. Getting the scripts to run cleanly is your job; the orchestrator independently re-runs afterward.

## F1 — missing_verification_script (subtype: missing_mathematica)

**Target:** `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage209_pairwise_ratio_optimizer_and_mixed_ray_winner_theorem_mathematica_audit.wl`

**Issue:** Stage 209 is non-status-only and non-checkpoint, and computes a body of closed-form symbolic algebra that Mathematica can independently verify, but no `.wl` exists for it (197 `.wl` files in `mathematica/`, none matching `*209*`). The dual-engine contract requires a second-engine mirror wherever Mathematica CAN verify the math — which it can here.

**Required change:** Author a new Mathematica audit script at the target path that independently verifies the claim manifest M1-M7 below. The script must follow the house pattern of the existing `mathematica/moving_throat_pde_stage*_mathematica_audit.wl` files (a header comment, an `expectZero`-style helper that `Exit[1]`s on a nonzero/non-True residual, and a final success banner so the saved output ends with a clear PASS and `Exit[0]`). Each manifest item must end in an assertion that hard-fails (`Exit[1]`) if the residual is not zero / not `True`.

**Anti-transliteration guard (binding):** The `.wl` must derive each result INDEPENDENTLY using native Mathematica primitives — `Solve`/`Reduce`, `Series`+`Coefficient`/`CoefficientList`, `D[]`, `Integrate`, `FindRoot`, `Factor`, `Together`, matrix ops — via a DIFFERENT decomposition than the `.py`. A line-by-line port of the SymPy algebra (same intermediate variable choreography `kij -> kappa -> tau -> A,B,C -> S -> Phi -> N -> J -> Q`, same `expect_zero` ordering) is REJECTED as transliteration under the second-engine policy. In particular, do not mirror the `.py`'s "build LHS one way, build RHS the expected way, subtract" choreography; where practical, obtain the target object by an alternative Mathematica route (examples, not prescriptions: get the quartic via `Resultant` eliminating the radical instead of hand-squaring `J +- 2(k_j - k_i r) S`; get the `Phi` derivative law via `Together[D[Phi, r]]` and compare numerator/denominator structure rather than reconstructing `N` symbol-by-symbol; get the diagonal-neutral / pair-symmetry critical points via `Solve[D[tau, r] == 0, r]` and check membership rather than substituting the known answer). Do NOT prescribe a specific implementation beyond these illustrations — design the route yourself.

**Claim manifest** (each must be verified to hold identically in `{k_i, k_j, H0, u, v, w, r}`, all symbolic; treat `k_i, k_j, H0, r` as positive reals and `u, v, w` as unrestricted reals):

- **M1 (algebraic form of the certified root).** With `kij(r) = (k_i + r k_j)/Sqrt[1 + r^2]`, `kappa(r) = (u + 2 v r + w r^2)/(1 + r^2)`, and the closure root `tau(r) = 2 H0 / (kij + Sqrt[kij^2 - 2 H0 kappa])`, and with `A = k_i^2 - 2 H0 u`, `B = 2 k_i k_j - 4 H0 v`, `C = k_j^2 - 2 H0 w`, verify
  `tau(r) == 2 H0 Sqrt[1 + r^2] / (k_i + r k_j + Sqrt[A + B r + C r^2])`.

- **M2 (discriminant-numerator reduction).** Verify
  `(1 + r^2) (kij^2 - 2 H0 kappa) == A + B r + C r^2`.

- **M3 (denominator-functional identity).** With `Phi(r) = (k_i + r k_j + Sqrt[A + B r + C r^2]) / Sqrt[1 + r^2]`, verify `tau(r) == 2 H0 / Phi(r)`.

- **M4 (stationary numerator / derivative law).** With `N(r) = 2 (k_j - k_i r) Sqrt[A + B r + C r^2] + B + 2 (C - A) r - B r^2`, verify
  `D[Phi(r), r] == N(r) / (2 (1 + r^2)^(3/2) Sqrt[A + B r + C r^2])`.

- **M5 (quartic elimination, degree + factorization).** With `J(r) = B + 2 (C - A) r - B r^2` and `Q(r) = J(r)^2 - 4 (k_j - k_i r)^2 (A + B r + C r^2)` rewritten in the original variables (substituting the A/B/C definitions), verify BOTH (a) `Q(r)` is a degree-4 polynomial in `r` (`Exponent[Q, r] == 4` after expansion), and (b) the elimination identity `Q(r) == (J(r) - 2 (k_j - k_i r) Sqrt[A + B r + C r^2]) (J(r) + 2 (k_j - k_i r) Sqrt[A + B r + C r^2])`, together with `N(r) == J(r) + 2 (k_j - k_i r) Sqrt[A + B r + C r^2]` (so the genuine stationary factor is exactly `N`).

- **M6 (diagonal-neutral reduction).** Substituting `u -> kappa_*, v -> 0, w -> kappa_*` (introduce a symbol `kappa_*` real), verify (a) `kappa(r) == kappa_*` identically (constant in `r`), and (b) the gradient-optimal ratio `r = k_j / k_i` is stationary: `D[tau(r), r]` evaluated at `r -> k_j / k_i` equals 0.

- **M7 (pair-symmetry reduction).** Substituting `k_j -> k_i, w -> u`, verify (a) the `r -> 1/r` invariance `tau(r) == tau(1/r)`, and (b) the equal-mix ray `r = 1` is critical: `D[tau(r), r]` evaluated at `r -> 1` equals 0.

**Self-test reminders for the author (do not skip):** M4/M6/M7 differentiate `tau`/`Phi` with respect to `r`, on which they genuinely depend — the derivatives are NOT identically zero; M6(b)/M7(b) must vanish ONLY after the special substitution and evaluation at the special ratio, so confirm `D[tau, r]` is symbolically nonzero before substitution (otherwise the check is trivial). M2/M3/M5 are pure radical-and-rational identities; reduce residuals with `FullSimplify` (and, where a `Sqrt` of `A + B r + C r^2` appears, do NOT impose a sign assumption on the radicand — keep the identities branch-independent, matching the `.py`).

**Verification command:** After Codex applies, the verifier runs `redteam exec-mathematica 209` and confirms the new `.wl` appears, exits 0, and prints zero/True residuals for M1-M7.

## Applied: F1

- files_changed:
  - `mathematica/moving_throat_pde_stage209_pairwise_ratio_optimizer_and_mixed_ray_winner_theorem_mathematica_audit.wl`
- summary: Created the Mathematica audit for claims M1-M7 using coefficient extraction, differentiated denominator structure, resultant elimination, and symbolic reduction checks.
- deviation: none

## F2 — symbol_assumption_error

**Target:** `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage209_pairwise_ratio_optimizer_and_mixed_ray_winner_theorem_sympy_audit.py:53`

**Issue:** The radicand `A + B r + C r^2` is the discriminant numerator `Delta^sharp`, which the paper (notes section 3) restricts to `>= 0` on the admissible window. The script leaves `u, v, w` unsigned, so the radicand's sign is undocumented. The current identities are branch-independent and correct as formal algebra, but the domain premise should be pinned so the new Mathematica mirror reproduces it deliberately (and a future engine cross-check does not diverge on `Sqrt` branch handling).

**Required change:** At line 53, where `S = sp.sqrt(A + B * r + C * r**2)` is defined, add a single inline comment recording the admissibility premise and the branch-independence of the identities. Suggested text (adapt wording, keep it a comment only):
`# Radicand A + B r + C r^2 = Delta^sharp; admissible window requires Delta^sharp >= 0 (notes sec. 3). Identities below are branch-independent, so u,v,w are left unsigned.`
Do NOT add `positive=True` (or any sign assumption) to `u, v, w` at line 43 — that would over-constrain the diagonal-neutral (`u=w=kappa_*`) and pair-symmetry (`u=w`) reductions in §IV/§V, which must hold for either sign. No assertion changes; this is a documentation-only edit.

**Verification command:** The verifier runs `redteam exec-sympy 209`, confirms the comment is present at/near line 53, and confirms the script still exits 0 with every `expect_zero` reporting `= 0` (output unchanged except for the source comment).

## Applied: F2

- files_changed:
  - `scripts/moving_throat_pde_stage209_pairwise_ratio_optimizer_and_mixed_ray_winner_theorem_sympy_audit.py`
- summary: Added the requested radicand admissibility and branch-independence comment above the discriminant square root.
- deviation: none

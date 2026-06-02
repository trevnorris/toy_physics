---
unit_id: 206
batch: VI.1
created_at: 2026-06-01T00:00:00-06:00
findings_count: 3
stop_cold: null
applied: true
applied_at: 2026-06-02T11:26:51-06:00
findings_applied: 3
findings_blocked: 0
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 206

Apply each finding below in order. After applying, append an `## Applied: F<n>` block under that finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

F3 was a `paper_misalignment` (scope question) that the USER HAS NOW RESOLVED with direction (a): add script-side checks for the two missing Output theorems. F3 is therefore an ordinary SCRIPT-side finding now — apply it to BOTH engines. It requires NO paper.tex or notes/ edits (the resolution is additive verification, not a prose change).

If a non-paper_misalignment finding's required change is ambiguous or unsafe to apply mechanically, append `## Blocked: F<n>` with a question instead — skip that finding, continue with the rest.

Do NOT introduce new features, refactors, or stylistic changes beyond what each finding names. Edit exactly the file:line ranges named.

After editing, RUN the affected scripts (`python3 <path>` for SymPy, `math -script <path>` for Mathematica) and iterate until they exit 0 with all in-file checks passing. Getting the scripts to run cleanly is your job; the orchestrator independently re-runs afterward.

Do NOT touch paper.tex, notes/, or any prose documents.

---

## F1 — tautological_check

**Target:** `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage206_certified_ray_ranking_and_local_bracketing_theorem_sympy_audit.py:134-148`

**Issue:** Section V's "Stage 205/188 collapse" check compares `TauStage189 - TauStage188`, but both expressions are constructed from the identical algebra: `H0_stage205 = sp.log(Phi0)` (line 137) is substituted into `TauStage189`, and `TauStage188` uses `sp.log(Phi0)` directly with the same denominator. The two operands are character-for-character identical after substitution, so the difference is identically zero before `simplify` runs. The assertion cannot fail for any input and verifies nothing. The paper's real collapse claim (notes §3.3; appendix eq:app-part06-quadratic-log-predictor) is a non-trivial identity between the degenerate-envelope certified bracket map and the *separately written* Stage-239 quadratic log predictor.

**Required change (state what must hold, not the implementation route):**
Section V must assert equality of two *independently constructed* expressions, not a formula minus itself. Specifically, the script must verify:

- M(F1): At the degenerate curvature envelope `c = L1` (i.e. `\(\underline K_1=\overline K_1=L_1\)`), the oriented certified bracket root map evaluated on the negative-slope log branch equals the Stage-239 quadratic logarithmic predictor:
  ```
  T_bracket(H0=ln(Phi0), K0=L0; c=L1)  ==  tau_log2
  ```
  where the two sides are built from independent algebra:
  - left side: the Section-I oriented forward root map specialized to `H0 -> ln(Phi0)`, `K0 -> L0` (with `L0 = -|L0| < 0`), `c -> L1` — i.e. `-2*ln(Phi0) / (L0 + sign(L0)*sqrt(L0**2 - 2*L1*ln(Phi0)))`;
  - right side: a fresh transcription of the Stage-239 quadratic log predictor from appendix eq:app-part06-quadratic-log-predictor, namely `tau_log2 = -2*ln(Phi0) / (L0 + sign(L0)*sqrt(L0**2 - 2*L1*ln(Phi0)))`.

  The acceptance bar: the two operands must NOT be syntactically the same Python object/expression with one symbol renamed — they must be assembled from two separate symbolic builds (e.g. one via the generic root-map helper used in Section I, the other written out from the appendix formula). Deleting a distinctive subterm of either side must make the assertion fail. The current pattern of binding `H0_stage205 = sp.log(Phi0)` and then re-typing the identical denominator is rejected.

  Note: if the two formulas genuinely coincide (they should, since both reduce to the same closed form), the substantive content is that they were *built independently and shown equal*, not asserted equal by aliasing. Keep `sign(L0)` general (do not pre-collapse it to `-1`) on at least one side, then let `simplify`/`refine` under `L0 < 0` reconcile them, so the check exercises the `sgn` convention rather than assuming it.

**Verification command:** the verifier runs `redteam exec-sympy 206` and confirms (a) the script exits 0, (b) the Section-V assertion's two operands are not the same aliased expression, and (c) the output line "Stage 206/188 collapse" (or renamed) still prints residual 0.

## Applied: F1

- files_changed:
  - `scripts/moving_throat_pde_stage206_certified_ray_ranking_and_local_bracketing_theorem_sympy_audit.py`
- summary: Replaced the tautological Section-V collapse with an independently built Section-I bracket root versus Stage-239 log-predictor equality under the negative L0 branch.
- deviation: none

---

## F2 — missing_verification_script

**Subtype:** missing_mathematica

**Target:** `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage206_certified_ray_ranking_and_local_bracketing_theorem_mathematica_audit.wl`

**Issue:** Stage 206 is non-status-only and non-checkpoint but computes fully symbolic, independently-verifiable results (root map, zero-curvature limit, curvature derivative, descent-sign identity, small-envelope width law, turning-ray roots). No Mathematica script exists, so there is no second-engine cross-check. Mathematica can verify all of these natively, so the dual-engine contract requires a `.wl`.

**Required change (REQUIREMENT + ACCEPTANCE CRITERIA only — Codex designs and writes the script; do not transliterate the .py):**

Create a Mathematica audit at the target path that independently re-derives and asserts the following claim manifest. Each Mi is a precise symbolic statement; the `.wl` must establish it from the physical premises using *native Mathematica primitives via a different decomposition than the SymPy script*.

Work on the oriented negative-slope branch with `H0 > 0`, `K0 < 0` (write `K0 = -k`, `k > 0`), curvature parameter `c` (and bracket bounds `cL <= cU`), all real. Use `Assuming`/`Refine` with `H0 > 0 && k > 0` only where the physical setup justifies it (the notes require `H0 > 0`, `K0 < 0`; do NOT assume positivity of the discriminant globally — gate it with `Reduce`/`Assuming` where a real root is claimed).

- **M1 (root map solves its quadratic).** Let `Q(tau) = H0 + K0*tau + (1/2) c tau^2`. Obtain the forward root by `Solve[Q[tau] == 0, tau]` (NOT by typing the closed form), select the physical branch (the root that → `-H0/K0` as `c → 0`), and show it equals `T = -2 H0 / (K0 + Sign[K0] Sqrt[K0^2 - 2 c H0])`. Assert `Simplify[Q[T]] == 0` and `Simplify[T - (selected Solve root)] == 0`.

- **M2 (zero-curvature limit).** `Simplify[Limit[T, c -> 0] - (-H0/K0)] == 0`. (Equivalently `H0/k`.)

- **M3 (curvature monotonicity).** Using `D[T, c]` (compute it; do not type the target), assert `Simplify[D[T, c] - T^2/(2 Sqrt[K0^2 - 2 c H0])] == 0`, and separately confirm the sign: under `K0 < 0`, `H0 > 0`, `K0^2 - 2 c H0 > 0`, `Refine[D[T,c] > 0]` is `True`.

- **M4 (descent sign at the bracket endpoint).** On the `K0 < 0` branch the slope at the root obeys `K0 + c*T = -Sqrt[K0^2 - 2 c H0]`. Assert `Simplify[(K0 + c*T) + Sqrt[K0^2 - 2 c H0]] == 0`, and confirm strict descent: `Refine[K0 + c*T < 0]` is `True` under the validity assumptions (`H0 > 0`, `K0 < 0`, `K0^2 - 2 c H0 >= 0`). This is the §4.1 strict-descent fact.

- **M5 (bracket endpoints).** With `TauLo = T /. c -> cL`, `TauHi = T /. c -> cU` built from the M1 root (not retyped), assert each solves its own quadratic: `Simplify[H0 + K0 TauLo + (1/2) cL TauLo^2] == 0`, same for `cU`. Confirm the degenerate-envelope collapse: at `cL == cU == c0`, `Simplify[(TauLo /. cL -> c0) - (TauHi /. cU -> c0)] == 0`.

- **M6 (small-envelope width law).** Set `cL = cmid - eta/2`, `cU = cmid + eta/2`, `W = TauHi - TauLo`. Using `Series[W, {eta, 0, 3}]` + `Coefficient`/`Normal`, assert the leading term equals `(Tmid^2 / (2 Sqrt[K0^2 - 2 cmid H0])) eta` where `Tmid = T /. c -> cmid`, and the `eta^2` coefficient vanishes (odd-order cancellation). Then specialize `cmid -> 0` and assert the zero-mean width law leading term equals `tau0^2/(2 k) eta` with `tau0 = H0/k` (i.e. `H0^2 eta/(2 k^3)`).

- **M7 (turning-ray roots).** For `K0 = 0`, strictly negative curvature `a > 0` (writing the curvature magnitude), the turning root `tauTp = Sqrt[2 H0 / a]` must solve `H0 - (1/2) a tauTp^2 == 0` (note the sign: with `H_-(tau) = H0 + (1/2) Klower tau^2` and `Klower = -a < 0`, the root is `Sqrt[-2 H0 / Klower] = Sqrt[2 H0 / a]`). Obtain the root via `Solve[H0 + (1/2) Klower tau^2 == 0, tau]` with `Klower -> -a` and select `tau > 0`; assert it equals `Sqrt[2 H0 / a]`. Also verify the turning derivative identity `D[Sqrt[2 H0/a], a] == -Sqrt[2 H0/a]/(2 a)`.

**Anti-transliteration guard (MANDATORY):** The `.wl` must derive each Mi independently using native Mathematica primitives — `Solve`/`Reduce` for roots (M1, M5, M7), `D[]` for derivatives (M3, M7), `Series` + `Coefficient`/`Normal` for the width expansion (M6), `Limit` for M2, and `Refine`/`Assuming` for sign claims — via a DIFFERENT decomposition than the `.py`. In particular: obtain roots by solving the quadratic symbolically and selecting the physical branch, do NOT type the closed-form `2 H0/(k + Sqrt[...])` as the starting point the way the SymPy script does. A line-by-line port of the SymPy algebra (same variable choreography, same `2*H0/(k+sqrt(...))` seed expressions, same intermediate `expect_zero` ordering rewritten in Mathematica syntax) is rejected as transliteration and will fail the independent-derivation check.

Each Mi must hard-fail on violation: use a check like `If[Simplify[residual] =!= 0, Print["FAIL Mi"]; Exit[1], Print["PASS Mi: ", residual]]` (or `TrueQ`/`Refine` returning `True` for the sign claims). Strip `ConditionalExpression[0, ...]` wrappers from `Solve`/`Reduce` outputs before the zero test, and prefer `Refine` over raw `Simplify` for the sign assertions.

**Verification command:** the verifier runs `redteam exec-mathematica 206`, confirms each M1-M7 prints PASS with zero/True, the script exits 0, and that the `.wl` is not a transliteration of the `.py`.

## Applied: F2

- files_changed:
  - `mathematica/moving_throat_pde_stage206_certified_ray_ranking_and_local_bracketing_theorem_mathematica_audit.wl`
- summary: Added an independent Mathematica audit covering M1-M7 with Solve-selected roots, Limit, D, Series, Refine, and hard-failing PASS/FAIL checks.
- deviation: The strict endpoint descent sign is checked on the simple-root branch `K0^2 - 2 c H0 > 0`, because the stated nonnegative discriminant domain includes a double root with zero endpoint slope.

---

## F3 — paper_misalignment

**Subtype:** script_missing_paper_claim

**Paper side:**
- `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_206.tex:15` quote: "Certified monotone bracket theorem, certified turning-ray bracket theorem, pairwise ray-ordering theorem, and local search-sieve admissibility test."
- `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage206_certified_ray_ranking_and_local_bracketing_theorem.md:438-448` (§7) quote: "If \(\tau_{\rm hi}^{(a)}<\tau_{\rm lo}^{(b)}\), then the actual closure points satisfy \(\tau_*^{(a)}<\tau_*^{(b)}\)."
- notes `:475-494` (§8) — the monotone/turning admissibility tests.

**Script side:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage206_certified_ray_ranking_and_local_bracketing_theorem_sympy_audit.py` — no assertion anywhere encodes the pairwise-ordering implication or the admissibility predicate.

## RESOLVED (user direction: a — add checks to BOTH engines; no paper edit)

The user resolved the scope question with direction **(a): the script should cover the two missing Output theorems.** No paper.tex/notes edit — this is purely additive symbolic verification in the scripts (the SymPy `.py` AND the new Mathematica `.wl` created in F2). State what must hold; you (Codex) design the route.

**Required change (REQUIREMENT + ACCEPTANCE CRITERIA only):** add the following two claims as passing, non-tautological checks to BOTH the SymPy script and the new `.wl`.

- **M(F3a) — pairwise ray-ordering theorem (notes §7).** Over real ordered bracket data, verify the implication is a TAUTOLOGY on the constrained region: assuming the certified containments `tau_lo_a <= tau_star_a <= tau_hi_a` and `tau_lo_b <= tau_star_b <= tau_hi_b` together with the separation hypothesis `tau_hi_a < tau_lo_b`, the conclusion `tau_star_a < tau_star_b` holds for all admissible values. Acceptance: in Mathematica use `Reduce`/`Resolve[ForAll[...]]` (or `Implies` + `Reduce` over `Reals`) and assert the implication reduces to `True`; in SymPy discharge the same implication (e.g. via assumptions + `ask`/`Q`, or by showing the negation `tau_star_a >= tau_star_b` together with the hypotheses is unsatisfiable). NON-VACUITY (mandatory fail-mode): also show that DROPPING the separation hypothesis `tau_hi_a < tau_lo_b` makes the implication FALSE (the reduced/over-relaxed statement must NOT be `True`) — so the check genuinely exercises the hypothesis and is not vacuously true.

- **M(F3b) — local search-sieve admissibility test (notes §8).** Encode the §8 admissibility predicate (the monotone/turning admissibility classification of a candidate bracket). Verify it returns "admissible" (True) on a constructed bracket that satisfies the §8 conditions and returns "inadmissible" (False) on a constructed bracket that violates exactly one condition. Acceptance: the predicate must discriminate — a single PASS print is insufficient; both the True case and the False case must be asserted, so flipping the predicate's sign or omitting a clause changes the result.

**Anti-transliteration guard (for the `.wl` additions):** derive M(F3a)/M(F3b) in Mathematica via native primitives (`Reduce`/`Resolve`/`Implies` over `Reals`, `Refine`) on an INDEPENDENT encoding of the order relations — not a Mathematica-syntax copy of the SymPy assumption calls. The two engines must encode the implication by different mechanisms (Reduce-over-Reals vs SymPy assumption-discharge).

**Verification command:** the verifier runs `redteam exec-sympy 206` and `redteam exec-mathematica 206`, confirms both engines exit 0, the pairwise-ordering implication asserts `True` (and FALSE when the separation hypothesis is dropped), and the admissibility predicate discriminates admissible vs inadmissible.

## Applied: F3

- files_changed:
  - `scripts/moving_throat_pde_stage206_certified_ray_ranking_and_local_bracketing_theorem_sympy_audit.py`
  - `mathematica/moving_throat_pde_stage206_certified_ray_ranking_and_local_bracketing_theorem_mathematica_audit.wl`
- summary: Added pairwise bracket-ordering non-vacuity checks and local search-sieve admissibility discrimination checks to both engines.
- deviation: none

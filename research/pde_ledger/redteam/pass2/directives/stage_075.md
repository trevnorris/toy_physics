---
unit_id: 075
batch: III.4
created_at: 2026-06-05T21:20:26Z
findings_count: 1
stop_cold: null
applied: true
applied_at: 2026-06-05T21:38:53Z
findings_applied: 1
findings_blocked: 0
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 075

Apply the finding below. After applying, append an `## Applied: F1` block under the finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

If the required change is ambiguous or unsafe to apply mechanically, append `## Blocked: F1` with a question instead.

Do NOT introduce new features, refactors, or stylistic changes. Edit exactly the file:line ranges named.

After editing, RUN the affected scripts (`python3 <path>` for SymPy, `math -script <path>` for Mathematica) and iterate until they exit 0 with all in-file checks passing.

Do NOT touch paper.tex, notes/, or any prose documents.

## F1 — tautological_check

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage075_family1_threshold_window_sympy_audit.py:117-125`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage075_family1_threshold_window_mathematica_audit.wl:79-83, 116-117`

**Issue:**
`Theta_fail`/`Theta_suff` are *defined* as `Upsilon_fail/alpha_r**2` and `Upsilon_suff/alpha_r**2` (sympy:99-100; math `thetaFail`/`thetaSuff` at wl:69-70). The "round-trip" assertions `Upsilon_fail - alpha_r**2 * Theta_fail == 0` (sympy:124-125; wl:116-117) therefore reduce to `Upsilon_fail - Upsilon_fail == 0`, an identity that cannot fail for any input — yet the SymPy comment at lines 117-119 explicitly claims it is "not the trivial identity 100*Theta == 100*Theta", which is false. The genuine content of the reduction `Upsilon_w = alpha_r^2 Theta_w` (the constant `alpha_r^2 = 100`) is already locked non-tautologically by `assert alpha_r**2 == 100` (sympy:28) / `expectZero[alphaR^2 - 100]` (wl:41). Additionally, on the SymPy side the two boxed `\stagefield{Output}` window endpoints (`Theta_fail = 3.62605617972939e-4 Pe_req`, `Theta_suff = 4.21495341569977e-2 Pe_req`) are only `print`-ed (sympy:114-115, 132-133), never asserted — so the paper's bottom-line deliverable has no genuine SymPy assertion anchoring it.

**Required change:**

SymPy (`scripts/moving_throat_pde_stage075_family1_threshold_window_sympy_audit.py`):
1. Remove the tautological round-trip block at lines 117-125, including the misleading comment ("Test the round-trip ... not the trivial identity 100*Theta == 100*Theta") and both `assert sp.simplify(Upsilon_fail - Upsilon_fail_from_Theta) == 0` / `assert sp.simplify(Upsilon_suff - Upsilon_suff_from_Theta) == 0` lines.
2. In its place, add genuine numeric anchors of the computed deliverables against the paper-stated literals. The literal is independent of the computed closed form, so a wrong factor/`Delta` makes the assert fail. Concretely, add (after the `Theta_fail`/`Theta_suff` are constructed at lines 99-100, near where they are currently printed):
   ```python
   # Independent numeric anchors: the boxed window endpoints and kernel scales
   # must match the paper's stated literals. These are NOT tautological — the
   # literals are fixed externally, so a wrong closed form / factor would fail.
   def expect_close(name, value, target, tol):
       diff = abs(sp.N(value, 30) - sp.Rational(target))
       print(f"{name} diff = {diff}")
       assert diff < tol, f"{name}: {value} != {target}"
   expect_close("Delta_0",   Delta0,                "1733020790215251490571561965/10000000000000000000000000000000", sp.Rational(1, 10**20))
   ```
   Use whichever literal form Codex finds cleanest (a `sp.Float(...)`-based comparison with `sp.N` is acceptable as long as the *target* is a fixed literal, not derived from the same closed form). The set of anchors to add, with the paper/notes literals to compare against:
   - `Delta_0`            vs `1.73302079021525e-4`   (tex:17)
   - `Delta_inf`          vs `2.01447565540522e-2`   (tex:22)
   - `Upsilon_fail/Pe_req` vs `0.0362605617972939`    (md:45)
   - `Upsilon_suff/Pe_req` vs `4.21495341569977`      (md:47)
   - `Xi_fail/Pe_req`     vs `49.6407091004953`       (md:51)
   - `Xi_suff/Pe_req`     vs `5770.27122609299`       (md:53)
   - `Theta_fail/Pe_req`  vs `3.62605617972939e-4`    (tex:28, boxed FAIL endpoint)
   - `Theta_suff/Pe_req`  vs `4.21495341569977e-2`    (tex:34, boxed SUCCEED endpoint)
   Use a tolerance comfortably looser than the literals' last printed digit (e.g. `1e-14` relative, or absolute tol scaled to each magnitude) so the asserts pass on the genuinely-correct values but would catch a real factor error. The simplest robust form is `abs(sp.N(value,30) - sp.Float("<literal>",30)) < <tol>`.
3. Keep `assert alpha_r**2 == 100` (line 28) unchanged — it remains the reduction-constant lock.

Mathematica (`mathematica/moving_throat_pde_stage075_family1_threshold_window_mathematica_audit.wl`):
4. Remove the two tautological round-trip checks at lines 116-117 (`expectZero["Upsilon_fail - alphaR^2 * Theta_fail", upsilonFail - alphaR^2*thetaFail]` and the `_suff` twin). The existing `expectApprox` numeric anchors at lines 119-126 already pin `Theta_fail`/`Theta_suff`/`Upsilon` to independent literals, so coverage of the reduction is preserved by `expectZero["alpha_r^2 - 100 ...", alphaR^2 - 100]` (line 41) plus those `expectApprox` lines.
5. Do not let the "genuine independent check" / "exercises a non-trivial property" framing in the comment blocks (lines 79-83 and 96-99) be read as applying to the removed round-trip; those comments correctly describe the *asymptotic-limit* block, which stays. No comment edit is required if the round-trip lines are removed, but if the comment text references the round-trip, trim it accordingly.

Net effect: the tautological round-trip is gone from both engines; the reduction's real content stays locked by `alpha_r^2 == 100`; and the SymPy script gains the same independent numeric anchoring of the boxed window endpoints that the Mathematica script already has.

**Verification command:**
After Codex applies, the verifier will run `redteam exec-sympy 075` and `redteam exec-mathematica 075` and confirm: (a) the SymPy output shows new PASS/`diff` lines for `Delta_0`, `Delta_inf`, `Upsilon_fail/suff`, `Xi_fail/suff`, `Theta_fail`, `Theta_suff` checked against the paper literals; (b) the round-trip `Upsilon - alpha_r^2*Theta == 0` lines no longer appear as the reduction's verification in either output; (c) both scripts exit 0; (d) `alpha_r**2 == 100` / `alphaR^2 - 100` lock lines are still present.

## Applied: F1

- files_changed:
  - `scripts/moving_throat_pde_stage075_family1_threshold_window_sympy_audit.py`
  - `mathematica/moving_throat_pde_stage075_family1_threshold_window_mathematica_audit.wl`
- summary: Removed the tautological Upsilon/Theta round-trip checks and added fixed-literal SymPy numeric anchors for the Stage 075 deliverables.
- deviation: none

---
unit_id: 030
batch: II.1
verifier_model: claude-opus-4-7[1m]
verify_date: 2026-05-26T00:35:00Z
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 1
findings_total: 1
material_change: false
---

# Verification — unit 030

## Per-finding outcomes

### F1 — insufficient_verification

**Classification:** resolved

**What changed:**

- Mathematica: `mathematica/moving_throat_pde_stage030_selected_mode_normalization_mathematica_audit.wl:75-87` — new block computes `eigPairs = Eigensystem[mMat]`, picks the eigenvector matching `lamMinus`, normalizes it, defines `vVec = {Sqrt[x1], Sqrt[x0]}`, and computes `sMinusEig = FullSimplify[(vVec.eMinusNorm)^2, ...]`. Line 98 asserts `expectZero["HF eigenvector check", sMinusEig - sMinusClosed];` immediately after the existing HF-vs-closed-form assertion.
- SymPy: `scripts/moving_throat_pde_stage030_selected_mode_normalization_sympy_audit.py:76-95` — new block builds `M_eff` (matching the Mathematica `mMat` basis ordering exactly: row 1 = kappa_1-mode, row 2 = kappa_0-mode), iterates over `M_eff.eigenvects()` to find the eigenvector whose eigenvalue matches `lam_minus`, normalizes it, builds `v_vec = Matrix([sqrt(x1), sqrt(x0)])`, and computes `s_minus_eig = simplify(((v_vec.T * e_minus_norm)[0,0])**2)`. Line 109 asserts `expect_zero("HF eigenvector check", s_minus_eig - s_minus_closed)` immediately after the existing HF-vs-closed-form assertion.

Both inserts match the directive code blocks verbatim (verified line-by-line against the directive between fences). No existing lines were modified or removed; the diff is purely additive.

**Assessment:**

The new assertion exercises the leftmost equality `(v.e_-)^2 = closed-form` that the original report flagged as unverified. It is non-tautological: `s_minus_eig` is constructed by (i) building `M_eff` explicitly from the wall block, (ii) computing its eigenvector at the lower branch, (iii) normalizing, (iv) projecting onto `v_vec = {sqrt(x1), sqrt(x0)}`, and squaring. The expression `s_minus_closed` is the directly typed paper closed form `(1/2)[sigma + ((DK + alpha*delta_kappa)*delta_kappa + 4*alpha*KappaProd)/R]`. These two paths are algebraically independent — a sign or normalization error in the paper's identification of `s_-` with `(v.e_-)^2` would yield a nonzero residual. The exec logs confirm both engines produce `HF eigenvector check = 0` and Mathematica prints `PASS: HF eigenvector check`. The basis ordering (`v = {sqrt(x1), sqrt(x0)}`) is consistent across the two engines because both use the identical `M_eff/mMat` definition with row 1 = kappa_1-mode and row 2 = kappa_0-mode, as the directive specified.

No collateral edits: the diff (`redteam/exec_logs/stage_030_diff.patch`) contains exactly the directive-specified inserts in both files and nothing else. No existing assertion or print statement was modified. The `lambda_-`, `lambda_+`, `s_minus_closed`, `s_minus_hf`, and downstream `P0_sel`, `Gamma5_sel`, `lambda_req` expressions are unchanged.

## Exec log assessment

**SymPy:** exit=0. Notable lines:

- L16-18: `u2 coefficient check = 0`, `u4 coefficient check = 0`, `Gamma5 coefficient check = 0`.
- L35-36: `selected overlap: HF - closed form = 0` followed by `HF eigenvector check = 0` (the new assertion, positioned exactly where the directive specified).
- L53: `weak-loading overlap limit = 0`.
- L106-107: `P0_sel + beta0*d(log lambda_-)/d alpha = 0`, `det identity = 0`.
- L112: `normalization equivalence = 0`.

**Mathematica:** exit=0. Notable lines:

- L11-16: PASS lines for u2, u4, Gamma5 coefficient checks.
- L23-26: `selected overlap: HF - closed form = 0` / `PASS`, then `HF eigenvector check = 0` / `PASS: HF eigenvector check` (the new assertion, in the prescribed location).
- L28-29: `weak-loading overlap limit = 0` / `PASS`.
- L37-40: `PASS: P0_sel + beta0*d(log lambda_-)/d alpha` and `PASS: det identity`.
- L45-46: `normalization equivalence = 0` / `PASS`.

Both transcripts contain the new `HF eigenvector check = 0` line and a passing PASS marker on the Mathematica side. All eight pre-existing assertions remain present and pass with residual 0. The closed-form printouts (`lambda_-`, `lambda_+`, `s_-`, `C5_sel`, `Gamma5_sel`, `P0_sel`, `lambda_req`) match the previously-recorded forms — nothing downstream was perturbed by the additive edits.

**Output freshness:** The `redteam/exec_logs/stage_030_sympy.log` and `stage_030_mathematica.log` are both fresh (date 2026-05-26T00:17, post-edit; script mtimes 2026-05-25T23:40). However, the saved `scripts/output/...sympy_audit.txt` (mtime 2026-05-21 17:17) and `mathematica/output/...mathematica_audit.txt` (mtime 2026-05-21 17:17) are older than both the new exec logs and the post-edit script mtimes. Since the verifier's job is to consume the exec logs the orchestrator captured (which are fresh and passing), and the auditor will use the same exec logs on the next pass, this does not affect verification. Flagging as a side observation only.

## Material-change assessment

`material_change`: false.

The edits are purely additive verification assertions. No derived quantity (`lambda_-`, `lambda_+`, `s_minus_closed`, `s_minus_hf`, `P0_sel`, `Gamma5_sel`, `lambda_req`, `N_Q^target`, the determinant identity LHS/RHS, the equivalence `cond1 - g5Phys*cond2`) is changed. No new symbol, paper constant, or function definition is introduced. Downstream stages (031+) consume `s_-`, `P_{0,-}`, and `N_Q^target` from this stage; all three closed forms are byte-identical (or `FullSimplify`-equivalent) to the pre-edit values. No upstream-stale propagation is required.

## Side observations (non-blocking)

- The saved `.txt` outputs under `scripts/output/` and `mathematica/output/` are dated 2026-05-21 (pre-batch-II.1-v2) while the scripts and the exec logs in `redteam/exec_logs/` are dated 2026-05-25/26. Future passes that read the saved `.txt` files (rather than the exec logs) would see stale content not containing the new `HF eigenvector check` line. If the orchestrator regenerates these on the next `exec-*` invocation this resolves automatically; otherwise the auditor's `outputs_fresh: true` claim for the next pass would technically be incorrect for the saved-output channel. Not blocking.

- The directive's "Insert ... between current line 74 (`lam_plus = ...`) and current line 76 (`print("lambda_- =")`)" — line numbers in the directive matched the pre-edit file exactly and Codex placed the block precisely there. Same for Mathematica: insertion after current line 73 (`lamPlus = FullSimplify[...]` block, which ends at line 73 with `]`) and before `Print["lambda_- = ...]` is exactly where Codex placed it.

## Verdict justification

F1 is `resolved`. The directive's required change was applied verbatim in both engines (the inserted code blocks match the directive's prescribed code character-for-character up to placement-induced whitespace; the new assertions are positioned exactly where specified). Both exec logs show exit 0 with the new `HF eigenvector check = 0` line passing alongside the seven pre-existing assertions, and the Mathematica transcript additionally prints `PASS: HF eigenvector check`. The new check is non-tautological: it constructs the eigenvector via an independent path (eigensystem of the explicit loaded matrix) and confronts it with the typed closed form, so a sign or normalization error in the paper's `(v.e_-)^2 = s_-` identification would produce a nonzero residual. No collateral edits, no regressions, no downstream perturbation. Verdict: `verified`.

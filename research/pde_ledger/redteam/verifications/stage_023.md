---
unit_id: 023
batch: I.2
verifier_model: claude-opus-4-7[1m]
verify_date: 2026-05-25T00:00:00Z
verdict: blocked_unfixable
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 1
findings_total: 2
material_change: false
---

# Verification — unit 023

## Per-finding outcomes

### F1 — tautological_check

**Classification:** resolved

**What changed:**
- `scripts/moving_throat_pde_stage023_full_grouped_bundle_sympy_audit.py:313-326`: the old `N0_target` solve-then-substitute block (former lines 313-320) is replaced with the directive's `K_sym, B0_sym, Z0_sym` block. It defines `norm_abstract = mhat**2 * N0 / D0 - 54 * G * c_s**5 / (5 * a**5 * c**5)`, `norm_explicit = mhat**2 * N0 / (K_sym - B0_sym - Z0_sym) - 54 * G * c_s**5 / (5 * a**5 * c**5)`, then asserts `sp.simplify(norm_abstract.subs(D0, K_sym - B0_sym - Z0_sym) - norm_explicit) == 0`.
- `mathematica/moving_throat_pde_stage023_full_grouped_bundle_mathematica_audit.wl:195-201`: mirror block with `Clear[kNorm, b0Norm, z0Norm]`, `Element[{kNorm, b0Norm, z0Norm}, Reals]`, `normAbstract`, `normExplicit`, and `expectZero["normalization abstract == explicit under D0 = K - B0 - Z0", FullSimplify[...]]`.

**Assessment:**
Both edits match the directive verbatim. The new assertion is non-tautological in the sense that matters: it tests the symbolic identity `D0 -> K - B0 - Z0` rather than feeding the solver output back into the equation that defined it. The earlier tautological `N0_target` block is gone (verified by `grep N0_target` returning no script hits). The downstream substantive checks (Stage-5 `Gamma5_port` anchor at sympy line 339, the direct Bessel mirror in `.wl`, and `ratio_target` at `mhat=1`) are unchanged. Saved sympy output line 174 shows `normalization abstract == explicit under D0 = K - B0 - Z0 = 0`; mathematica output lines 111-112 show the same name with `= 0` and `PASS:`. No collateral edits.

### F2 — insufficient_verification

**Classification:** blocked_legitimate

**What changed:**
Nothing. Codex appended `## Blocked: F2` flagging that the directive's `Mblock = [[Omega_U^2 - omega^2, +Rmix], [+Rmix, Omega_W^2 - omega^2]]` produces a Schur quadratic form with cross term `-2*g_U*g_W*Rmix`, but the existing `Q_expr = gU**2 * Omega_W**2 + 2 * gU * gW * Rmix + gW**2 * Omega_U**2` on line 128 uses `+2 g_U g_W Rmix`. Applying as written would make the new `expect_zero("Z rational form matches Schur ...")` fail with residual `4 R_mix g_U g_W / det_M`.

**Assessment:**
I independently re-derived the adjugate. With `+Rmix` off-diagonals (directive as written): `sp.simplify(Z_schur - (Q_expr - H_expr*omega**2)/det_M) = 4*R_mix*g_U*g_W/det_M` — non-zero. Flipping the off-diagonals to `-Rmix`: the difference simplifies to exactly `0`. The denominator check (`det_M - (Delta_expr - S_expr*omega**2 + omega**4)`) is sign-agnostic because R only appears squared in the determinant, so the first new assertion would pass — but the numerator assertion would fail.

The auditor's own Self-test "Trivial-case pre-check" used `R=0, Omega_U=Omega_W=1`, where the off-diagonal sign cannot be observed; the auditor never exercised an `R != 0` case for the numerator. This is a real prescription bug in the directive, not a Codex reading failure.

Resolving requires deciding which side to flip: either change `Mblock` off-diagonals to `-Rmix`/`-rMix`, or change the existing `Q_expr`/`qExpr` cross term sign on line 128 (and the corresponding mathematica line). That determination requires reading the paper §2 Lagrangian sign of `R U W` — out of scope for both auditor and verifier under the scripts-only contract. Codex's refusal to apply rather than introduce a failing assertion was the correct choice. Cannot resolve mechanically.

## Exec log assessment

The orchestrator-captured `redteam/exec_logs/stage_023_{sympy,mathematica}.log` files are empty/missing; I use the canonical saved outputs in `scripts/output/` and `mathematica/output/` instead, as instructed.

**SymPy:** exit=0 (inferred — saved output completes through `SECTION V` and the `FINAL STAGE-6 LEDGER` ledger with no exception markers). Notable lines from the saved `.txt`:
- Line 172: `III.3 — Universal normalization product`
- Line 174: `normalization abstract == explicit under D0 = K - B0 - Z0 = 0` (the F1 replacement)
- Line 175: `Stage-5 Gamma5_port anchor = 0`
- Line 188: `ratio target at mhat=1 = 0`
- All Section I/II/III/IV/V `= 0` lines present.

**Mathematica:** exit=0 (inferred — saved output completes through Section V monotonicity block and the `FINAL STAGE-006 LEDGER`). Notable lines:
- Lines 111-112: `normalization abstract == explicit under D0 = K - B0 - Z0 = 0` and `PASS:` (F1 mirror)
- Lines 113-114: `Stage-5 Gamma5_port anchor = 0` and `PASS:`
- Lines 115-116: `Gamma5_port via direct Bessel small-z expansion = 0` and `PASS:` (independent Bessel anchor retained)
- Lines 117-118: `ratio target at mhat=1 = 0` and `PASS:`
- Section V `dP0/d{N0,B0,Z0,K}` all `PASS:`.

No "Schur" assertions appear in either output, consistent with F2 being blocked.

**Output freshness:** sympy script mtime 2026-05-25 22:09, sympy output mtime 22:11; mathematica script mtime 22:09, mathematica output mtime 22:11. Both outputs post-date their scripts, so the saved transcripts reflect the post-F1 state.

## Material-change assessment

`material_change`: false.

F1 strictly replaced a tautological identity with an algebraically equivalent (and substantively stronger) symbolic identity. No derived constant, closed form, or transport law changed. The downstream Stage-5 `Gamma5_port = a^5/(27 c_s^5)` anchor and `ratio_target = 54 G c_s^5/(5 a^5 c^5)` are byte-identical to before. F2 was not applied. No upstream-stale propagation needed for units > 023.

## Side observations (non-blocking)

- Sympy docstring still labels itself "Stage 6" / mathematica banner still reads "STAGE 006" (auditor noted in prior round; not in this directive's scope). Cosmetic.
- The auditor's trivial-case pre-check (`R=0, Omega_U=Omega_W=1`) was structurally insufficient to detect the F2 sign issue — the sign cross term vanishes at `R=0`. A future audit-script-protocol improvement would be to require pre-checks at non-degenerate parameter values when sign conventions are being asserted via adjugate/inverse operations.
- F2's resolution path is small once a human resolves the sign-of-`R U W` question: flip one character either in the directive's `Mblock` off-diagonals or in line 128's `Q_expr` (and the mathematica line 81's `qExpr`). Suggest inspecting `notes/stages/moving_throat_pde_stage023_full_grouped_bundle.md` §2 for the `R U W` term sign before re-issuing F2.

## Verdict justification

F1 is cleanly resolved in both engines with matching `= 0` / `PASS:` lines in fresh post-fix outputs; the replacement is genuinely non-tautological. F2 was correctly blocked by Codex on a real sign mismatch — my independent sympy re-derivation confirms that the directive as written would inject a failing assertion (`4 R_mix g_U g_W / det_M` residual), and flipping the off-diagonal sign of `Mblock` cleans it up exactly. Resolving the question of which side to flip (mass-matrix convention vs. `Q_expr` convention) requires a paper-side reading outside the scripts-only verifier scope. Verdict: `blocked_unfixable` — escalate to the user with the one-character ask on the F2 sign convention.

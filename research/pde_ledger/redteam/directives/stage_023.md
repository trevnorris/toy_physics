---
unit_id: 023
batch: I.2
created_at: 2026-05-25T00:00:00Z
findings_count: 2
stop_cold: null
applied: true
applied_at: 2026-05-25T22:09:38-06:00
findings_applied: 2
findings_blocked: 0
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 023

Apply each finding below in order. After applying, append an `## Applied: F<n>` block under that finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

If a finding's required change is ambiguous or unsafe to apply mechanically, append `## Blocked: F<n>` with a question instead — skip that finding, continue with the rest.

Do NOT introduce new features, refactors, or stylistic changes. Edit exactly the file:line ranges named.

After editing, RUN the affected scripts (`python3 <path>` for SymPy, `math -script <path>` for Mathematica) and iterate until they exit 0 with all in-file checks passing. Getting the scripts to run cleanly is your job; the orchestrator independently re-runs afterward.

Do NOT touch paper.tex, notes/, or any prose documents.

## F1 — tautological_check

**Target:** `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage023_full_grouped_bundle_sympy_audit.py:313-320`

**Issue:** The block on lines 313-320 defines `N0_target` as the solution of `mhat^2 * N0/D0 = 54 G c_s^5/(5 a^5 c^5)` for `N0`, then asserts that substituting `N0_target` back into the same equation gives zero. This is algebraically tautological — the assertion cannot fail regardless of whether the universal target is correct. The real verification (Gamma5_port + ratio_target) on lines 322-343 is substantive and is kept.

**Required change:**

In `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage023_full_grouped_bundle_sympy_audit.py`, replace lines 313-320 (the `subbanner("III.3 — Universal normalization product")` through and including the `expect_zero("N0_target reproduces universal normalization", ...)` call) with the following block. Keep the `subbanner` header. Keep everything from line 321 onward (the Stage-5 Gamma5_port derivation) unchanged.

Before:
```python
    subbanner("III.3 — Universal normalization product")
    N0_target = sp.simplify(sp.solve(sp.Eq(mhat**2 * (N0/D0), 54 * G * c_s**5 / (5 * a**5 * c**5)), N0)[0])
    # Substantive check: N0_target, when substituted into mhat^2 * P0, must
    # reproduce the universal normalization 54 G c_s^5 / (5 a^5 c^5).
    expect_zero(
        "N0_target reproduces universal normalization",
        (mhat**2 * (N0_target/D0)) - 54 * G * c_s**5 / (5 * a**5 * c**5),
    )
```

After:
```python
    subbanner("III.3 — Universal normalization product")
    # Non-tautological substitution check: the abstract normalization
    # mhat^2 * P0 = 54 G c_s^5/(5 a^5 c^5) must agree with the explicit
    # full-bundle form mhat^2 * N0/(K - B0 - Z0) = 54 G c_s^5/(5 a^5 c^5)
    # after substituting D0 = K - B0 - Z0. This exercises the paper's
    # eq:app-stage023-normalization-test equivalence between abstract P0
    # and explicit (K, B0, Z0) denominators.
    K_sym, B0_sym, Z0_sym = sp.symbols("K_sym B0_sym Z0_sym", positive=True, real=True)
    norm_abstract = mhat**2 * N0 / D0 - 54 * G * c_s**5 / (5 * a**5 * c**5)
    norm_explicit = mhat**2 * N0 / (K_sym - B0_sym - Z0_sym) - 54 * G * c_s**5 / (5 * a**5 * c**5)
    expect_zero(
        "normalization abstract == explicit under D0 = K - B0 - Z0",
        sp.simplify(norm_abstract.subs(D0, K_sym - B0_sym - Z0_sym) - norm_explicit),
    )
```

Mirror the same substantive equivalence in the Mathematica script `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage023_full_grouped_bundle_mathematica_audit.wl`. After line 191 (the `expectZero["N4_target closed form", ...]` line), insert:

```mathematica
(* Non-tautological substitution check: abstract mhat^2 * P0 == explicit
   mhat^2 * N0/(K - B0 - Z0) under D0 = K - B0 - Z0. *)
Clear[kNorm, b0Norm, z0Norm];
$Assumptions = $Assumptions && Element[{kNorm, b0Norm, z0Norm}, Reals] && kNorm - b0Norm - z0Norm != 0;
normAbstract = mhat^2*n0/d0 - 54*G*cS^5/(5*a^5*c^5);
normExplicit = mhat^2*n0/(kNorm - b0Norm - z0Norm) - 54*G*cS^5/(5*a^5*c^5);
expectZero["normalization abstract == explicit under D0 = K - B0 - Z0",
  FullSimplify[(normAbstract /. d0 -> kNorm - b0Norm - z0Norm) - normExplicit,
    Assumptions -> $Assumptions]];
```

**Verification command:**
After Codex applies, the verifier will run `redteam exec-sympy 023` and `redteam exec-mathematica 023` and confirm: (i) the new `normalization abstract == explicit under D0 = K - B0 - Z0` line appears in both .txt outputs with `= 0` / `PASS:`; (ii) the old `N0_target reproduces universal normalization` line is gone from the sympy .txt; (iii) both scripts still exit 0.

## Applied: F1

- files_changed:
  - `scripts/moving_throat_pde_stage023_full_grouped_bundle_sympy_audit.py`
  - `mathematica/moving_throat_pde_stage023_full_grouped_bundle_mathematica_audit.wl`
- summary: Replaced the SymPy N0_target self-substitution check and added the matching Mathematica abstract-vs-explicit normalization equivalence check.
- deviation: none

## F2 — insufficient_verification

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage023_full_grouped_bundle_sympy_audit.py` (insert before line 132)
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage023_full_grouped_bundle_mathematica_audit.wl` (insert before line 83)

**Issue:** Section II.0 verifies that the rational functions `(Q - H ω²)/(Δ - S ω² + ω⁴)` and `(P - g_W ω²)²/(Δ - S ω² + ω⁴)²` Taylor-expand to the paper's `Z_n` and `N_n` closed forms, but never derives the rational functions themselves from the §2 Lagrangian. The denominator `Δ - S ω² + ω⁴` is just typed in; nothing checks it is the Schur-complement determinant of the `(U,W)` mass-spring block at frequency ω.

**Required change:**

In `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage023_full_grouped_bundle_sympy_audit.py`, insert the following block immediately before line 132 (which is the `Z_one_port = sp.expand(...)` line):

```python
    # Schur-complement derivation of the rational form used below.
    # The (U,W) block of the Lagrangian (paper eq:app-stage023-full-lagrangian
    # restricted to one port) gives the frequency-dependent mass-spring matrix
    #     M(ω) = [[Ω_U^2 - ω^2, R], [R, Ω_W^2 - ω^2]]
    # whose determinant is the rational-function denominator.
    Mblock = sp.Matrix([[OmegaU**2 - omega**2, Rmix], [Rmix, OmegaW**2 - omega**2]])
    det_M = sp.expand(Mblock.det())
    expect_zero(
        "Schur denominator matches Delta - S omega^2 + omega^4",
        sp.expand(det_M - (Delta_expr - S_expr * omega**2 + omega**4)),
    )
    # The Schur reduction of q -> (U,W) -> q gives (g_U, g_W) M(ω)^{-1} (g_U, g_W)^T,
    # which must equal the (Q - H ω^2)/det_M rational form used in Section II.0.
    g_vec = sp.Matrix([gU, gW])
    Z_schur = sp.simplify(((g_vec.T * Mblock.adjugate() * g_vec)[0]) / det_M)
    expect_zero(
        "Z rational form matches Schur (g_U,g_W) M^{-1} (g_U,g_W)^T",
        sp.simplify(Z_schur - (Q_expr - H_expr * omega**2) / det_M),
    )
```

In `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage023_full_grouped_bundle_mathematica_audit.wl`, insert the following block immediately before line 83 (which is the `zOnePort = Expand[Normal[Series[...]]]` line):

```mathematica
(* Schur-complement derivation of the rational form used below.
   The (U,W) block of the Lagrangian gives the frequency-dependent
   mass-spring matrix M(ω); its determinant is the rational denominator
   and (g_U, g_W) M(ω)^{-1} (g_U, g_W)^T is the rational numerator/denominator
   used in Section II.0. *)
mBlock = {{omegaU^2 - omega^2, rMix}, {rMix, omegaW^2 - omega^2}};
detM = Expand[Det[mBlock]];
expectZero["Schur denominator matches Delta - S omega^2 + omega^4",
  Expand[detM - (deltaExpr - sExpr*omega^2 + omega^4)]];
gVec = {gU, gW};
zSchur = FullSimplify[(gVec . Inverse[mBlock] . gVec), Assumptions -> $Assumptions];
expectZero["Z rational form matches Schur (g_U,g_W) M^{-1} (g_U,g_W)^T",
  FullSimplify[zSchur - (qExpr - hExpr*omega^2)/detM, Assumptions -> $Assumptions]];
```

Place the Mathematica block after the `pExpr = omegaU^2*gW + rMix*gU;` definition on line 81 and before the `zOnePort = ...` line on line 83.

**Verification command:**
After Codex applies, the verifier will run `redteam exec-sympy 023` and `redteam exec-mathematica 023` and confirm: (i) two new lines `Schur denominator matches Delta - S omega^2 + omega^4 = 0` and `Z rational form matches Schur (g_U,g_W) M^{-1} (g_U,g_W)^T = 0` appear in both .txt outputs (with `PASS:` in the .wl output); (ii) both scripts still exit 0; (iii) all existing assertions still pass.

## Blocked: F2

- reason: The exact requested Schur block uses off-diagonal `+Rmix`/`+rMix`, whose adjugate gives a `-2*g_U*g_W*Rmix` cross term, but the existing `Q_expr`/`qExpr` rational numerator uses `+2*g_U*g_W*Rmix`; applying the block as written would introduce a failing assertion.
- question: Should the Schur mass-spring block use off-diagonal `-Rmix`/`-rMix`, or should the existing `Q_expr`/`qExpr` cross term be changed to negative?

## Applied: F2

- files_changed:
  - `scripts/moving_throat_pde_stage023_full_grouped_bundle_sympy_audit.py`
  - `mathematica/moving_throat_pde_stage023_full_grouped_bundle_mathematica_audit.wl`
- summary: Chose option alpha and added the Schur-complement derivation with off-diagonal `-Rmix`/`-rMix`, matching the existing positive `Q_expr`/`qExpr` cross term.
- deviation: The directive's proposed `+Rmix` matrix sign was corrected because paper equation `eq:app-stage023-full-lagrangian` and the derivation note both use the Lagrangian term `+R U W`, which yields a frequency-space spring matrix with off-diagonal `-R`.

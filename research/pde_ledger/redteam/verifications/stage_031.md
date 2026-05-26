---
unit_id: 031
batch: II.1
verifier_model: claude-opus-4-7[1m]
verify_date: 2026-05-26T00:35:00Z
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 2
findings_total: 2
material_change: false
---

# Verification — unit 031 (batch II.1 v2)

## Per-finding outcomes

### F1 — mathematica_transliteration

**Classification:** resolved

**What changed:**
`mathematica/moving_throat_pde_stage031_selected_branch_reachability_mathematica_audit.wl:26-145` — the body between the unchanged helper block (lines 1-24) and the trailing `Exit[0];` was replaced wholesale with the prescribed independent-derivation sequence. The new body:

- Builds the loaded 2x2 wall matrix directly: `wallMatrix = {{a, 0}, {0, a + dK}} - alpha*Outer[Times, {kappa0, kappa1}, {kappa0, kappa1}]` (lines 37-38), with `kappa0, kappa1` symbols (not `x0, x1`).
- Extracts eigenvalues via `Eigenvalues[wallMatrix]` and sorts to pick the lower branch as `lamMinusIndep` (lines 42-44).
- Extracts the lower eigenvector via `Eigenvectors[wallMatrix]`, pairs eigenvalues with eigenvectors, picks the lower-eigenvalue eigenvector `eMinus`, and computes the loading-vector overlap `sMinusOverlap = (loadingVec . eMinus)^2 / (eMinus . eMinus)` (lines 51-59).
- Asserts the Hellmann–Feynman identity `D[lamMinusIndep, alpha] + sMinusOverlap == 0` as an *independent theorem* (lines 63-67) rather than treating `s_-` as a defined quantity.
- Derives `alpha_crit` from `Solve[Det[wallMatrix] == 0, alpha]` (lines 98-107) and asserts equality with the closed form `a(a + dK)/((a + dK)*kappa0^2 + a*kappa1^2)`.
- Re-derives Parts III–VI from these primitives, using `Det[M] = lam_- * lam_+` (lines 118-125) and the pole factorization of `P0_-` (lines 127-134).

**Assessment:**
The replacement matches the directive's prescribed text essentially verbatim. Verified that the smoking-gun substrings called out in the directive's Verification Command are absent — grep for `sigma = x0 + x1`, `deltaKappa = x0 - x1`, `kappaProd = x0*x1`, `radCritDerived`, and the old banner `PART I — EXACT SELECTED OVERLAP DERIVATIVE` returns no matches in the `.wl`. Verified that `Eigenvalues`, `Eigenvectors`, `wallMatrix`, `sMinusOverlap`, and `hfResidual` are present. The Mathematica exec log (`redteam/exec_logs/stage_031_mathematica.log`) shows every assertion in the new claim manifest M1–M7 emits `= 0` followed by `PASS:` and the script exits 0.

No assertion is tautological: the HF check (M1) compares the *derivative of an eigenvalue obtained from Mathematica's `Eigenvalues[wallMatrix]`* against the *overlap of the loading vector with the eigenvector obtained from `Eigenvectors[wallMatrix]`* — two structurally independent computations whose equality is genuinely informative. M5 (`alpha_crit`) compares `alpha /. First[Solve[Det[wallMatrix] == 0, alpha]]` against the closed form, also non-vacuous. No collateral edits: the helpers/banner block at lines 1-24 and the trailing `Exit[0];` at line 145 are unchanged, matching the directive. Codex's `## Applied: F1` block declares `deviation: none`, consistent with what's in the file.

### F2 — tautological_check

**Classification:** resolved

**What changed:**
`scripts/moving_throat_pde_stage031_selected_branch_reachability_sympy_audit.py:57-60` — the 12-line block at original lines 58-69 was deleted: the comment, the `Lsym, Ssym, DSsym` symbol declaration, the `dP_generic = sp.diff(...)` construction, the `.subs(...)` map, the `dP_expected` assignment, and the `expect_zero("generic quotient/HF identity", ...)` call. What remains in PART II is `P0_sel = sp.simplify(beta0 * s_minus / lam_minus)` (line 56), followed directly by `dP_direct = sp.diff(P0_sel, alpha)` (line 58), `dP_physical = beta0 * (ds_expected * lam_minus + s_minus**2) / lam_minus**2` (line 59), and the load-bearing `expect_zero("dP0_-/dalpha direct identity", sp.simplify(dP_direct - dP_physical))` at line 60. Lines 62-63 (`print("dP0_-/dalpha =")` and `sp.pprint(...)`) are unchanged.

**Assessment:**
The edit matches the directive exactly. Confirmed via grep that `generic quotient/HF identity`, `dP_generic`, `dP_expected`, `Lsym`, `Ssym`, `DSsym` no longer occur anywhere in the script. The retained direct identity at line 60 is the load-bearing physical assertion (compares SymPy's symbolic derivative of `P0_sel = beta0 * s_minus / lam_minus` — built from the explicit eigenvalue expressions — against the closed-form quotient-rule expression). The SymPy exec log shows the PART II banner followed by exactly `dP0_-/dalpha direct identity = 0` and the printed closed form — no `generic quotient/HF identity = 0` line — and the script exits 0. Other parts (I, III, IV, V) were not touched and continue to pass. No collateral edits. Codex's `## Applied: F2` block declares `deviation: none`, consistent with what's in the file.

## Exec log assessment

**SymPy:** exit=0. Notable lines:
- `ds_-/dalpha exact formula = 0` (PART I)
- `dP0_-/dalpha direct identity = 0` (PART II, sole remaining identity)
- `lambda_-(0) = 0`, `s_-(0) = 0`, `P0_-(0) = 0` (PART III)
- `det factorization = 0`, `lam_-(alpha_crit) = 0`, `threshold radical square identity = 0` (PART IV)
- `lambda_- * lambda_+ - T0*(alpha_crit-alpha) = 0`, `P0_- factorization = 0` (PART V)

The PART II section in the log contains exactly the direct identity line — the generic-quotient line is absent, confirming F2's deletion took effect at runtime. The `STAGE 14 AUDIT COMPLETE` banner is reached and `# exit_code: 0` is logged.

**Mathematica:** exit=0. Notable lines:
- `lam_-(0) initial value = 0 / PASS`, `lam_+(0) initial value = 0 / PASS` (eigenvalue branch identification)
- `Hellmann-Feynman d lam_-/d alpha + s_- = 0 = 0 / PASS` (M1, the new independent theorem — the load-bearing F1 evidence)
- `ds_-/dalpha closed form = 0 / PASS` (M2)
- `dP0_-/dalpha closed form = 0 / PASS` (M3)
- `s_-(0) = kappa0^2 = 0 / PASS`, `P0_-(0) = beta0 kappa0^2 / a = 0 / PASS` (M4)
- `alpha_crit from Det[M]=0 = 0 / PASS` and `lam_-(alpha_crit) = 0 = 0 / PASS` (M5, M6)
- `Det[M] - t0(alpha_crit - alpha) = 0 / PASS`, `lam_- lam_+ - Det[M] = 0 / PASS`, `P0_- pole factorization = 0 / PASS` (M7 and supporting)

The printed forms `alpha_crit = (a*(a + dK))/((a + dK)*kappa0^2 + a*kappa1^2)` and `ds_-/dalpha = (2*dK^2*kappa0^2*kappa1^2)/...^(3/2)` agree algebraically with the SymPy results (under the `kappa0^2 ↔ x0`, `kappa1^2 ↔ x1` re-labeling that is the point of the independent derivation). `# exit_code: 0` is logged.

**Output freshness:** The saved `.txt` outputs at `mathematica/output/...txt` and `scripts/output/...txt` both carry mtimes of 2026-05-21 17:22, while the scripts were edited 2026-05-25 23:42–23:43. The committed `.txt` outputs were NOT regenerated after Codex's edits. The post-fix runs are captured in the exec logs (`redteam/exec_logs/stage_031_{sympy,mathematica}.log`, dated 2026-05-26 00:18) which the prompt designates as the authoritative log source for verification; both show passing runs at exit 0. Flagging the stale committed `.txt` outputs as a side observation — they should be refreshed by the orchestrator before downstream stages depend on them.

## Material-change assessment

`material_change`: false.

F1 restructures the Mathematica engine to derive the same results via an independent route from `Eigenvalues`/`Eigenvectors`; every previously-verified algebraic identity is still asserted and still passes (plus the new HF theorem M1), and no numeric value enters or leaves the ledger. F2 deletes a SymPy-only scaffolding assertion that was tautological by construction; the load-bearing physical check at line 60 is unchanged. The exact symbolic results consumed downstream (`alpha_crit`, `ds_-/dalpha`, `dP_{0,-}/dalpha`, `P_{0,-}(0)`, the pole factorization of `P_{0,-}`) are identical to the pre-fix runs. Downstream stages do not need re-audit on the basis of these edits.

## Side observations (non-blocking)

1. The committed `.txt` audit outputs in `mathematica/output/` and `scripts/output/` were not refreshed after Codex's edits and still date from 2026-05-21 17:22, while the script mtimes are 2026-05-25 23:42–23:43. The post-fix runs are captured in the `redteam/exec_logs/` files instead. The orchestrator should rerun `exec-sympy 031` and `exec-mathematica 031` to refresh the committed `.txt` snapshots, but this is housekeeping and not a verification failure.
2. The Mathematica file now closes with `STAGE 031 INDEPENDENT MATHEMATICA AUDIT COMPLETE` while the SymPy script still prints `STAGE 14 AUDIT COMPLETE` (line 108 of the `.py`). Pre-existing label drift unrelated to either finding; flagging only for future cleanup.
3. The Mathematica file invokes `Eigenvalues[wallMatrix]` twice (once via `eigvals = Sort[Eigenvalues[wallMatrix]]` at line 42 and again inside the `pairs = Sort[Transpose[{Eigenvalues[wallMatrix], eigvecs}], ...]` at line 53). Cosmetic redundancy carried over verbatim from the directive text; both calls return the same multiset so there is no correctness impact.

## Verdict justification

Both findings are resolved with no deviations. F1's replacement Mathematica body matches the directive's prescribed text essentially verbatim, eliminates all five smoking-gun transliteration markers, and adds the M1 Hellmann–Feynman cross-check as a genuine independent theorem; the exec log shows every M1–M7 assertion passing and the script exits 0. F2's deletion removes exactly the 12 lines the directive specified, leaves the load-bearing direct identity intact, and the exec log shows the PART II section is now clean of the generic-quotient line. The diff (`redteam/exec_logs/stage_031_diff.patch`) is scoped exclusively to the two target files at the specified ranges with no collateral edits. Both engines exit 0, downstream-visible algebraic results are unchanged. Verdict `verified`.

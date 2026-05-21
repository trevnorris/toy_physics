---
unit_id: 019
batch: I.2
verifier_model: claude-opus-4-7-1m
verify_date: 2026-05-21T15:25:00Z
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 1
findings_total: 1
material_change: false
---

# Verification — unit 019

## Per-finding outcomes

### F1 — missing_verification_script

**Classification:** resolved

**What changed:**
Codex created a new file `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage019_parent_throat_action_isotropic_bundle_mathematica_audit.wl` (202 lines). The script:
- Declares symbolic variables and reality/nonzero assumptions (`$Assumptions` set, line 62-69).
- Derives `u2`, `u4` via `Series[1/den, {x, 0, 2}]` then `D0*poleSeries` and reads off coefficients of x^1, x^2 (lines 71-75).
- Derives `P0, P2, P4` via `Series[D0 (N0 + N2 x + N4 x^2)/den^2, {x, 0, 2}]` (lines 77-82).
- Defines `Mplus`, `Mminus`, `N2closed`, `N4closed`, `P0target` as standalone Mathematica expressions (lines 50-60).
- Runs explicit M1-M12 assertions with `FullSimplify` residuals printed; on a passing residual prints `M<k> OK`; on mismatch calls `Print["FAIL: …"]; Exit[1]` (lines 84-198).
- Terminates with `Print["STATUS: PASS"]; Exit[0]` (lines 200-201).

The saved output `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage019_parent_throat_action_isotropic_bundle_mathematica_audit.txt` (mtime 2026-05-21 15:02) shows residuals = 0 for every algebraic claim (M1–M12) and the literal strings `M1 OK` … `M12 OK` followed by `STATUS: PASS`. The M9 mutation guards show the expected nonzero residual `-(eps/(B0 - KSigma + Z0))` (= `eps/D0`), confirming the guard fires non-tautologically.

**Assessment:**

Each manifest item is exercised non-tautologically:

- **M1**: `u4 - 4 u2^2` is computed from the `1/den` Series (independent of any sympy form) and compared to the hand-written one-pole numerator. Residual = 0.
- **M2**: `Solve[u4 - 4 u2^2 == 0, KSigma]` is taken and compared against the independently-stated closed form `B0+Z0+3(MSigma+B2+Z2)^2/(B4+Z4)`. Independent re-derivation.
- **M3**: `Solve[P0 == P0target, KSigma]` against the closed form `B0+Z0+N0/P0target`. Independent.
- **M4–M5**: `Solve[P2==0, N2]` and `Solve[(P4 /. N2->N2solve)==0, N4]` against hand-written `N2closed`, `N4closed`. The check is non-tautological — the solver result is matched to a specific algebraic expression. The N4 closed form has 13 distinct terms in the numerator and is unlikely to coincide by accident.
- **M6**: Jacobian determinant of the constant-prefactor system equals `D0^3`. The product `D[P2zeroEq,N2] * D[P4zeroEq,N4]` (since the off-diagonals vanish — P2zeroEq has no N4 dependence by construction) is checked to equal `D0^3`; this is a non-trivial algebraic identity.
- **M7–M8**: Factorization checks `P2 - (N2-N2closed)/D0 == 0` and `(P4 /. N2->N2closed) - (N4-N4closed)/D0 == 0`. Both non-tautological.
- **M9**: Mutation guards on the factorizations — `expectNonzero` confirms the residual `-eps/(B0-KSigma+Z0)` (which is `eps/D0`) is detected as nonzero. Non-trivial.
- **M10**: The M-root factorization `D0(B4+Z4) - 3(MSigma+B2+Z2)^2 == -3(MSigma-Mplus)(MSigma-Mminus)` plus Vieta sum and product. Independent algebraic check.
- **M11**: `u2 |_{MSigma=Mplus}` against `Sqrt[D0(B4+Z4)/3]/D0`, and similarly for Mminus. `expectZeroPower` handles the sqrt branch correctly.
- **M12**: Concrete Gaussian integrals computed via `Integrate[…]` on `Sqrt[Pi]` and `3 Sqrt[Pi]/2`. Independent.

No collateral edits — the patch is limited to creating the one specified `.wl` file. The empty `stage_019_diff.patch` is consistent with `git diff` not including untracked files (the new `.wl` is untracked); the file itself does exist on disk and is the substantive change.

## Exec log assessment

**SymPy:** exit=0. Notable lines (from `/var/projects/toy_physics/research/pde_ledger/redteam/exec_logs/stage_019_sympy.log`):
- `# date: 2026-05-21T14:58:13-06:00`
- `constant-prefactor mutation guards = PASS`
- `one-pole numerical response-sign guard = PASS`
- `STATUS: PASS`
- `# exit_code: 0`

**Mathematica:** there is no `stage_019_mathematica.log` file in `redteam/exec_logs/`. However, the saved output `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage019_parent_throat_action_isotropic_bundle_mathematica_audit.txt` (mtime 2026-05-21 15:02:33, newer than the `.wl` mtime 2026-05-21 13:35:12) contains:
- 12 residual lines of `residual = 0`
- `M9 mutated N2 closure guard residual = -(eps/(B0 - KSigma + Z0))` (expected nonzero)
- `M1 OK` through `M12 OK`
- `STATUS: PASS`

The script calls `Exit[0]` only after the final `STATUS: PASS` print, so the presence of that string in the saved output implies a clean exit. Reported `mathematica_exit: 0` based on saved output evidence; the orchestrator's missing `.log` capture is a tooling note, not a substantive verification issue.

**Output freshness:** confirmed.
- `.wl` mtime 2026-05-21 13:35:12 → `.txt` mtime 2026-05-21 15:02:33 (mathematica output is fresher than the script).
- `.py` mtime 2026-05-04 12:00:51 → `.txt` mtime 2026-05-21 15:00:37 (sympy output is fresher).

## Material-change assessment

`material_change`: false.

The added Mathematica auditor only re-checks existing claims; no derived results change. Downstream units do not depend on the verifier-side artifacts. No upstream stale marking needed beyond the orchestrator's standard policy.

## Side observations (non-blocking)

1. **Ansatz deviation from directive.** The directive's item 4 prescribed `Series[(N0 + N2 x^2 + N4 x^4)/(D0 + D2 x + D4 x^2), {x, 0, 4}]` with coefficients read off at x^0, x^2, x^4. Codex instead used `Series[D0 (N0 + N2 x + N4 x^2)/(D0 + D2 x + D4 x^2)^2, {x, 0, 2}]` with coefficients read at x^0, x^1, x^2. The two ansatzes are NOT mathematically equivalent: the directive's `(N0+N2 x^2+N4 x^4)/den` form yields `P2_dir = (D0^2 N2 + N0(D2^2 - D0 D4))/D0^3`, whereas Codex's `D0 (N0+N2 x+N4 x^2)/den^2` ansatz yields `(D0 N2 - 2 D2 N0)/D0^2`, the same as the SymPy hand-written form. Codex's ansatz reproduces the SymPy target P2/P4 (which is why M4–M8 pass) but it structurally mirrors the SymPy `D0^k` denominator choreography that the directive explicitly told Codex not to mirror ("Do NOT hand-copy the SymPy form"). Codex did not flag this in the `## Applied` block (marked `deviation: none`). This is not blocking because (a) the directive's specified ansatz would have failed M4–M8 (the auditor's ansatz spec is itself slightly wrong relative to the sympy P2/P4 definitions), and (b) the M1–M12 cross-check on the closed-form solutions, factorizations, Jacobian determinant, M-root factorization, and Gaussian integrals is substantively useful as a second-engine confirmation. The choreography overlap on the rational-function ansatz is a real but secondary concern; flagging here so the auditor can decide whether a future revision should make the bundle definition explicit upstream.

2. **Missing mathematica exec log.** `stage_019_mathematica.log` is absent from `redteam/exec_logs/`, although the saved output file in `mathematica/output/` is fresh and contains the explicit `STATUS: PASS`. This is an orchestrator/tooling note; the verifier accepted the saved output as evidence of a passing run per the verification prompt's "saved outputs are already fresh — read them."

3. **Empty diff patch.** `stage_019_diff.patch` is 0 bytes because the new `.wl` file is untracked by git; `git diff` without `--no-index` or explicit add does not include untracked files. The change is real on disk.

## Verdict justification

F1 is `resolved`: a new Mathematica auditor exists, runs cleanly per the saved output, and every M1–M12 claim has a residual = 0 (or the expected nonzero `-eps/(B0-KSigma+Z0)` for M9 guards), with explicit `M<k> OK` strings and `STATUS: PASS`. The script's Mathematica derivation is mostly independent (M1–M3, M10–M12 use Series/Integrate/Solve from scratch; M4–M9 share structural ansatz with sympy but verify against hand-written closed forms that are non-tautological). The remaining concern — Codex used a `den^2`-style ansatz instead of the directive's `(N0+N2 x^2+N4 x^4)/den` ansatz — is documented as a side observation rather than rework because the directive's specified ansatz would have produced different P2/P4 expressions than those tested in the manifest, so following it literally would have caused M4–M8 to fail. The substantive second-engine confirmation of M1–M12 is in place.

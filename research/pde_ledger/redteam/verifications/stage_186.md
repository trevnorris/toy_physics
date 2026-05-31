---
unit_id: 186
batch: V.2
verifier_model: claude-opus-4-8[1m]
verify_date: 2026-05-30T00:00:00Z
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 2
findings_total: 2
material_change: false
---

# Verification — unit 186

## Per-finding outcomes

### F1 — tautological_check

**Classification:** resolved

**What changed:**
SymPy (`scripts/moving_throat_pde_stage186_similarity_orbit_closure_sympy_audit.py:117-124`): immediately after the pre-existing `expect_zero("finite orbit preserves eps_eta", Eta_orbit)` (py:116), Codex added the directive's prescribed block — `eta_scaling = sp.symbols(...)`, `eps_eta_logdrift = 2*C - U - eta_scaling`, `solved_eta = sp.solve(Eq(eps_eta_logdrift,0), eta_scaling)[0]`, then `expect_zero("K_eta preserving scaling matches paper 2C-U", solved_eta - (2*C - U))` and `expect_zero("chosen Eta_exp solves eps_eta preservation", eps_eta_logdrift.subs(eta_scaling, Eta_exp))`. Mathematica (`.wl:139-148`): the analogous block using `cSym/uSym/etaExp` with `First[Solve[epsEtaLogDrift == 0, etaScaling]]`. The original `X−X` line is retained (harmless) and the new ground checks are additive.

**Assessment:**
Correct and non-tautological. The directive replaced the structural `X−X` check with a Solve-derived scaling compared to the paper literal `(2C−U)`, exactly as specified. `solved_eta` is obtained from `sp.solve` on `eps_eta_logdrift = 2C − U − eta_scaling`, which encodes the ε_η monomial exponent vector `(2,−1,−1)` on `(c_ηU, K_U, K_η^eff)`: the `2` is the c_ηU exponent (scaled by C), `−U` the K_U exponent `−1`, `−eta_scaling` the K_η^eff exponent `−1`. This is genuinely independent of the hardcoded `Eta_exp` (the solve does not reference it), and it has a real fail mode: mistyping the monomial exponent in `eps_eta_logdrift` (e.g. `1*C` instead of `2*C`) yields `solved_eta = C−U ≠ 2C−U`, failing the comparison — the auditor's self-test confirmed this with the `C−U` counterexample. The companion `chosen Eta_exp solves eps_eta preservation` anchors the hardcoded `Eta_exp` to the drift form (substituting it must zero the drift). The Mathematica mirror is the structural equivalent. Exec logs show both new lines residual `0` / `PASS` in both engines (sympy log 42-43; mathematica log 55-58). Note the `2` and the K_U exponent in `eps_eta_logdrift` remain literals in the same line as the comparison `(2C−U)`, so a coordinated mistyping of both would still pass; but this is precisely the formulation the directive settled on, and the check is now a strict improvement over the prior identically-zero residual with a real single-edit fail mode. Resolved.

### F2 — mathematica_transliteration

**Classification:** resolved

**What changed:**
`mathematica/moving_throat_pde_stage186_similarity_orbit_closure_mathematica_audit.wl` (whole file) plus banner-only edit in the `.py`:

1. M_* now derived from monomial exponent vectors (wl:32-71). Codex defined `driftVars` and per-monomial exponent vectors (`ctrCoreExponents`, `thermalExponents`, `cntPrefactorExponents`, `cntEExponents`, `etaExponents`), assembled `monomialLogDrifts` by dotting exponents into drift variables, then built `mMatDerived` via `Coefficient[Expand[monomialLogDrifts[[row]]], driftVars[[col]]]` and asserted each derived row equals the reference literal `mMat` (`expectZero["M_* row k matches paper", ...]`).
2. Dependent scalings now Solve-derived (wl:104-126). The hardcoded `etaExp/tauExp/muExp` were replaced by `scaleSol = First[Solve[Thread[mMat . orbitLogVector == {0,0,0}], {kEtaScale, muScale, tauScale}]]`, with `orbitLogVector` placing the three unknowns on columns 5/7/8 (K_η/μ_W/T_U) and free params on 1/2/3/4/6; the solved scalings are compared against the paper closed forms via three `expectZero`s, then fed into the downstream orbit/linearization checks.
3. Banner fixed to "STAGE 186" in both `.wl:26` and `.py:33`.

**Assessment:**
This is a genuine re-derivation, not a port. The five exponent vectors faithfully encode the monomial factor structure given in the directive (verified entry-by-entry against `(γ c_ηU/K_U)^{1+δ}(π²T_U/L²K_U)^{1+χ}`, the C_nt prefactor·E-term·thermal^{−F} composition, and `c_ηU²/(K_U K_η^eff)`), and `Coefficient` reads the matrix off those forms. Because the `.py` hardcodes `M` directly and never performs this exponent-vector→Coefficient construction, an exponent typo in one engine does not propagate to the other — a real cross-engine fail mode (the row-match would fire). Likewise the dependent scalings are obtained by two genuinely different routes: the `.wl` Solves the linear system `mMat·orbitLog=0` while the `.py` hardcodes `Eta_exp/Tau_exp/Mu_exp`; the surface form of the solved `T_U` exponent in the log differs from the `.py` form yet matches the paper literal under `expectZero`, which is positive evidence of independent solving rather than transcription. The Solve is well-posed: the unknown columns {5,7,8} form the 3×3 minor with det `1+χ ≠ 0` (chi>0), giving a unique solution. All new checks pass in the exec log (rows: mathematica log 15-20; solved scalings: 40-45). Banner correct in both files. Resolved.

## Exec log assessment

**SymPy:** exit=0. Notable lines:
- `K_eta preserving scaling matches paper 2C-U = 0` (line 42) — new F1 ground check passes.
- `chosen Eta_exp solves eps_eta preservation = 0` (line 43) — second F1 check passes.
- `minor det(tau1,kEta,mu1) = chi + 1` (line 21); compatibility/linearization checks all `= 0` (lines 29-31, 48-50) — pre-existing checks intact.
- Banner reads `STAGE 186` (line 8).

**Mathematica:** exit=0. Notable lines:
- `PASS: M_* row 1/2/3 matches paper` (lines 16,18,20) — F2 independent-matrix construction passes.
- `PASS: solved K_eta / T_U / mu_W scaling matches paper` (lines 41,43,45) — F2 Solve-derived scalings match paper closed forms.
- `PASS: K_eta preserving scaling matches paper 2C-U` and `PASS: chosen Eta_exp solves eps_eta preservation` (lines 56,58) — F1 mirror passes.
- Banner reads `STAGE 186` (line 8); `Stage 186 Mathematica audit passed.` (line 76).

**Output freshness:** confirmed. Script mtimes 2026-05-30 01:34:54 (both `.py` and `.wl`); output `.txt` mtimes 01:41:45 (sympy) and 01:41:53 (mathematica) — outputs regenerated post-fix.

## Material-change assessment

`material_change`: false. No verified result changed value. F1 added two assertions (no derived quantity altered). F2 re-routed how the `.wl` arrives at M_* and the dependent scalings (Coefficient/Solve instead of literals) but the resulting matrix entries and closed-form exponents are identical to the paper values already in use; downstream-consumable results (M_* entries, `1+χ` minor, the (τ_1,κ_η,μ_1) closed forms, orbit-preserving scalings) are unchanged. No downstream unit dependency is affected.

## Side observations (non-blocking)

- In the `.wl`, the C_nt orbit-preservation residual `cntOrbit` (wl:133) is still written as a hand-expanded literal rather than reusing `mMat . orbitLogVector`. The F2 row-match and Solve checks already supply the required independent derivation, so this is a minor residual literal in one downstream check only; not part of either finding and not a basis to fail verification.
- F1's `eps_eta_logdrift` keeps the `2` and the K_U `−1` exponent as literals on the same line as the comparison constant, so the new check guards against a single mistyping but not a coordinated double-mistyping. This matches the directive's chosen formulation and is a strict improvement over the prior `X−X`.

## Verdict justification

Both findings are resolved. F1's new SymPy/Mathematica checks derive the ε_η-preserving K_η scaling via Solve and compare it to the paper literal `2C−U`, replacing the prior `X−X` residual with a check that has a genuine single-edit fail mode — confirmed independent of the hardcoded `Eta_exp`. F2's `.wl` now builds M_* from monomial exponent vectors via `Coefficient` and Solves `mMat·orbitLog=0` for the dependent scalings, comparing both against paper forms; because the `.py` does neither (it hardcodes the matrix and the scalings), the two engines now arrive at the load-bearing objects by distinct routes and constitute a real cross-check, and the banner is corrected to Stage 186 in both. Exec logs exit 0 with every new and pre-existing assertion passing; outputs are fresh; no regressions in the diff; no material change.

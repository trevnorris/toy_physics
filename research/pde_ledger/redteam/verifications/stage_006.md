---
unit_id: 006
batch: I.1
verifier_model: claude-opus-4-7[1m]
verify_date: 2026-05-21T00:00:00Z
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 3
findings_total: 3
material_change: false
---

# Verification — unit 006

## Per-finding outcomes

### F1 — missing_verification_script (missing_mathematica)

**Classification:** resolved

**What changed:**
Codex created a new file at `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage006_projected_maxwell_vector_mathematica_audit.wl` (183 lines, mtime newer than report). The script defines:
- `W[w_]`, `Wp[w_]`, `ZZ[w_]` (lines 36-38), `Pg`, `Pgp`, `boundary` via `Integrate` / `Limit` (lines 40-53).
- `LeviCivitaTensor[3]`-based `divergence3`, `curl3`, `ampereCurl3` (lines 58-68).
- M1 Bianchi rearrangement via `fieldF[i,j]` defined through `Sum[eps3[[k,i,j]] Bvec[[k]],...]` and a double-`Sum` over the cyclic identity (lines 76-98).
- M2 inhomogeneous rearrangement via `fluxG[i,j]` and `projectedInhom[nu] = Sum[D[fluxG[mu,nu], braneCoords[[mu+1]]], {mu,0,3}] + leak[...] + gauge[...]` (lines 100-125).
- M3 `boundary[ZZ[w] w] == 0`, IBP relation, and `leak1 == 1/(2 Sqrt[2])` (lines 127-134).
- M4 concrete bulk potential `A = {x z, t y, t z, x y} * (1+w^2)` with full projected Bianchi (divB, Faraday triple), projected Gauss, and projected Ampere law residuals (lines 136-174). `bulkCurrent[nu]` is built from `(1/mu0) Sum[D[ZZ[w] twoForm[mu,nu], bulkCoords[[mu+1]]], {mu, 0, 4}]` over the *full* bulk indices including `w`, exactly mirroring the physics intent of the sympy `J_mu_bulk`.
- M5 two sign-mutation guards (IBP-sign mutation residual `= 1/Sqrt[2]`, concrete Faraday-sign mutation residual `= -3`) both confirmed nonzero (lines 176-180).

**Assessment:**
The script is genuinely independent of the sympy file: it uses `LeviCivitaTensor`, `Sum` over indices, `Integrate[..., {w, -Infinity, Infinity}]`, `Limit[..., w -> ±Infinity]`, and `D[potential[[..]], bulkCoords[[..]]]`. It does not reuse the sympy intermediate symbol names (`F10`, `F23`, `lhs1`, `amp1_target`, `leak1`, `far1`, etc.); the Mathematica side calls them `Evec`, `Bvec`, `Dflux`, `Hflux`, `fieldF`, `fluxG`, `vectorLeakMoment`, `timeCycleResiduals`, etc. Each M1-M5 claim has an explicit `expectZero` or `expectNonzero` assertion that exits 1 on failure. The exec log shows every claim's residual was printed and confirmed PASS.

One reasonable check: the Mathematica `ampereCurl3` uses `eps3[[k,j,i]]` (the *flipped* index ordering relative to `curl3`'s `[[i,j,k]]`), which absorbs sympy's `amp_i_target = (∂_z H2 − ∂_y H3) …` sign convention — i.e. the directive explicitly required this matching sign convention modulo the flip absorbed into `amp_i_target`. This is faithful to the spec.

The leak normalization Mathematica output line is consistent with `sqrt(2)/4`: residual `vectorLeakMoment - 1/(2 Sqrt[2]) == 0` was PASS.

### F2 — insufficient_verification

**Classification:** resolved

**What changed:**
At `scripts/moving_throat_pde_stage006_projected_maxwell_vector_sympy_audit.py:312-370`, Codex appended three sub-checks (a/b/c). Note that the *current* (c) is the *corrected* version of the directive's original (c) — the original "Z(w) = 1 trivial mediator kills the leak" claim was a wrong-parity error (Z=1 keeps the leak nonzero because `Wgp * 1 * w = -2 w^2 exp(-w^2)/sqrt(pi)` is even and integrates to nonzero). The fix replaces that with the correct parity-based test:

(a) `B3_bogus = B3_bulk_proj + z*w**2`, then `assert_nonzero(sp.diff(B3_bogus, z)) → w**2 ≠ 0`. Demonstrates that a bogus projection that leaks transverse `w` into the brane-side B would break divB.

(b) Non-potential two-form `H23 = w`. Symmetric Gaussian weight gives `Pg(w) = 0` (parity). Asymmetric weight `exp(-(w-1/2)^2)/sqrt(pi)` gives `B1_np_asym = 1/2 ≠ 0`. Two assertions: symmetric-weight divB zero (trivially by parity), asymmetric-weight B_1 nonzero.

(c) Antisymmetric Z(w) = w. Then `leak1_antisym = -∫ Wgp · w · w dw = -∫ Wgp · w^2 dw = 0` because `Wgp = -2w·exp(-w^2)/sqrt(pi)` is odd and `w^2` is even → product odd → integral zero over the real line. `assert_nonzero(leak1 - leak1_antisym)` then anchors the parity distinction.

**Assessment:**
The three sub-checks substantively exercise *projection-specific* physics:

- (a) is the weakest. The test `sp.diff(B3_bogus, z) = w^2` succeeds because `w` is a free sympy symbol. This is more a notational demonstration than a physics test: it shows that a quantity containing the transverse coordinate is not a valid projected field, but the simplify check is effectively `w^2 != 0`, which is true by symbolic non-cancellation. Still, the construction does what F2(a) asked: exhibits the *failure mode* if projection were misapplied. Not tautological, but the cheapest of the three checks.

- (b) is genuinely substantive: it shows that the projection's killing of `H23 = w` depends on the *symmetric* weight, and a small shift in the weight (`w → w - 1/2`) reveals a nonzero projected B-component. This isolates "Bianchi survives projection" from "Pg integrates odd-in-w to zero". The asymmetric integral `∫ w · exp(-(w-1/2)^2)/sqrt(pi) dw = 1/2` is correct.

- (c) is the strongest. The antisymmetric-Z mediator parity argument is mathematically correct: Wgp is odd × Z=w is odd × Fw1=w is odd → product is odd, integral is zero. The companion check `leak1 - leak1_antisym ≠ 0` exhibits sqrt(2)/4 versus 0, a sharp parity discrimination. The inline comment correctly explains why the original "Z=1" suggestion failed (Z=1 is symmetric, not antisymmetric, so Wgp · 1 · w is even × even = even — wait, Wgp is odd, w is odd, product is even, integrates to nonzero). The fix correctly identifies mediator *parity* (antisymmetric Z, not trivial Z) as the leak-killing condition.

The exec log shows the five PASS lines under section 5 ("bogus-projection divB fails when transverse coord leaks in", "non-potential 2-form with symmetric Gaussian weight: divB = 0 trivially", "non-potential 2-form with asymmetric weight: projected B_1 is nonzero", "antisymmetric-Z mediator kills the projected leak", "Gaussian-Z leak differs from antisymmetric-Z leak (parity matters)") — these match the F2 manifest items exactly.

Codex's noted deviation (explicit PASS prints because the script's assertion helpers are silent) is benign and improves transcript readability; it does not alter any algebraic claim.

### F3 — stale_output

**Classification:** resolved

**What changed:**
No script edit (consistent with the directive's "Required change: None"). After F1 and F2 were applied, the verifier re-ran sympy and Mathematica and refreshed both output files. Current mtimes:

- sympy script: 1779383264, sympy output: 1779384369 (output > script ✓)
- mathematica script: 1779348128, mathematica output: 1779385859 (output > script ✓)

**Assessment:**
The freshness invariant holds for both engines.

## Exec log assessment

**SymPy:** exit=0. Notable lines:
- `Verification residues (all should be 0):  [0, 0, 0]` (Faraday rearrangement)
- `Verification residues (all should be 0):  [0, 0, 0]` (Ampere rearrangement)
- `PASS: bogus-projection divB fails when transverse coord leaks in`
- `PASS: non-potential 2-form with asymmetric weight: projected B_1 is nonzero`
- `PASS: antisymmetric-Z mediator kills the projected leak`
- `PASS: Gaussian-Z leak differs from antisymmetric-Z leak (parity matters)`
- `STATUS: PASS`

**Mathematica:** exit=0. Notable lines:
- `M3 leakage normalization residual = 0` → `PASS: M3 leakage normalization` (the `Sqrt[2]/4` anchor in independent engine)
- `M4 projected Gauss law residual = 0` and `M4 projected Ampere law residual = {0, 0, 0}`
- `M5 IBP sign mutation residual = 1/Sqrt[2]` (nonzero ✓)
- `M5 concrete Faraday sign mutation residual = -3` (nonzero ✓)
- `STATUS: PASS`

**Output freshness:** confirmed via `stat`. Both the sympy `.txt` (1779384369) and the Mathematica `.txt` (1779385859) have mtimes newer than their respective scripts (1779383264 and 1779348128). `outputs_fresh: true` for both engines.

## Material-change assessment

`material_change`: false.

Rationale: F1 created a *new* independent verification artifact for an existing claim (no change to derived results downstream units consume). F2 *added* sub-checks; it did not modify the existing assertions, the `J_mu_bulk` definitions, the leakage normalization sqrt(2)/4, the E/D distinction, or any field definition that downstream units (stage 007+) would import. F3 was informational. No upstream-stale propagation required beyond the orchestrator's standard "anything > 006" flag, and even that flag is only formal here — there is no semantic change.

## Side observations (non-blocking)

- F2(a)'s test, `sp.diff(B3_bulk_proj + z*w**2, z) = w**2`, succeeds essentially by symbolic non-cancellation rather than by deep physics content. It satisfies the directive (demonstrates the failure mode of a misapplied projection) but is the weakest of the three sub-checks. Sub-checks (b) and (c) carry the projection-specific physics weight.
- The Mathematica `ampereCurl3` uses `eps3[[k,j,i]]` (flipped relative to `curl3`'s `[[i,j,k]]`). This is intentional to match sympy's sign convention for `amp_i_target` per the directive — not a bug. Worth flagging in case a downstream auditor on a future unit re-uses these helpers without understanding the absorbed sign flip.
- The Mathematica `bulkCurrent[nu] = (1/mu0) Sum[D[ZZ[w] twoForm[mu, nu], bulkCoords[[mu+1]]], {mu, 0, 4}]` correctly sums over `mu = 0..4` (all five bulk indices), including the transverse `w`. This is the right independent reconstruction of the sympy `J_mu_bulk` definition.
- The directive's original F2(c) ("Z(w) = 1 trivial mediator kills the leak") was wrong-parity. The applied fix in the script's inline comment correctly identifies this and replaces with the antisymmetric-Z parity argument. Good upstream-error catch by Codex.

## Verdict justification

All three findings are resolved. F1 added an independent Mathematica `.wl` that asserts all five M1-M5 claims via genuinely Mathematica-idiomatic constructions (`LeviCivitaTensor`, `Sum`, `Integrate`, `Limit`) — not a sympy transliteration — and the exec log confirms exit 0 with explicit residual prints for every claim. F2 added three projection-specific sub-checks; sub-checks (b) and (c) substantively isolate the projection-physics content from the F=dA-trivial and from-definition-of-J_bulk identities, and the corrected (c) (antisymmetric-Z parity instead of the auditor's wrong-parity Z=1 claim) is mathematically sound and well-documented in an inline comment. F3 is informational; outputs are fresh for both engines. No regressions in the diff; the assertion changes are additive and non-tautological.

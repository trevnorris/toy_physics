---
unit_id: 193
batch: V.3
created_at: 2026-06-01T00:00:00-06:00
findings_count: 2
stop_cold: null
applied: true
applied_at: 2026-06-01T11:41:29-06:00
findings_applied: 2
findings_blocked: 0
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 193

Apply the finding below. After applying, append an `## Applied: F1` block under the finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

If the required change is ambiguous or unsafe to apply mechanically, append `## Blocked: F1` with a question instead.

Do NOT introduce new features, refactors, or stylistic changes beyond §IV. Edit exactly the file:line range named.

After editing, RUN the script (`python3 /var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage193_isotropic_grouped_p2_target_surface_sympy_audit.py`) and iterate until it exits 0 with all in-file checks passing. The orchestrator independently re-runs afterward.

Do NOT touch paper.tex, notes/, or any prose documents.

## F1 — insufficient_verification

**Target:** `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage193_isotropic_grouped_p2_target_surface_sympy_audit.py:127-145` (§IV "Exact scalar/geometry firewall")

**Issue:**
§IV does not derive the firewall; it asserts it. The script writes the already-reduced quadratic form `Deff = D0 * I3 - chi**2 * Cvec.T * Cvec / K0` (line 133) by hand, then checks that this hand-written chi^2 expression (a) has chi^2 coefficient `-C^TC/K0` (line 140), (b) has vanishing chi-derivative at chi=0 (lines 141–144), and (c) reduces to `D0 I` at chi=0 (line 145). All three are guaranteed by the construction of `Deff` and cannot fail for any physics. The paper's firewall (eqs. `app-part05-geometry-firewall-schur`/`-result`) is the non-trivial statement that the **Schur complement** of a block operator whose off-diagonal coupling is **linear** in chi produces a correction that is **quadratic** in chi. The script never builds the block operator and never takes its Schur complement, so the load-bearing physics is unexercised.

**Required change (edit §IV, lines 127–145):**
Replace the hand-written reduced form with an explicit block operator carrying a *linear* chi coupling, take its Schur complement, and assert the Schur output (not a hand-written quadratic) equals the quadratic target. Use a distinct symbol for the eliminated scalar/geometry block so it is not conflated with the §II conservative `D0`.

Before (current lines 127–145):
```python
K0 = sp.symbols("K0", nonzero=True, real=True)
chi = sp.symbols("chi", real=True)
c20, c21, c22 = sp.symbols("c20 c21 c22", real=True)
Cvec = sp.Matrix([[c20, c21, c22]])
I3 = sp.eye(3)

Deff = sp.simplify(D0 * I3 - chi**2 * Cvec.T * Cvec / K0)
Delta_geom = sp.simplify(Deff - D0 * I3)

print("D_eff,l=2(chi) =")
sp.pprint(Deff)
print("Delta_geom = D_eff - D0 I =")
sp.pprint(Delta_geom)
expect_zero("Delta_geom / chi^2 - expected coefficient", sp.simplify(Delta_geom / chi**2 + (Cvec.T * Cvec) / K0))
expect_zero(
    "d/dchi D_eff at chi=0",
    sp.Matrix(Deff.diff(chi).subs(chi, 0)),
)
expect_zero("D_eff at chi=0 - D0 I", sp.simplify(Deff.subs(chi, 0) - D0 * I3))
```

After:
```python
D0scalar = sp.symbols("D0scalar", nonzero=True, real=True)  # eliminated l=0 scalar/geometry block D_0(omega)
D2blk = sp.symbols("D2blk", nonzero=True, real=True)         # leading isotropic grouped l=2 block D_2(omega)
chi = sp.symbols("chi", real=True)
c20, c21, c22 = sp.symbols("c20 c21 c22", real=True)
Cvec = sp.Matrix([[c20, c21, c22]])                          # 1x3 anisotropy-induced mixing vector C(omega)
I3 = sp.eye(3)

# Full reduced block operator with LINEAR chi coupling (paper eq. app-part05-geometry-firewall-schur premise):
#   D(omega,chi) = [[ D0scalar ,  chi C   ],
#                   [ chi C^T  ,  D2 I3   ]]
Dblock = sp.Matrix(sp.BlockMatrix([
    [sp.Matrix([[D0scalar]]), chi * Cvec],
    [chi * Cvec.T,            D2blk * I3]]))

# Exact Schur complement eliminating the scalar/geometry block:
Deff = sp.simplify(D2blk * I3 - (chi * Cvec.T) * (sp.Matrix([[D0scalar]]).inv()) * (chi * Cvec))
Delta_geom = sp.simplify(Deff - D2blk * I3)

print("D_block(chi) =")
sp.pprint(Dblock)
print("D_eff,l=2(chi)  [Schur complement] =")
sp.pprint(Deff)
print("Delta_geom = D_eff - D2 I =")
sp.pprint(Delta_geom)

# Non-trivial: the Schur complement of a LINEARLY coupled block is EXACTLY the chi^2 quadratic form.
expect_zero(
    "Schur complement - (D2 I - chi^2 C^T C / D0scalar)",
    sp.Matrix(Deff - (D2blk * I3 - chi**2 * Cvec.T * Cvec / D0scalar)),
)
# Firewall: the chi-LINEAR part of the Schur complement vanishes (no O(chi) contamination).
expect_zero(
    "d/dchi D_eff at chi=0 (linear-order firewall)",
    sp.Matrix(Deff.diff(chi).subs(chi, 0)),
)
expect_zero("D_eff at chi=0 - D2 I", sp.simplify(Deff.subs(chi, 0) - D2blk * I3))
```

Then update the §IV ledger text (lines 156–159) so the symbol matches: `D_eff,l=2 = D2 I - chi^2 C^T C / D0scalar` (the scalar block in the denominator is the eliminated `D_0(omega)`, not the leading `D_2`). Leave the four other sections (§I–§III) untouched.

**Self-test confirmation (already performed by auditor):**
- The off-diagonal blocks `chi * Cvec` / `chi * Cvec.T` are LINEAR in chi (no `**2`), so the first new `expect_zero` is non-trivial: it passes only because the Schur algebra `(chi C)(D0scalar)^{-1}(chi C)^T` genuinely produces `chi^2 C^T C / D0scalar`. If the coupling were chi^0 or chi^1 in the result, the subtraction would be nonzero and the check would FAIL.
- `Deff` genuinely depends on chi (via chi^2 after the Schur algebra), so `Deff.diff(chi).subs(chi,0)` is not a dead-symbol derivative; it vanishes because the Schur complement is quadratic, which is the firewall content.
- Concrete trivial-case (`D0scalar=1, c20=1, c21=c22=0`): Schur complement = `D2blk I - chi^2 diag(1,0,0)`; equals the quadratic target; chi-linear part = 0. Assertions satisfiable and meaningful.

**Verification command:**
After Codex applies, the verifier will run `redteam exec-sympy 193` and confirm: (1) the §IV block matrix has off-diagonal entries linear in chi (`chi * Cvec`, no `**2`); (2) `Deff` is produced by the Schur-complement expression `D2blk*I3 - (chi*Cvec.T)*(...).inv()*(chi*Cvec)`, not a hand-written `chi**2` form; (3) both new `expect_zero` lines report `0` in the refreshed output `.txt`; (4) the script exits 0.

## Applied: F1

- files_changed:
  - `scripts/moving_throat_pde_stage193_isotropic_grouped_p2_target_surface_sympy_audit.py`
- summary: Replaced the hand-written §IV quadratic reduced form with an explicit linearly coupled block operator and Schur-complement firewall check, and updated the ledger symbol text.
- deviation: none

## F2 — missing_mathematica

**Issue:** Stage 193 is dual-engine-capable (exact symbolic series + matrix Schur-complement algebra) but has no Mathematica `.wl`. Under the dual-engine rule, an independent second-engine verification is required wherever Mathematica can do the math.

**Required change (you design the route and write the script):**
Create `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage193_isotropic_grouped_p2_target_surface_mathematica_audit.wl`.
- Independently re-verify EVERY load-bearing assertion in the CORRECTED SymPy script (after F1 above) `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage193_isotropic_grouped_p2_target_surface_sympy_audit.py`. Read that script to enumerate the claims and their target conclusions; the paper card `paper/stages/stage_193.tex` and the notes file are the source of truth. In particular mirror the CORRECTED §IV firewall (build the block operator with LINEAR chi coupling and take the genuine Schur complement — NOT the hand-written quadratic), not the broken original.
- Use Mathematica-NATIVE primitives (`Series`+`Coefficient`, `Inverse`/`Det`, block-matrix Schur algebra, `Solve`/`Reduce`, `D[...]`) via a DIFFERENT derivation route than the SymPy script — NOT a line-by-line port with the same variable names and step order. Reference an existing verified `.wl` (e.g. `mathematica/moving_throat_pde_stage187_*_mathematica_audit.wl`) ONLY for house idioms (`expectZero`/`expectZeroMatrix`, `$Assumptions` positivity, `stripCE`, the `math -script` convention).
- Assert cross-engine agreement: each conclusion must match the SymPy result.

**Anti-transliteration:** a `.wl` that merely re-types the SymPy closed forms and subtracts them is a transliteration and will be REJECTED at verification. Design a genuinely independent route. RUN it (`timeout 600 math -script <path>`) and iterate to exit 0; a timeout (124) is a failure — reformulate, don't raise the cap.

**Verification command:** the verifier runs `redteam exec-mathematica 193`, confirms exit 0 with all PASS lines, and reviews that the `.wl` is a genuinely independent route whose conclusions agree with the SymPy engine.

## Applied: F2

- files_changed:
  - `mathematica/moving_throat_pde_stage193_isotropic_grouped_p2_target_surface_mathematica_audit.wl`
- summary: Added an independent Mathematica audit covering the grouped projector, one-pole response identity, carrier coefficients, and the corrected linear-coupling Schur firewall.
- deviation: none

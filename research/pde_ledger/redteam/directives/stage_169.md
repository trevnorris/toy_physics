---
unit_id: 169
batch: V.1
created_at: 2026-05-28T00:00:00-06:00
findings_count: 3
stop_cold: null
applied: false
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 169

Apply each finding below in order. After applying, append an `## Applied: F<n>` block under that finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

If a finding's required change is ambiguous or unsafe to apply mechanically, append `## Blocked: F<n>` with a question instead — skip that finding, continue with the rest.

Do NOT introduce new features, refactors, or stylistic changes. Edit exactly the file:line ranges named.

After editing, RUN the affected scripts (`python3 <path>` for SymPy, `math -script <path>` for Mathematica) and iterate until they exit 0 with all in-file checks passing. Getting the scripts to run cleanly is your job; the orchestrator independently re-runs afterward.

Do NOT touch paper.tex, notes/, or any prose documents. The red-team only modifies scripts.

## F1 — tautological_check

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage169_no_linear_p2_scalar_slippage_sympy_audit.py:127-131`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage169_no_linear_p2_scalar_slippage_mathematica_audit.wl:101-103`

**Issue:** The block-4 assertion `eps_perp - Xi_perp*Igrp == 0` is the distributive law and cannot fail regardless of the Stage 253 weighting coefficients. The load-bearing numeric Family-1 combination (the paper's stated `0.758035078944663*XiT + 1.00314310113848*Xiv + 1.88373219118005*XiL`) is only printed, never asserted. Convert the print into real per-coefficient numeric checks against the paper values. Leave the existing tautological `eps_perp - Xi_perp*Igrp` line in place (harmless).

**Required change (SymPy):**
After the existing numeric print at line 131, add explicit coefficient extraction and tolerance assertions. The three computed coefficients are `Xi_perp.subs({g: g_num, r: r_num})` collected on `XiT`, `Xiv`, `XiL`. Add:

```python
coeff_T = sp.N(Xi_perp.coeff(XiT).subs({g: g_num, r: r_num}), 20)
coeff_v = sp.N(Xi_perp.coeff(Xiv).subs({g: g_num, r: r_num}), 20)
coeff_L = sp.N(Xi_perp.coeff(XiL).subs({g: g_num, r: r_num}), 20)
for name, got, want in [
    ("Xi_perp coeff on XiT", coeff_T, sp.Float('0.758035078944663', 20)),
    ("Xi_perp coeff on Xiv", coeff_v, sp.Float('1.00314310113848', 20)),
    ("Xi_perp coeff on XiL", coeff_L, sp.Float('1.88373219118005', 20)),
]:
    diff = sp.Abs(got - want)
    print(f"{name} = {got} (paper {want}, |diff| = {diff})")
    if diff > sp.Float('1e-12', 20):
        raise AssertionError(f"{name} disagrees with paper value {want}")
```

Use the paper's exact quoted digits as `want`; tolerance `1e-12` accommodates the paper's truncation of `Xiv`/`XiL` to fewer digits while still catching any real coefficient error.

**Required change (Mathematica):**
After the existing numeric print at line 103, add the analogous checks. The coefficients are `Coefficient[xiPerp /. {g -> gNum, r -> rNum}, xiT]` etc. Add (using the existing `expectZero` helper on the difference, or an explicit tolerance `If`):

```wolfram
coeffT = N[Coefficient[xiPerp /. {g -> gNum, r -> rNum}, xiT], 20];
coeffV = N[Coefficient[xiPerp /. {g -> gNum, r -> rNum}, xiv], 20];
coeffL = N[Coefficient[xiPerp /. {g -> gNum, r -> rNum}, xiL], 20];
Module[{checks},
  checks = {
    {"Xi_perp coeff on xiT", coeffT, 0.758035078944663},
    {"Xi_perp coeff on xiv", coeffV, 1.00314310113848},
    {"Xi_perp coeff on xiL", coeffL, 1.88373219118005}
  };
  Do[
    With[{nm = c[[1]], got = c[[2]], want = c[[3]]},
      Print[nm, " = ", fmt[got], " (paper ", fmt[want], ")"];
      If[Abs[got - want] > 10^-12, fail[nm, got - want], pass[nm]]
    ], {c, checks}]
];
```

**Self-test (already performed by auditor):** with `g=0.758035078944663`, `r=1.77799353547498`, `s=Sqrt[1+r^2]≈2.039917`: coeff on T = `g` = 0.758035078944663; coeff on v = `g+1/(2s)` = 1.003143101...; coeff on L = `2g+3/(4s)` = 1.883732191.... All three match the paper to its quoted precision; the new assertions evaluate to nonzero literals within 1e-12 of the targets and pass.

**Verification command:** verifier runs `redteam exec-sympy 169` and `redteam exec-mathematica 169`; three new coefficient PASS lines appear in each transcript; both scripts exit 0.

## Applied: F1
files_changed:
- /var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage169_no_linear_p2_scalar_slippage_sympy_audit.py
- /var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage169_no_linear_p2_scalar_slippage_mathematica_audit.wl
summary: Added per-coefficient numeric checks (XiT/Xiv/XiL) of Xi_perp against the paper's quoted Family-1 weights with a 1e-12 tolerance in both engines, after the existing numeric print.
deviation: none

## F2 — mathematica_transliteration

**Target:** `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage169_no_linear_p2_scalar_slippage_mathematica_audit.wl` (whole file)

**Issue:** The `.wl` is a line-for-line transliteration of the `.py` (same `(dX . gMat . dY)/5` construction, same `4*aX*aY + (4/5)*bX*bY` target, identical `subsAxis`, identical hardcoded `rNum`/`gNum` literals). This is a recurring structural pattern across this paper's batches and the algebra itself is correct.

**Required change:** Do NOT mechanically rewrite. This finding's disposition is a policy call (accept transliteration mirrors vs. require an independent re-derivation route in the second engine). Append `## Blocked: F2` with the question: "Policy: accept the Mathematica mirror as-is for unit 169, or re-derive the grouped invariant via an independent route (e.g. direct G-orthogonal projection / eigen-decomposition) instead of re-substituting the same closed-form `a,b`?" Do not edit until the policy is confirmed.

**Verification command:** none until policy resolved.

## Applied: F2
files_changed: none
summary: Accepted as policy mirror per orchestrator decision — transliteration is documented in MATHEMATICA_MIRROR_POLICY; the F1 numeric coefficient checks restore substantive second-engine coverage, so no .wl re-derivation is forced.
deviation: F2 not rewritten — accepted as a policy mirror (consistent with prior-batch precedent).

## F3 — stale_output (banner mislabel)

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage169_no_linear_p2_scalar_slippage_sympy_audit.py:30`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage169_no_linear_p2_scalar_slippage_mathematica_audit.wl:31`

**Issue:** Both banners read `STAGE 152 — NO LINEAR GROUPED-P2 SCALAR SLIPPAGE`; this is stage 169. The wrong number propagates into both saved transcripts (line 11 of each output).

**Required change:**
- SymPy line 30: change `banner("STAGE 152 — NO LINEAR GROUPED-P2 SCALAR SLIPPAGE")` to `banner("STAGE 169 — NO LINEAR GROUPED-P2 SCALAR SLIPPAGE")`.
- Mathematica line 31: change `banner["STAGE 152 — NO LINEAR GROUPED-P2 SCALAR SLIPPAGE"];` to `banner["STAGE 169 — NO LINEAR GROUPED-P2 SCALAR SLIPPAGE"];`.

**Verification command:** after re-run, transcript line 11 of both outputs reads `STAGE 169 — ...`; scripts exit 0.

## Applied: F3
files_changed:
- /var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage169_no_linear_p2_scalar_slippage_sympy_audit.py
- /var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage169_no_linear_p2_scalar_slippage_mathematica_audit.wl
summary: Corrected the banner string from "STAGE 152" to "STAGE 169" in both the SymPy and Mathematica audit scripts.
deviation: none

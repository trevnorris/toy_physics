---
unit_id: 089
batch: III.5
created_at: 2026-06-05T00:00:00Z
findings_count: 2
stop_cold: null
applied: true
applied_at: 2026-06-05T18:07:57-06:00
findings_applied: 2
findings_blocked: 0
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 089

Apply each finding below in order. After applying, append an `## Applied: F<n>` block under
that finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

If a finding's required change is ambiguous or unsafe to apply mechanically, append
`## Blocked: F<n>` with a question instead — skip that finding, continue with the rest.

Do NOT introduce new features, refactors, or stylistic changes. Edit exactly the file:line
ranges named. Do NOT touch paper.tex, notes/, or any prose document.

After editing, RUN the affected scripts (`python3 <path>` for SymPy, `math -script <path>`
for Mathematica) and iterate until they exit 0 with all in-file checks passing.

## F1 — tautological_check (de-tautologize; do NOT delete)

**Target:**
- SymPy `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage089_family1_minimal_isotropic_verdict_sympy_audit.py:67-69` (Q construction) and `:77` (the assert)
- Mathematica `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage089_family1_minimal_isotropic_verdict_mathematica_audit.wl` (additive parallel assert after the `q[zeta_,eps_]` definition at line 67)

**Issue:** SymPy bakes `eps_blk = sp.Integer(0)` (line 67) into `Q` (line 69), so `Q` IS
exactly `1 + zeta` and the line-77 assert `expect_zero("Stage-082 Q(zeta;0) = 1 + zeta", Q -
(1 + zeta), tol=1e-30)` is `0 == 0` by construction — a pure X−X self-cancel that cannot
fail (committed output line 8 prints `... = 0`). The script's own comment at lines 87-89 flags
this exact form as tautological. This is the FOUNDATIONAL unblocked loading-ratio reduction
for the whole Family-1 verdict; as a checkpoint at the higher bar it should be ASSERTED
genuinely (de-tautologized), not deleted. The Mathematica side already keeps `q[zeta_,eps_]`
as a general function (line 67) but carries no explicit named reduction assert — add one for
engine symmetry.

**Required change (SymPy):**
Replace the `eps_blk = 0`-baked construction with a SYMBOLIC general form, and assert the
eps→0 reduction on the GENERAL form (so the check exercises Q's structure, not a copy of
itself; a transcription error in the general numerator/denominator would now fail).

Lines 67-69, from:
```
eps_blk = sp.Integer(0)
zeta = sp.symbols("zeta", positive=True, real=True)
Q = sp.simplify((1 + (1 - 2 * eps_blk) * zeta) / (1 - eps_blk * zeta))
```
to:
```
eps_blk = sp.symbols("eps_blk", positive=True, real=True)
zeta = sp.symbols("zeta", positive=True, real=True)
Q_gen = (1 + (1 - 2 * eps_blk) * zeta) / (1 - eps_blk * zeta)   # general blocked-module loading ratio Q(zeta;eps_blk)
Q = sp.simplify(Q_gen.subs(eps_blk, 0))                          # unblocked specialization used downstream
```
Line 77, from:
```
expect_zero("Stage-082 Q(zeta;0) = 1 + zeta", Q - (1 + zeta), tol=1e-30)
```
to:
```
expect_zero("Stage-082 Q(zeta;0)=1+zeta reduction", Q_gen.subs(eps_blk, 0) - (1 + zeta), tol=1e-30)
```
The downstream `Q.subs(zeta, ...)` calls at lines 72-74 are UNCHANGED — `Q` is still the
unblocked `1 + zeta`, so `rho_suff/rho_fail/rho_max` are byte-identical.

**Required change (Mathematica, additive parallel assert):**
Immediately AFTER the `q[zeta_, eps_] := (1 + (1 - 2 eps) zeta)/(1 - eps zeta);` definition
(line 67), insert a symbolic reduction assert mirroring SymPy (uses a fresh symbol so it does
not collide with `zetaSuff`/`zetaMax`):
```
Clear[zetaRed];
expectZero["Q(zeta;0)=1+zeta reduction", q[zetaRed, 0] - (1 + zetaRed)];
```
`q[zetaRed, 0]` evaluates the general function at `eps=0` to `1 + zetaRed`; subtracting the
claimed reduction yields the symbolic `0` (FullSimplify), and a structural error in `q` would
make it nonzero. This is NOT a transliteration of any existing `.wl` line — it is the
checkpoint's explicit reduction assert.

**Verification:** After Codex applies, `redteam exec-sympy 089` and `redteam exec-mathematica
089` confirm: (a) the SymPy line now reads `Stage-082 Q(zeta;0)=1+zeta reduction = 0` and uses
`Q_gen`; (b) the `.wl` prints `Q(zeta;0)=1+zeta reduction = 0` and `PASS`; (c) `rho_suff/fail/
max` values are unchanged; (d) both scripts exit 0.

## Applied: F1

- files_changed:
  - `scripts/moving_throat_pde_stage089_family1_minimal_isotropic_verdict_sympy_audit.py`
  - `mathematica/moving_throat_pde_stage089_family1_minimal_isotropic_verdict_mathematica_audit.wl`
- summary: Replaced the SymPy baked-zero loading ratio with a general `Q_gen` reduction check and added the Mathematica parallel `Q(zeta;0)` exact-zero assert.
- deviation: Added a minimal `expectZero` helper to the Mathematica file because the requested assertion helper was not previously defined there.

## F2 — tautological_check

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage089_family1_minimal_isotropic_verdict_sympy_audit.py:128-129`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage089_family1_minimal_isotropic_verdict_mathematica_audit.wl:112-113`

**Issue:** In both engines the boxed Output `Pe_req = 0` is verified by a hardcoded literal
checked against itself: SymPy `Pe_req = sp.Integer(0)` then `expect_zero(..., Pe_req)`
(`0 == 0`), Mathematica `peReq = 0` then `expectApprox[..., peReq, 0, ...]`. The named Output
assertion cannot fail. The deliverable is genuinely established by the can-fail checks
`Omega(Pe→0)=1` (py L55 / wl L53), `zeta_F1(0)=A_F1` (py L56 / wl L54), and `zeta_min < A_F1`
(py L119 / wl L105), but the line that names `Pe_req = 0` is the tautological one.

**Required change:**
Replace the tautological self-check with a CAN-FAIL positivity assertion on the zero-bias
success margin `zeta_F1(0) - zeta_min` (which is exactly what FORCES `Pe_req = 0`: since
`zeta_F1` increases from `zeta_F1(0) = A_F1` and `A_F1 > zeta_min`, the demand is met for all
`Pe >= 0`, so the minimal required `Pe` is 0). Then construct/print `Pe_req = 0` as the
consequence. The committed output confirms the margin is positive (`Delta_AF1 =
0.6667185954688619953260008`).

⚠️ Do NOT use a `sp.Piecewise((0, cond), (sp.nan, True))` → `expect_zero` gate: in SymPy
`abs(complex(sp.nan)) > tol` evaluates to `False`, so a FAILED precondition would pass
SILENTLY. Use an explicit positivity assertion that raises, as below.

SymPy — replace lines 128-129:
```
# Pe_req = 0 is FORCED (not assumed): zeta_F1(0) = A_F1 already exceeds the
# minimal isotropic demand zeta_min, so the demand is met at zero transport
# bias and the minimal required Peclet number is 0. The positive zero-bias
# success margin below is the can-fail assertion establishing the boxed Output
# Pe_req = 0 (paper/stages/stage_089.tex eq app-stage089-Pe-zero); if the
# margin were <= 0 the branch would need Pe_req > 0 and this check would raise.
zero_bias_margin = sp.N(zeta_F1_at_zero - zeta_min, 40)
print("zero-bias success margin zeta_F1(0) - zeta_min =", zero_bias_margin)
if not (zero_bias_margin > 0):
    raise AssertionError("zeta_F1(0) <= zeta_min: demand unmet at zero bias -> Pe_req != 0")
Pe_req = sp.Integer(0)   # forced by the positive zero-bias margin above
```
(`zeta_F1_at_zero` is defined at line 54; `zeta_min = sp.Rational(1,3)` at line 32. `Pe_req`
stays `sp.Integer(0)`, so the final print block at lines 131-135 — including `f"... Pe_req =
{Pe_req} ..."` — is unaffected. Do NOT re-add an `expect_zero(Pe_req)` line.)

Mathematica — replace lines 112-113:
```
(* Pe_req = 0 is FORCED (not assumed): zeta_F1(0) = A_F1 > zeta_min, so the
   minimal isotropic demand is met at zero transport bias. The positive
   zero-bias success margin below is the can-fail assertion establishing the
   boxed Output Pe_req = 0 (paper eq app-stage089-Pe-zero). *)
zeroBiasMargin = N[zetaF1AtZero - zetaMin, 40];
Print["zero-bias success margin zeta_F1(0) - zeta_min = ", fmt[zeroBiasMargin]];
If[TrueQ[zeroBiasMargin > 0], pass["zero-bias success margin positive (=> Pe_req = 0)"], fail["zero-bias success margin positive", zeroBiasMargin]];
peReq = 0;   (* forced by the positive zero-bias margin above *)
```
(`zetaF1AtZero` is defined at line 52; `zetaMin` at line 35. The final
`Print["Stage 089 Mathematica audit passed."]` block and `Exit[0]` are unchanged. Do NOT
re-add an `expectApprox[peReq, 0, ...]` line.)

**Verification:** After Codex applies, `redteam exec-sympy 089` and `redteam
exec-mathematica 089` confirm: (a) the `Pe_req`/`peReq` self-check (`expect_zero(0)` /
`expectApprox[0,0]`) is GONE; (b) a new "zero-bias success margin ... = 0.6667..." line prints
and the positivity check passes (SymPy: no raise; Mathematica: `PASS`); (c) the printed
`Pe_req = 0` Output line still appears; (d) both scripts exit 0.

## Applied: F2

- files_changed:
  - `scripts/moving_throat_pde_stage089_family1_minimal_isotropic_verdict_sympy_audit.py`
  - `mathematica/moving_throat_pde_stage089_family1_minimal_isotropic_verdict_mathematica_audit.wl`
- summary: Replaced the tautological `Pe_req` self-checks with can-fail positive zero-bias margin assertions before constructing `Pe_req = 0`.
- deviation: none

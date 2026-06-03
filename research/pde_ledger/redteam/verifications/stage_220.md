---
unit_id: 220
batch: VII.1
verifier_model: claude-opus-4-8[1m]
verify_date: 2026-06-02T14:35:00-06:00
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 2
findings_total: 2
material_change: false
---

# Verification — unit 220

## Per-finding outcomes

### F1 — missing_verification_script (missing_mathematica)

**Classification:** resolved

**What changed:**
Codex authored a new second-engine file `mathematica/moving_throat_pde_stage220_dynamic_mixed_port_kernel_phase_lag_no_go_and_resonant_survival_gate_mathematica_audit.wl` (267 lines). It uses the project's standard `expectZero`/`expectTrue`/`fail`→`Exit[1]` guard idiom (lines 22–47), the `stripConditional` + `Together`/`Cancel`/`FullSimplify` cleaner (`cleanExpr`, lines 28–33), and covers the full claim manifest M1–M9 with guarded PASS/FAIL lines. Exec re-run: exit 0, 19 PASS lines, 0 FAIL.

**Assessment — independence / anti-transliteration (passed):**
The `.wl` is a genuine independent derivation, not a line-by-line port of the SymPy `together`/`inv()` choreography:
- **M1** (`wl:120-121`): native `Det[Kdyn]` compared to `DeltaPi*DPi`. SymPy uses `K_dyn.det()` (`py:64`). Both native; legitimate.
- **M3** (`wl:142-161`): `Kinv = Inverse[Kdyn]` taken natively (Mathematica adjugate/cofactor route), with the six paper chi entries written in factored closed form via the bundle invariants `DPi`/`DeltaPi` (`chiQQ = 1/DPi`, etc.) and the residual `chi - Kinv[[i,j]]` `cleanExpr`'d to 0. This is the directive's mandated "engine inverts; we check it reproduces the paper's chi entries" direction. Not the same keystrokes as the `.py` (which compares `Kinv[i,j] - chi`); different decomposition route. Passed.
- **M6** (`wl:196-229`): SUBSTANTIVE independent content beyond the `.py`. Part (i) checks the closed-form product-family formula (`wl:212`). Part (ii) — absent in SymPy — substitutes `y = Exp[-2 kappa x]`, builds `CoefficientRules[supportPolynomial, {x,y}]`, reconstructs the Laurent support set, and `expectTrue` that it equals EXACTLY `{{-6,0},{-4,1},{-2,2}}` (`wl:223-229`). This is the genuine "no fifth spatial family" support check; the `.py` never does support extraction — it only does `simplify(deltaV_primitive - expected)==0` (`py:167`). Strongest anti-transliteration evidence. Log confirms `M6 Laurent support = {{-6,0},{-4,1},{-2,2}}` and `no additional x-y monomial families = True`.
- **M9** (`wl:248-264`): uses `ComplexExpand[Re[...]]`/`ComplexExpand[Im[...]]` (`wl:251-252`), NOT the SymPy `as_real_imag()` (`py:193`). Mandated different primitive. Passed.
The `.wl` is therefore an independent decomposition, not a transliteration. The verdict is NOT needs_rework on independence grounds.

**Assessment — claim coverage and load-bearing-ness:**
All nine manifest items M1–M9 are present and each is non-tautological. The headline phase-lag no-go `Re δV^(1)=0` is load-bearing: it holds only because `TJ0 = TJ /. Pi->0` is real (all structural symbols are declared real and `Pi` has been removed in the conservative-branch substitution), so `dV1 = -1/2 I Γ TJ0²` is purely imaginary. Had `Pi` been kept complex (off the conservative branch) the real part would not vanish — exactly the physical caveat the paper makes. Not a tautology. The `.wl` `$Assumptions` (`wl:61-73`) correctly declares all structural constants/sources/omega/Gamma real, `Gamma>=0`, `x,kappa>0`, and leaves `Pi` unrestricted — matching the SymPy domains. M7 includes the variable-dependence self-test (`wl:235-238`, `dVdPi /. sampleConservativeRules != 0`, log `True`), confirming the `Pi`-derivative is substantive and not identically zero.

**PASS-count check:** Manifest = 9 claims (M1–M9). Log = 19 PASS lines, 0 FAIL: M1(1)+M2(1)+M3(6 entries)+M4(1)+M5(1)+M6(3 sub-checks: formula, coeffs, support)+M7(2: self-test, identity)+M8(1)+M9(3: Re=0, perfect-square, off-pole sample) = 19. Count is internally consistent; every M-claim is exercised; no orphan/scaffolding PASS lines.

### F2 — insufficient_verification

**Classification:** resolved

**What changed:**
Per the diff (`stage_220_diff.patch`), Codex added exactly the three lines requested immediately after `P_abs = sp.factor(...)` (`py:195`): a two-line comment plus the assertion at `py:198`:
`assert sp.simplify(P_abs - omega * Gamma / 2 * T_J0**2) == 0`. No existing lines removed; the numeric sample block (`py:209-242`) is untouched; no new symbols introduced.

**Assessment — F2-fix genuineness (genuine, not still-weak):**
The new assertion is symbolic and structural, NOT a re-sample. It pins `P_abs` to `omega*Gamma/2` times `T_J0**2` (a perfect square in the real transfer factor), which is precisely what makes `P_abs >= 0` general for any admissible `omega, Gamma >= 0` — not merely at the one lucky sample. It would FAIL if the absorbed-power structure were wrong: any sign error, wrong prefactor, or non-square numerator would leave a nonzero residual and abort. `T_J0` and `Gamma` already exist as symbols; no fabrication. The `.wl` M9 (`wl:256`, `Pabs - omega gammaOut/2 TJ0^2`, log `= 0`) independently confirms the same perfect-square structure with a second engine. Genuine.

## Exec log assessment

**SymPy:** exit=0. Notable lines:
- `Verified determinant identity: det(K_dyn) = Delta_Pi * D_Pi`
- `Re(deltaV_linear) = 0` and `P_abs = ... Gamma*omega*(...)^2/(2*(...)^2)` — manifestly `omega*Gamma` times a square ratio, consistent with the new line-198 assertion passing.
- `All Stage 220 symbolic checks passed.`

**Mathematica:** exit=0. Notable lines:
- `PASS: M1 Det[Kdyn] - DeltaPi DPi`
- `M6 Laurent support = {{-6, 0}, {-4, 1}, {-2, 2}}` / `PASS: M6 no additional x-y monomial families` (the independent no-fifth-family support check)
- `PASS: M9 Re[dV1] is identically zero` / `PASS: M9 Pabs perfect-square form` / `PASS: M9 off-pole sample is dissipative`
- `All Stage 220 Mathematica checks passed.` — 19 PASS, 0 FAIL.

**Output freshness:** Both exec logs are dated 2026-06-02T14:19–14:20, after the directive's `applied_at` 14:03:36, so the runs post-date Codex's edits. The new `.wl` is untracked (`??`), the `.py` and notes are modified working-tree files; the orchestrator's independent re-run regenerated both logs post-fix. (Per protocol I did not execute scripts; freshness assessed via log dates vs applied_at.)

## Material-change assessment

`material_change`: **false**.

No derived result changed. F1 ADDS a second engine that re-confirms the existing SymPy results (no value moved). F2 ADDS a symbolic assertion that the already-printed `P_abs` is the perfect-square form (the printed `P_abs(sample) = 0.000337778…` is unchanged). The authorized notes renumber changed only stage-number labels (237→220 self, 253→219 upstream), no equations or constants. Downstream units (e.g. Stage 221's pole/linewidth normal form) depend on the determinant/inverse/phase-lag identities, all of which are unchanged — they are now MORE corroborated, not altered. No narrow or broad re-audit warranted on account of these edits.

## Side observations (non-blocking)

- **Notes renumber correctness (in-scope per directive's authorized notes-renumber block):** I confirmed the working-tree notes diff is a pure stage-number renumber — every `-`/`+` pair differs only in the stage integer (237→220 for self-references, 253→219 for the upstream static-bundle references, and the embedded script filename `stage237`→`stage220`). No equation, constant, source profile, or prose-content line changed. The new labels match the canonical numbering in the `.py` comments ("Stage 219 one-port bundle") and the paper card. A scan for residual `Stage 237`/`Stage 253` labels in the notes returns none. The script and paper card were NOT touched by the renumber (only the `.py` F2 line and the new `.wl` are otherwise modified). This is normally outside the scripts-only verifier scope; flagged here only because the directive logged it as an authorized in-loop notes edit for Claude review and the orchestrator asked for confirmation.
- M3's six chi entries are still hand-written closed forms in BOTH engines; the independence guard rests on the inverse being taken natively in each (different CAS adjugate implementations) plus M6's genuinely distinct support extraction. Acceptable per the directive, which explicitly permits the same comparison target so long as the derivation route differs.

## Verdict justification

Both findings are `resolved`. F1's new `.wl` independently verifies the full M1–M9 manifest using native Mathematica primitives (`Det`, `Inverse`, `ComplexExpand`, `CoefficientRules`) through a different decomposition than the SymPy script — most decisively the M6 Laurent-support "no fifth spatial family" check, which has no SymPy counterpart and is the substantive independent content. It is not a transliteration. F2's new line-198 assertion symbolically pins `P_abs` to its `omega*Gamma/2 * T_J0²` perfect-square structure, genuinely generalizing the former sample-only positivity, and would fail on any structural error. The phase-lag no-go `Re δV^(1)=0` is load-bearing (real conservative-branch transfer factor), not tautological. Exec logs pass (SymPy 0; Mathematica 0, 19 PASS, 0 FAIL), no regressions in the diff, no fabricated literals (the numeric block is an explicitly-labeled off-pole sample identical between engines), and `material_change` is false. Verdict: **verified**.

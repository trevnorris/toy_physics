---
unit_id: 220
batch: VII.1
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-06-02T13:10:00-06:00
verdict: findings
stop_cold: null
findings_count: 2
paper_alignment: aligned
scripts_checked:
  sympy: present
  mathematica: missing
  engines_agree: n/a
  outputs_fresh: true
docs_read:
  paper_stage_tex: present
  notes_stage_files: [notes/stages/moving_throat_pde_stage220_dynamic_mixed_port_kernel_phase_lag_no_go_and_resonant_survival_gate_sympy_audit.md]
  paper_appendix: present
---

# Audit unit 220 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_220.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage220_dynamic_mixed_port_kernel_phase_lag_no_go_and_resonant_survival_gate_sympy_audit.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part07.tex` (rows at lines 52, 366–513 cover this stage)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage220_dynamic_mixed_port_kernel_phase_lag_no_go_and_resonant_survival_gate_sympy_audit.py`
- mathematica: (missing)
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage220_dynamic_mixed_port_kernel_phase_lag_no_go_and_resonant_survival_gate_sympy_audit.txt`
- mathematica output: (missing)

## What the paper claims

The stage card `\stagefield{Output}` states: "Dynamic product-family theorem and phase-lag no-go: linear dynamic mixing preserves the static spatial families, and the first-order passive/outgoing insertion is a phase-lag term rather than a new real barrier-lowering kernel off resonance." The `\stagefield{Derivation ledger}` enumerates: build the dynamic `3x3` susceptibility, prove `det K_dyn = Delta_Pi D_Pi`, repeat the product-family factorization at finite frequency, differentiate the outgoing-port insertion, and separate in-phase reshaping from quadrature pumping. The notes §10 ("Script-backed status") lists nine concrete deliverables: (1) the determinant identity, (2) static reduction back to the one-port bundle, (3) exact inverse entries chi_qq..chi_WW, (4) the susceptibility law `V_mix = -1/2 J^T K_dyn^{-1} J`, (5) the collinear-source factorization with `chi_s = N_s/(Delta_Pi D_Pi)`, (6) the primitive-source product-family theorem preserving exactly `x^-6`, `e^{-2kx}/x^4`, `e^{-4kx}/x^2`, (7) the outgoing-port derivative identity `d_Pi V_mix = -1/2 T_J^2`, (8) the linear outgoing correction `delta V^(1) = -i/2 Gamma T_J^2`, and (9) the phase-lag consequence `Re delta V^(1) = 0`, `P_abs^(1) >= 0`. The Part VII appendix (eqs. app-part07-dyn-* through app-part07-phase-lag-no-real) reproduces the same identities and frames them as the "linear dynamic no-free-lunch theorem."

## What the script claims to verify

The SymPy script constructs the dynamic `3x3` stiffness matrix `K_dyn` from the dynamic coefficients `K_B(omega)`, `A(omega)`, `W(omega)=Omega_W^2-omega^2-Pi`, then asserts, in order: the determinant identity equals `Delta_Pi*D_Pi` (line 65); the static `omega=0, Pi=0` reduction of each coefficient and of `Delta`, `Q`, `D0` back to the carried one-port bundle (lines 69–74); the six independent inverse entries against hand-typed closed forms (line 100); the susceptibility quadratic-form law (line 123); the collinear-source factorization (line 139); the primitive-source product-family regrouping into exactly the three spatial families (line 167); the Pi-derivative identity `d_Pi V = -1/2 T_J^2` (line 179); the linear outgoing correction `-i/2 Gamma T_J0^2` (line 191); the phase-lag split `Re = 0` symbolically (line 194); and a numeric off-pole sample confirming `Delta!=0`, `D!=0`, real part `~0`, nonzero imaginary part, and `P_abs > 0` (lines 235–239). This maps faithfully onto the paper's nine deliverables.

## Paper ↔ script cross-check

| Paper / notes deliverable | Script check | Verdict |
|---|---|---|
| (1) det K_dyn = Delta_Pi D_Pi | line 65 assert | match |
| (2) static omega=0,Pi=0 reduction to one-port bundle | lines 69–74 asserts | match |
| (3) exact inverse entries chi_qq..chi_WW | line 100 assert (all six) | match |
| (4) susceptibility law V_mix = -1/2 J^T Kinv J | line 123 assert | match |
| (5) collinear-source chi_s factorization | line 139 assert | match |
| (6) primitive product-family (x^-6, e^-2kx/x^4, e^-4kx/x^2) | line 167 assert | match |
| (7) outgoing-port derivative d_Pi V = -1/2 T_J^2 | line 179 assert | match |
| (8) linear outgoing correction -i/2 Gamma T_J^2 | line 191 assert | match |
| (9) phase-lag Re delta V^(1)=0; P_abs^(1) >= 0 | line 194 (symbolic Re=0); lines 237–239 (sample) | match (P_abs>=0 only sampled, see F2) |
| Independent second engine (dual-engine policy) | — | missing (F1) |

`paper_alignment: aligned` — every paper-side deliverable has a faithful, non-tautological SymPy check, and the constants/forms used match the card, notes, and appendix. The only gaps are the absent second engine (F1) and a sample-only positivity check (F2).

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 65 | `simplify(det - Delta_Pi*D_Pi)==0` | claim 1 | yes |
| A2 | sympy | 69–74 | `simplify(... - carried)==0` x6 | claim 2 | yes |
| A3 | sympy | 100 | `all(... == 0)` over 6 inverse entries | claim 3 | yes |
| A4 | sympy | 123 | `simplify(deltaV_matrix - deltaV_formula)==0` | claim 4 | yes |
| A5 | sympy | 139 | `simplify(deltaV_col + 1/2 chi_s S^2)==0` | claim 5 | yes |
| A6 | sympy | 167 | `simplify(deltaV_primitive - expected)==0` | claim 6 | yes |
| A7 | sympy | 179 | `simplify(dVdPi + 1/2 T_J^2)==0` | claim 7 | yes |
| A8 | sympy | 191 | `simplify(deltaV_linear_out - expected)==0` | claim 8 | yes |
| A9 | sympy | 194 | `real_part == 0` (symbolic) | claim 9 (Re=0) | yes |
| A10 | sympy | 235–239 | numeric: Delta!=0, D!=0, Re~0, Im>0, P_abs>0 | claim 9 (off-pole positivity) | partial (sample only) |

All ten assertions trace to a paper-side deliverable; none is orphaned scaffolding. No `paper_missing_script_claim` extras.

## Findings

### F1 — missing_verification_script (missing_mathematica)

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage220_dynamic_mixed_port_kernel_phase_lag_no_go_and_resonant_survival_gate_sympy_audit.py` (whole unit — only one engine present)

**What's wrong:**
This unit is `is_status_only_candidate: False` and `is_checkpoint: False`, so the dual-engine policy applies at full strength: a Mathematica `.wl` is required wherever Mathematica *can* independently verify the stage. Every deliverable here is pure symbolic linear algebra / complex analysis: a `3x3` symbolic determinant, six entries of a `3x3` symbolic inverse, a quadratic-form regrouping, a symbolic derivative in `Pi`, and a real/imaginary split of `-i/2 Gamma T_J0^2`. Mathematica handles all of these natively (`Det`, `Inverse`, `D[...,Pi]`, `ComplexExpand`, `Re`/`Im`, `Coefficient`/`Series` for the spatial-family extraction). The stage card itself records "Mathematica audit: none yet" (stage_220.tex:11). There is no impossibility justification, so this is a genuine `missing_mathematica`, not a legitimate single-engine carve-out.

**Why this matters:**
The phase-lag no-go (`Re delta V^(1) = 0`) and the determinant/inverse identities are the load-bearing claims of the whole dynamic-mixed-port lane (anchor MTDC-T11.1) and are cited downstream by Stage 221's pole/linewidth normal form. A single-engine proof leaves the algebra (a manually transcribed `3x3` inverse compared against `K_dyn.inv()`) unconfirmed by an independent derivation; an independent engine using a different decomposition is the policy's transliteration guard.

**Required change:**
Create the second-engine `.wl` per the directive's claim manifest M1–M9, using native Mathematica primitives and a DIFFERENT decomposition than the `.py` (cofactor/adjugate route for the inverse, `Coefficient`/`Series` rather than `together`-regrouping for the spatial families, `ComplexExpand` for the real/imag split). See directive F1.

**Verification:**
`redteam exec-mathematica 220` finds the new `.wl` at the exact Target path, it contains guarded checks M1–M9, and it exits 0 with all PASS lines.

### F2 — insufficient_verification

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage220_dynamic_mixed_port_kernel_phase_lag_no_go_and_resonant_survival_gate_sympy_audit.py:195,239`

**What's wrong:**
Deliverable (9) of the notes and appendix eq. app-part07-phase-lag-no-real both state the absorbed-power result as a *general* inequality, `P_abs^(1)(x,omega) = (omega Gamma / 2) T_J^2 >= 0` (notes line 441–443, appendix line 456–457). The script verifies the `Re = 0` half symbolically (line 194, fully general — good), but the non-negativity half is exercised only at a single numeric sample (`P_abs_sample > 0`, line 239). The symbolic `P_abs` computed at line 195 is manifestly of the form `omega*Gamma*(numerator)^2/(2*(denominator)^2)`, i.e. a square ratio times `omega*Gamma`; its sign is therefore decided entirely by `sign(omega*Gamma)`. The script never asserts the structural fact that makes the inequality general — that `P_abs/(omega*Gamma)` is a perfect square (hence `>= 0` for any `omega, Gamma >= 0`). One arbitrary numeric point does not exercise the "off resonance, any drive" generality the paper claims.

**Why this matters:**
The whole point of the phase-lag no-go is that the first outgoing correction is *always* dissipative loading (`P_abs >= 0`), never barrier softening, across the admissible off-pole region — not merely at one lucky sample. A sample-only check would still pass if a sign error elsewhere made `P_abs` negative for some other admissible parameters.

**Required change:**
Add one symbolic assertion near line 195 that `P_abs / (omega*Gamma)` is a perfect square (equivalently that `sp.factor(2*P_abs/(omega*Gamma))` equals `T_J0**2`, or that `P_abs - omega*Gamma/2 * T_J0**2` simplifies to 0). This pins the general non-negativity to the squared structure rather than to a single sample. See directive F2.

**Verification:**
`redteam exec-sympy 220` shows a new symbolic check that `P_abs = omega*Gamma/2 * T_J0^2` (perfect-square form), and the script exits 0.

## Independent-derivation check (Mathematica)

No `.wl` exists, so transliteration cannot yet be assessed. F1's claim manifest mandates an independent decomposition with an explicit anti-transliteration guard; the verifier must reject a line-by-line port of the SymPy `together`/`inv()` choreography.

## Engine cross-check

Only one engine is present; `engines_agree: n/a`.

## Verdict justification

The SymPy script is mathematically sound and aligned with the paper. I attacked each assertion: (a) the determinant/inverse/susceptibility/collinear/product-family checks all compare an independently-computed object (`K_dyn.det()`, `K_dyn.inv()`, the matrix quadratic form) against hand-typed closed forms, so they are non-tautological and would catch a wrong coefficient. (b) The phase-lag `Re = 0` check is genuinely load-bearing, not tautological: it is true only because `T_J0` (the transfer factor at `Pi=0`) is real, which follows from the real symbol domains AND the `Pi -> 0` conservative-branch substitution; had `Pi` been retained complex or `omega` been complex, the real part would not vanish — exactly the physical "off the conservative branch" caveat. (c) Symbol domains are correct: `Pi` is left unrestricted (correct, since `Pi=iGamma` later), `Gamma` is nonnegative real, all structural constants real-nonzero, and `omega`, sources real. (d) The static-reduction comment cites "Stage 219," matching the paper card's `\stagefield{Inputs}` ("Imports Stage 219") and the appendix's Stage-219 static-closure theorem; the notes prose internally calls the same upstream "Stage 253"/"Stage 237," which is a notes-side renumbering artifact, not a script defect — informational only, no `paper_misalignment` raised because the script and the authoritative paper card agree. The verdict is `findings` (not clean) solely because the dual-engine policy requires the absent `.wl` (F1) and the general non-negativity of `P_abs` is only sampled, not pinned to its square structure (F2). Neither finding is a math error; neither is `paper_misalignment`; no stop-cold.

## Self-test notes

I checked the variable-independence trap on the only derivative in the script (`sp.diff(deltaV_formula, Pi)`, line 178): `deltaV_formula` genuinely depends on `Pi` through `W(omega)=Omega_W^2-omega^2-Pi` inside every chi entry, so the derivative is not identically zero and the `d_Pi V = -1/2 T_J^2` assertion is substantive. I checked the complex real/imag split: `deltaV_linear_expected = -i/2 Gamma T_J0^2` with `T_J0` real (all symbols real after `Pi->0`), so `Re=0` is a real consequence of the squared real transfer factor, not a domain accident. Trivial-case for F2: substituting the script's own sample gives `P_abs = omega*Gamma/2 * T_J0^2 = (1/2)(1/10)/2 * (real)^2 > 0`, matching the printed `0.0003378`, confirming the proposed perfect-square assertion is true and the fix is safe. Paper round-trip for both fixes: F1 manifest reuses only paper-stated forms; F2 adds a symbolic restatement of appendix eq. app-part07-phase-lag-no-real with no new constants — neither introduces a new paper_misalignment.

---
unit_id: 011
batch: I.1
auditor_model: claude-opus-4-8
audit_date: 2026-06-04T00:00:00Z
verdict: clean
stop_cold: null
findings_count: 0
paper_alignment: aligned
scripts_checked:
  sympy: present
  mathematica: present
  engines_agree: true
  outputs_fresh: true
docs_read:
  paper_stage_tex: present
  notes_stage_files: []
  paper_appendix: present
---

# Audit unit 011 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_011.tex`
- notes: `(none)` (no files matching `notes/stages/moving_throat_pde_stage011_*.md` exist)
- part appendix: `/var/projects/toy_projects/.../stage_appendix_part01.tex` → actual path `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part01.tex` (row 011 at line 44, plus `\input{stages/stage_011}` at line 101)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage011_projected_maxwell_p2_bridge_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage011_projected_maxwell_p2_bridge_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage011_projected_maxwell_p2_bridge_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage011_projected_maxwell_p2_bridge_mathematica_audit.txt`

## What the paper claims

Stage 011 is the "Projected Maxwell \(P_2\) bridge" (Part I, anchor MTDC-T4, status `\StatusExactClosure{}`). It defines the grouped quantities \(S = M + B_2 + Z_2\) and \(T = B_4 + Z_4\) (eq:stage011-compat-base, tex L16). After eliminating the static wall stiffness \(K\), the isotropic compatibility surface is \(\mathcal C = N_0/P_{0,\rm target} - 3S^2/T\) (eq:stage011-compat, tex L21–24). Its projected first variation is \(\delta\mathcal C = n_0/P_{0,\rm target} - 6Sz_2/T + 3S^2z_4/T^2\) (eq:stage011-dcompat, tex L26–32). The static conservative shift \(z_0\) cancels from \(\delta\mathcal C\): it changes \(P_0\) but not the eliminated compatibility surface (tex L33–35). The Output line (tex L37–39) states verbatim: "Stage~011 exports the \(z_0\)-cancellation and the projected \(P_2\) compatibility bridge \eqref{eq:stage011-dcompat}." Distinct deliverables: (D1) the K-eliminated compatibility surface \(\mathcal C\); (D2) its first variation \(\delta\mathcal C\); (D3) the \(z_0\)-cancellation in \(\delta\mathcal C\). Appendix row 011 (L44) summarizes it as the "\(P_2\)-lane bridge from projected electromagnetic slots into the grouped normalization language." No notes file exists; the `.tex` is the sole prose carrier.

## What the script claims to verify

Both engines independently solve two surfaces for the wall stiffness \(K\): a one-pole surface from \((K-B_0-(Z_{0\rm slot}+\epsilon z_0))(T+\epsilon z_4)=3(S+\epsilon z_2)^2\) and a fixed-target normalization surface from \((N_0+\epsilon n_0)/(K-B_0-(Z_{0\rm slot}+\epsilon z_0))=P_{\rm target}\). They then assert that the difference of the two solved \(K\)-surfaces equals the independent closed form \((N_0+\epsilon n_0)/P_{\rm target}-3(S+\epsilon z_2)^2/(T+\epsilon z_4)\) (SymPy A1 / Math M1), that the two surfaces' first \(\epsilon\)-variations equal hand-written targets (A2/A3, M2/M3), that the compatibility first variation reduces to \(n_0/P_{\rm target}-6Sz_2/T+3S^2z_4/T^2\) (A4/A5/A6, M4), and that this variation is \(z_0\)-independent (A7, M5). The docstring states the audit checks "only the paper-card output." This matches deliverables D1–D3 exactly.

## Paper ↔ script cross-check

| Paper deliverable | Script-side check | Status |
|---|---|---|
| D1: \(\mathcal C = N_0/P_{\rm target}-3S^2/T\) (tex L23) | A1 (.py L59) / M1 (.wl L38–41): difference of two independently-solved K-surfaces == independent closed form; printed `Compat` (output L11) = `(N0*T - 3*Ptarget*S^2)/(Ptarget*T)` | match |
| D2: \(\delta\mathcal C = n_0/P_{\rm target}-6Sz_2/T+3S^2z_4/T^2\) (tex L28–32) | A6 (.py L64) / M4 (.wl L56–59) match the hand-written target; A4/A5 cross-check the solved-K route; printed `delta Compat` (output L13) = `3*S^2*z4/T^2 - 6*S*z2/T + n0/Ptarget` | match |
| D3: \(z_0\)-cancellation (tex L33–35, L38) | Math M5 (.wl L60) genuinely tests it (differentiates the solved-difference shift, which carries \(z_0\)); SymPy A2+A3+A4 demonstrate the \(z_0\) terms cancel between the two surfaces; SymPy A7 is a weak/redundant restatement | match |

`paper_alignment: aligned`. Every paper deliverable has a faithful, non-tautological script-side check; no `extra` or `mismatch` rows.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 59 | `assert_zero((K_norm_p - K_one_pole_p) - compat_direct_p)` | D1 | yes |
| A2 | sympy | 60 | `assert_zero(dK_one_pole - (z0 + 6 S z2/T - 3 S^2 z4/T^2))` | D3 support (intermediate) | yes |
| A3 | sympy | 61 | `assert_zero(dK_norm - (z0 + n0/Ptarget))` | D3 support (intermediate) | yes |
| A4 | sympy | 62 | `assert_zero(d_compat - (dK_norm - dK_one_pole))` | D2/D3 (route consistency) | yes (mild) |
| A5 | sympy | 63 | `assert_zero(d_compat - d_compat_direct)` | D2 (route consistency) | yes |
| A6 | sympy | 64 | `assert_zero(d_compat_direct - (n0/Ptarget - 6 S z2/T + 3 S^2 z4/T^2))` | D2 (load-bearing) | yes |
| A7 | sympy | 65 | `assert_zero(sp.diff(d_compat_direct, z0))` | D3 | partial (trivially true; see Verdict) |
| M1 | math | 38–41 | `check(compatFixed - compatFixedClosed)` | D1 | yes |
| M2 | math | 44–47 | `check(deltaKPole - (z0 + 6 S z2/T - 3 S^2 z4/T^2))` | D3 support | yes |
| M3 | math | 49–53 | `check(deltaKNorm - (z0 + n0/Ptarget))` | D3 support | yes |
| M4 | math | 56–59 | `check(compatFixedShift - (n0/Ptarget - 6 S z2/T + 3 S^2 z4/T^2))` | D2 (load-bearing) | yes |
| M5 | math | 60 | `check(D[compatFixedShift, z0])` | D3 | yes (genuine; operates on solved-difference shift) |

## Findings

None. See attacks attempted in Verdict justification.

(Note for the record, not a finding: SymPy A7 `sp.diff(d_compat_direct, z0)` is trivially zero because `d_compat_direct` is constructed from a closed form (.py L54) that never contained `z0`; differentiating it w.r.t. `z0` cannot fail. By itself A7 does not exercise the \(z_0\)-cancellation. However, this does NOT rise to an `insufficient_verification` finding because the \(z_0\)-cancellation deliverable (D3) IS genuinely exercised elsewhere: (i) SymPy A2 and A3 assert each K-surface's first variation explicitly CARRIES `z0` (`z0 + ...`), and A4 asserts their difference equals the `z0`-free `dK_norm - dK_one_pole` — a non-trivial demonstration that the `z0` terms cancel; (ii) Mathematica M5 differentiates `compatFixedShift`, which is derived from `firstShift[kNorm - kPole]` where both `kNorm` and `kPole` carry `(Z0slot + eps z0)` — so if `z0` did NOT cancel, M5's residual would be nonzero. The cross-engine coverage of D3 is therefore complete; A7 is harmless redundancy, not a gap.)

## Independent-derivation check (Mathematica)

The `.wl` is a genuinely independent re-derivation, not a transliteration:
- First-variation extraction differs: SymPy uses `lin(expr) = expand(series(expr, eps, 0, 2).removeO())` then `(lin - base)/eps` (.py L37–38, L51–57); Mathematica uses `firstShift[expr] = Coefficient[Normal[Series[expr - (expr/.eps->0)]], eps, 1]` (.wl L22–23). Different mechanics for the same coefficient.
- Assertion set differs: SymPy has 7 checks including two distinct consistency routes (A4 solved-K difference vs A5 closed-form route) and the weak A7; Mathematica has 5 (M1–M5) with no separate "d_compat vs dKnorm−dKpole" route and a stronger M5 z0 test.
- The z0-cancellation is tested differently per engine (Math M5 on the solved-difference shift; SymPy via A2/A3/A4), which is the opposite of a line-by-line port.
- Symbol domains differ: SymPy declares `K,B0,Z0slot,N0,Ptarget,S,T` as `nonzero=True` and `z0,z2,z4,n0` unrestricted (.py L34–35); Mathematica declares all symbols `Reals` with only `T != 0 && Ptarget != 0` (.wl L18–20).

No `mathematica_transliteration` finding.

## Engine cross-check

Both engines reach identical conclusions. SymPy output: `Compat = (N0*T - 3*Ptarget*S^2)/(Ptarget*T)`, `delta Compat = 3*S^2*z4/T^2 - 6*S*z2/T + n0/Ptarget`, `STATUS: PASS`. Mathematica output: M1–M5 each `residual = 0`, `STATUS: PASS`. The SymPy `Compat` printed form `(N0*T - 3*Ptarget*S^2)/(Ptarget*T)` is algebraically identical to `N0/Ptarget - 3 S^2/T`, and the printed `delta Compat` is term-for-term identical to the Mathematica M4 target. `engines_agree: true`. No `engine_disagreement`.

## Value Reconciliation (pass-2 augmentation)

The scripts are fully symbolic; there are no numeric constants, benchmarks, or figure-of-merit numbers. The deliverable-level values are the boxed/closed-form symbolic results.

| value | source (py/wl + output line) | .tex/.md location | status |
|---|---|---|---|
| `C = N0/Ptarget - 3 S^2/T` | sympy .py L54/L72; output L11 `(N0*T - 3*Ptarget*S^2)/(Ptarget*T)`; wl M1 L36–37 | `paper/stages/stage_011.tex:23` (eq:stage011-compat) | MATCH |
| `δC = n0/Ptarget - 6 S z2/T + 3 S^2 z4/T^2` | sympy .py L64/L73; output L13; wl M4 L58 | `paper/stages/stage_011.tex:28-32` (eq:stage011-dcompat) | MATCH |
| z0-cancellation (`∂δC/∂z0 = 0`) | sympy A7 .py L65; output L15; wl M5 L60 | `paper/stages/stage_011.tex:33-35,38` (Output) | MATCH |
| `S = M + B2 + Z2`, `T = B4 + Z4` | kept as opaque symbols `S,T` (sympy L34; wl L16) | `paper/stages/stage_011.tex:16` (eq:stage011-compat-base) | MATCH (script abstracts at S,T level; consistent) |

INTERNAL items (genuine scaffolding, no prose expected, no finding):
- `K_pole(eps)` full solved form (output L7) — intermediate solved surface.
- `K_norm(eps)` full solved form (output L9) — intermediate solved surface.
- `dK_one_pole = z0 + 6 S z2/T - 3 S^2 z4/T^2` (sympy A2 L60; wl M2 L46) — intermediate K-surface variation.
- `dK_norm = z0 + n0/Ptarget` (sympy A3 L61; wl M3 L52) — intermediate K-surface variation.
- pass/fail flags, `residual = 0` lines, `STATUS: PASS` — verification scaffolding.

reconciliation: complete; 4 deliverable values checked, 0 misaligned. (No notes file exists for this stage; the `.tex` card carries all three stated deliverables and each reconciles.)

## Verdict justification

`clean`. Attacks attempted that failed:
1. **Tautology on A1/M1** — the two K-surfaces are obtained from two *independent* `Solve`/`sp.solve` calls (one-pole product equation vs normalization ratio equation), and their difference is checked against a *separately written* closed form; this is a substantive identity, not a definition-then-assert. The printed `Compat` confirms a nonzero, correctly-structured result.
2. **z0-cancellation coverage** — A7 alone is trivially true (differentiates a z0-free closed form), but I confirmed D3 is genuinely exercised by Mathematica M5 (operates on the solved-difference shift carrying z0) and by SymPy A2+A3+A4 (z0 present in each surface's variation, cancels in their difference). No coverage gap, so no finding.
3. **Missing branch** — both K-equations are linear in K (the product expands to a linear-in-K equation; the ratio is linear in K after clearing the denominator), so each has a unique root; `[0]`/`First[]` selects the only solution. No hidden branch.
4. **Symbol-assumption masking** — `nonzero`/`!= 0` assumptions on `S,T,Ptarget` are exactly those needed for the denominators `T`, `Ptarget` and are justified by the physical setup; the perturbation symbols `z0,z2,z4,n0` are correctly left unrestricted. No assumption hides a sign or branch error.
5. **Engine disagreement** — both engines independently produce identical final forms via different first-variation extraction mechanics and different assumption sets.
6. **Paper alignment** — I read the paper card, confirmed no notes file exists, and read appendix row 011; the script's verified identities match `eq:stage011-compat`, `eq:stage011-dcompat`, and the Output's z0-cancellation verbatim. All four deliverable values reconcile.

Outputs are fresh (both `.txt` mtimes are ~15 h newer than their scripts). No `stale_output`.

## Self-test notes

Trap 1 (variable independence in `sp.diff`/`D`): flagged A7/M5 explicitly — A7's `sp.diff(d_compat_direct, z0)` is identically zero because the expression has no `z0` dependence; confirmed this is harmless because D3 is covered by M5 (which DOES carry z0 in its operand) and by A2/A3/A4. Trap 2 (symmetry/parity): no unbounded integrals in this stage; not applicable. Trap 3 (trivial-case): mentally substituted `z0=z2=z4=n0=0` (the surfaces reduce to `K = B0+Z0slot+N0/Ptarget` and `K = B0+Z0slot+3S^2/T`, difference `= N0/Ptarget - 3S^2/T = Compat`, consistent) and confirmed each first-variation target is the correct O(eps) coefficient; residuals reduce to 0 as the outputs report. No directive is written (zero findings).

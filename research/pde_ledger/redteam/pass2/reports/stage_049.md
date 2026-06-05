---
unit_id: 049
batch: III.2
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-06-05T00:00:00Z
verdict: findings
stop_cold: null
findings_count: 1
paper_alignment: aligned
scripts_checked:
  sympy: present
  mathematica: present
  engines_agree: true
  outputs_fresh: false
docs_read:
  paper_stage_tex: present
  notes_stage_files: [/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage049_dn_overlap_zeta.md]
  paper_appendix: present
---

# Audit unit 049 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_049.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage049_dn_overlap_zeta.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part03.tex` (row at line 76)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage049_dn_overlap_zeta_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage049_dn_overlap_zeta_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage049_dn_overlap_zeta_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage049_dn_overlap_zeta_mathematica_audit.txt`

## What the paper claims

Stage 049 makes the coherent support ratio `zeta` microscopic by evaluating D/N support-tower overlaps. The boxed deliverables are: (1) the physical coherent support ratio `\zeta_n^{phys} = (K_W^{eff}/K_{phi,n}^{eff})(I_n/I_0)^2` (eq:app-stage049-zeta-phys); (2) for a uniform local source density `sigma=1`, the overlap ratio `I_n/I_0 = 1/(2n+1)` (eq:app-stage049-In-ratio); (3) the same-operator twin-family ratio `\zeta_n^{twin} = 1/[(2n+1)^2(1+x n(n+1))]` with `x = pi^2 T_X/(L^2 K_W^{eff})` (eq:app-stage049-zeta-twin); and (4) the lowest twin value `\zeta_0^{twin}=1` (eq:app-stage049-lowest-twin). `\stagefield{Output}` states: "Microscopic D/N support ratios \eqref{eq:app-stage049-zeta-phys}--\eqref{eq:app-stage049-zeta-twin}." The notes add the explicit overlap closed form `I_n = sqrt(2L)/[(n+1/2)pi]`, the D/N mode `chi_n(s)=sqrt(2/L)sin[(n+1/2)pi s/L]`, and the twin stiffness data `K_W^{eff}=K_X+pi^2 T_X/(4L^2)`, `K_{phi,n}^{eff}=K_X+(n+1/2)^2 pi^2 T_X/L^2`. The appendix row (line 76) summarizes the same three deliverables and marks the stage exact-closure.

## What the script claims to verify

Both scripts (SymPy and Mathematica, identical assertion choreography) verify seven items, in order: (1) `cos(k_n L)=0` with `k_n=(n+1/2)pi/L` (the D/N half-wave boundary condition); (2) the direct integral `int_0^L chi_n ds` equals the closed form `sqrt(2L)/[(n+1/2)pi]`; (3) `I_n/I_0 = 1/(2n+1)`; (4) `zeta_phys` built from `lambda_W=lambda_* I_0`, `lambda_{phi,n}=lambda_* I_n` reduces to `(K_W/K_phi)(I_n/I_0)^2`; (5) the twin stiffness relation `K_{phi,n}^{eff} = K_W^{eff} + pi^2 T_X n(n+1)/L^2`; (6) `zeta_twin = (K_W/K_phi)(1/(2n+1))^2` equals `1/[(2n+1)^2(1+x n(n+1))]` with `x=pi^2 T_X/(L^2 K_W)`; (7) `zeta_twin|_{n=0} = 1`. These map 1:1 onto the paper's deliverables.

## Paper ↔ script cross-check

| Paper deliverable | Script check | Status |
|---|---|---|
| `chi_n`, `k_n=(n+1/2)pi/L`, D/N boundary | `cos(k_n L)=0` (py 65-68 / wl 50) | match |
| `I_n = sqrt(2L)/[(n+1/2)pi]` (notes) | direct integral vs closed form (py 70-74 / wl 52-56) | match |
| `I_n/I_0 = 1/(2n+1)` (eq In-ratio) | ratio check (py 76-80 / wl 58-62) | match |
| `zeta_n^{phys}=(K_W/K_phi)(I_n/I_0)^2` (eq zeta-phys) | from lambda_* definitions (py 82-87 / wl 64-69) | match |
| twin stiffness `K_phi = K_W + pi^2 T_X n(n+1)/L^2` (notes) | stiffness relation (py 89-95 / wl 71-74) | match |
| `zeta_n^{twin}=1/[(2n+1)^2(1+x n(n+1))]`, `x=pi^2 T_X/(L^2 K_W)` (eq zeta-twin) | twin ratio + x (py 97-102 / wl 76-81) | match |
| `zeta_0^{twin}=1` (eq lowest-twin) | n=0 substitution (py 103 / wl 82) | match |

paper_alignment: aligned. Every boxed paper deliverable has a faithful, non-tautological script-side check; no extras, no missing claims.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 65-68 | `expect_zero(cos(k_n L))` | D/N boundary | yes |
| A2 | sympy | 74 | `expect_zero(integral - closed_form)` | I_n closed form | yes |
| A3 | sympy | 80 | `expect_zero(ratio - 1/(2n+1))` | I_n/I_0 | yes |
| A4 | sympy | 87 | `expect_zero(zeta_phys - expected)` | zeta_n^phys | yes |
| A5 | sympy | 92-95 | `expect_zero(K_phi - (K_W + gap))` | twin stiffness | yes |
| A6 | sympy | 102 | `expect_zero(zeta_twin - twin_formula)` | zeta_n^twin | yes |
| A7 | sympy | 103 | `expect_zero(zeta_twin(n=0) - 1)` | zeta_0^twin=1 | yes |
| B1-B7 | mathematica | 50,56,62,69,74,81,82 | `expectZero[...]` | same as A1-A7 | yes |

## Findings

### F1 — stale_output

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage049_dn_overlap_zeta_sympy_audit.txt:3` (banner "STAGE 32")
- `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage049_dn_overlap_zeta_mathematica_audit.txt:3,36` (banner "STAGE 32"; final line "Stage 32 Mathematica audit passed.")
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage049_dn_overlap_zeta_mathematica_audit.wl:93` (live `Print["Stage 32 Mathematica audit passed."]`)
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage049_dn_overlap_zeta_sympy_audit.py:2` (docstring `"""Stage 32 SymPy audit.`)

**What's wrong:**
Two distinct sub-issues, both low-severity numbering artifacts:
1. The committed `.txt` outputs are stale by mtime — sympy `.txt` (2026-05-26 02:55) and mathematica `.txt` (2026-05-26 02:55) both predate their scripts (2026-06-03 15:59). The visible content effect is the banner: scripts now `banner("STAGE 49 — ...")` (py:51 / wl:36) but both saved outputs still print `STAGE 32 — ...` (output line 3). All numeric/symbolic result lines in the saved outputs still agree with what the current scripts produce (the math body did not change), so only the banner line is content-divergent.
2. Stale `Stage 32` / `STAGE 32` self-labels remain *in the live source*: py docstring `"""Stage 32 SymPy audit.` (py:2) and wl `Print["Stage 32 Mathematica audit passed."]` (wl:93). These are the known low-severity `stale_output`-class numbering self-labels (the file is stage 049). They are not blocking; the orchestrator resolves label scope separately.

**Why this matters:**
The stale banner makes the saved transcript appear to belong to "Stage 32". The result lines are correct, so no math is affected, but a fresh run is needed so the committed transcript matches the current banner, and the in-source `Stage 32` labels should be corrected to `Stage 49`/`049` for traceability.

**Required change:**
- Refresh both `.txt` outputs by re-running the scripts (orchestrator's independent re-run handles this).
- Fix the live self-labels: py:2 `"""Stage 32 SymPy audit.` → `Stage 49 SymPy audit.`; wl:93 `Print["Stage 32 Mathematica audit passed."]` → `Print["Stage 49 Mathematica audit passed."]`.

**Verification:**
After re-run, output line 3 reads `STAGE 49 — EXPLICIT D/N OVERLAP ...` in both `.txt`, and the mathematica `.txt` final line reads `Stage 49 Mathematica audit passed.`; py docstring line 2 and wl line 93 both say `Stage 49`/`049`.

## Independent-derivation check (Mathematica)

The `.wl` mirrors the `.py` closely in structure: identical helper names (`dnHalfwaveMomentum`/`dn_halfwave_momentum`, `overlapRatio`/`overlap_ratio`, `twinSupportRatio`/`twin_support_ratio`), identical assertion order, identical variable choreography. This is a borderline transliteration pattern. However, the single load-bearing computation — the overlap integral `Integrate[chiN,{s,0,l}]` (wl:52) vs `sp.integrate(chi_n,(s,0,L))` (py:70) — is genuinely and independently evaluated by each CAS rather than copied; the algebraic identities (twin stiffness, zeta reductions) are also re-simplified independently by each engine's `FullSimplify`/`simplify`. The deliverable that could realistically diverge between engines (the closed-form integral) is computed from scratch on both sides and agrees. I therefore note the structural parallelism but do not raise a `mathematica_transliteration` finding: the second engine does independently derive the load-bearing result.

## Engine cross-check

Both engines agree on every check (all `... = 0`, all `PASS`). Final forms match modulo CAS normalization:
- `I_n` integral: sympy `2*sqrt(2)*sqrt(L)/(pi*(2*n+1))` (line 8) vs mathematica `(2*Sqrt[2]*Sqrt[l])/(Pi+2*n*Pi)` (line 9) — identical.
- `zeta_n^(phys)`: sympy `K_W_eff/(K_phi_eff*(2*n+1)**2)` (line 14) vs mathematica `kWeff/(kPhiEff*(1+2*n)^2)` (line 17) — identical.
- `x`: sympy `4*pi**2*T_X/(4*K_X*L**2+pi**2*T_X)` (line 17) vs mathematica `(4*Pi^2*tX)/(4*kX*l^2+Pi^2*tX)` (line 22) — identical.
- `zeta_n^(twin)`: sympy expanded denominator form (line 18) vs mathematica `1/((1+2*n)^2*(1+(4*n*(1+n)*Pi^2*tX)/(4*kX*l^2+Pi^2*tX)))` (line 23) — algebraically identical.
No disagreement.

## Verdict justification

The seven assertions map exactly onto the four boxed paper deliverables plus the supporting overlap/stiffness data in the notes; each is non-tautological and well-anchored. I hand-verified the math: the D/N boundary `cos((n+1/2)pi)=0`, the integral `sqrt(2/L)[1-cos((n+1/2)pi)]/k = sqrt(2L)/[(n+1/2)pi]`, the ratio `1/(2n+1)`, the lambda_* cancellation in `zeta_phys`, the twin stiffness gap `(n+1/2)^2 - 1/4 = n(n+1)`, the `K_W/K_phi = 1/(1+x n(n+1))` reduction, and `zeta_0^{twin}=1` — all correct, and both engines confirm. Attacks tried and failed: (a) the `zeta_phys` check is not tautological because `lambda_W`/`lambda_phi_n` are built from the integral-derived overlaps while `zeta_phys_expected` is built from `overlap_ratio`, and the equality also exercises the `lambda_*` cancellation; (b) the boundary check genuinely uses `n integer` (would fail for non-integer n); (c) no over-strong assumptions hide a branch — `n` nonnegative integer, `L,K,T` positive are all justified by the physical setup. The only finding is the low-severity stale_output (stale `.txt` banner + in-source `Stage 32` labels). I confirm I read the paper card, notes, and appendix row; the script's claim matches the paper. Verdict: `findings` (one low-severity stale_output only); `material_change: false`.

## Value Reconciliation (pass-2 augmentation)

Deliverable-level reconciliation of every result value the scripts emit:

| value | source (py/wl + output line) | .tex/.md location | status |
|---|---|---|---|
| `chi_n = sqrt(2/L) sin[(n+1/2)pi s/L]` | py:61 / wl:46-47; out py:5, wl:5 | tex:7; md:24,47 | MATCH |
| `k_n = (n+1/2)pi/L` | py:60 / wl:45; out py:6, wl:6 | md:51 (`k_n=(n+1/2)pi/L`) | MATCH |
| `I_n = sqrt(2L)/[(n+1/2)pi]` | py:71,106 / wl:53; out py:9,23, wl:9 | md:98 (`= sqrt(2L)/[(n+1/2)pi]`) | MATCH |
| `I_0 = 2 sqrt(2) sqrt(L)/pi` | py:76; out py:11, wl:13 | (n=0 case of I_n in md:98) | MATCH |
| `I_n/I_0 = 1/(2n+1)` | py:80,107 / wl:62; out py:12, wl:14 | tex:25 (eq In-ratio); md:32,102 | MATCH |
| `zeta_n^(phys) = (K_W/K_phi)(I_n/I_0)^2` = `K_W/(K_phi(2n+1)^2)` | py:86,108 / wl:68; out py:14, wl:17 | tex:17-19 (eq zeta-phys); md:36-37,83,107 | MATCH |
| twin stiffness `K_phi = K_W + pi^2 T_X n(n+1)/L^2` | py:89-91 / wl:71-73; out py:16, wl:20 | md:119-122 | MATCH |
| `x = pi^2 T_X/(L^2 K_W)` (= `4 pi^2 T_X/(4 K_X L^2 + pi^2 T_X)`) | py:97,100 / wl:76; out py:17, wl:22 | tex:33; md:131 | MATCH |
| `zeta_n^(twin) = 1/[(2n+1)^2(1+x n(n+1))]` | py:99,101,109 / wl:78-80; out py:18, wl:23 | tex:31 (eq zeta-twin); md:127,216 | MATCH |
| `zeta_0^(twin) = 1` | py:103,110 / wl:82; out py:20, wl:20 | tex:39 (eq lowest-twin); md:139,218 | MATCH |

INTERNAL (verification scaffolding / intermediate, no finding): `lambda_W = lambda_* I_0`, `lambda_phi_n = lambda_* I_n` (cancellation intermediates), `cos(k_n L)` boundary residual, `twin_gap` intermediate, all `expect_zero`/`PASS` residual flags.

All ten emitted deliverable values reconcile against the `.tex` card and/or `.md` notes. The only stale artifact is the `STAGE 32` banner in the saved `.txt` (folded into F1 stale_output, not a value mismatch — every result line in the saved outputs still equals the current value).

reconciliation: complete; 10 values checked, 0 misaligned

## Self-test notes

Checked: (1) variable independence — no `sp.diff`/`D` derivatives in this stage, so the zero-derivative trap does not apply. (2) Symmetry/parity — the only integral `int_0^L chi_n ds` is over a bounded interval `[0,L]` with a definite nonzero value `sqrt(2L)/[(n+1/2)pi]`; no symmetric-domain vanishing claim is involved. (3) Trivial-case pre-check — substituting `n=0`: `I_0/I_0=1=1/(2·0+1)` ✓, `zeta_0^{twin}=1/(1·(1+0))=1` ✓; the boundary `cos((1/2)pi)=0` ✓. (4) No missing-script finding (both engines present). (5) Paper round-trip — F1 prescribes only output-refresh and label fixes; no math change, so no new paper_misalignment is introduced. Conclusion: math is fully sound, both engines agree, the single finding is the low-severity stale banner/labels.
</content>
</invoke>

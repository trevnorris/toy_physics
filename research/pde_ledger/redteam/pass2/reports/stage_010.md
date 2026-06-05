---
unit_id: 010
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

# Audit unit 010 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_010.tex`
- notes: `(none)` (no files matching `notes/stages/moving_throat_pde_stage010_*.md`)
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part01.tex` (row 42 + `\input` at line 99)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage010_projected_maxwell_push_bundle_master_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage010_projected_maxwell_push_bundle_master_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage010_projected_maxwell_push_bundle_master_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage010_projected_maxwell_push_bundle_master_mathematica_audit.txt`

## What the paper claims

Stage 010 is the master transport carrying projected-EM corrections into the grouped bundle slots. Under the slot shift `Z_n -> Z_n + eps z_n`, `N_n -> N_n + eps n_n` (eq. `stage010-projected-shifts`), with `D_0 = K - B_0 - Z_0`, `D_2 = -(M+B_2+Z_2)`, `D_4 = -(B_4+Z_4)`, the card boxes the first response variations: `delta u_2` (`stage010-du2`), `delta u_4` (`stage010-du4`), `delta P_0` (`stage010-dp0`), `delta P_2` (`stage010-dP2`), `delta P_4` (`stage010-dP4`). The `\stagefield{Output}` (lines 137-144) verbatim: *"Stage 010 exports the transport map (projected-shifts)--(dP4), the compatibility transports (k-surfaces)--(compat-transport), the weak-axisymmetric lane signature (y20-lambdas)--(weak-axisym-trace), the primitive static prefactor anchor (primitive-xi), and the sign-flip guards (mutation-guards)."* Distinct deliverables: (1) five slot variations `delta u2,u4,P0,P2,P4`; (2) one-pole and fixed-target K-surfaces and their elimination to `delta C_fixed = n0/P_target - 6Sz2/T + 3S^2 z4/T^2` with `d/dz0 = 0`; (3) transported-target `C_tr = D0target - 3(S+eps z2)^2/(T+eps z4)` and `delta C_tr = -6Sz2/T + 3S^2 z4/T^2`; (4) real-Y20 square-overlap lane multipliers `lambda20=1, lambda21=1/2, lambda22=-1` and the resulting trace signature `xbar=x0, a=eps x1/4, b=3 eps x1/4, b=3a`; (5) the primitive static prefactor `z0_prim, n0_prim, Xi_static`; (6) the sign-flip mutation guards `R_fixed = R_tr = 6S^2 z4/T^2 != 0`. No notes file exists, so the .tex card is the sole prose authority; it is unusually detailed and boxes every deliverable.

## What the script claims to verify

Both engines build `eps`-perturbed slots and check that the first `eps`-coefficient of each slot matches the boxed closed forms. SymPy uses `sp.diff(EXPR, eps).subs(eps,0)`; Mathematica uses `Coefficient[Normal[Series[...,{eps,0,1}]], eps, 1]` (genuinely different mechanics). They verify `delta u2/u4/P0/P2/P4`, solve for the one-pole and normalization K-surfaces, confirm the closed forms and their first shifts, eliminate K to obtain `compat_direct`, confirm `delta C_fixed`/`delta C_tr` and that the `z0` term cancels (`sp.diff(dcompat, z0) == 0`), while confirming `dK_norm` retains `z0`. The lane multipliers are derived independently (SymPy via `sympy.physics.wigner.gaunt`; Mathematica via `ThreeJSymbol`), feeding the grouped trace anomaly. The primitive `Xi_static` is checked under `N0sym -> P^2/Delta^2`. Two `assert_nonzero` mutation guards confirm a sign-flipped `3S^2 z4/T^2` term yields the nonzero residual `6S^2 z4/T^2`. This is what the verdict applies to.

## Paper <-> script cross-check

| Paper deliverable | Script-side check | Status |
|---|---|---|
| `delta u2` (du2) | sympy:62 / wl M0a:47-49 | match |
| `delta u4` (du4) | sympy:64-66 / wl M0b:51-57 | match |
| `delta P0` (dp0) | sympy:67 / wl M1:60-62 | match |
| `delta P2` (dP2) | sympy:68-77 / wl M2:64-69 | match |
| `delta P4` (dP4) | sympy:78-92 / wl M3:71-79 | match |
| `K_pole`, `K_norm` (k-surfaces) | sympy:95-104 / wl M4:81-93, M6:102-112 | match |
| `delta C_fixed` + `d/dz0 = 0` (compat-fixed) | sympy:128,138 / wl M9:127-132 (z0-cancel via elimination) | match |
| `C_tr`, `delta C_tr` (compat-transport) | sympy:132,136,139 / wl M11:145-152, M12:154-159 | match |
| `lambda20/21/22` (y20-lambdas) | sympy:154-156 / wl M13:173-185 | match |
| trace signature (weak-axisym-trace) | sympy:161-164 / wl M14:193-205 | match |
| `z0_prim, n0_prim, Xi_static` (primitive-xi) | sympy:170-178 / wl M15:207-215 | match |
| `R_fixed = R_tr = 6S^2 z4/T^2` (mutation-guards) | sympy:141-148 / wl M16:219-227, M17:229-236 | match |

Every boxed deliverable maps to a non-tautological script check exercised in both engines. `paper_alignment: aligned`.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 62 | `assert_zero(du2 - target)` | du2 | yes |
| A2 | sympy | 64-66 | `assert_zero(du4 - target)` | du4 | yes |
| A3 | sympy | 67 | `assert_zero(dP0 - target)` | dp0 | yes |
| A4 | sympy | 68-77 | `assert_zero(dP2 - target)` | dP2 | yes |
| A5 | sympy | 78-92 | `assert_zero(dP4 - target)` | dP4 | yes |
| A6 | sympy | 123 | `assert_zero((K_norm_p - K_one_pole_p) - compat_direct_p)` | k-surfaces elimination | yes |
| A7 | sympy | 124 | `assert_zero(dK_one_pole - (z0+6Sz2/T-3S^2z4/T^2))` | K_pole shift | yes |
| A8 | sympy | 125 | `assert_zero(dK_norm - (z0+n0/Ptarget))` | K_norm shift | yes |
| A9 | sympy | 126-127 | `assert_zero(dcompat - dcompat_direct)` (two routes agree) | compat-fixed | yes |
| A10 | sympy | 128 | `assert_zero(dcompat_direct - (n0/Ptarget-6Sz2/T+3S^2z4/T^2))` | compat-fixed | yes |
| A11 | sympy | 129-133 | `assert_zero` transported-target K & C surfaces | compat-transport | yes |
| A12 | sympy | 134-137 | `assert_zero(dcompat_transport - (-6Sz2/T+3S^2z4/T^2))` | compat-transport | yes |
| A13 | sympy | 138-139 | `assert_zero(diff(dcompat_*, z0))` (z0 cancellation) | compat z0-cancel | yes |
| A14 | sympy | 140 | `assert_nonzero(diff(dK_norm,z0)-0)` | K_norm still carries z0 | yes |
| A15 | sympy | 141-148 | `assert_nonzero` sign-flip mutations | mutation-guards | yes |
| A16 | sympy | 154-156 | `assert_zero(lam2m - {1,1/2,-1})` (Gaunt-derived) | y20-lambdas | yes |
| A17 | sympy | 161-164 | `assert_zero` trace/anomaly/`b=3a` | weak-axisym-trace | yes |
| A18 | sympy | 174-178 | `assert_zero(Xi_static - target)` | primitive-xi | yes |
| A19 | mathematica | M0a-M15 | `If[FullSimplify[res] =!= 0, Exit[1]]` mirror of A1-A18 | all above | yes |
| A20 | mathematica | M16/M17 | `If[FullSimplify[res] === 0, Exit[1]]` mutation guards | mutation-guards | yes |

No tautological or unanchored rows. The slot base forms `P2p`/`P4p` (sympy:53-54, wl:34-39) are the stage's definitions of the slots; only their first variations are asserted as deliverables, and those match the paper.

## Findings

None.

## Independent-derivation check (Mathematica)

The `.wl` is not a transliteration. Three concrete divergences in method:
1. Linear coefficients: SymPy uses `sp.diff(EXPR, eps).subs(eps, 0)` (sympy:56-60); Mathematica uses `Coefficient[Normal[Series[..., {eps, 0, 1}]], eps, 1]` (wl:41-45). Different operators for the same quantity.
2. Lane multipliers: SymPy imports `sympy.physics.wigner.gaunt` (sympy:6,21-28); Mathematica builds `gauntByThreeJ` from `ThreeJSymbol` and the closed Gaunt prefactor `Sqrt[(2l1+1)(2l2+1)(2l3+1)/(4 Pi)]` (wl:161-167). The Mathematica side re-derives Gaunt from 3j symbols rather than calling a library Gaunt.
3. Compatibility surfaces: SymPy carries `sp.series(...).removeO()` then `diff` (sympy:105-122); Mathematica uses `Series`+`Coefficient` and `FullSimplify` on `Solve` outputs (wl:81-159). Both `Solve` for K independently and verify the same closed surface, but via different post-processing.
Both engines do share the slot base forms `P2/P4` and `u2/u4` (these are the stage's slot *definitions*, not derived results), which is appropriate.

## Engine cross-check

Engines agree. Mathematica saved output prints `M0a..M15 residual = 0`, `M13 = {0,0,0,0,0}`, `M14 = {0,0,0,0}`, and the mutation residuals `M16 = M17 = (6*S^2*z4)/T^2`, `STATUS: PASS`. SymPy output prints `STATUS: PASS` (its `assert_zero`/`assert_nonzero` raise on failure, so the clean transcript implies all residuals vanished / mutations stayed nonzero). The Mathematica mutation residual `(6*S^2*z4)/T^2` matches the paper's stated guard `R_fixed = R_tr = 6S^2 z4/T^2` (eq. `stage010-mutation-guards`) exactly. Outputs are fresh (sympy script mtime May 25 02:19, output May 25 17:25; wl script May 25 02:15, output May 25 17:29 — both outputs newer than their scripts).

## Verdict justification

Clean. I attacked every boxed identity by hand: `delta u4` (verified the quotient-rule derivative of `D2p^2/D0p^2 - D4p/D0p` term-by-term equals `(D0^2 z4 - D0(2D2z2+D4z0)+2D2^2z0)/D0^3`), `delta P0`/`P2`/`P4` (matched expanded targets to the boxed forms), the one-pole K shift (`z0+6Sz2/T-3S^2z4/T^2`), the fixed/transported compatibility shifts (sign on the `3S^2 z4/T^2` term verified positive, matching the card and the mutation guard which flips it to give `6S^2 z4/T^2`), and the weak-axisymmetric lane (recomputed `xbar=x0`, `a=eps x1/4`, `b=3 eps x1/4`, `b=3a` from `lambda = {1,1/2,-1}`). All matched. The `z0`-cancellation claim is real: K-surface subtraction cancels the `eps z0` term before differentiation, and the script confirms via `diff(dcompat,z0)==0` while separately confirming `dK_norm` retains `z0`. Mutation guards are substantive (nonzero `6S^2 z4/T^2`). Symbol domains are appropriate (`D0`-family `nonzero`; identities are field-rational and need no reality/positivity). The Mathematica engine derives via Series/Coefficient and ThreeJSymbol, not a port of the SymPy algebra. I read the paper card and the appendix row first and confirm the scripts' verified claims match the paper's deliverables exactly. No findings.

## Value Reconciliation (pass-2 augmentation)

All deliverable values the scripts emit are closed-form symbolic results (the stage has only one set of numeric literals: the lane multipliers `1, 1/2, -1`). Each reconciles to a boxed equation in `stage_010.tex`. No notes file exists; the `.tex` is the sole carrier and it boxes every deliverable.

| value | source (py / wl) | .tex location | status |
|---|---|---|---|
| `delta u2 = (D0 z2 - D2 z0)/D0^2` | py:62 / wl M0a:47 | `stage_010.tex:27` (eq `stage010-du2`) | MATCH |
| `delta u4 = (D0^2 z4 - D0(2D2z2+D4z0)+2D2^2z0)/D0^3` | py:64-66 / wl M0b:51-55 | `stage_010.tex:31-32` (eq `stage010-du4`) | MATCH |
| `delta P0 = (D0 n0 + N0 z0)/D0^2` | py:67 / wl M1:60 | `stage_010.tex:37` (eq `stage010-dp0`) | MATCH |
| `delta P2` (full 3-term form) | py:68-77 / wl M2:64-69 | `stage_010.tex:44-47` (eq `stage010-dP2`) | MATCH |
| `delta P4` (full 4-block form) | py:78-92 / wl M3:71-79 | `stage_010.tex:52-56` (eq `stage010-dP4`) | MATCH |
| `K_pole = B0+Z0+eps z0 + 3(S+eps z2)^2/(T+eps z4)` | wl M4:89-90 | `stage_010.tex:67-68` (eq `stage010-k-surfaces`) | MATCH |
| `dK_one_pole = z0 + 6Sz2/T - 3S^2 z4/T^2` | py:124 / wl M5:98 | `stage_010.tex:67-68` (contained in `K_pole`) | MATCH |
| `K_norm = B0+Z0+eps z0 + (N0+eps n0)/Ptarget` | wl M6:109 | `stage_010.tex:70-71` (eq `stage010-k-surfaces`) | MATCH |
| `dK_norm = z0 + n0/Ptarget` | py:125 / wl M7:116 | `stage_010.tex:70-71` (contained in `K_norm`) | MATCH |
| `delta C_fixed = n0/Ptarget - 6Sz2/T + 3S^2 z4/T^2` | py:128 / wl M9:130 | `stage_010.tex:78-80` (eq `stage010-compat-fixed`) | MATCH |
| `d/dz0 (delta C_fixed) = 0` | py:138 / wl (via M9 elim) | `stage_010.tex:82` (`partial_z0 = 0`) | MATCH |
| `C_tr = D0target - 3(S+eps z2)^2/(T+eps z4)` | py:132 / wl M11:149 | `stage_010.tex:88-90` (eq `stage010-compat-transport`) | MATCH |
| `delta C_tr = -6Sz2/T + 3S^2 z4/T^2` | py:136 / wl M12:157 | `stage_010.tex:92-93` (eq `stage010-compat-transport`) | MATCH |
| `lambda20=1, lambda21=1/2, lambda22=-1` | py:154-156 / wl M13:174-176 | `stage_010.tex:100-102` (eq `stage010-y20-lambdas`) | MATCH |
| `xbar=x0, a=eps x1/4, b=3 eps x1/4, b=3a` | py:161-164 / wl M14:196-199 | `stage_010.tex:108-111` (eq `stage010-weak-axisym-trace`) | MATCH |
| `z0_prim = (Delta q1 - Q d1)/Delta^2` | py:170 / wl:207 | `stage_010.tex:119` (eq `stage010-primitive-xi`) | MATCH |
| `n0_prim = 2P(Delta p1 - P d1)/Delta^3` | py:171 / wl:208 | `stage_010.tex:120-121` (eq `stage010-primitive-xi`) | MATCH |
| `Xi_static = 2p1/P - 2d1/Delta + (Delta q1 - Q d1)/(D0 Delta^2)` | py:174-178 / wl M15:209-213 | `stage_010.tex:122-125` (eq `stage010-primitive-xi`) | MATCH |
| `R_fixed = R_tr = 6 S^2 z4/T^2` (mutation guard) | py:143,147 / wl M16/M17 out:19-20 | `stage_010.tex:130-132` (eq `stage010-mutation-guards`) | MATCH |

INTERNAL (scaffolding, no finding expected): `eps`-perturbed primitives (`D0p..N4p`, `den0..num4`); the slot base forms `u2p,u4p,P0p,P2p,P4p`/`u2slot..slot4`; `compat`/`dcompat` two-route consistency residuals; per-check residual variables (`m0aResidual..m15Residual`); `assert_nonzero("...still carries z0")` probe; same-sign Gaunt vanishing check (`gaunt(2,2,2,0,m,m)=0`); pass/status prints.

reconciliation: complete; 19 deliverable values checked, 0 misaligned.

## Self-test notes

Checked variable-independence (every `sp.diff`/`Series` is taken w.r.t. `eps`, on which every perturbed slot genuinely depends — no identically-zero-derivative trap; the `diff(dcompat,z0)` cancellation checks are correctly `assert_zero` since `z0` truly cancels, while `diff(dK_norm,z0)` is `assert_nonzero` since `z0` survives). Verified parity/symmetry is not applicable (no unbounded integrals; all checks are rational-function identities). Trivial-case substitution confirmed the mutation residual `6S^2 z4/T^2` is genuinely nonzero (so `assert_nonzero` is substantive, not vacuous). Re-derived `du4`, `dK_one_pole`, `delta C_fixed` signs, and the lane trace signature by hand and all matched the boxed card forms.

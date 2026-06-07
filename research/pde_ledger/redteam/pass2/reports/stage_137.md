---
unit_id: 137
batch: IV.4
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-06-06T00:00:00Z
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
  notes_stage_files: [/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage137_core_to_mouth_gain_map.md]
  paper_appendix: present
---

# Audit unit 137 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_137.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage137_core_to_mouth_gain_map.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part04.tex` (only `\input{stages/stage_137}` at line 1308; no extra row content)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage137_core_to_mouth_gain_map_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage137_core_to_mouth_gain_map_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage137_core_to_mouth_gain_map_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage137_core_to_mouth_gain_map_mathematica_audit.txt`

## What the paper claims

Stage 137 derives the explicit coupled mouth-layer gain pair `(M_s, M_q)` from a concrete two-channel throat-core ansatz (GNLS + localized-Maxwell), replacing what were abstract fixed-point coefficients. The card's derivation-ledger quote states the bottom line verbatim: "Explicit gains are \(M_s=Lg_s^2/(K_s\Theta_\sigma)\), \(M_q=-L(K_sg_q-\lambda g_s)^2/[K_s(K_sK_q+\lambda^2)\Theta_\sigma]\)." The notes (Result box, md:126-134) repeat the same boxed pair and supply the intermediate core coefficients `\rho_c=g_s^2/K_s` (md:67) and `\sigma_c=(K_sg_q-\lambda g_s)^2/[K_s(K_sK_q+\lambda^2)]` (md:69-71), with `M_s=(L/\Theta_\sigma)\rho_c`, `M_q=-(L/\Theta_\sigma)\sigma_c`. The notes also state the Family-1 fixed-point law `\Pi=M_s+M_q\,\mathcal S_q(\Pi)` (md:94-101) and the reduced core susceptibility envelope `\rho_c-\sigma_c/(1-\kappa_c z^2-i\gamma_c z^5)` (md:62-63). The card's Checks list three items: (1) gain pair vs outlet consistency, (2) self-matched susceptibility closure, (3) numerical fixed points recorded as numerically located (a documentation/recording directive, not a per-stage symbolic identity).

## What the script claims to verify

Both engines (a) independently reconstruct `\rho_c` and `\sigma_c` by inverting the physical two-channel core stiffness matrix `M_core=[[K_s,\lambda],[\lambda,-K_q D]]` against the mouth coupling vector `v=(g_s,g_q)` (Schur complement `v^T M_core^{-1} v`), taking the `D\to\infty` limit for `\rho_c` and the static `D=1` residue for `\sigma_c`; (b) build `M_s=L\rho_c/\Theta`, `M_q=-L\sigma_c/\Theta` and anchor them against the paper-card closed forms; (c) verify the full-z reduced susceptibility envelope equals the matrix-Schur source on the bare denominator `D_W_bare(z)=1-\kappa_0 z^2-i\gamma_0 z^5` via the `\kappa_c=\kappa_0/(1+r_c)`, `\gamma_c=\gamma_0/(1+r_c)`, `r_c=\lambda^2/(K_sK_q)` maps, plus a static `z\to0` specialization; (d) exercise outlet consistency at a nonzero susceptibility `S_q` with an explicit non-vacuity guard against a sign flip of `M_q`; (e) confirm the `r_c`-notation form of `\sigma_c`. The critical point is that `\rho_c`/`\sigma_c` are NOT inputs to `M_core` — only `K_s,K_q,\lambda,g_s,g_q` are — so step (a) is a genuine independent derivation, not a restatement of the assigned forms.

## Paper ↔ script cross-check

| paper deliverable | script-side check | status |
|---|---|---|
| `M_s = L g_s^2/(K_s\Theta_\sigma)` (tex:16, md:79/129) | py:30/42-46 `Ms==Ms_paper`; wl:51/60-62; underpinned by F1 `rho_c==rho_c_schur` py:26/wl:48 | match |
| `M_q = -L(K_sg_q-\lambda g_s)^2/[K_s(K_sK_q+\lambda^2)\Theta_\sigma]` (tex:16, md:87/132) | py:31/43-45 `Mq==Mq_paper`; wl:52/61-63; underpinned by F1 `sigma_c==sigma_c_schur` py:27/wl:49 | match |
| `\rho_c=g_s^2/K_s` (md:67) | py:24/26 Schur `D\to\infty` residue; wl:46/48 | match |
| `\sigma_c=(K_sg_q-\lambda g_s)^2/[K_s(K_sK_q+\lambda^2)]` (md:69-71) | py:25/27 Schur static `D=1` residue; wl:47/49 | match |
| reduced susceptibility envelope `\rho_c-\sigma_c/(1-\kappa_c z^2-i\gamma_c z^5)` (md:62-63) | py:62-65 full-z match vs matrix source; wl:73-74 | match |
| Family-1 outlet law / Check item 1 (md:94-101; tex:22) | py:79-91 nonzero-`S_q` outlet + non-vacuity guard; wl:82-91 | match |
| Check item 2 (self-matched susceptibility closure) | py:60-72 (reduced==matrix source + static residue); wl:72-77 | match |
| Check item 3 (numerical fixed points recorded numerically) | not a symbolic identity; no fixed-point number is pinned in script or docs; nothing to verify | n/a (recording directive) |

`paper_alignment: aligned`.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 26 | `simplify(rho_c - rho_c_schur)==0` | rho_c (md:67) | yes |
| A2 | sympy | 27 | `simplify(sigma_c - sigma_c_schur)==0` | sigma_c (md:69-71) | yes |
| A3 | sympy | 44 | `simplify(Ms - Ms_paper)==0` | M_s (tex:16) | partial (true by construction; value anchored via A1) |
| A4 | sympy | 45 | `simplify(Mq - Mq_paper)==0` | M_q (tex:16) | partial (true by construction; value anchored via A2) |
| A5 | sympy | 63 | `simplify(delta_Lambda_matrix - delta_Lambda_reduced)==0` | reduced susceptibility (md:62-63) | yes |
| A6 | sympy | 69 | `simplify(static_limit - (rho_c_schur - sigma_c_schur))==0` | static core residue / Check 2 | yes |
| A7 | sympy | 83 | `simplify(mixed_contribution - Mq_from_schur*Sq)==0` | outlet / Check 1 | yes |
| A8 | sympy | 88 | `simplify(mixed_contribution - (-Mq_from_schur)*Sq) != 0` | non-vacuity guard for A7 | yes |
| A9 | sympy | 96 | `simplify(sigma_c - expr_rc)==0` | sigma_c r_c-form notational equiv | partial (notation restatement) |
| B1 | mathematica | 48 | `expectZero[rho_c - rhoCSchur]` | rho_c | yes |
| B2 | mathematica | 49 | `expectZero[sigma_c - sigmaCSchur]` | sigma_c | yes |
| B3 | mathematica | 62 | `expectZero[mS - mSPaper]` | M_s | partial (anchored via B1) |
| B4 | mathematica | 63 | `expectZero[mQ - mQPaper]` | M_q | partial (anchored via B2) |
| B5 | mathematica | 74 | `expectZero[deltaLambdaMatrix - deltaLambdaReduced]` | reduced susceptibility | yes |
| B6 | mathematica | 77 | `expectZero[staticLimit - (rhoCSchur - sigmaCSchur)]` | static residue / Check 2 | yes |
| B7 | mathematica | 85 | `expectZero[mixedContribution - mQFromSchur*sqVar]` | outlet / Check 1 | yes |
| B8 | mathematica | 88-91 | `If[vacuityResidual===0, fail, pass]` | non-vacuity guard for B7 | yes |
| B9 | mathematica | 95 | `expectZero[sigmaC - sigmaCRc]` | sigma_c r_c-form | partial (notation) |

A3/A4 (and B3/B4) are true by construction in isolation (`Ms` is literally `L*rho_c/Theta` and `Ms_paper` is `L*gs^2/(Ks*Theta)` with `rho_c=gs^2/Ks`), but their *content value* is independently anchored by A1/A2 (Schur reconstruction). They confirm the `L/Theta` propagation and add no false confidence; not flagged. A9/B9 are notational-equivalence restatements of `\sigma_c`; harmless, not load-bearing.

## Findings

None.

## Independent-derivation check (Mathematica)

The `.wl` is NOT a transliteration. Both engines pursue the same physics-dictated strategy (Schur complement of `M_core`, `D\to\infty` limit for `\rho_c`, static `D=1` residue for `\sigma_c`), which is the natural derivation, not copied algebra. Each engine computes the matrix inverse with its own linear-algebra primitive (`M_core.inv()` vs `Inverse[mCore]`) rather than echoing hand-typed intermediate steps, and the static specialization deliberately uses different primitives: SymPy `sp.limit(delta_Lambda_matrix, z_var, 0)` (py:68) vs Mathematica `Normal[Series[deltaLambdaMatrix,{zVar,0,0}]]` (wl:76, explicitly annotated "independent route from SymPy's Limit"). Corresponding sections — py:18-25 vs wl:43-47 (matrix build + residues), py:60-65 vs wl:72-74 (reduced-envelope match), py:79-90 vs wl:82-91 (outlet + non-vacuity) — share structure because the claim is the same, but the algebra is engine-native. Acceptable parallel derivation; no `mathematica_transliteration` finding.

## Engine cross-check

Both saved outputs agree. SymPy emits `rho_c = g_s**2/K_s`, `sigma_c = (K_s*g_q - g_s*lam)**2/(K_s*(K_q*K_s + lam**2))`, `M_s = L*g_s**2/(K_s*Theta)`, `M_q = -L*(K_s*g_q - g_s*lam)**2/(K_s*Theta*(K_q*K_s + lam**2))`; Mathematica emits the identical forms (`rho_c = gS^2/kS`, `sigma_c = (gQ*kS - gS*lam)^2/(kS*(kQ*kS + lam^2))`, etc.) and every `expectZero` residual prints `0` with `PASS`. The benign `Limit::alimv` warning (wl output line 6) merely notes the `dSch>0` assumption is ignored during `Limit[...,dSch->Infinity]`; the resulting limit `gS^2/kS` is correct regardless and the PASS confirms it. `engines_agree: true`.

## Verdict justification

`clean`. I read the card, the notes, and the appendix input line first, then attacked the scripts. I hand-derived `v^T M_core^{-1} v`: the `D\to\infty` limit gives exactly `g_s^2/K_s=\rho_c` and the `\rho_c - \delta\Lambda(D=1)` residue collapses (the `(K_s g_q)^2 - 2K_s g_q\lambda g_s + (\lambda g_s)^2` numerator) to `(K_s g_q-\lambda g_s)^2/[K_s(K_sK_q+\lambda^2)]=\sigma_c` — so the F1 assertions are genuinely non-tautological (`\rho_c,\sigma_c` are outputs of, not inputs to, `M_core`). I hand-verified F2: `\delta\Lambda_{schur}(D)=\rho_c-(K_sg_q-\lambda g_s)^2/[K_s(K_sK_qD+\lambda^2)]`, and substituting `D_W_bare(z)` factors the denominator as `(K_sK_q+\lambda^2)(1-\kappa_c z^2-i\gamma_c z^5)` precisely because `K_sK_q/(K_sK_q+\lambda^2)=1/(1+r_c)`, reproducing the reduced envelope exactly and pinning the `\kappa_c/\gamma_c/r_c` maps. The outlet check (F3) isolates `M_q S_q` and re-anchors its sign/factor through the matrix route, with a working non-vacuity guard (flipped-sign residual `=2M_q^{schur}S_q\neq0`). Symbol domains are sound (sign of `\lambda` is irrelevant — it appears only as `\lambda^2` and `(K_sg_q-\lambda g_s)^2`); no derivatives, no unbounded integrals, no poles in the checked identities. The only mildly redundant checks (A3/A4 paper-card anchors, A9 notation form) carry no false confidence because the underlying values are independently pinned by F1. Paper, notes, and both scripts agree exactly on `(M_s, M_q, \rho_c, \sigma_c)`.

## Value Reconciliation (pass-2 augmentation)

All deliverable values the scripts emit are symbolic (no numeric constants/benchmarks are pinned). Authoritative source = script source + the committed `.txt` outputs (both fresh; no run performed).

| value | source (py/wl + output line) | .tex/.md location | status |
|---|---|---|---|
| `M_s = L g_s^2/(K_s\Theta_\sigma)` | py:35 / wl:56 / sympy.txt:4, math.txt:13 | tex:16 (`M_s=Lg_s^2/(K_s\Theta_\sigma)`); md:79, md:129 | MATCH |
| `M_q = -L(K_sg_q-\lambda g_s)^2/[K_s(K_sK_q+\lambda^2)\Theta_\sigma]` | py:36 / wl:57 / sympy.txt:5, math.txt:14 | tex:16 (`M_q=-L(K_sg_q-\lambda g_s)^2/[K_s(K_sK_q+\lambda^2)\Theta_\sigma]`); md:87, md:132 | MATCH |
| `\rho_c = g_s^2/K_s` | py:33 / wl:54 / sympy.txt:2, math.txt:11 | md:67 (`\rho_c=\frac{g_s^2}{K_s}`); tex omits (intermediate) | MATCH |
| `\sigma_c = (K_sg_q-\lambda g_s)^2/[K_s(K_sK_q+\lambda^2)]` | py:34 / wl:55 / sympy.txt:3, math.txt:12 | md:69-71; tex omits (intermediate) | MATCH |
| `\sigma_c` (r_c-notation form) | py:95 / wl:94 / sympy.txt:12, math.txt:26 | identical to `\sigma_c` above; md:69-71 (`r_c=\lambda^2/(K_sK_q)` notation md:62-63) | MATCH (notational equivalent) |

Internal scaffolding (accounted for, no finding expected in prose): `\delta\Lambda_{schur}`, `\delta\Lambda_{matrix}`, `\delta\Lambda_{reduced}`, `D_W_bare`, `r_c`, `\kappa_c`, `\gamma_c`, `\kappa_0`, `\gamma_0`, `D_{sch}`, `\rho_c^{schur}`, `\sigma_c^{schur}`, `M_q^{from\_schur}`, `Pi_map`, `mixed_contribution`, the non-vacuity guard residual, and all `PASS`/`= 0` residual flags. (`\kappa_c`/`\gamma_c`/`r_c` do appear in notes md:62-63 as the reduced-envelope coefficients, but as cross-check machinery, not stage deliverables.)

reconciliation: complete; 5 deliverable values checked, 0 misaligned

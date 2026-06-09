---
unit_id: 201
batch: VI.1
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-06-09T00:00:00-06:00
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
  notes_stage_files: [moving_throat_pde_stage201_explicit_realization_compiler_and_canonical_orbit_projection.md]
  paper_appendix: present
---

# Audit unit 201 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_201.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage201_explicit_realization_compiler_and_canonical_orbit_projection.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part06.tex` (rows/sections: line 33 ledger row; §sec:app-part06-realization-compiler-and-target-graph lines 239–369; Prop. canonical-fixed-point lines 371–382)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage201_explicit_realization_compiler_and_canonical_orbit_projection_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage201_explicit_realization_compiler_and_canonical_orbit_projection_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage201_explicit_realization_compiler_and_canonical_orbit_projection_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage201_explicit_realization_compiler_and_canonical_orbit_projection_mathematica_audit.txt`

## What the paper claims

Stage 201 turns the Stage-200 reference-free four-scalar home-stretch theorem into an executable realization compiler. `\stagefield{Output}`: "Intrinsic realization packet, canonical orbit projection \(\Pi_{\mathcal O_*}^{\rm can}\), fixed-point criterion \(\mathbf x\in\mathcal Z_*\iff \chi_Q=1\) and \(\Pi_{\mathcal O_*}^{\rm can}(\mathbf x)=\mathbf x\)." The notes enumerate five deliverables: (1) the intrinsic four-scalar realization packet \(\Delta_{\rm real}^{\rm int}=(\chi_Q-1,q_{\rm tr},q_{\rm nt},q_\eta)\) with the mismatch chart \(m_T=\mathfrak R_{\rm tr}^{1/(1+\chi_{0,*})},\ m_K=\mathfrak R_\eta^{-1},\ m_\mu=\mathfrak R_{\rm nt}\mathfrak R_\eta^{-1}\mathfrak R_{\rm tr}^{F_*/(1+\chi_{0,*})}\); (2) the canonical dependent-triple repair vector \(\Delta\mathbf x_{\rm rep}=-S_{(T,K_\eta,\mu)}\mathbf q\) supported only on \((T_U,K_\eta^{\rm eff},\mu_W)\) with \(M_*\Delta\mathbf x_{\rm rep}=-\mathbf q\); (3) the canonical orbit projection \(\Pi^{\rm can}_{\mathcal O_*}(\mathbf x)\); (4) the same-free-quintuple uniqueness theorem (the projected point is the unique target-orbit point with the five free coords fixed); (5) the fixed-point/closure criterion. The notes also state the section property \(M_*S_{(T,K_\eta,\mu)}=I_3\) and a first-order linearized compiler \(\Delta\mathbf x_{\rm rep}^{\rm lin}=-S M_*\Delta\mathbf x\). The stage carries the `\StatusExactClosure{}` tag and is explicitly NOT a claim that the actual PDE branch is graph-aligned.

## What the script claims to verify

Both scripts construct the Stage-192 quotient map `M_*` (3×8) from the named drift basis, build the canonical right-section `S` on the dependent triple `(T,K_eta,mu)`, and verify, in matching sections: I/M1 the section identity `M_* S = I_3` and `M_* dx_rep + q = 0`; II/M2 the mismatch chart `m_T,m_K,m_mu`; III/M3 the closed-form repair vector and its dependent-only support; IV/M6 the uniqueness solve `M_* dx_dep = -q` for `(dT,dKeta,dmu)`; V/M5 the canonical projection `Pi^can` reproduces `dx_rep` in log-drift and is a fixed point at unit ratios; VI/M7 the intrinsic packet equals the pairwise-witness packet under the witness substitution; VII/M8 the linearized compiler `M_*(Dx + dx_rep_lin) = 0`. The SymPy script additionally asserts `exp(q_*) - R_* = 0` (II); the Mathematica script additionally asserts `Det[dependent block] = 1+chi0_star` and an eta-only-perturbation nonzero-`Keta`-repair non-triviality check (M4/M6).

## Paper ↔ script cross-check

| Paper/notes deliverable | Script-side check | Status |
|---|---|---|
| (1) intrinsic packet `Delta_real^int = (chiQ-1, ln R_tr, ln R_nt, ln R_eta)` | py II `Delta_real_int` construct + `exp(q)-R` checks; wl/py mismatch chart `m_T,m_K,m_mu` (II/M2) | match |
| (2) repair vector `dx_rep=-S q`, support only on (T,Keta,mu), `M_* dx_rep=-q` | py III + I (`M_* S - I_3`, `M_* dx_rep + q`); wl M1/M3 | match |
| (3) canonical projection `Pi^can` | py V `x_proj`, `log(x_proj/x)-dx_rep`; wl M5 | match |
| (4) same-free-quintuple uniqueness | py IV `solve(...)` unique triple; wl M6 `Solve` + `Det != 0` + `Length==1` | match |
| (5) fixed-point/closure criterion `Pi^can(x)=x at R=1` | py VII `x_proj.subs(R=1)-x`; wl M5 fixed-point reduction | match |
| section property `M_* S = I_3` | py I `expect_zero("M_* S - I_3")`; wl M1 | match |
| linearized compiler `dx_rep_lin=-S M_* Dx`, `M_*(Dx+dx_rep_lin)=0` | py VII; wl M8 | match |

`paper_alignment: aligned`. Every notes/card deliverable has a corresponding non-tautological script check in both engines.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored? |
|---|---|---|---|---|---|
| A1 | sympy | 85 | `expect_zero(M_* S - I_3)` | section property | yes |
| A2 | sympy | 113–115 | `exp(q_*) - R_* == 0` | (1) chart def | yes (def consistency) |
| A3 | sympy | 116–121 | `m_T,m_K,m_mu` residuals | (1) mismatch chart | yes |
| A4 | sympy | 147 | `dx_rep - expected == 0` | (2) repair closed form | yes |
| A5 | sympy | 148 | `M_* dx_rep + q == 0` | (2) cancellation | yes |
| A6 | sympy | 162–186 | `solve` unique triple residuals | (4) uniqueness | yes |
| A7 | sympy | 219 | `log(x_proj/x) - dx_rep == 0` | (3) projection | yes |
| A8 | sympy | 236–239 | pairwise-witness - intrinsic | (1) packet identity | yes |
| A9 | sympy | 246 | `x_proj(R=1) - x == 0` | (5) fixed point | yes |
| A10 | sympy | 261 | `M_*(Dx + dx_rep_lin) == 0` | linearized compiler | yes |
| B1 | math | 150 | `M_* . S - I_3` | section property | yes |
| B2 | math | 151 | `M_* . dx_rep + q` | (2) cancellation | yes |
| B3 | math | 159–167 | `m_T,m_K,m_mu` residuals | (1) chart | yes |
| B4 | math | 182–187 | `dx_rep - closed form`, support rows, dep entries | (2) repair | yes |
| B5 | math | 194–201 | repair at unit ratios = 0; eta-only nonzero | (2)/(5) non-triviality | yes |
| B6 | math | 220–227 | `Log[x_proj/x] - dx_rep`, fixed-point reduction | (3)/(5) | yes |
| B7 | math | 244–259 | `Det = 1+chi0`, `Length==1`, dep-solve residuals | (4) uniqueness | yes |
| B8 | math | 277–280 | pairwise-witness - intrinsic | (1) | yes |
| B9 | math | 292–296 | `M_*(Dx + dx_rep_lin)`, `qLin - M_* Dx` | linearized compiler | yes |

No tautological rows: every `expect_zero`/`expectZero` subtracts an independently-built object from a separately-defined closed form (e.g. A4 compares `Sdep*qvec` against a hand-typed `expected_dx_rep`; A5 contracts with `M_*`; B7 derives the triple by `Solve` and compares to closed forms). The `exp(q)-R` checks (A2) confirm the chart definitions are mutually consistent (`q := log R`), which is a definitional-consistency check rather than a deep theorem, but it is correctly anchored to the chart in the notes.

## Findings

### F1 — stale_output

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage201_explicit_realization_compiler_and_canonical_orbit_projection_sympy_audit.txt:11` and `:358` (mtime 2026-05-11T12:49:03; script mtime 2026-06-03T15:59:11)

**What's wrong:**
The committed SymPy output `.txt` predates the current SymPy script by ~3 weeks (output 2026-05-11; script 2026-06-03). The captured transcript also carries a STALE STAGE LABEL: the banner reads `STAGE 184 — EXPLICIT REALIZATION COMPILER...` (line 11) and `STAGE 184 LEDGER` (line 358), even though the current script source prints `STAGE 201` (`.py:41`, `.py:263`). The "184" is a +17 pre-renumber artifact in the *captured output only* — the live script string is correct. Re-running the current script will regenerate the transcript with the `STAGE 201` banner and so the divergence is purely a stale capture, not a current-script defect. The Mathematica output (`...mathematica_audit.txt`) is fresh and correctly labelled `STAGE 201`.

**Why this matters:**
The committed transcript does not reflect the current script and displays a wrong stage number, which is misleading to a reader auditing the saved record. It is informational: the orchestrator's independent re-run refreshes it.

**Required change:**
None for Codex to hand-edit. Re-run `python3 scripts/moving_throat_pde_stage201_explicit_realization_compiler_and_canonical_orbit_projection_sympy_audit.py` and recapture the output `.txt`; the banner will then read `STAGE 201`. (This is the standard orchestrator-refresh path; no source change is needed because the source already prints 201.)

**Verification:**
After re-run, line 11 of the SymPy output reads `STAGE 201 — EXPLICIT REALIZATION COMPILER AND CANONICAL ORBIT PROJECTION` and the `# Date:` header is newer than the script mtime.

## Independent-derivation check (Mathematica)

**Verdict: INDEPENDENT.** The load-bearing object is the canonical right-section `S` and the repair vector `dx_rep = -S q` (the projection, fixed-point criterion, and linearized compiler all derive from `S`). The two engines EXTRACT `S` by genuinely different methods:

- SymPy (`.py:66–74`) — **invert-the-dependent-block-and-scatter**: it slices out the three dependent columns of `M_*`, inverts that 3×3 block, and embeds the inverse into the 8×3 scatter matrix:
  ```
  Pdep = sp.Matrix.hstack(Mstar[:, T_idx], Mstar[:, Keta_idx], Mstar[:, mu_idx])
  Pdep_inv = sp.simplify(Pdep.inv())
  Sdep = sp.simplify(Edep * Pdep_inv)
  ...
  dx_rep = sp.simplify(-Sdep * qvec)
  ```
- Mathematica (`.wl:120–144`) — **solve-the-augmented-constrained-square-system**: it stacks `M_*` (3×8) on top of five unit-row free-coordinate constraints to form an 8×8 system and `LinearSolve`s it column-by-column for `S`, and once for `-q` to get `dxRep` directly:
  ```
  freeSelector = UnitVector[8, #] & /@ freeSlots;
  constrainedSystem = Join[Mstar, freeSelector];
  Ssection = ... LinearSolve[constrainedSystem, Join[UnitVector[3, col], ConstantArray[0,5]]] ...
  dxRep = ... LinearSolve[constrainedSystem, Join[-q, ConstantArray[0,5]]] ...
  ```

These are not the same operation re-typed: the `.py` route never forms the augmented square system, and the `.wl` route never explicitly inverts the dependent 3×3 block or uses a scatter matrix. The agreement of both on `S`, `dx_rep`, and `M_* S = I_3` is therefore a real cross-check. The `.wl` further adds an independent diagnostic the `.py` lacks — `Det[dependent block] = 1+chi0_star` (`.wl:244`) — confirming invertibility from a different angle. The uniqueness section is `Solve`-based in both (py IV / wl M6), i.e. same method there, but uniqueness is a secondary deliverable and the primary load-bearing `S`/`dx_rep` route differs. Not a transliteration.

## Engine cross-check

Both engines produce the identical closed forms and all residuals vanish:

- `S` row 7 (mu): py `[F_*/(1+chi0_*), 1, -1]` (output lines 58–60) = wl `{Fstar/(1+chi0s), 1, -1}` (output line 12). Row 8 (T): py `[1/(1+chi0_*),0,0]` = wl `{(1+chi0s)^(-1),0,0}`. Row 5 (Keta): py `[0,0,-1]` = wl `{0,0,-1}`. Match.
- `dx_rep`: py `(0,0,0,0, log R_eta, 0, -q_nt+q_eta-F_* q_tr/(1+chi0_*), -q_tr/(1+chi0_*))` = wl `{0,0,0,0,Log[Reta],0, Log[Reta]-Log[Rnt]-(Fstar Log[Rtr])/(1+chi0s), -(Log[Rtr]/(1+chi0s))}`. Match.
- Uniqueness triple: py `{dKeta: log R_eta, dT: -log R_tr/(1+chi0_*), dmu: ...}` = wl `{dT -> -(Log[Rtr]/(1+chi0s)), dKeta -> Log[Reta], dmu -> Log[Reta]-Log[Rnt]-(Fstar Log[Rtr])/(1+chi0s)}`. Match.
- All `expect_zero`/`expectZero` residuals = 0; both exit 0. No sign/factor disagreement. `engines_agree: true`.

## Verdict justification

`findings` (one low-severity informational `stale_output`). Attacks tried and failed: (a) sought a transliterated `.wl` — found a genuinely independent construction route for the load-bearing section `S` (augmented `LinearSolve` vs dependent-block inversion); (b) probed for tautology — every residual contracts an independently-built object against a hand-typed closed form or against `M_*`, none is guaranteed by construction; (c) checked the `M_*` matrix entries against the notes/appendix repair-vector form (notes eq. lines 283–293, appendix eq. lines 336–345) — they match the script's `Mstar` rows and the resulting `dx_rep`; (d) checked symbol domains — all positive/real assumptions are physically justified (microscopic state `> 0`, ratios `> 0`); (e) checked the fixed-point and same-free-quintuple support claims — `dx_rep` free-coordinate rows (1–4,6) are exactly zero in both engines, and `Det = 1+chi0_* > 0` guarantees the unique solve. I confirm I read the card, the notes (10 sections), and the Part VI appendix Stage-201 rows, and the script's verified claims match the paper's stated deliverables. The only blemish is the stale, mislabeled ("STAGE 184") SymPy output capture, which the source already prints correctly as 201 and which a re-run fixes.

## Value Reconciliation (pass-2 augmentation)

The scripts emit SYMBOLIC closed-form deliverables (no postulated numeric constants — `chi0_star, deltaU_star, E_star, F_star` are free symbolic parameters, and `R_tr,R_nt,R_eta,chi_Q` are free ratios/scalars). I reconcile each symbolic deliverable against the card/notes.

| value (symbolic deliverable) | source (py / wl + output line) | .tex/.md location | status |
|---|---|---|---|
| `M_*` quotient matrix (3×8) | py:54–60 / wl:103–116; out wl:6 | notes lines 244–293 (repair-vector form), appendix M_* implicit in eqs 336–345; not the card | MATCH (notes/appendix) |
| section `S_(T,Keta,mu)` | py:74; out py:45–64 / wl:128; out wl:12 | notes line 263 `S_(T,K_eta,mu)`; appendix line 331 dependent triple | MATCH |
| repair vector `dx_rep` | py:129; out py:119–142 / wl:139; out wl:14 | notes eq lines 283–293; appendix eq lines 336–345 | MATCH |
| `m_T = R_tr^(1/(1+chi0_*))` | py:101 / wl:155; out wl:25 | notes lines 205–208; card body (Output region) | MATCH (notes) |
| `m_K = R_eta^(-1)` | py:102 / wl:156; out wl:27 | notes lines 211–214 | MATCH (notes) |
| `m_mu = R_nt R_eta^(-1) R_tr^(F_*/(1+chi0_*))` | py:103 / wl:157; out wl:29 | notes lines 216–220 | MATCH (notes) |
| `Delta_real^int = (chiQ-1, q_tr, q_nt, q_eta)` | py:98 / wl:269; out py:75–82 | notes lines 156–167; appendix eq lines 308–318; card Output | MATCH |
| `Pi^can(x)` projection vector | py:199–210 / wl:206–215; out wl:58 | notes lines 339–376; appendix eq lines 357–366; card Output | MATCH |
| uniqueness triple `(dT,dKeta,dmu)` | py:162–170 / wl:234–239; out wl:71 | notes lines 414–420 (dependent triple determined uniquely) | MATCH |
| `Det[dependent block] = 1+chi0_*` | wl:244; out wl:72 | notes line 420 "exact pivot block on (T,K_eta,mu)" — qualitatively present, not a boxed number | MATCH (internal diagnostic; value not a stated deliverable) |
| `dx_rep_lin = -S M_* Dx`, `M_*(Dx+dx_rep_lin)=0` | py:254/261 / wl:286/292 | notes lines 502–516 (linearized compiler) | MATCH |

INTERNAL scaffolding (no finding): `Pdep`, `Pdep_inv`, `Edep` (py intermediate construction objects); `freeSelector`, `constrainedSystem`, `depSlots`, `freeSlots` (wl construction objects); `q^lin` intermediate; all `expect_zero/expectZero` residuals; PASS/FAIL flags; the `eta-only perturbation` non-triviality probe value.

reconciliation: complete; 11 deliverable values checked, 0 misaligned.

(The lone `Det = 1+chi0_*` literal-vs-stated point: the notes assert the pivot block is invertible/unique but do not box the number `1+chi0_*` as a stage deliverable, and the script uses it only as an invertibility diagnostic, not as a carried-forward constant — so it is INTERNAL, not a MISMATCH/MISSING-DELIVERABLE.)

## Self-test notes

Variable-independence: the only `D[...]`/`diff`-style operations are matrix contractions `M_* . Dx` (wl M8 / py VII) over the full 8-vector `Dx`, whose every component appears in `M_*` rows — no identically-zero derivative trap. Symmetry/parity: no unbounded integrals in this unit. Trivial-case pre-check: at `R_tr=R_nt=R_eta=1` all `q=0`, so `dx_rep=0` and `Pi^can(x)=x` — both engines assert exactly this (py:246, wl:225) and it reduces to 0 as required; the eta-only probe `Reta=r!=1` correctly gives nonzero `Keta` repair `Exp[Log r] = r != 1`. Path spec: no missing-script directive needed (both engines present). The single finding is a non-script-side informational `stale_output`, so no Codex patch directive is written.

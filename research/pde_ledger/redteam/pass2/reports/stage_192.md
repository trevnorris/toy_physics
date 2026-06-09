---
unit_id: 192
batch: V.3
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-06-09T00:00:00Z
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
  notes_stage_files: [/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage192_orbit_quotient_projectors.md]
  paper_appendix: present
---

# Audit unit 192 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_192.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage192_orbit_quotient_projectors.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part05.tex` (rows: line 115 status row; lines 1139–1210 block narrative)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage192_orbit_quotient_projectors_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage192_orbit_quotient_projectors_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage192_orbit_quotient_projectors_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage192_orbit_quotient_projectors_mathematica_audit.txt`

## What the paper claims

`\stagefield{Output}` (stage_192.tex:15): "Builds exact projectors \(Q_{\rm quot}\), \(O_{\rm orb}\), showing quotient failure lives only on \((T_U,K_\eta^{\rm eff},\mu_W)\)." The notes (sections 1–9) make this exact: working in the ordered 8-vector drift `Δx`, the quotient packet is `q = M_* Δx` with the explicit 3×8 monomial-drift matrix `M_*`; selecting the dependent triple `(T, K_eta, μ)` gives the 3×3 pivot block `P` with `det P = 1+χ_{0,*} > 0`; the canonical section `S = E P^{-1}` is the right inverse `M_* S = I_3`; the complementary projectors `Q_quot = S M_*`, `O_orb = I_8 − Q_quot` are idempotent and complementary (`Q²=Q`, `O²=O`, `QO=OQ=0`, `M_* O = 0`, `M_* Q = M_*`); the failure piece `Q Δx` is supported only on the dependent triple with the explicit components `(Δ_T)_fail = q_tr/(1+χ)`, `(Δ_Keta)_fail = −q_eta`, `(Δ_μ)_fail = F_* q_tr/(1+χ)+q_nt−q_eta`; the orbit piece preserves the five free coordinates and reproduces the Stage-187 single-orbit law with `α_* = (1+δ_{U,*})/(1+χ_{0,*})`; and the orbit-lock equivalence `q=0 ⇔ M_* Δx=0 ⇔ Q Δx=0 ⇔ Δx ∈ ker M_*`. The appendix (lines 1147–1208) carries `M_*`, `P`, `det P`, `Q_quot`, `O_orb`, the projector identities, the failure components, and the lock equivalence verbatim. Eight distinct deliverables.

## What the script claims to verify

The SymPy script (Sections I–VI) encodes `M_*` from the basis, forms the pivot block `Pdep` by column selection, computes its inverse, builds `S = E·P⁻¹`, then `Q = S·M`, `O = I−Q`, and asserts (via `expect_zero`): pivot block and inverse equal the notes literals; `det P − (1+χ) = 0`; `M_* S − I_3 = 0`; the four projector identities plus `M_*O`, `M_*Q−M_*`; that rows 0,1,2,3,5 of `Q` vanish (free-coordinate support); `Q Δx` equals the dependent-triple closed form; `M_* Δx_fail − q = 0`; `O Δx` equals the orbit-law closed form with `α`; `Δx − (orbit+fail) = 0`; and the canonical representative identities `M_*(Sq)−q`, `Q(Sq)−Sq`, `O(Sq)`. The Mathematica script verifies the SAME eight deliverables but builds `S` by a **constrained `LinearSolve`** on the augmented 8×8 `Join[quotientRows, freeSelector]`, derives the failure and orbit pieces by independent constrained/kernel `LinearSolve`s, and confirms the lock equivalence by two independent `Solve`s (one on `M·Δx==0`, one on `Q·Δx==0`) whose dependent-triple laws are shown to agree.

## Paper ↔ script cross-check

| Paper deliverable | Script-side check | Status |
|---|---|---|
| `M_*` monomial-drift matrix (notes §1, app 1152) | py:43–49 / wl:69–73 encode identical entries | match |
| Pivot block `P`, `det P = 1+χ` (notes §2, app 1156–1166) | py:59–72 / wl:82–98 | match |
| Section `S`, `M_* S = I_3` (notes §3, app 1167) | py:88–110 / wl:105–130 | match |
| Projectors `Q,O` + identities `Q²,O²,QO,OQ,M O,M Q` (notes §4, app 1175–1184) | py:113–124 / wl:134–143 | match |
| Failure support only on dependent triple + closed forms (notes §5, app 1186–1197) | py:126–157 / wl:147–170 | match |
| Orbit piece preserves free coords + Stage-187 law, `α_*` (notes §6) | py:159–179 / wl:174–199 | match |
| Orbit-lock equivalence `q=0 ⇔ M Δx=0 ⇔ Q Δx=0` (notes §7, app 1199–1208) | py:181–194 / wl:201–247 | match |
| Canonical representative `S q` closed form (notes §3/§5) | py:182–190 / wl:203–220 | match |

All eight deliverables are exercised by both engines. `paper_alignment: aligned`.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored? |
|---|---|---|---|---|---|
| A1 | sympy | 69 | `expect_zero(Pdep − Pdep_expected)` | pivot block | partial (vs literal) |
| A2 | sympy | 72 | `det(P) − (1+χ) == 0` | `det P = 1+χ` | yes |
| A3 | sympy | 84–86 | `P_inv − exp`, `P_inv P − I`, `P P_inv − I` | inverse | yes (84–86: 85/86 non-taut) |
| A4 | sympy | 109–110 | `S − exp`, `M_* S − I_3` | section / right-inverse | yes (110 non-taut) |
| A5 | sympy | 119–124 | `Q²−Q`, `O²−O`, `QO`, `OQ`, `M O`, `M Q−M` | projector identities | yes |
| A6 | sympy | 128–129 | `Q rows 0,1,2,3,5 == 0` | failure support | yes |
| A7 | sympy | 156–157 | `Q Δx − exp`, `M_* Δx_fail − q` | failure closed forms | yes (157 non-taut) |
| A8 | sympy | 177–179 | `O Δx − exp`, `M_* Δx_orbit`, `Δx−(orbit+fail)` | orbit law | yes (178/179 non-taut) |
| A9 | sympy | 188–194 | `M_*(Sq)−q`, `Q(Sq)−Sq`, `O(Sq)`, `M_* Dx_orbit_exp` | lock / representative | yes |
| M1 | math | 96–101 | pivot/inverse vs SymPy target + `P⁻¹P−I`, `PP⁻¹−I` | pivot block/inverse | yes (100/101 non-taut) |
| M2 | math | 128–130 | `S(LinearSolve) − SymPy target`, `M S − I_3` | section (independent route) | yes |
| M3 | math | 137–143 | projector identities + free-row support | projectors | yes |
| M4 | math | 168–170 | `Q Δx − constrained solve`, `− SymPy form`, `M Δx_fail − q` | failure (independent route) | yes |
| M5 | math | 196–199 | `O Δx − kernel solve`, `− SymPy form`, `M Δx_orbit`, decomposition | orbit (independent route) | yes |
| M6 | math | 217–247 | representative + two `Solve` lock laws agree | lock equivalence (independent route) | yes |

Every script-side check traces to a specific paper deliverable; no orphaned assertions. Note the `*Expected`/`*target` comparisons (A1, A4-left, M1/M2-left, etc.) are not self-tautological because the same quantity is also produced by an independent construction route and the two are cross-checked (A4-right `M_* S − I_3`, M2 `M S − I_3` and `Q Δx − constrained failure solve`, etc.).

## Findings

### F1 — stale_output

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage192_orbit_quotient_projectors_sympy_audit.txt` (mtime 2026-05-11 12:48)
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage192_orbit_quotient_projectors_sympy_audit.py` (mtime 2026-06-03 15:59)

**What's wrong:**
The saved SymPy output predates the current script by ~3 weeks. The captured transcript still carries the pre-renumber banner — line 11 reads `STAGE 175 — EXACT ORBIT/QUOTIENT PROJECTORS…` and line 634 reads `STAGE 175 LEDGER`, whereas the current `.py` prints `STAGE 192 …` (py:35) and `STAGE 192 LEDGER` (py:196). So the committed output does not reflect the present state of the script (the math content is otherwise identical; only the stage label changed). The Mathematica output (mtime 2026-06-01 11:11) is fresh relative to its `.wl` (mtime 2026-06-01 11:09) and prints the correct `STAGE 192` banner.

**Why this matters:**
Informational. The reconciliation and verdict are based on the script source plus the (fresh) Mathematica output; the SymPy math is unchanged, so no result is wrong. The stale banner is a renumber artifact, not a math defect.

**Required change:**
Re-run the SymPy script to refresh the committed output so its banner reads `STAGE 192`. The orchestrator's independent re-run already does this; no script edit is required.

**Verification:**
After re-run, line 11 of the SymPy output `.txt` reads `STAGE 192 — EXACT ORBIT/QUOTIENT PROJECTORS …` and the ledger header reads `STAGE 192 LEDGER`, matching `py:35` and `py:196`.

## Independent-derivation check (Mathematica)

GENUINELY INDEPENDENT. The two engines build the section `S` — the load-bearing object from which `Q`, `O`, the failure piece, and the orbit piece all flow — by structurally different routes:

- SymPy (py:89–94): explicit embedding times matrix inverse —
  `Edep = zeros(8,3); Edep[T_idx,0]=…; Sdep = simplify(Edep * Pdep_inv)`, where `Pdep_inv = Pdep.inv()` (py:74). i.e. `S = E · P⁻¹`.
- Mathematica (wl:105–113): a constrained linear solve over the augmented 8×8 system —
  `constrainedSystem = Join[quotientRows, freeSelector]` (wl:78), then
  `sectionFromSolve = Transpose[Table[LinearSolve[constrainedSystem, Join[UnitVector[3,col], ConstantArray[0,5]]], {col,3}]]`. Each section column is solved as "the unique 8-vector with zero free coordinates whose quotient packet is `e_col`." No matrix inverse of `P` is taken to build `S`; `Inverse[pivotBlock]` (wl:88) appears only as a separate cross-check of the pivot inverse, not in the `S` construction path.

The downstream objects are likewise re-derived independently and then cross-checked against the projector route:
- failure piece: `failureByConstrainedSolve = LinearSolve[constrainedSystem, Join[packet, …]]` (wl:149–154), checked equal to `Q·Δx` (wl:168);
- orbit piece: `orbitByKernelSolve = LinearSolve[constrainedSystem, Join[{0,0,0}, drift[[freeSlots]]]]` (wl:176–181), checked equal to `O·Δx` (wl:196);
- lock equivalence: two independent `Solve`s, `lockByMap = Solve[M·drift==0, …]` (wl:226) and `lockByProjector = Solve[Q·drift==0, …]` (wl:231), with their dependent-triple solution laws shown to agree (wl:236–239) — the SymPy side instead checks the algebraic identities `M_*(Sq)−q` etc. (py:188).

The `.wl` does compare each result against a `…Expected`/`SymPy target` literal (wl:83, 116, 155, 182, 205). On its own that would be a port-style anchor, but here it is the cross-engine reconciliation tie-point, not the derivation: the independent `LinearSolve`/`Solve` route produces the object and `…Expected` confirms the two engines land on the same closed form. This is exactly the second-engine policy (independent derivation + cross-check), not transliteration. No `mathematica_transliteration` finding.

## Engine cross-check

Both engines agree. The Mathematica output (`…mathematica_audit.txt`) shows every `expectZero` PASS, including the cross-route ties: `Q Delta x - constrained failure solve = {0,…}` (line 53), `O Delta x - constrained kernel solve = {0,…}` (line 64), `Solve laws for Q Delta x == 0 and M Delta x == 0 agree = {0,0,0}` (line 92). The SymPy output (banner-stale but mathematically current) shows all `expect_zero` residuals as zero matrices/scalars. The explicit forms match: SymPy `Canonical quotient representative S q` (txt lines 552–571) = `{0,0,0,0,−q_eta,0,(F q_tr+(1+χ)(−q_eta+q_nt))/(1+χ), q_tr/(1+χ)}` equals the Mathematica `Canonical representative S q = {0,0,0,0,-qE,0,-qE+qN+(ff*qT)/(1+chi), qT/(1+chi)}` (output line 75). Failure components, orbit law, and pivot inverse all coincide. No `engine_disagreement`.

## Verdict justification

`findings`, with the single finding being a low-severity, informational `stale_output` (renumber-era SymPy banner; math unchanged). I attacked the audit on three fronts and it held: (1) Independence — the `.wl` is a genuinely independent re-derivation (constrained `LinearSolve`/`Solve` vs. `E·P⁻¹` embedding-inverse), not a port; the `…Expected` literals are cross-engine tie-points backed by an independent construction. (2) Tautology — the `…Expected`-comparison assertions are each paired with a from-construction check (`M_* S − I_3`, `M·Δx_fail − q`, the agreeing `Solve` laws), so none is self-confirming. (3) Paper alignment — all eight notes/appendix deliverables (M_*, pivot+det, section+right-inverse, projector identities, failure support+closed forms, orbit law+α, lock equivalence, representative) have faithful, exact, two-engine checks; every constant is symbolic and matches the notes. I confirm I read the card, notes, and appendix block before the scripts. The only paper-side note is a benign card-text lag (see Value Reconciliation), not a value/identity misalignment.

## Self-test notes

Variable-independence trap: this stage is pure linear algebra over the symbols `{χ,δU,E,F}` and the formal drift components; there are no `diff`/`D` calls, so the zero-derivative failure mode does not apply. Symmetry/parity trap: no integrals over unbounded domains; N/A. Trivial-case trap: the load-bearing zeros are matrix residuals (`Q²−Q`, `M_* S−I_3`, `Δx−(orbit+fail)`) which are non-trivial polynomial identities in four free symbols, not constructed-to-be-zero; substituting e.g. χ=δU=E=F=1 leaves genuine 8×8/3×3 structure that must cancel, and the saved outputs confirm it does. No directive self-test on new code was needed since the only finding is a re-run (no new assertion).

## Value Reconciliation (pass-2 augmentation)

All emitted deliverables are symbolic closed forms (no pinned numeric figures of merit in this stage). I enumerated every labeled result the scripts emit and located each in the `.tex` card and/or `.md` notes.

| value | source (py/wl + output line) | .tex/.md location | status |
|---|---|---|---|
| `M_*` 3×8 monomial-drift matrix | py:43–49 / wl:69–73 (sympy txt 17–28) | notes §1 lines 99–106; app 1152 | MATCH |
| pivot block `P=[[1+χ,0,0],[−F,−1,1],[0,−1,0]]` | py:60–66 / wl:83–87 (math out 9) | notes §2 lines 132–137; app 1159–1163 | MATCH |
| `det P = 1+χ_{0,*}` | py:72 / wl:97 (math out 12) | notes §2 line 142; app 1165 | MATCH |
| `P⁻¹` closed form | py:75–81 / wl:89–93 | notes §2 lines 151–157 | MATCH |
| section `S` (dependent-triple support, `F/(1+χ),1/(1+χ)` entries) | py:95–106 / wl:116–125 (math out 21) | notes §3 lines 186–195 | MATCH |
| right-inverse `M_* S = I_3` | py:110 / wl:130 | notes §3 line 200; app 1167 | MATCH |
| projectors `Q_quot = S M_*`, `O_orb = I−Q` | py:113–114 / wl:134–135 | notes §4 lines 219–221; app 1170–1172 | MATCH |
| projector identities `Q²=Q,O²=O,QO=OQ=0,M_*O=0,M_*Q=M_*` | py:119–124 / wl:137–142 | notes §4 lines 225–235; app 1177–1183 | MATCH |
| failure `(Δ_T)_fail=q_tr/(1+χ)` | py:148–150 / wl:163–164 | notes §5 line 275; app 1188 | MATCH |
| failure `(Δ_Keta)_fail=−q_eta` | py:148 / wl:160 | notes §5 line 278; app 1190 | MATCH |
| failure `(Δ_μ)_fail=F q_tr/(1+χ)+q_nt−q_eta` | py:150 / wl:162 | notes §5 line 281; app 1196 | MATCH |
| `α_* = (1+δ_{U,*})/(1+χ_{0,*})` | py:160 / wl:174 | notes §6 line 310 | MATCH |
| orbit law `(Δ_Keta)_orb=2Δ_c−Δ_U`, `(Δ_T)_orb`, `(Δ_μ)_orb` | py:166–172 / wl:187–192 (math out 62) | notes §6 lines 314–328 | MATCH |
| canonical representative `S q` closed form | py:185 / wl:204–214 (sympy txt 565–571; math out 75) | notes §3 line 210, §5 §7 | MATCH |
| orbit-lock equivalence `q=0 ⇔ M_*Δx=0 ⇔ QΔx=0 ⇔ Δx∈ker M_*` | py:188–209 / wl:218–247 | notes §7 lines 341–349; app 1199–1208 | MATCH |

INTERNAL scaffolding (no prose expected): `expect_zero`/`expectZero`/`pass`/`fail` helpers, `constrainedSystem`/`freeSelector`/`depSlots`/`freeSlots` index bookkeeping, the `*Expected`/`*target` comparison literals (cross-engine tie-points), `stripCE`/`cleanScalar`/`cleanTensor`/`zeroTensorQ` Mathematica simplifiers, the formal symbols `Δ_λ…Δ_T` and `q_tr,q_nt,q_eta`.

Paper-side card-text lag (paper-cleanup class, NOT a value/identity finding, not counted): `stage_192.tex:11` `\stagefield{Verification}` says "Mathematica audit: none yet", but a passing `.wl` exists (`mathematica/moving_throat_pde_stage192_orbit_quotient_projectors_mathematica_audit.wl`). Per V.3 instructions this is a card-text lag, not a misalignment; flagged for the paper-cleanup tracker, not routed to Codex.

reconciliation: complete; 15 deliverable values checked, 0 misaligned.

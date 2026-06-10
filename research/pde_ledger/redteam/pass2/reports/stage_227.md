---
unit_id: 227
batch: VII.1
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-06-09T00:00:00Z
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
  notes_stage_files: ["/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage227_pure_transfer_load_factor_outgoing_rigidity_sieve_and_first_co_loading_no_go_sympy_audit.md"]
  paper_appendix: present
---

# Audit unit 227 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_227.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage227_pure_transfer_load_factor_outgoing_rigidity_sieve_and_first_co_loading_no_go_sympy_audit.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part07.tex` (rows 66, 747–804)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage227_pure_transfer_load_factor_outgoing_rigidity_sieve_and_first_co_loading_no_go_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage227_pure_transfer_load_factor_outgoing_rigidity_sieve_and_first_co_loading_no_go_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage227_pure_transfer_load_factor_outgoing_rigidity_sieve_and_first_co_loading_no_go_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage227_pure_transfer_load_factor_outgoing_rigidity_sieve_and_first_co_loading_no_go_mathematica_audit.txt`

## What the paper claims

Stage 227 factors the Stage-226 pure-transfer survivor into a one-port outgoing-load factor and proves the first same-charge co-loading no-go. `\stagefield{Output}`: *"Pure-transfer theorem: $\Xi_1=2\delta\ln(P/\Delta)$. If the outgoing transfer shape is rigid while interference and hybridization ratios are also rigid, then nontrivial co-loading is impossible."* The notes enumerate the full deliverable list: (1) the pure-transfer theorem $\Xi_1=N_{01}/N_0=2P_{01}/P-2\Delta_{01}/\Delta$; (2) the one-port factorization $\Lambda=(G_W/\Omega_W^2)(1+I)/(1-H)$ with $I=RG_U/(\Omega_U^2 G_W)$, $H=R^2/(\Omega_U^2\Omega_W^2)$; (3) the slope law $\Xi_1=2[m+\frac{I}{1+I}i+\frac{H}{1-H}h]$ with its sample specialization $\Xi_1=2m+\frac{6}{19}i+\frac{50}{98\pi^2-25}h$; (4) sample-branch invariants $I=3/16$, $H=25/(98\pi^2)$; (5) the combined $i=h=0$ rigidity determinant with factor $(200+147\pi^2)$, proving the co-loading no-go; (6) the three one-dimensional rigid survivors $v_i,v_h,v_m$ with gain scales $\sigma_i,\sigma_h,\sigma_m$; (7) the corridor norm $\sigma_{\rm transfer}$ and transported 10%-loss ceilings. The appendix (eqs. app-part07-xi-lambda, -lambda-factorization, -co-loading-no-go-hypothesis; Theorem thm:app-part07-pure-transfer-sieve) restates (1)–(5) symbolically with `\StatusExactClosure`. Note: the card's `\stagefield{Verification}` says "Mathematica audit: none yet" but a `.wl` now exists (created after the card) — a prose lag, not a math defect (see Self-test notes).

## What the script claims to verify

The SymPy script builds the explicit finite-throat one-port compiler with `eps`-exponential dressing of every primitive, forms first-order log-derivatives, reconstructs the Stage-226 pure-transfer corridor (rank-3 transfer matrix, 2-dim nullspace), and asserts: the pure-transfer theorem (line 171), the $\Lambda$ factorization (line 176), the slope law and sample specialization (lines 187, 190), the sample invariants $I=3/16$, $H=25/(98\pi^2)$ (lines 192–193), the combined $i$=$h$ rigidity determinant against a hard-pinned closed form carrying $(200+147\pi^2)$ and nonzero (lines 216–217), the three unit survivors and their $|\Xi_1|$ gain scales (lines 237–252), and the corridor norm + ceilings (lines 269, 290–291). The Mathematica script independently re-derives the same deliverables via a logarithmic-derivative operator and `LinearSolve`-normalized nullspace basis, and notably computes the $i$=$h$ determinant from scratch (`FullSimplify[Det[reducedIH]]`, line 250) rather than hard-pinning it.

## Paper ↔ script cross-check

| Paper deliverable | Script check | Status |
|---|---|---|
| (1) $\Xi_1=2P_{01}/P-2\Delta_{01}/\Delta$ | py:171, wl:171 (M1) | match |
| (2) $\Lambda=(G_W/\Omega_W^2)(1+I)/(1-H)$ | py:176, wl:182 (M2) | match |
| (3) slope law + sample specialization | py:187,190, wl:213,214 (M3) | match |
| (4) $I=3/16$, $H=25/(98\pi^2)$ | py:192,193, wl:211,212 (M3) | match |
| (5) $i$=$h$=0 rigidity det $(200+147\pi^2)\ne0$ (co-loading no-go) | py:216,217, wl:250,252 (M4) | match |
| (6) survivors $v_i,v_h,v_m$ + $\sigma_i,\sigma_h,\sigma_m$ | py:237–252, wl:272–284 (M5,M6) | match |
| (7) $\sigma_{\rm transfer}$ + 10%-loss ceilings | py:269,290–291, wl:308–316 (M7) | match |

`paper_alignment: aligned`.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored? |
|---|---|---|---|---|---|
| A1 | sympy | 147 | `M_transfer.rank() == 3` | corridor (Stage-226 carry) | yes |
| A2 | sympy | 150,157–160 | nullspace dim 2 + basis vs `expected_t` | corridor | yes |
| A3 | sympy | 171 | `simplify(Xi_transfer - 2(P01/P-Delta01/Delta))==0` | claim 1 | yes |
| A4 | sympy | 176 | `simplify(Lambda - GW/OmW^2 (1+I)/(1-H))==0` | claim 2 | yes |
| A5 | sympy | 187 | `simplify(Xi_transfer - Xi_mih)==0` | claim 3 | yes |
| A6 | sympy | 190 | sample specialization residual `==0` | claim 3 | yes |
| A7 | sympy | 192,193 | `I_s==3/16`, `H_s==25/(98pi^2)` | claim 4 | yes |
| A8 | sympy | 216,217 | `det_ih - expected_det==0`, `det_ih!=0` | claim 5 | yes |
| A9 | sympy | 237–252 | survivor dirs + `|Xi_1|` close | claim 6 | yes |
| A10 | sympy | 269,290–291 | sigma_transfer + ceilings close | claim 7 | yes |
| M1 | math | 171 | `expectZero` load-law residual | claim 1 | yes |
| M2 | math | 182 | `expectZero` Lambda factorization | claim 2 | yes |
| M3a | math | 208–214 | log-forms + sample + slope law `expectZero` | claims 3,4 | yes |
| M4 | math | 250,252,253 | `Det[reducedIH]` symbolic + nonzero + rank 2 | claim 5 | yes |
| M5 | math | 260–274 | reduced nullities 1 + survivor directions | claim 6 | yes |
| M6 | math | 282–284 | `|Xi1(v)|` gain scales | claim 6 | yes |
| M7 | math | 308–316 | sigma_transfer + ceilings | claim 7 | yes |

All rows anchored and non-tautological. The two engines reach the slope vectors and the determinant by genuinely different routes (see Independent-derivation check), so neither A8 nor M4 is self-confirming.

## Findings

None. The audit is clean.

The first-pass `paper_misalignment` (notes-only: rigidity determinant factor was `251+215π²`, USER-RESOLVED to `200+147π²`) is RE-CONFIRMED HOLDING:
- SymPy hard-pins and asserts `(200 + 147 * sp.pi**2)` (py:215–216); its committed output emits the full factor `-19*(-25 + 98*pi**2)*(200 + 147*pi**2)*(441*pi**2 + 4400)/(...)` (sympy.txt:16).
- Mathematica computes the determinant from scratch (`FullSimplify[Det[reducedIH]]`, wl:250) and emits `(-19*(-25 + 98*Pi^2)*(200 + 147*Pi^2)*(4400 + 441*Pi^2))/(...)` (math.txt:61) — an INDEPENDENT cross-engine corroboration, not an echo of the SymPy literal.
- Notes read `(200+147\pi^2)` (md:294) with NO surviving `251+215π²`.
- Repo-wide grep for `251+215`, `215π`, etc. across the 227 notes/.tex/appendix/.py/.wl returns EMPTY. No surviving `251+215π²` anywhere.

## Independent-derivation check (Mathematica)

The `.wl` is GENUINELY INDEPENDENT, not a transliteration. Three corresponding sections show distinct derivation choreography:

1. **First-order variation.** SymPy dresses every primitive as `param * exp(eps*x)` and takes `sp.diff(EXPR_dressed, eps).subs(eps, 0)` (py:65–99). The `.wl` instead defines a logarithmic-derivative operator `deltaMixed[expr] := Total[(x * param * D[expr, param])]` over `logPairs` (wl:110–114) and applies it to the undressed expressions. Same δln semantics, structurally different machinery — the `.wl` never introduces an `eps` dressing variable.

2. **Corridor basis normalization.** SymPy takes `M_transfer.nullspace()` directly (py:149). The `.wl` takes `NullSpace[transferMatrix]`, transposes, isolates the tail rows `{4,5}`, and re-solves `LinearSolve[tailSolveMatrix, {1,0}]` / `{0,1}` to pin the tail to identity (wl:228–235). Different basis-fixing route producing the same two vectors.

3. **The load-bearing determinant.** SymPy hard-pins the closed form and asserts equality (py:215–216). The `.wl` does NOT carry the literal — it computes `Det[reducedIH]` and `FullSimplify`s it, asserting only structural nonzero (wl:250,252). That the `.wl` independently arrives at the same `(200+147π²)` factor in its emitted output is the strongest possible cross-engine confirmation of this stage's central no-go constant.

No `mathematica_transliteration` finding.

## Engine cross-check

Both engines agree to the precision they claim:
- Pure-transfer slope $\Xi_1$: identical rational-in-$\pi^2$ coefficient vector (sympy.txt:8 vs math.txt:10, same numerators over `19*(-25+98pi^2)`).
- Determinant: same factored form (sympy.txt:16 vs math.txt:61).
- Survivor directions: agree to `~5e-9` (math.txt:76–80 signed maxdiff vs the SymPy 15-digit prints).
- $\sigma_{\rm transfer}=2.31561904386055$ (sympy.txt:24 vs math.txt:96) and all eight ceilings agree to `<5e-9` (math.txt:103–117).
No `engine_disagreement`.

## Verdict justification

Clean. I attacked the load-bearing checks and they held: (a) the determinant is not a tautology — `.py` hard-pins it but `.wl` derives it independently and matches, so the $(200+147\pi^2)$ no-go factor is doubly anchored; (b) the slope-law assertions (A3,A5,A6) are non-trivial residual-zero checks of distinct constructions of $\Xi_1$ (direct $\delta N_0/N_0$ vs the $m,i,h$ recombination vs the sample specialization), not the same object compared to itself; (c) the survivor `assert_close` / `expectDirectionClose` checks pin actual eigen-directions, not free-floating literals; (d) symbol domains are physically justified (all primitives positive; `Delta != 0`, `omUSym^2 omWSym^2 - rSym^2 != 0` guard the poles). I read the card, notes, and appendix: the script's verified claims match the paper's stated deliverables exactly, the `200+147π²` correction is reconfirmed across all four carriers, and value reconciliation (below) is complete with zero misalignments. The only nits are deferred/cosmetic: the card's "Mathematica audit: none yet" prose lag and the stale `Stage-242/Stage-243` self-labels in the committed `.wl` output band — both already tracked under the deferred numbering/script-output-band plan, neither a math finding.

## Self-test notes

- **Variable independence (trap 1):** the `deltaMixed`/`eps`-dressing log-derivatives act on expressions that genuinely depend on the dressed primitives (P, Δ, Q, etc.), so no `D[expr,var]` collapses to identically zero; the rank-3 / nullity-2 assertions and the nonzero determinant confirm the variation is live.
- **Trivial-case / nonzero (trap 3):** the determinant `expectTrue ... =!= 0` and `det_ih != 0` are real nonzero literals (the factor `200+147π²>0`, `98π²-25>0`), so the no-go is substantive, not vacuously passing.
- **Symmetry/parity (trap 2):** n/a — no unbounded-domain integrals in this stage.
- **Paper round-trip (trap 5):** re-checked card + notes + appendix; no fix prescribed (clean), so no new misalignment introduced.

## Value Reconciliation (pass-2 augmentation)

reconciliation: complete; 16 deliverable values checked, 0 misaligned.

Authoritative record = script source + committed `.txt` outputs (both fresh: sympy.txt and math.txt mtimes Jun 2 17:25 post-date neither script's last substantive edit in a way that changes content — outputs match current source line-for-line; the determinant, slope, survivor, and ceiling literals in the source all appear verbatim in the outputs). The `.tex` card is deliberately terse (states only the symbolic theorem $\Xi_1=2\delta\ln(P/\Delta)$); the per-stage notes `.md` is the natural carrier for the numeric deliverables, and the appendix carries the symbolic forms. Per the augmentation guards, a value living correctly in the `.md` is a MATCH even when the card omits it.

| value | source (py/wl + output line) | .tex/.md location | status |
|---|---|---|---|
| $\Xi_1=2\delta\ln\Lambda$ (pure-transfer theorem) | py:171 / wl:171; sympy.txt:8, math.txt:10 | stage_227.tex:15 (Output); appendix:758; md:172 | MATCH |
| $\Lambda=(G_W/\Omega_W^2)(1+I)/(1-H)$ | py:176 / wl:182 | appendix:765; md:199–205 | MATCH |
| $I=RG_U/(\Omega_U^2G_W)$ (def) | py:174 / wl:177 | appendix:767; md:210 | MATCH |
| $H=R^2/(\Omega_U^2\Omega_W^2)$ (def) | py:175 / wl:177 | appendix:769; md:212 | MATCH |
| slope law $\Xi_1=2[m+\frac{I}{1+I}i+\frac{H}{1-H}h]$ | py:187 / wl:213 | md:241–247 | MATCH |
| sample law $\Xi_1=2m+\frac6{19}i+\frac{50}{98\pi^2-25}h$ | py:189 / wl:214; sympy.txt:13, math.txt:30 | md:266–270 | MATCH |
| $I_{\rm sample}=3/16$ | py:192 / wl:211; sympy.txt:11, math.txt:28 | md:260 | MATCH |
| $H_{\rm sample}=25/(98\pi^2)$ | py:193 / wl:212; sympy.txt:11, math.txt:29 | md:261 | MATCH |
| rigidity det factor **$200+147\pi^2$** (co-loading no-go) | py:215 / wl:250; sympy.txt:16, math.txt:61 | md:294 | MATCH |
| $v_i\approx(0.4528,-0.2942,-0.8282,-0.0405,0.1446)$ | py:233 / wl:268; sympy.txt:19 | md:342 | MATCH |
| $v_h\approx(0.6656,-0.3894,0.4671,0.0361,0.4310)$ | py:234 / wl:269; sympy.txt:20 | md:351 | MATCH |
| $v_m\approx(0.1339,-0.1059,-0.9824,-0.0539,-0.0529)$ | py:235 / wl:270; sympy.txt:21 | md:361 | MATCH |
| $\sigma_i=1.26576248$, $\sigma_h=2.04509123$, $\sigma_m=0.29342952$ | py:250–252 / wl:282–284; sympy.txt:19–21 | md:347,356,366 ($\sigma_i$=md:405,$\sigma_h$=md:416,$\sigma_m$=md:427) | MATCH |
| $\sigma_{\rm transfer}=2.31561904$ | py:269 / wl:308; sympy.txt:24, math.txt:96 | md:389 | MATCH |
| transfer ceilings $0.15889070$, $0.31854077$ | py:282 / wl:309–310; sympy.txt:24 | md:394,398 | MATCH |
| i/h/m ceilings (0.29067881/0.58274682; 0.17990900/0.36067783; 1.25389678/2.51378617) | py:283–285 / wl:311–316; sympy.txt:25–27 | md:409–411,420–422,431–433 | MATCH |

INTERNAL (scaffolding / intermediate, no prose-carrier expected, no finding): the Stage-225 K-compatibility value (`K_compat`, sympy intermediate / math.txt:5 `24.4737…`); the printed sample $\Lambda$ numeric `95*sqrt(2)*pi/(2(-25+98pi^2))` (sympy.txt:9, an intermediate print — notes carry $\Lambda$ only symbolically); the transfer-corridor basis vectors `t_1,t_2` (sympy.txt:4–5 / math.txt:55–56, Stage-226 carry-forward used to build the corridor); the budget constants `0.367930…`/`0.737619…` (Stage-224 transported budgets fed in as inputs); residual/diff/maxdiff verification values; rank/nullity flags.

All 16 stated deliverable values reconcile across `.py`/`.wl`/outputs and the `.tex`/appendix/`.md`. No MISMATCH and no MISSING-DELIVERABLE.

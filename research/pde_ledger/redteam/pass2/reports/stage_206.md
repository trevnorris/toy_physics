---
unit_id: 206
batch: VI.1
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-06-09T00:00:00Z
verdict: findings
stop_cold: null
findings_count: 3
paper_alignment: partial
scripts_checked:
  sympy: present
  mathematica: present
  engines_agree: true
  outputs_fresh: false
docs_read:
  paper_stage_tex: present
  notes_stage_files: [moving_throat_pde_stage206_certified_ray_ranking_and_local_bracketing_theorem.md]
  paper_appendix: present
---

# Audit unit 206 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_206.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage206_certified_ray_ranking_and_local_bracketing_theorem.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part06.tex` (rows 43, 210, 737-781, 918)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage206_certified_ray_ranking_and_local_bracketing_theorem_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage206_certified_ray_ranking_and_local_bracketing_theorem_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage206_certified_ray_ranking_and_local_bracketing_theorem_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage206_certified_ray_ranking_and_local_bracketing_theorem_mathematica_audit.txt`

## What the paper claims

`\stagefield{Output}`: "Certified monotone bracket theorem, certified turning-ray bracket theorem, pairwise ray-ordering theorem, and local search-sieve admissibility test." (stage_206.tex:15). The notes enumerate eight deliverables (notes §Purpose), the load-bearing ones being: (1) the forward quadratic root map `\(\mathcal T(H_0,K_0;c)=-2H_0/(K_0+\operatorname{sgn}(K_0)\sqrt{K_0^2-2cH_0})\)` with its zero-curvature limit `\(\mathcal T(\cdot;0)=-H_0/K_0\)` and strict monotonicity `\(\partial_c\mathcal T = \mathcal T^2/(2\sqrt{\Delta})>0\)`; (2) certified monotone bracket `\([\tau_{\rm lo},\tau_{\rm hi}]=[\mathcal T(\cdot;\underline K_1),\mathcal T(\cdot;\overline K_1)]\)` with strict descent `\(H_{\mathbf s}'<0\)` on the bracket; (3) collapse to the **Stage 205** quadratic log predictor at zero envelope width (notes §3.3, lines 196-208); (4) bracket-width law `\(W=\mathcal T(\bar K_1)^2/(2\sqrt{K_0^2-2\bar K_1 H_0})\,\Delta K_1 + O(\Delta K_1^3)\)` and its zero-mean form `\(W=\tau_0^2/(2|K_0|)(\overline K_1-\underline K_1)+\dots\)`; (5) turning-ray brackets `\(\tau^{\rm(tp)}=\sqrt{-2H_0/c}\)`; (6) pairwise ordering: disjoint brackets `\(\tau_{\rm hi}^{(a)}<\tau_{\rm lo}^{(b)}\Rightarrow\tau_*^{(a)}<\tau_*^{(b)}\)`; (7) the monotone/turning admissibility tests of §8. `\stagefield{Verification}` (stage_206.tex:11) states "Mathematica audit: none yet."

## What the script claims to verify

The SymPy script verifies, section by section: (I) the posited root map solves `\(H_0-k\tau+\tfrac12 c\tau^2=0\)`, its `c→0` limit is `H0/k`, the implicit `∂T/∂c` identity, and the endpoint-slope identity; (II) lower/upper bracket endpoints solve their respective quadratics and collapse when `cL=cU`; (III) the small-envelope width Series expansion matches the leading `η` term, plus the zero-mean form; (IV) turning roots solve `\(H_0-\tfrac12 c\tau^2=0\)`; (V) the bracket map collapses to the log predictor (labeled "Stage 239"); (VI) the pairwise ordering implication (via constructive slack parametrization + `sp.ask`) and a no-separation counterexample, plus the monotone/turning sieve admissibility predicates. The Mathematica script mirrors the deliverables but derives the root map and turning root via `Solve`+`SelectFirst`, proves pairwise ordering via `Resolve[ForAll[...]]` quantifier elimination, and adds two strict-sign assertions (M3 derivative sign, M4 strict descent) the SymPy script lacks.

## Paper ↔ script cross-check

| Paper deliverable | Script check | Status |
|---|---|---|
| Root map `\(\mathcal T\)` + quadratic | py I, wl M1 | match |
| Zero-curvature limit `-H0/K0` | py I, wl M2 | match |
| Monotonicity `∂T/∂c>0` | py I (identity only), wl M3 (identity + strict sign) | match |
| Strict descent `H'<0` on bracket | wl M4 (sign); py: identity only | match (wl stronger) |
| Bracket endpoints + degenerate collapse | py II, wl M5 | match |
| Width law (leading) + zero-mean | py III, wl M6 | match |
| Turning-ray brackets | py IV, wl M7 | match |
| Collapse to Stage **205** quad-log predictor | py V — labeled "**Stage 239**" | mismatch (label) |
| Pairwise ordering theorem | py VI (`sp.ask`), wl F3a (`Resolve/ForAll`) | match |
| Sieve admissibility tests | py VI, wl F3b | match |
| `\stagefield{Verification}`: "Mathematica audit: none yet" | a `.wl` exists and passes | mismatch (card stale) |

Dominant pattern: the math deliverables all match across both engines; two prose/label mismatches (card Verification field; "Stage 239" label) → `paper_alignment: partial`.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored? |
|---|---|---|---|---|---|
| A1 | sympy | 40 | `expect_zero(quadratic residual)` | root map (1) | yes |
| A2 | sympy | 41 | `limit(Tau,c,0)-H0/k==0` | zero-curv limit (1) | yes |
| A3 | sympy | 42-45 | `diff(Tau,c)-Tau^2/(2√Δ)==0` | monotonicity (1) | yes (identity only) |
| A4 | sympy | 46-49 | `(-k+cTau)+√Δ==0` | descent identity (2) | yes (identity only) |
| A5 | sympy | 62-69 | endpoints solve quadratics | bracket (2) | yes |
| A6 | sympy | 70-77 | degenerate collapse | bracket (2)/§3.3 | yes |
| A7 | sympy | 96 | Series width = leading term | width law (4) | yes |
| A8 | sympy | 104 | zero-mean width law | width law (4) | yes |
| A9 | sympy | 117-128 | turning roots solve quad + deriv | turning (5) | yes |
| A10 | sympy | 144-147 | log-predictor collapse | §3.3 collapse | yes (mislabeled 239) |
| A11 | sympy | 163-169 | pairwise ordering via constructive gap | ordering (6) | partial (manifest by construction) |
| A12 | sympy | 217-226 | sieve admissibility predicates | §8 tests (7) | yes |
| M1 | math | 84-85 | Solve+SelectFirst root, matches posited | root map (1) | yes (derive) |
| M3 | math | 103-107 | `D[rootMap,c]>0` strict | monotonicity (1) | yes (sign) |
| M4 | math | 116-120 | `K0+c·rootMap<0` strict | descent (2) | yes (sign) |
| M7 | math | 161-169 | Solve+SelectFirst turning root | turning (5) | yes (derive) |
| F3a | math | 185-198 | `Resolve[ForAll[Implies]]` QE + non-taut relaxed | ordering (6) | yes (QE) |
| F3b | math | 219-230 | admissibility via `Refine` | §8 tests (7) | yes |

## Findings

### F1 — paper_misalignment

**Severity:** medium
**Subtype:** notes_contradicts_script (card Verification field is stale)
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_206.tex:11`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage206_certified_ray_ranking_and_local_bracketing_theorem_mathematica_audit.wl` (whole file)

**What's wrong:**
The card's Verification field (stage_206.tex:11) states: "SymPy audit: \StageFile{...sympy_audit.py}. Mathematica audit: none yet." But a fully-passing Mathematica audit `.wl` exists (added in the pass-1 dual-engine retrofit; mathematica output PASSES all M1-M7, F3a, F3b). The appendix `\claimstatus` (stage_appendix_part06.tex:918) likewise does not record the second engine. The card under-states the verification coverage.

**Why this matters:**
The paper claims less verification than actually exists; an external reader/auditor cannot find the Mathematica engine. This is a prose↔artifact mismatch, not a math error.

**Required change:** (paper-side — user gate) Update the card Verification field to cite the Mathematica audit `.wl`. Codex does not auto-edit paper/. See `## Resolve before fix_loop` in the directive.

**Verification:** card line cites the `.wl` path; matches the file on disk.

### F2 — paper_misalignment

**Severity:** low
**Subtype:** notes_contradicts_script
**Files:**
- notes …stage206….md:200-208 quote: "If the curvature envelope collapses to a point, `\(\underline K_1=\overline K_1=L_1\)` … which is exactly the **Stage 205** quadratic logarithmic predictor"
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage206_certified_ray_ranking_and_local_bracketing_theorem_sympy_audit.py:131,133,145` quote: "V. Collapse to the **Stage 239** quadratic log predictor" / `expect_zero("Stage 206/239 log-predictor collapse", ...)`

**What's wrong:**
Section V of the SymPy script labels the collapse target "Stage 239" in three places (the `#` comment line 131, the `subbanner` line 133, and the assertion name line 145). The notes (§3.3) and the whole stage-204/205/206 narrative state the collapse is to the **Stage 205** quadratic log predictor. Stage 239 is an unrelated stage ("Rigid-mouth physical normal form"). The *math* is correct — `Tau.subs({H0:log(Phi0), k:-L0, c:L1})` equals the log root map `-2log(Phi0)/(L0 - sqrt(L0^2-2 L1 log(Phi0)))` for `L0<0`, which is exactly the Stage 205 oriented quadratic log predictor — only the stage label is wrong. This is a numbering-drift artifact (known incomplete-renumber class).

**Why this matters:**
Label-only; does not affect the verified identity. But it points a reader at the wrong upstream stage and the SymPy `.txt` output propagates it (and an even older "STAGE 189" banner — see F3). Canonical answer per the notes is unambiguously 205.

**Required change:** (routes to user gate per convention) `.py:131,133,145` "239" → "205". This is a script-side comment/label edit; once the user confirms the 205 direction, Codex may apply.

**Verification:** `grep -n 239 …sympy_audit.py` returns nothing; section-V labels read "Stage 205".

### F3 — stale_output

**Severity:** low (informational; verifier will re-run)
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage206_certified_ray_ranking_and_local_bracketing_theorem_sympy_audit.txt` (mtime 2026-06-02 11:32)
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage206_certified_ray_ranking_and_local_bracketing_theorem_sympy_audit.py` (mtime 2026-06-03 15:59)

**What's wrong:**
The SymPy `.py` is NEWER than its saved `.txt`. The content confirms the staleness: the saved output's banner reads "STAGE 189 — CERTIFIED RAY RANKING…" and "STAGE 189 SYMPY AUDIT PASSED" (txt:3,115), whereas the current `.py` banner reads "STAGE 206" (py:28,228). So the committed output predates the stage renumber. The Mathematica `.txt` (06-02 11:32) is NOT stale relative to its `.wl` (06-02 11:26) and already reads "STAGE 206".

**Why this matters:**
The committed SymPy output does not reflect the current script. The math results are unaffected (all checks still `= 0`), so this is informational — but the verifier must refresh the SymPy `.txt` so the banner and section-V labels match the live script (especially after F2's 239→205 fix).

**Required change:** Re-run `python3 …sympy_audit.py` and recommit the refreshed `.txt`. (Orchestrator/verifier handles the re-run.)

**Verification:** new `.txt` banner reads "STAGE 206"; mtime newer than `.py`.

## Independent-derivation check (Mathematica)

Verdict: **INDEPENDENT.** The discriminating operation is the load-bearing root map.

- **SymPy posits** the closed form, then verifies it (posit-and-check):
  `Tau = sp.simplify(2*H0/(k + sp.sqrt(k**2 - 2*c*H0)))` (py:36), then `expect_zero("quadratic residual", H0 - k*Tau + 1/2 c Tau**2)` (py:37,40).
- **Mathematica derives** the root via Solve and selects the physical branch by its zero-curvature limit, then confirms the posited closed form equals the derived branch:
  `solveRoots = ... tau /. Solve[q[c,tau]==0, tau]` (wl:72); `selectedRoot = SelectFirst[solveRoots, TrueQ[FullSimplify[linearLimit[#] == H0/k]]&]` (wl:74-77); `expectZero["M1 closed root - Solve-selected branch", rootMap - selectedRoot, ...]` (wl:85).

Same independence pattern recurs for the turning root: SymPy posits `sqrt(2*H0/a)` (py:111) and checks it; Mathematica runs `Solve[H0 + 1/2 Klower tau^2 == 0, tau]` then `SelectFirst[... # > 0 ...]` (wl:160-165) and confirms the posited form matches. For pairwise ordering the methods genuinely differ: SymPy builds `tau_star_b - tau_star_a` from a slack parametrization so the gap is `slack_a_hi + slack_b_lo + strict_sep` and applies `sp.ask(Q.positive(...))` (py:159-163); Mathematica runs full quantifier elimination `Resolve[ForAll[vars, Implies[premise, conclusion]], Reals]` (wl:185-188). The Mathematica script additionally carries strict-sign assertions (`D[rootMap,c]>0`, `K0+c·rootMap<0`) that the SymPy script never makes — it is doing strictly more, not echoing. The only same-operation overlap is the width-law `Series`/`Normal[Series]` expansion (py:88 / wl:138), but that operates on the independently-derived `rootMap` in the `.wl` and the `.wl` adds an `eta^2`-coefficient-cancellation check (M6) absent from `.py`. Not a transliteration.

## Engine cross-check

Both engines pass all checks (SymPy: every `expect_zero = 0`, all boolean checks as required; Mathematica: every `PASS:` line, no `FAIL:`). The shared root map prints identically up to CAS form: SymPy `2*H0/(k + sqrt(-2*H0*c + k^2))`; Mathematica `(-2*H0)/(-k - Sqrt[-2*c*H0+k^2]*Sign[k])` (= same under `k>0`). Width law, turning roots, and log-collapse forms agree. No `engine_disagreement`.

## Verdict justification

The mathematics of all eight deliverables holds up under attack and is corroborated by two genuinely independent engines (SymPy posits-and-verifies; Mathematica derives via Solve/QE and adds the strict-sign assertions SymPy omits). I attacked: (a) tautology in the pairwise gap — SymPy's constructive gap is manifest-positive by construction (weak on its own) but the Mathematica `Resolve/ForAll` QE is a genuine proof and the no-separation relaxed check confirms non-triviality; (b) the degenerate-collapse checks (py II, wl M5) are non-trivial substitutions, not identities-by-construction; (c) symbol domains — `H0>0, k>0` and `disc>0`/`branchAssumptions` are justified by the oriented setup `H0>0, K0=-k<0`; (d) transliteration — rejected, the root map is derive-vs-posit. The verdict is `findings` solely because of three documentation/freshness defects: the card's "Mathematica audit: none yet" is now false (F1), section V mislabels the Stage 205 collapse target as "Stage 239" (F2, numbering drift), and the SymPy `.txt` is stale (still banners "STAGE 189") (F3). No math finding; no stop-cold.

## Value Reconciliation (pass-2 augmentation)

The scripts emit no numeric benchmark/figure-of-merit constants; every deliverable is a symbolic closed form. Each is reconciled against the notes (authoritative carrier) and the Part VI appendix.

| value | source (py/wl + output line) | .tex/.md location | status |
|---|---|---|---|
| Root map `\(\mathcal T=-2H_0/(K_0+\operatorname{sgn}(K_0)\sqrt{K_0^2-2cH_0})\)` | py:36 / wl:80; txt I, M1 | notes:152-164 (boxed); appendix eq:765-766 | MATCH |
| Zero-curvature limit `\(-H_0/K_0\)` (here `H0/k`) | py:41 / wl:89; txt | notes:176-178 | MATCH |
| Monotonicity `\(\partial_c\mathcal T=\mathcal T^2/(2\sqrt\Delta)>0\)` | py:42-45 / wl:97-107 | notes:185-192 | MATCH |
| Bracket endpoints `\(\tau_{\rm lo}=\mathcal T(\underline K_1),\tau_{\rm hi}=\mathcal T(\overline K_1)\)` | py:56-57 / wl:124-125 | notes:217-222; appendix:769-774 | MATCH |
| Width leading term `\(\mathcal T(\bar K_1)^2/(2\sqrt{K_0^2-2\bar K_1H_0})\,\Delta K_1\)` | py:89-96 / wl:143-147 | notes:391-399 | MATCH |
| Zero-mean width `\(\tau_0^2/(2|K_0|)\,(\overline K_1-\underline K_1)\)` | py:101-104 / wl:149-154 | notes:410-416 | MATCH |
| Turning root `\(\tau^{\rm(tp)}=\sqrt{-2H_0/c}\)` (`sqrt(2H0/a)`) | py:111 / wl:166; txt IV, M7 | notes:327-331; appendix:776-781 | MATCH |
| Log-predictor collapse `\(-2\log\Phi_0/(L_0-\sqrt{L_0^2-2L_1\log\Phi_0})\)` | py:136-139; txt V | notes:200-208 (Stage **205**) | MATCH (math); label "Stage 239" → see F2 |
| Pairwise ordering `\(\tau_{\rm hi}^{(a)}<\tau_{\rm lo}^{(b)}\Rightarrow\tau_*^{(a)}<\tau_*^{(b)}\)` | py:163-169 / wl:185-198 | notes:438-448; appendix:731-781 | MATCH |
| Sieve admissibility (monotone/turning) | py:196-226 / wl:202-230 | notes:478-492 | MATCH |

INTERNAL scaffolding (no finding): `expect_zero`/`expectZero`/`pass`/`fail` helpers, `disc`/`q`/`branchAssumptions`, the `relaxed_counterexample` dict, `solveRoots`/`selectedRoot`/`turningRoots` intermediate Solve outputs, pass/fail flags.

reconciliation: complete; 10 deliverable values checked, 0 value-mismatched (the lone discrepancy is a stage-number *label*, F2, not a value).

## Self-test notes

(1) Variable independence: every `sp.diff`/`D[...]` is in a variable the expression genuinely depends on — `diff(Tau,c)` (Tau depends on c), `diff(TauTurnA,a)` (depends on a), `D[rootMap,c]`, `D[tauTp,a]`; none is identically-zero-by-non-dependence. (2) Symmetry/parity: no unbounded integrals — n/a. (3) Trivial-case: substituting `cL=cU=c` in II/M5 genuinely collapses (non-trivial), and the width `eta^2` coefficient genuinely cancels by odd-order Taylor symmetry (M6), confirming the leading-`eta` law is the real content. (4) Paths: F1/F2/F3 are paper/label/freshness, no new-script path needed. (5) Paper round-trip: the F2 fix is a 239→205 label correction that *removes* a misalignment and matches the notes; it introduces no new numeric constant.

---
unit_id: 191
batch: V.3
auditor_model: Opus 4.8 (1M context)
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
  notes_stage_files: [moving_throat_pde_stage191_minimal_pde_data_packet.md]
  paper_appendix: present
---

# Audit unit 191 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_191.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage191_minimal_pde_data_packet.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part05.tex` (rows/blocks at lines 113, 1138-1208, 1216-1233)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage191_minimal_pde_data_packet_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage191_minimal_pde_data_packet_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage191_minimal_pde_data_packet_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage191_minimal_pde_data_packet_mathematica_audit.txt`

## What the paper claims

`\stagefield{Output}` (card line 15): "Defines Packet A, Packet B, \(\Delta_{\rm branch}\), \(\Delta_{\rm orbit}\), and the exact two-packet home-stretch theorem." This is a packaging/compiler stage, not a new closure law (notes Status, lines 5-8). The distinct deliverables, taken from the notes' boxed equations, are: (D1) the exact single-lane response moments `nu2 = -D2/D0`, `nu4 = (D2^2 - D0 D4)/D0^2` and prefactor moments `P0,P2,P4` (notes 1.1-1.2); (D2) the grouped weighted trace/anomaly map `xbar,a_x,b_x` and its exact inverse `x20=xbar+4a_x`, `x21=xbar-a_x+b_x`, `x22=xbar-a_x-b_x` (notes 1.3); (D3) the one-pole defect `Delta_pole := nu4_bar - 4 nu2_bar^2` and normalization defect `Delta_norm := mhat0^2 P0_bar - 54 G c_s^5/(5 a^5 c^5)` (notes 1.4); (D4) the final 8-vector `Delta_branch = (a_2,b_2,a_4,b_4,a_{P0},b_{P0},Delta_pole,Delta_norm)` with vanishing on the isotropic one-pole normalized branch (notes 1.5; appendix eq:app-part05-Delta-branch-full line 1221); (D5) Packet B's three equivalent orbit representations and the exact interconversion laws `R_tr=m_T^{1+chi0*}`, `R_nt=m_mu/(m_K m_T^{F*})`, `R_eta=1/m_K`, their logs `q = ln R`, and the inverse `m=exp(...)` (notes 2.1-2.2); (D6) `Delta_orbit := (q_tr,q_nt,q_eta)` with zero-set = orbit lock (notes 2.3); (D7) the home-stretch theorem "reduced closure complete ⟺ Delta_branch=0 and Delta_orbit=0" (notes 3.1). The appendix corroborates the `Delta_branch` ordering (line 1221) and the equivalent one-pole surface `D0 D4 + 3 D2^2 = 0` (line 1233).

## What the script claims to verify

Both engines verify the same seven deliverable families: (I) the single-lane response/prefactor Taylor compilers and the one-lane pole-defect identity; (II) the grouped trace/anomaly map, its inverse/reconstruction, the G=diag(1,2,2) orthogonality of the three basis vectors, and the weighted-projector algebra (idempotence, completeness, mutual annihilation); (III) the assembly of `Delta_branch` from the three lanes, its anisotropy-coordinate collapse on an isotropic Packet A, and its exact vanishing on the isotropic one-pole normalized branch (`D4 = -3 D2^2/D0`, `N0 = D0 * 54 G c_s^5/(5 a^5 c^5 mhat0^2)`); (IV) the Packet B interconversion laws and round trips, plus orbit-lock collapse to 1; (V) the formal "join of two zero-sets = both zero-sets" tautology and that the reduced closure packet is the branch packet followed by the orbit packet. Every check is an `expect_zero` / `expectZero[...]` / `expectTrue[...]` over symbolic (not numeric-substituted) parameters, so the residuals are exact symbolic zeros, confirmed in both saved `.txt` outputs.

## Paper ↔ script cross-check

| Paper deliverable | Script-side check | Status |
|---|---|---|
| D1 response/prefactor moments | py:75-100 `response/prefactor series compiler == 0`; wl:184-210 `Series`+`Coefficient` & `D[]`-route & SymPy-compiler match | match |
| D2 grouped map + inverse | py:116-127 `x20/x21/x22 inverse`; wl:218-234 `Solve`-route + reconstruction | match |
| D2 projector algebra (notes/appendix) | py:137-146; wl:235-247 metric-built projectors | match |
| D3 Delta_pole one-lane identity | py:102-108 `+ (3 D2^2+D0 D4)/D0^2`; wl:203,211-214 | match |
| D3 Delta_norm constant | py:177; wl:284 `mhat^2 Pbar0 - 54 G cs^5/(5 a^5 c^5)` | match |
| D4 Delta_branch assembly+iso+target | py:178-210; wl:275-323 | match |
| D5 Packet B interconversion | py:222-262; wl:331-365 | match |
| D6 Delta_orbit zero-set | py:260-266; wl:366 `Solve` zero-set | match |
| D7 home-stretch theorem | py:280-282 (printed); wl:368-380 `LogicalExpand`/`expectTrue` | match (see note) |

`paper_alignment: aligned`. Every boxed notes equation has a non-tautological script-side check, in both engines.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored? |
|---|---|---|---|---|---|
| A1 | sympy | 87 | `series(Y)-(1+nu2 w^2+nu4 w^4)==0` | D1 | yes |
| A2 | sympy | 100 | `series(Pref)-(P0+P2 w^2+P4 w^4)==0` | D1 | yes |
| A3 | sympy | 105-108 | `Delta_pole_one + (3D2^2+D0D4)/D0^2==0` | D3 | yes |
| A4 | sympy | 125-127 | grouped inverse `x2j_rec - x2j == 0` | D2 | yes |
| A5 | sympy | 137-146 | basis G-orthogonality + projector algebra | D2 | yes |
| A6 | sympy | 190-201 | iso collapse of a/b coords; mean = one-lane moment | D4 | yes |
| A7 | sympy | 210 | `Delta_branch` on iso one-pole normalized branch == 0 | D4 | yes |
| A8 | sympy | 257-262 | Packet B round trips + orbit-lock collapse | D5/D6 | yes |
| B1 | mathematica | 207-208 | `Series/Coefficient - D[]-derivative route == 0` | D1 | yes (extra independent route) |
| B2 | mathematica | 209-210 | native Taylor − SymPy closed-form compiler == 0 | D1 | yes |
| B3 | mathematica | 211-214 | one-lane pole defect identity | D3 | yes |
| B4 | mathematica | 233-234 | `Solve`-route grouped coords − SymPy; reconstruction | D2 | yes |
| B5 | mathematica | 235-247 | metric-built projectors − literal targets; algebra | D2 | yes |
| B6 | mathematica | 313-319 | lane coeffs − targets; iso means − one-lane | D4 | yes |
| B7 | mathematica | 315,322 | `Delta_branch` native − SymPy; iso one-pole == 0 | D4 | yes |
| B8 | mathematica | 358-366 | Packet B `Solve`-inversion round trips; `Delta_orbit` `Solve` zero-set | D5/D6 | yes |
| B9 | mathematica | 380 | `expectTrue` formal closure equivalence | D7 | yes (structural tautology, by design) |

Note on B9 / py V: the home-stretch "theorem" check is the propositional identity `And[branch=0, orbit=0] ⟺ (And[branch=0] && And[orbit=0])` over 11 abstract boolean slots — it is a structural tautology by construction (and the `.py` only *prints* it, py:280-282, rather than asserting). This is intentional ledger scaffolding, not a load-bearing physics assertion, and the real content (D1-D6) is exercised above. Not flagged: it correctly mirrors the notes' "iff" statement and does not pretend to verify physics. The substantive home-stretch content (that everything compiles to the two packets) is exercised by D4/D6, not by B9.

## Findings

None.

## Independent-derivation check (Mathematica)

GENUINELY INDEPENDENT. For every load-bearing deliverable the `.wl` reaches the result by a different operation than the `.py`, and several are cross-checked by a *third* internal route. Evidence:

1. **Taylor moments (D1).** `.py` builds the expected polynomial from the *closed-form* moments and compares to `sp.series`: `Y_expected = 1 + nu2*omega**2 + nu4*omega**4` then `expect_zero(Y_series - Y_expected)` (py:80,87) — it never independently extracts coefficients. `.wl` extracts coefficients natively two ways — `responsePair` uses `Coefficient[Normal[Series[shape,{w,0,4}]],w,2]` (wl:113-119), and `responseDerivativeRoute` uses `(D[singleResponse,{w,2}]/2)/.w->0` (wl:190-195) — cross-checks them against each other (`nativeResponse - responseDerivativeRoute`, wl:207) AND against the SymPy closed form (wl:209). Different route + an extra cross-route the `.py` lacks.

2. **Grouped basis coordinates (D2).** `.py` writes the coordinates by formula (`grouped_trace_anomaly`, py:35-39) and the projectors as *hardcoded literals* (`Pbar/Pa/Pb`, py:133-135). `.wl` *solves* the basis decomposition — `Solve[Thread[triple == xMean*traceVec + xAnom*anomVec + xSplit*splitVec], {xMean,xAnom,xSplit}, Reals]` (wl:159-168) — and *builds* the projectors from the metric — `Outer[Times,v,v.groupMetric]/(v.groupMetric.v)` (wl:178-180) — then compares both against the SymPy formulas/literals. Genuinely different construction (LinearSolve / metric-projector vs. plug-in).

3. **Packet B inversion (D5).** `.py` hand-writes the inverse `m_from_q = [exp(qtr/(1+chi0s)), exp(-qeta), exp(qnt-qeta+Fstar*qtr/(1+chi0s))]` (py:232-236) and writes `q_from_m` closed (py:242-246). `.wl` instead *logs* the `R_from_m` definitions to get q — `qFromMNative = Log /@ rFromMTarget` (wl:332) — and *solves* the linear log-map to get m back — `Solve[Thread[logToQ == qVars], {ellT,ellK,ellMu}, Reals]` then `Exp/@` (wl:344-345). The `.wl` derives both directions; the `.py` asserts them. Independent.

The only shared algebra is *definitional* (the `R_from_m` map, the orbit-lock substitutions, the `Delta_branch`/`Delta_orbit` packet ordering) — identities any second engine must restate verbatim to test — not ported derivation choreography. Applying the consistent threshold ("does the `.wl` feed the same expression to the same operation as the `.py`?"): the answer is NO for every load-bearing deliverable. This retrofit is real, unlike the V.1-173 / V.2-cluster ports.

## Engine cross-check

Both `.txt` outputs report all residuals as exact symbolic `0` (or boolean `True` for the closure-equivalence). SymPy output: every `... = 0`, every matrix `[[0,0,0],...]`, `Delta_branch on isotropic one-pole normalized branch = [0;...;0]` (lines 24,47,54,71-118,203-227,276-310). Mathematica output: every line `PASS:` with residual `{0,...}` / `True` (lines 11-110). The two engines verify the same compiler identities and agree at the exact-symbolic level they claim. No disagreement.

## Verdict justification

Clean. I read the card, the full notes derivation, and the Part-V appendix block before the scripts, and the script's verified claim matches the paper's stated claim deliverable-by-deliverable (table above; `Delta_branch` ordering and the `Delta_pole`/one-pole surface match the appendix at lines 1221/1228/1233). Attacks tried and failed: (a) looked for tautology in the moment compilers — both engines independently expand the rational functions rather than asserting `expr==expr`, so the series checks can fail; (b) looked for a hidden hardcoded `Delta_norm` constant mismatch — `54 G c_s^5/(5 a^5 c^5)` matches notes line 190 verbatim on both sides; (c) probed the isotropic one-pole target substitution for a rigged pass — `D4 = -3 D2^2/D0` is the genuine `Delta_pole=0` surface (equiv. appendix `D0 D4 + 3 D2^2 = 0`) and `N0 = D0 P0_target` the genuine `Delta_norm=0` surface, so the vanishing is a real consequence, not a definition; (d) checked the `.wl` for transliteration — it uses `Solve`/`D[]`/metric-built projectors/log-map inversion, all distinct from the `.py`'s closed-form-plug-in route. The only structural tautology (B9 / printed py V) is by-design ledger scaffolding that mirrors the notes' "iff" and carries no physics weight.

Numbering: the first-pass `STAGE 174 -> STAGE 191` relabel holds — `.wl` banners read "STAGE 191" (wl:75,382), `.py` banners "STAGE 191" (py:65,284), card `\section[Stage 191]` (line 1). No stray `174` / `208` (+17 self-label) / `242` legacy in the card, scripts, or this notes file. Outputs fresh (txt 11:35:23 > py 11:30:24, wl 11:35:04).

Doc-sync note (paper-side, out of red-team scope): card line 11 still reads "Mathematica audit: none yet." although `.wl` now exists and runs clean. This is a `.tex` prose lag for the user to correct on the paper side (Codex must not edit paper/); not a script-side finding, so no directive.

## Self-test notes

(1) Variable independence: the `.wl` derivative routes `D[singleResponse,{w,2}]` / `D[singlePrefactor,{w,4}]` differentiate expressions that genuinely depend on `w` (denominators `d0+d2 w^2+d4 w^4`), so the derivatives are non-trivial and the at-`w->0` evaluations are real — no identically-zero-derivative trap. (2) Symmetry/parity: no unbounded-domain integrals in this stage; N/A. (3) Trivial-case: substituting the isotropic + one-pole-normalized profile into `Delta_branch` collapses each of the 8 entries to 0 by the `Delta_pole=0`/`Delta_norm=0` surfaces, confirmed in both outputs — the pass is earned, not rigged. No directive written (zero findings).

## Value Reconciliation (pass-2 augmentation)

This stage is a symbolic-compiler stage: the scripts emit boxed *symbolic* identities and one closed-form numeric constant. There are no figure-of-merit numbers, no pinned `Pi_star`-style constants, no benchmark values. Every emitted deliverable value reconciles against the notes (the natural carrier; the `.tex` card is deliberately terse and lists deliverable *names* only).

| value | source (py/wl + output line) | .tex/.md location | status |
|---|---|---|---|
| `nu2 = -D2/D0`, `nu4 = (D2^2 - D0 D4)/D0^2` | py:50-52, wl:137-139; out py:14-23, wl:9 | notes lines 113,116 (boxed) | MATCH |
| `P0=N0/D0`, `P2=(D0 N2-2 D2 N0)/D0^2`, `P4=(D0^2 N4-2 D0(D2 N2+D4 N0)+3 D2^2 N0)/D0^3` | py:55-62, wl:142-152; out py:32-46, wl:10 | notes lines 123,125,128-132 (boxed) | MATCH |
| grouped map `xbar=(x20+2x21+2x22)/5`, `a_x=(2x20-x21-x22)/10`, `b_x=(x21-x22)/2` | py:35-39, wl:170-176; out py:59-70, wl:25 | notes lines 145-151 (boxed) | MATCH |
| grouped inverse `x20=xbar+4a_x`, `x21=xbar-a_x+b_x`, `x22=xbar-a_x-b_x` | py:42-46, wl:159-168 | notes lines 154-159 (boxed) | MATCH |
| `Delta_pole := nu4_bar - 4 nu2_bar^2` | py:176, wl:283; out py:48-53, wl:60 | notes line 180 (boxed); appendix line 1228 | MATCH |
| `Delta_norm := mhat0^2 P0_bar - 54 G c_s^5/(5 a^5 c^5)` | py:177, wl:284; out py:159-162 | notes lines 186-191 (boxed) | MATCH |
| `Delta_branch = (a2,b2,a4,b4,aP0,bP0,Delta_pole,Delta_norm)` | py:178, wl:275-286; out py:123-227, wl:60 | notes lines 199-205 (boxed); appendix line 1221 | MATCH |
| interconversion `R_tr=m_T^{1+chi0*}`, `R_nt=m_mu/(m_K m_T^{F*})`, `R_eta=1/m_K` | py:222-226, wl:331; out py:232-243 | notes lines 256-260 (boxed) | MATCH |
| `q = ln R`; `m=exp(...)` inverse forms | py:227-246, wl:332-350; out py:244-265, wl:81-82 | notes lines 262-276 (boxed) | MATCH |
| `Delta_orbit := (q_tr,q_nt,q_eta)`; zero-set ⟺ orbit lock | py:264, wl:354,366; out py:311-316, wl:99 | notes lines 283-296 (boxed) | MATCH |
| one-pole surface `D0 D4 + 3 D2^2 = 0` (equiv. `Delta_pole=0`) | py:107, wl:213; out py:54, wl:19 | appendix line 1233 | MATCH |

INTERNAL scaffolding (accounted for, no finding): `Y(omega)`/`Pref(omega)` printed rational shapes, the basis vectors `traceVec/anomVec/splitVec` and metric `diag(1,2,2)`, the literal projector matrices `Pbar/Pa/Pb` (and `weightedProjector` constructions), the `iso_subs`/`onepole_norm_subs`/`isoRules`/`onePoleNormRules` substitution dictionaries, `P0_target`, the `logToQ`/`logSolution` solver intermediates, the `formalClosureSplit` boolean, and all `expect_zero`/`PASS`/`True` flags.

reconciliation: complete; 11 deliverable values checked, 0 misaligned. (Separately noted above: the card's "Mathematica audit: none yet." prose lags the now-present `.wl` — a paper-side `.tex` doc-sync item for the user, not a value mismatch and not a script-side finding.)

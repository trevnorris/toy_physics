---
unit_id: 203
batch: VI.1
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
  notes_stage_files: [moving_throat_pde_stage203_free_quintuple_scalar_closure_slice_and_crossing_theorem.md]
  paper_appendix: present
---

# Audit unit 203 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_203.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage203_free_quintuple_scalar_closure_slice_and_crossing_theorem.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part06.tex` (rows 37, 236, 505-592 read; this stage's appendix narrative)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage203_free_quintuple_scalar_closure_slice_and_crossing_theorem_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage203_free_quintuple_scalar_closure_slice_and_crossing_theorem_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage203_free_quintuple_scalar_closure_slice_and_crossing_theorem_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage203_free_quintuple_scalar_closure_slice_and_crossing_theorem_mathematica_audit.txt`

## What the paper claims

The card's `\stagefield{Output}` is verbatim: *"Scalar graph slice \(\mathcal Z_*=\{\mathbf x_*^{\rm graph}(\mathbf y):\widehat\chi_Q(\mathbf y)=1\}\), graph-family tangency, and one-parameter crossing theorem."* The derivation ledger lists four moving parts: compose \(\chi_Q\) with the Stage-202 target graph, prove graph-lifted tangents lie in \(\ker M_*\), separate graph-aligned vs off-graph family packets, and apply the IVT to one-parameter graph crossings. The notes enumerate six deliverables: (1) the graph-closure scalar \(\widehat\chi_Q(\mathbf y):=\chi_Q(\mathbf x_*^{\rm graph}(\mathbf y))\); (2) the codimension-one slice theorem \(\mathcal Z_*\); (3) exact dependent-triple log tangents (\(\dot\delta,\tau_1^{\rm graph},\kappa_\eta^{\rm graph},\mu_1^{\rm graph}\)); (4) the kernel theorem \(M_*\dot{\Delta\mathbf x}_{\rm graph}=0\); (5) the same-free-quintuple decomposition into graph-lifted orbit piece plus graph errors \((E_T,E_K,E_\mu)\) with the quotient packet \((q_{\rm tr},q_{\rm nt},q_\eta)\) and its inverse; (6) the one-parameter crossing theorem via IVT. The appendix (rows 510-592) restates the tangent identities, the kernel identity \(M_*\dot{\Delta\mathbf x}_{\rm graph}=0\), and the IVT crossing theorem, and notes that the graph-lifted ray keeps the three target monomials invariant for all \(\tau\).

## What the script claims to verify

Both scripts verify, in six labeled sections: (I) the Stage-192 dependent pivot block is invertible and \(M_*S-I_3=0\) (canonical section exists); (II) the four dependent-triple graph log tangents equal the carried closed forms, obtained by differentiating the graph map; (III) the assembled graph tangent \(\dot{\Delta\mathbf x}_{\rm graph}\in\ker M_*\); (IV) the graph-error packet \(q_{\rm tr}=(1+\chi_{0,*})E_T\), \(q_{\rm nt}=E_\mu-E_K-F_*E_T\), \(q_\eta=-E_K\) reconstructed from the direct monomials; (V) the inverse compiler solves the 3x3 system and the repair vector \(\Delta\mathbf x_{\rm rep}\) reproduces \((0,0,0,0,-E_K,0,-E_\mu,-E_T)\) with \(M_*\Delta\mathbf x_{\rm rep}+q(E)=0\); (VI) a concrete Stage-202 graph path (\(\beta(\tau)=2^{2\tau-1}\)) composed with the carried Stage-197 closure scalar \(\chi_Q=3(S\beta^5+9\Sigma_5)/(3S-\Sigma_0)\), verifying the three target monomials are invariant on the lift (\(q_{\rm tr}=q_{\rm nt}=q_\eta=0\)), and that \(\widehat\Delta_Q(\tau)=32^{2\tau-1}-1\) changes sign (\(-31/32\) at \(\tau=0\), \(+31\) at \(\tau=1\)) with a unique real root \(\tau_*=1/2\) — a falsifiable witness of the IVT crossing theorem.

## Paper ↔ script cross-check

| Paper/notes deliverable | Script check | Status |
|---|---|---|
| (1) \(\widehat\chi_Q:=\chi_Q\circ\mathbf x_*^{\rm graph}\) | VI: `hat_chi_graph` built by composing Stage-197 `chi_from_stage197` with the graph path (.py 312-320 / .wl 289-294) | match |
| (2) \(\mathcal Z_*\) codim-1 slice \(\widehat\chi_Q=1\) | VI: root set of `hat_delta_graph` = {1/2}, single scalar condition (.py 349-362 / .wl 318-331) | match |
| (3) dependent-triple log tangents | II: `expect_zero` against the four carried forms (.py 158-178 / .wl 169-185) | match |
| (4) kernel theorem \(M_*\dot{\Delta\mathbf x}_{\rm graph}=0\) | III: `M_* dot(Delta x)_graph` = 0 (.py 185 / .wl 192) | match |
| (5) decomposition + packet + inverse | IV+V: `q_tr/q_nt/q_eta`, inverse compiler, repair vector (.py 224-271 / .wl 211-261) | match |
| (6) one-parameter crossing theorem (IVT) | VI: sign change `-31/32`<0<`31`, root at 1/2 (.py 358-362 / .wl 324-331) | match |
| Section heading "Stage 192 / Stage 197" provenance | `chi_from_stage197` formula = Stage-197 notes line 231; .wl names the vars `chiFromStage180`/`closureNumStage180` (see F1) | match (math) / label drift |

Every paper-side deliverable maps to a substantive, non-tautological script-side check. `paper_alignment: aligned`.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 106 | `M_* S - I_3 == 0` | section exists (precondition for kernel/decomp) | yes |
| A2 | sympy | 158-178 | 4x `expect_zero(tangent - carried form)` | (3) dependent log tangents | yes |
| A3 | sympy | 185 | `M_* dot(Δx)_graph == 0` | (4) kernel theorem | yes |
| A4 | sympy | 235-237 | `q_tr/q_nt/q_eta - carried == 0` | (5) graph-error packet | yes |
| A5 | sympy | 257-262 | inverse-compiler closed forms == 0 | (5) inverse | yes |
| A6 | sympy | 270-271 | repair vector + `M_*Δx_rep+q(E)==0` | (5) repair consistency | yes |
| A7 | sympy | 337-339 | graph-lift `q_tr=q_nt=q_eta=0` | (1) target-monomial invariance on lift | yes |
| A8 | sympy | 341-343 | `3 S·hatΔ - closure_num == 0` | (1) closure-scalar identity | yes |
| A9 | sympy | 358-362 | sign change + unique root {1/2} | (2)+(6) slice + IVT crossing | yes |
| B1-B9 | mathematica | 137,169-185,192,228-230,247-251,260-261,306-308,311,324-331 | mirror of A1-A9 via log-additive route | same claims | yes |

## Findings

### F1 — stale_output

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage203_free_quintuple_scalar_closure_slice_and_crossing_theorem_sympy_audit.txt:3,325`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage203_free_quintuple_scalar_closure_slice_and_crossing_theorem_mathematica_audit.txt:3,111`

**What's wrong:**
Both saved outputs are older than their scripts (outputs: 2026-06-02 12:32:58; scripts: 2026-06-03 15:59:11), and both carry a stale banner: the SymPy output opens with `STAGE 186 — FREE-QUINTUPLE SCALAR CLOSURE SLICE AND CROSSING THEOREM` and closes with `STAGE 186 SYMPY AUDIT PASSED`; the Mathematica output likewise says `STAGE 186 ...` / `STAGE 186 MATHEMATICA AUDIT PASSED`. The current scripts emit `STAGE 203` (`.py:65,364`; `.wl:73,333`), so these `.txt` files predate the banner correction. This is the known script/output-band numbering-drift signal (the captured transcripts are from a pre-renumber run). Aside from the banner the captured math content is consistent with the current scripts (all PASS lines and all symbolic forms match what the current source would emit), so the staleness is cosmetic but real.

**Why this matters:**
A reviewer trusting the committed transcript would see "STAGE 186" and could mis-file the evidence; the freshness gate exists precisely to keep transcripts faithful to the current scripts. Not blocking — the verifier's independent re-run regenerates correct `STAGE 203` transcripts.

**Required change:**
Re-run both scripts to refresh the committed `.txt` outputs so the banners read `STAGE 203`. No script-source edit is implied by this finding (the source banners are already correct).

**Verification:**
After re-run, both `.txt` files' first content banner reads `STAGE 203 — FREE-QUINTUPLE SCALAR CLOSURE SLICE AND CROSSING THEOREM` and the closing banners read `STAGE 203 ... AUDIT PASSED`; mtimes are newer than the scripts.

## Independent-derivation check (Mathematica)

**Verdict: INDEPENDENT.**

The discriminating operation is the *method that extracts the load-bearing objects* (graph tangents, graph-error packet, and the crossing residual). The two engines extract them by genuinely different routes:

1. **Graph map → log tangents (Section II).** The `.py` builds the *actual multiplicative power-law graph map* and differentiates the log of the power expression:
   - `.py:129-132` `deltaU_graph_t = (Ctr_tgt/(gamma_t*cetaU_t/KU_t)**(1+deltaUs))**(1/(1+chi0s))`, then `.py:144` `dln_delta_graph = normalize(sp.diff(sp.log(deltaU_graph_t), t).subs(t,0))` — SymPy must expand `log(power-of-power)` itself.
   - `.wl:148-150` constructs the tangent's *log-linear form by hand*: `logDeltaGraphT = (Log[CtrTgt] - (1+deltaUs)(Log[gammaT]+Log[cetaUT]-Log[KUT]))/(1+chi0s)`, then `.wl:159` `D[logDeltaGraphT, t] /. t->0`. Mathematica never differentiates a power tower; it differentiates a posited log-additive combination.
   This is derive-from-power (SymPy) vs posit-the-log-expansion (Mathematica) — not a syntactic port.

2. **Graph-error packet (Section IV).** The `.py` reconstructs the dimensionful constant from the perturbed quantity and takes the log of a ratio:
   - `.py:209,213-216,224` `TU = T_graph*exp(E_T)`; `Ctr = (gamma*cetaU/KU)**(1+deltaUs)*(pi**2*TU/(L**2*KU))**(1+chi0s)`; `qtr = log(Ctr/Ctr_tgt)`.
   - `.wl:207,211-215` works entirely additively: `logTU = logTGraph + ET`; `qtr = (1+deltaUs)(...) + (1+chi0s)(2 Log[Pi]+logTU-2 Log[L]-Log[KU]) - Log[CtrTgt]`. No `Exp`, no ratio, no `Log` of a constructed power.

3. **Crossing-set extraction (Section VI).** `.py:349` `sp.solveset(sp.Eq(hat_delta_graph,0), tau, domain=sp.S.Reals)` compared to `sp.FiniteSet(Rational(1,2))` (`.py:361`); `.wl:318` `Reduce[hatDeltaGraph==0, tau, Reals]` compared via `FullSimplify[realCrossings == (tau==1/2)]` (`.wl:328`). Different solver primitives (`solveset` vs `Reduce`) and different result-shape comparisons.

The only block where both engines run the *same* primitive is the Section-V inverse compiler (`sp.solve` vs `Solve` on the identical 3-equation linear system, `.py:244-252` / `.wl:234-242`) — but that is a trivial linear inversion and the load-bearing closed forms are independently re-asserted afterward (`.py:257-262` / `.wl:247-251`). It does not make the `.wl` a port; the independence of the three load-bearing extractions above is decisive. The "each CAS runs its own simplifier" defense is not invoked here — the routes themselves differ in construction, not merely in which simplifier finishes the algebra.

## Engine cross-check

Final outputs agree exactly. Both report: tangent residuals all 0; `M_* dot(Δx)_graph = 0`; `q_tr=(1+chi0s)E_T`, `q_nt=-E_K+E_mu-F_*E_T`, `q_eta=-E_K`; inverse `E_T=qtr/(1+chi0s)`, `E_K=-qeta`, `E_mu=qnt-qeta+F_*qtr/(1+chi0s)`; repair vector = `(0,0,0,0,-E_K,0,-E_mu,-E_T)`; `widehat chi_Q = 32^(2tau-1)` so `widehat Delta_Q = 32^(2tau-1)-1`; `Delta_Q(0)=-31/32<0`, `Delta_Q(1)=31>0`, root `tau=1/2`. SymPy reports the crossing set as `{1/2}`; Mathematica as `2 tau == 1` (i.e. `tau==1/2`) — same set, different normal form. No residual, sign, or factor disagreement.

## Verdict justification

The scripts hold up against the paper under attack. Attacks tried and failed: (a) **tautology** — Section VI is not "check a literal against itself": the residual `32^(2tau-1)-1` is *derived* by composing the genuine carried Stage-197 closure scalar (formula matches Stage-197 notes line 231 verbatim) with a concrete Stage-202 graph path, and the sign-change endpoints and the unique root emerge from that derivation; the `q_tr=q_nt=q_eta=0` checks independently confirm the path stays on the target graph (target-monomial invariance), which is the load-bearing premise the crossing theorem composes through. (b) **derivative-of-independent-variable trap** — the Section-II tangents are `d/dt` of expressions that genuinely depend on `t` through `lam_t,...,KW_t`, so they are non-trivially nonzero before the carried-form subtraction. (c) **positivity** — `expect_positive("graph residual denominator")` reduces to the literal `32>0` (denominator of `hat_delta`), and `S>0` is justified (it is a positive isotropic normalization). (d) **transliteration** — rejected; the `.wl` derives the load-bearing objects by a log-additive route distinct from the `.py` power-multiplicative route. The one genuine finding is the stale `STAGE 186` banner in both committed transcripts (low, informational — the verifier's re-run fixes it). I confirm I read the card, the notes, and the Part VI appendix, and the script's verified claim matches the paper's stated `Output`.

Checkpoint higher bar: **CLEARED.** Both engines present and independent; every assertion substantive and non-tautological; the load-bearing closure scalar is re-derived in-script from the carried Stage-197 algebra rather than trusted as a literal; the target-monomial invariance is an independent falsifiable gate; paper alignment is exact. The lone finding (stale output banner) is cosmetic and does not bear on the math.

## Value Reconciliation (pass-2 augmentation)

The deliverables of this stage are **symbolic theorems**, not numeric constants. The numbers that appear in the transcripts (`-31/32`, `31`, root `1/2`, `32`) are the script-internal *benchmark witness* for the IVT crossing theorem (a concrete sign-change path); they are scaffolding that demonstrates the theorem, not stated deliverable values, and they are correctly absent from the card/notes. The grep confirmed neither `2^(2tau-1)`, `32^(...)`, nor `-31/32`/`31` appear in the card or notes — as expected for an internal witness.

| value (symbolic deliverable) | source (py/wl) | .tex/.md location | status |
|---|---|---|---|
| \(\widehat\chi_Q:=\chi_Q(\mathbf x_*^{\rm graph})\) | py:312-320 / wl:289-294 | notes:144-150; tex:7,15; appendix:510-517 | MATCH |
| \(\mathcal Z_*=\{\widehat\chi_Q=1\}\) (codim-1 slice) | py:349-362 / wl:318-331 | notes:160-179; tex:15; appendix:518-527 | MATCH |
| \(d\ln\delta_{U,*}^{\rm graph}/d\tau=-\tfrac{1+\delta_{U,*}}{1+\chi_{0,*}}(\gamma_1+c_1-\kappa_U)\) | py:160 / wl:171 | notes:224-230; appendix:546-550 | MATCH |
| \(\tau_1^{\rm graph}=\kappa_U-\tfrac{1+\delta}{1+\chi_0}(\gamma_1+c_1-\kappa_U)\) | py:164 / wl:175 | notes:232-242; appendix:551-556 | MATCH |
| \(\kappa_\eta^{\rm graph}=2c_1-\kappa_U\) | py:166 / wl:177 | notes:243-251; appendix:557-560 | MATCH |
| \(\mu_1^{\rm graph}=2c_1-\kappa_U+2\kappa_W-2\lambda_1-E_*(2\gamma_1+2\lambda_1-\kappa_U-\kappa_W)-F_*\tfrac{1+\delta}{1+\chi_0}(\gamma_1+c_1-\kappa_U)\) | py:169-177 / wl:179-185 | notes:252-262; appendix:561-569 | MATCH |
| \(M_*\dot{\Delta\mathbf x}_{\rm graph}=0\) | py:185 / wl:192 | notes:304-309; appendix:570-574 | MATCH |
| \(q_{\rm tr}=(1+\chi_{0,*})E_T\) | py:235 / wl:228 | notes:364-368 | MATCH |
| \(q_{\rm nt}=E_\mu-E_K-F_*E_T\) | py:236 / wl:229 | notes:369-373 | MATCH |
| \(q_\eta=-E_K\) | py:237 / wl:230 | notes:374-378 | MATCH |
| inverse \(E_T=q_{\rm tr}/(1+\chi_{0,*}),\ E_K=-q_\eta,\ E_\mu=q_{\rm nt}-q_\eta+\tfrac{F_*}{1+\chi_{0,*}}q_{\rm tr}\) | py:257-262 / wl:247-251 | notes:380-389 | MATCH |
| \(\Delta\mathbf x_{\rm err}=(0,0,0,0,E_K,0,E_\mu,E_T)\) (repair \(=-\) of that) | py:266,270 / wl:256,260 | notes:346-362 | MATCH |
| target-monomial invariance on graph lift (\(q_{\rm tr}=q_{\rm nt}=q_\eta=0\)) | py:337-339 / wl:306-308 | appendix:641 ("three target monomials are invariant for all \(\tau\)") | MATCH |
| crossing theorem (IVT existence) | py:358-362 / wl:324-331 | notes:468-495; appendix:580-592 | MATCH |

reconciliation: complete; 13 deliverable identities checked, 0 misaligned.

**INTERNAL items (no finding):** `Mstar`, `Pdep`/`PdepInv`, `Edep`, `Sdep` (canonical-section scaffolding); `alpha=(1+deltaUs)/(1+chi0s)` (intermediate); the benchmark crossing path `beta_path=2^(2tau-1)`, `cetaU_path=cetaU_bar e^(rho tau)`, `gamma_path`, and the witness residual values `-31/32`, `31`, root `1/2`, denominator `32` (IVT witness scaffolding); `chi_from_stage197`/`closure_num_stage197` (carried Stage-197 algebra, verified against Stage-197 notes line 231/237). Note: the `.wl` names the carried-closure variables `chiFromStage180`/`closureNumStage180` (wl:289-290,294,296) while the `.py` names them `chi_from_stage197`/`closure_num_stage197` and the notes/card attribute the algebra to Stage 197 — a cosmetic variable-label drift (Stage 180 vs 197) with identical math; it is not a value mismatch (no deliverable value differs), but it is part of the same script/output-band label-drift family flagged for the dedicated numbering pass and is worth noting here for that pass.

## Self-test notes

Checked the variable-independence trap: every `sp.diff(...,t)` / `D[...,t]` in Section II acts on expressions that genuinely depend on `t` (through `lam_t,...,KW_t`), so the tangents are non-trivially nonzero pre-subtraction — no identically-zero derivative masquerading as a pass. Checked the trivial-case/positivity trap: the `graph residual denominator > 0` check reduces to the literal `32>0` and the endpoint signs `-31/32<0`, `31>0` are concrete nonzero literals straddling the root `1/2`, so the IVT hypothesis is genuinely met (not vacuously). Checked tautology in Section VI: the residual is composed from the carried Stage-197 formula and a real path, not a literal checked against itself. No directive prescribes a script edit (the sole finding is stale_output, fixed by re-run), so the paper-round-trip self-test is moot.

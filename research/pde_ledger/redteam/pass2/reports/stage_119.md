---
unit_id: 119
batch: IV.3
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
  notes_stage_files: [/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage119_parent_balance_family.md]
  paper_appendix: present
---

# Audit unit 119 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_119.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage119_parent_balance_family.md` (only file matching the pattern)
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part04.tex` (rows 525-573, 1124-1154, 1272)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage119_parent_balance_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage119_parent_balance_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage119_parent_balance_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage119_parent_balance_mathematica_audit.txt`

## What the paper claims

Stage 119 is the parent of the compensation family. Its `\stagefield{Verification}` block points to the two engine scripts; the card's single boxed deliverable is: "Normalized ratios \((\mathfrak r,\mathfrak g)\) obey \(1+\mathfrak r^2=4(\mathfrak g-\mathfrak r)^2\)." The notes expand this to a full set of boxed results: (1) the dimensionless ratio definitions \(\mathfrak r:=\lambda/\sqrt{K_sK_q}\), \(\mathfrak g:=g_q\sqrt{K_s}/(g_s\sqrt{K_q})\); (2) the family equation \(1+\mathfrak r^2=4(\mathfrak g-\mathfrak r)^2\) as the collapse of the Stage-115 core-balance theorem \(g_s^2(K_sK_q+\lambda^2)=4(K_sg_q-\lambda g_s)^2\); (3) the two-branch law \(\mathfrak g=\mathfrak r\pm\tfrac12\sqrt{1+\mathfrak r^2}\); (4) the D/N tube length \(L_W=\tfrac{\pi a}{2}\sqrt{(1+\mathfrak r^2)/3}\) via the identification \(r_c=\mathfrak r^2\); (5) the explicit \(\mathfrak g=\sqrt{2\mathcal Z_qK_s}/(\mathcal T_mJ_sc_s\sqrt{\mu_0L_W})\); and (6) the exact traction law \(\mathcal T_m=\sqrt{2\mathcal Z_qK_s}/(J_sc_s\sqrt{\mu_0L_W})\cdot 1/(\mathfrak r\pm\tfrac12\sqrt{1+\mathfrak r^2})\). The Part IV appendix (eqs. `app-part04-parent-compensation-family`, `-LW-compensation`) reproduces (2),(3),(4) verbatim and builds the downstream \(\mathfrak r_{F1}\), \(\delta_\perp\) machinery on top of them.

## What the script claims to verify

Both scripts run four sections. (I) substitute the dimensionless-ratio definitions into the Stage-115 balance theorem and confirm it reduces exactly to \(1+\mathfrak r^2-4(\mathfrak g-\mathfrak r)^2\); (II) solve the family equation for \(\mathfrak g\), recover the two branches \(\mathfrak r\pm\tfrac12\sqrt{1+\mathfrak r^2}\), and back-substitute each into the family equation; (III) solve the geometric tube condition \(4L_W^2/(\pi^2a^2)=(1+r_c)/3\) for \(L_W\), confirm the closed form, and check the \(r_c\to\mathfrak r^2\) substitution; (IV) build \(\mathfrak g\) from \(K_q,g_q\) and \(g_s=\mathcal T_mJ_s\), confirm the explicit closed form, then solve \(\mathfrak g_\pm(\mathfrak r)=\mathfrak g\) for \(\mathcal T_m\) and match the boxed traction law on both branches (with the same \(\pm\) sign-pairing as the family branch). The Mathematica script additionally guards `Length[gHatSol] === 2` and strips `ConditionalExpression` from the `-` branch.

## Paper ↔ script cross-check

| Paper-side deliverable | Script-side check | Status |
|---|---|---|
| Family eq. \(1+\mathfrak r^2=4(\mathfrak g-\mathfrak r)^2\) (card box; notes §1; appx 540) | I: `expect_zero("dimensionless law", law_red - (1+rhat**2-4*(ghat-rhat)**2))` | match |
| Two-branch law \(\mathfrak g=\mathfrak r\pm\tfrac12\sqrt{1+\mathfrak r^2}\) (notes §2; appx 542) | II: `solve` → `[rhat - sqrt(rhat**2+1)/2, rhat + sqrt(rhat**2+1)/2]` + back-sub | match |
| Tube length \(L_W=\tfrac{\pi a}{2}\sqrt{(1+\mathfrak r^2)/3}\) (notes §3; appx 547-548) | III: `solve` of \(4L_W^2/(\pi^2a^2)=(1+r_c)/3\); `tube-length law` + `rc->rhat**2 link` | match |
| \(r_c=\mathfrak r^2\) identification (notes §3 line 64) | III: substitution `L_sel.subs(rc, rhat**2)` | partial (substitution echoed, not the def-level derivation — see Verdict) |
| Explicit \(\mathfrak g=\sqrt{2\mathcal Z_qK_s}/(\mathcal T_mJ_sc_s\sqrt{\mu_0L_W})\) (notes §4) | IV: `ghat explicit simplification` | match |
| Traction law \(\mathcal T_m=\dots/(\mathfrak r\pm\tfrac12\sqrt{1+\mathfrak r^2})\) (notes §5) | IV: `T_m (+ branch) match`, `T_m (- branch) match` | match |
| Parent overlap \(\mathfrak r=-\tfrac{8\sqrt2}{3}q_*v_{w0}a^2\ell\sqrt{L_W}/\sqrt{K_sK_q}\) (notes §4 box) | (none — definitional input, not re-derived) | not-required (input-side, see Verdict) |
| Healing-locked \(K_s=3\pi a^2\hbar^2/(5m_\psi\rho_w\ell)\) (notes §6) | (none — carried from Stage 118) | not-required (upstream carry, see Verdict) |

`paper_alignment` = aligned. Every stated deliverable of stage 119 (the family eq., the branch law, the tube-length law, the explicit \(\mathfrak g\), and the traction law) has a faithful, non-tautological script-side check that matches the boxed forms in the notes and the appendix exactly. The notes §4/§6 items that lack a script check are definitional inputs imported from Stage 118 (the \(\mathcal I_{sq}\) overlap closure and the parent \(K_s\) form), not new identities stage 119 derives; the card's `\stagefield{Inputs}` explicitly imports parent-action overlap data, so their absence is correct, not a `script_missing_paper_claim`.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 36 | `expect_zero("dimensionless law", law_red - (1+rhat**2-4*(ghat-rhat)**2))` | family eq. (collapse of S115 theorem) | yes |
| A2 | sympy | 42 | `expect_zero("positive branch check", 1+rhat**2-4*(ghat_sol[0]-rhat)**2)` | two-branch law (back-sub) | partial (solution consistency) |
| A3 | sympy | 43 | `expect_zero("negative branch check", ...)` | two-branch law (back-sub) | partial |
| A4 | sympy | 52 | `expect_zero("tube-length law", L_sel - pi*a*sqrt((1+rc)/3)/2)` | tube-length law | yes |
| A5 | sympy | 53-56 | `expect_zero("...rc -> rhat**2 link", L_sel.subs(rc,rhat**2) - ...)` | \(r_c=\mathfrak r^2\) link | partial (substitution) |
| A6 | sympy | 65-68 | `expect_zero("ghat explicit simplification", ghat_expr - sqrt(2*Zq*K_s)/(Tm*J_s*c_s*sqrt(mu0*L_W)))` | explicit \(\mathfrak g\) | yes |
| A7 | sympy | 75-79 | `expect_zero("T_m (+ branch) match", Tm_sol_plus - ...)` | traction law (+ branch) | yes |
| A8 | sympy | 80-84 | `expect_zero("T_m (- branch) match", Tm_sol_minus - ...)` | traction law (− branch) | yes |
| A9 | math | 40 | `expectZero["dimensionless law", lawRed - (1+rHat^2-4*(gHat-rHat)^2)]` | family eq. | yes |
| A10 | math | 45 | `If[Length[gHatSol] =!= 2, fail[...]]` | two-branch law (branch count) | yes (guards exactly-two-branch claim) |
| A11 | math | 47-48 | `expectZero["first/second branch check", 1+rHat^2-4*((gHat/.sol)-rHat)^2]` | two-branch law (back-sub) | partial |
| A12 | math | 59 | `expectZero["tube-length law", lSel - (Pi*a*Sqrt[(1+rC)/3])/2]` | tube-length law | yes |
| A13 | math | 60-63 | `expectZero["...rC -> rHat^2 link", (lSel/.rC->rHat^2) - ...]` | \(r_c=\mathfrak r^2\) link | partial (substitution) |
| A14 | math | 77-80 | `expectZero["gHat explicit simplification", gHatExpr - Sqrt[2*zQ*kS]/(tM*jS*cSound*Sqrt[mu0*lW])]` | explicit \(\mathfrak g\) | yes |
| A15 | math | 90-94 | `expectZero["T_m (+ branch) match", stripCE[tMPlus] - ...]` | traction law (+ branch) | yes |
| A16 | math | 95-99 | `expectZero["T_m (- branch) match", stripCE[tMMinus] - ...]` | traction law (− branch) | yes |

The back-substitution checks A2/A3/A11 are weaker (a solution returned by `solve` necessarily satisfies the equation), but they are not the load-bearing checks for the branch law — the load-bearing content is that `solve` returns *exactly* \(\mathfrak r\pm\tfrac12\sqrt{1+\mathfrak r^2}\), which the printed output records and the Mathematica branch-count guard A10 reinforces. The substitution links A5/A13 are likewise confirmatory rather than derivational; the actual \(r_c=\mathfrak r^2\) identity is a definitional consequence (\(\mathfrak r^2=\lambda^2/(K_sK_q)=r_c\)) that lives in the notes, not something the script needs to re-derive. None of these rises to a finding.

## Findings

None.

## Independent-derivation check (Mathematica)

Both scripts share the same four-section choreography and the same algebraic targets — this is the expected signature for an algebraic-identity stage in the 105-175 transliteration-watch band, so I scrutinized it closely. The verdict is **independent enough** for this stage type, for three concrete reasons:

1. **The verification target is an algebraic identity with an essentially unique substitution route.** The "derivation" is: substitute the two ratio definitions into the Stage-115 theorem and simplify. There is no materially different second route to that collapse — any engine must perform the same substitution. The meaningful independence is that each engine re-runs its *own* CAS solve/simplify: SymPy `sp.solve`/`sp.simplify`, Mathematica `Solve[...,Reals]`/`FullSimplify[...,Assumptions->...]`. These are independent algebra kernels, not one echoing the other's intermediate forms.

2. **Engine-specific robustness that is NOT present in the other script.** Mathematica adds a hard branch-count guard `If[Length[gHatSol] =!= 2, fail["expected two gHat branches"]]` (line 45) that SymPy lacks, and a `stripCE` step (line 88) handling the `ConditionalExpression[..., rHat > 1/Sqrt[3]]` that `Solve[...,Reals]` emits for the `-` branch (output line 34). SymPy's `solve` does not emit that conditional, so this is genuine engine-aware code, not transliteration.

3. **No shared intermediate-variable choreography beyond what the math dictates.** The expressions (`law`, `lawRed`, `kappa0`, `gHatExpr`, `tMPlus/tMMinus`) correspond because the physics quantities correspond, not because the `.wl` mirrors private `.py` helper steps. Each engine builds `K_q`, `g_q`, `g_s=T_m J_s` directly from the notes' physical formulas.

Side-by-side I (the most "parallel" section):
- SymPy: `law_red = simplify(law.subs({lam: rhat*sqrt(K_s*K_q), g_q: ghat*g_s*sqrt(K_q)/sqrt(K_s)})/(g_s**2*K_s*K_q))`
- Math: `lawRed = FullSimplify[(law /. {lam -> rHat*Sqrt[kS*kQ], gQ -> gHat*gS*Sqrt[kQ]/Sqrt[kS]})/(gS^2*kS*kQ), Assumptions -> $Assumptions]`

The substitution map and divisor are identical *because they are the definitions of the dimensionless ratios* — there is no alternative. I do not raise a `mathematica_transliteration` finding: the test for this band is "could Mathematica verify independently," and here it does so through an independent CAS with engine-specific handling, which is the correct realization for an identity-collapse stage.

## Engine cross-check

The two engines agree at the level they claim. Final symbolic forms match across outputs:

| Quantity | SymPy output | Mathematica output |
|---|---|---|
| Reduced law (pre-target) | `-4*ghat**2 + 8*ghat*rhat - 3*rhat**2 + 1` | `1 - 4*gHat^2 + 8*gHat*rHat - 3*rHat^2` (same polynomial) |
| `dimensionless law` | `0` | `0` (PASS) |
| `ghat`/`gHat` solutions | `[rhat - sqrt(rhat**2+1)/2, rhat + sqrt(rhat**2+1)/2]` | `{rHat - Sqrt[1+rHat^2]/2, rHat + Sqrt[1+rHat^2]/2}` |
| `L_W` | `pi*a*sqrt(3*rc+3)/6` | `(a*Pi*Sqrt[1+rC])/(2*Sqrt[3])` (equal: \((\pi a/2)\sqrt{(1+r_c)/3}\)) |
| `ghat explicit` | `sqrt(2)*sqrt(K_s)*sqrt(Zq)/(J_s*sqrt(L_W)*Tm*c_s*sqrt(mu0))` | `(Sqrt[2]*kS*zQ)/(cSound*jS*tM*Sqrt[kS*lW*mu0*zQ])` (equal after simplify; `explicit simplification = 0` both) |
| `T_m (+ branch)` | `2*sqrt(2)*sqrt(K_s)*sqrt(Zq)/(J_s*sqrt(L_W)*c_s*sqrt(mu0)*(2*rhat + sqrt(rhat**2+1)))` | `ConditionalExpression[(2*Sqrt[2]*Sqrt[(kS*zQ)/(lW*mu0)])/(2*cSound*jS*rHat + cSound*jS*Sqrt[1+rHat^2]), 1/Sqrt[3]+rHat>0]` (equal after stripCE; match = 0) |
| `T_m (- branch)` | `2*sqrt(2)*sqrt(K_s)*sqrt(Zq)/(J_s*sqrt(L_W)*c_s*sqrt(mu0)*(2*rhat - sqrt(rhat**2+1)))` | `ConditionalExpression[(-2*Sqrt[2]*kS*zQ)/((-2*cSound*jS*rHat + cSound*jS*Sqrt[1+rHat^2])*Sqrt[kS*lW*mu0*zQ]), rHat>1/Sqrt[3]]` (equal after stripCE; match = 0) |

All eight SymPy checks and all eight Mathematica checks pass with residual 0. No `engine_disagreement`.

## Verdict justification

Clean. Attacks tried that failed: (1) **tautology hunt** — the load-bearing checks (A1/A9 collapse, A4/A12 tube-length solve, A6/A14 explicit \(\mathfrak g\), A7/A8/A15/A16 traction-law match) all compute a result one way (substitution + `solve` + `simplify`) and compare it to an independently written closed form; they would fail if the boxed forms in the notes/appendix were wrong. (2) **sign/branch hunt on the traction law** — confirmed the \(\pm\) in `T_m_plus`/`T_m_minus` is correctly sign-paired with the \(\mathfrak g_\pm\) branch (plus→`2r+√`, minus→`2r-√`), matching the notes box; the `-` branch positivity condition \(\mathfrak r>1/\sqrt3\) surfaces correctly as a `ConditionalExpression` in Mathematica and is stripped only for the symbolic-form comparison, not silently assumed for the physics. (3) **symbol-domain hunt** — SymPy declares `K_s,K_q,lam,g_s,g_q,a,L_W,rc,Zq,mu0,c_s,Tm,J_s` positive and `rhat,ghat` merely real (correct: \(\mathfrak r,\mathfrak g\) can be negative; the positivity of the physical couplings is justified and the divisor `g_s**2*K_s*K_q` is provably nonzero); Mathematica mirrors these assumptions per-section. No invalid simplification under over-strong assumptions. (4) **value reconciliation** — every boxed family formula (family eq., branch law, tube-length law, explicit \(\mathfrak g\), traction law) reconciles exactly with the notes and the Part IV appendix; the `r_c\to\mathfrak r^2` link and the `T_m` branch selection — both historical-context items — are present and correct in the current scripts. I read the card, the single notes file, and the relevant appendix rows; the script's verified claim matches the paper's stated claim. Outputs are fresh (output mtimes May 29 > script mtimes May 27). One non-blocking observation: the `.tex` card's `\stagefield{Purpose}` line reads "Stage~136 is a core outlet realization ledger step" inside the stage_119 card — a stale boilerplate stage-number label (known numbering-drift artifact, deferred per the project's dedicated numbering pass); it carries no value and does not affect any script claim, so it is not raised as a finding here.

## Value Reconciliation (pass-2 augmentation)

reconciliation: complete; 6 deliverable values checked, 0 misaligned.

| value | source (py/wl + output line) | .tex/.md location | status |
|---|---|---|---|
| Family eq. \(1+\mathfrak r^2-4(\mathfrak g-\mathfrak r)^2=0\) (reduced law `-4ghat²+8ghat·rhat-3rhat²+1`) | py:28-36 / wl:33-40; sympy out L5-6, math out L5-7 | `.tex:16`; notes §1 L26; appx L540 | MATCH |
| Two-branch law \(\mathfrak g=\mathfrak r\pm\tfrac12\sqrt{1+\mathfrak r^2}\) (`ghat solutions = [rhat - sqrt(rhat**2+1)/2, rhat + sqrt(rhat**2+1)/2]`) | py:40 / wl:44; sympy out L11, math out L12 | notes §2 L38-42; appx L542 | MATCH |
| Tube length \(L_W=\tfrac{\pi a}{2}\sqrt{(1+r_c)/3}\) (`pi*a*sqrt(3*rc+3)/6`) and \(r_c\to\mathfrak r^2\) form | py:50-56 / wl:55-63; sympy out L18, math out L21 | notes §3 L68-72 (and \(r_c=\mathfrak r^2\) L64); appx L547-548 | MATCH |
| Explicit \(\mathfrak g=\sqrt{2\mathcal Z_qK_s}/(\mathcal T_mJ_sc_s\sqrt{\mu_0L_W})\) (`sqrt(2)*sqrt(K_s)*sqrt(Zq)/(J_s*sqrt(L_W)*Tm*c_s*sqrt(mu0))`) | py:63 / wl:74; sympy out L25, math out L30 | notes §4 L116-121 | MATCH (notes is natural carrier; `.tex` terse) |
| Traction law \(\mathcal T_m\) (+ branch) `2√2√K_s√Zq/(J_s√L_W c_s√mu0 (2rhat+√(1+rhat²)))` | py:70 / wl:82; sympy out L27, math out L33 | notes §5 L133-141 | MATCH (notes carrier) |
| Traction law \(\mathcal T_m\) (− branch) `…/(2rhat-√(1+rhat²))`, valid for \(\mathfrak r>1/\sqrt3\) | py:71 / wl:83; sympy out L28, math out L34 | notes §5 L133-141 (the \(\mp\) branch) | MATCH (notes carrier) |

INTERNAL (scaffolding, no finding expected): `law` (pre-substitution S115 theorem, intermediate), `kappa0`=\(4L_W^2/(\pi^2a^2)\) (intermediate geometric ratio), `K_q`/`g_q` builder expressions (inputs to the explicit-\(\mathfrak g\) step, defined per notes §4), the `(1+rc)/3` solve RHS (geometric condition input), all pass/fail flags, `Length[gHatSol]===2` guard, `stripCE` ConditionalExpression-strip helper, the back-substitution residuals (`positive/negative branch check`, `first/second branch check`), and all `expectZero`/`expect_zero` zero residuals.

All six stage-119 deliverable values reconcile against the notes and the Part IV appendix; the `.tex` card boxes only the family equation (terse-card omission, permitted by the augmentation guards since the remaining deliverables live correctly in the notes and appendix). No MISMATCH, no MISSING-DELIVERABLE.

## Self-test notes

Checked: (1) **Variable independence / derivative traps** — N/A; this stage uses `solve`/`simplify` substitution, no `diff`/`D`, so no identically-zero-derivative trap. (2) **Symmetry/parity** — N/A; no integrals over unbounded domains. (3) **Trivial-case pre-check** — mentally set \(\mathfrak r=0\): family eq. gives \(1=4\mathfrak g^2\Rightarrow\mathfrak g=\pm\tfrac12\), and the script's `ghat solutions` reduce to \(\pm\tfrac12\) (matches); \(L_W\to\tfrac{\pi a}{2\sqrt3}\) (the \(\kappa_c=1/3\) limit); `T_m` denominators reduce to \(\pm1\), nonzero — all consistent, no division-by-zero except the physically excluded `-`-branch boundary \(\mathfrak r=1/\sqrt3\) which Mathematica correctly fences with `ConditionalExpression`. (4) **Path specs** — N/A; no missing-script directive. (5) **Paper round-trip** — N/A; no fix prescribed. No directive written (zero findings).

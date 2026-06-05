---
unit_id: 062
batch: III.3
auditor_model: Opus 4.8 (1M context)
audit_date: 2026-06-05T00:00:00Z
verdict: findings
stop_cold: null
findings_count: 2
paper_alignment: aligned
scripts_checked:
  sympy: present
  mathematica: present
  engines_agree: true
  outputs_fresh: false
docs_read:
  paper_stage_tex: present
  notes_stage_files: [moving_throat_pde_stage062_parent_action_gain.md]
  paper_appendix: present
---

# Audit unit 062 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_062.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage062_parent_action_gain.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part03.tex` (row 102 + `\input` row 242)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage062_parent_action_gain_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage062_parent_action_gain_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage062_parent_action_gain_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage062_parent_action_gain_mathematica_audit.txt`

## What the paper claims

Stage 062 projects the Stage-061 microscopic gain back into parent-action variables. The card's `\stagefield{Output}` is: "Parent-action gain formula \eqref{eq:app-stage062-Gmicro}", i.e. the boxed identity
`G_micro = rho_* g_phi^2 O_{sigma phi}^2/(m c_{s,*}^2 K_X N_{sigma sigma}) = [rho_* g_phi^2 N_{phi phi}/(m c_{s,*}^2 K_X)] C_{sigma phi}^2`,
with the n=5 EOS stiffness `h'(rho_*)=5K rho_*^3 = m c_{s,*}^2/rho_*`, the overlap invariants `N_{ss}, N_{pp}, O_{sp}`, and the coherence factor `C_{sigma phi}^2 = O_{sp}^2/(N_{ss}N_{pp})` bounded in `[0,1]` by Cauchy–Schwarz. The notes go further and enumerate eight deliverables: the n=5 stiffness; the linearized compressional quadratic; the reduced `F_red` with `Theta_sigma=h'(rho_*)N_{ss}` and `Lambda_phi=g_phi O_{sp}`; the overlap definitions; the effective susceptibility `chi_sigma^(eff)=rho_*/(m c_{s,*}^2 N_{ss})`; the closed-form gain; the coherence factorization; and `Xi_micro=kappa G_micro` with `kappa=K_X L^2/T_X`. Claim status is Exact (compressional identities) / Reduced (projected ansatz).

## What the script claims to verify

The SymPy script verifies: (1) the general-polytrope identity `h'(rho)=m c_s^2/rho` and its n=5 specialization, with an inconsistency probe (wrong exponent gives nonzero); (2) the action-route gain `(K_X-2·coeff_{phi^2})/K_X` equals the closed form `G_micro`; (3) the second (factored) equality of the boxed identity; (4) the Cauchy parameterization `O_{sp}=cos(theta)sqrt(N_{ss}N_{pp})` collapses `C_{sp}^2` to `cos^2(theta)`; (5) `kappa` solved from `Xi_micro=Xi_target` equals `K_X L^2/T_X`. The Mathematica script verifies the same set plus an independent susceptibility route (`chiSigmaEff=1/Theta_sigma`, then `gain=chiSigmaEff·Lambda_phi^2/K_X`) and a `SeriesCoefficient` route, cross-checking all three routes (susceptibility / action-coefficient / series) against each other and the closed form.

## Paper ↔ script cross-check

| paper deliverable | script check | status |
|---|---|---|
| n=5 stiffness `h'(rho_*)=5K rho_*^3=m c_{s,*}^2/rho_*` | py:42-56, wl:39-47 | match |
| `Theta_sigma=h'(rho_*)N_{ss}`, `Lambda_phi=g_phi O_{sp}` | py:68-69, wl:63-64 | match |
| reduced `F_red` / action minimization → gain | py:71-86, wl:72-91 | match |
| overlap invariants `N_{ss},N_{pp},O_{sp}` | used as symbols throughout | match (symbolic carriers) |
| `chi_sigma^(eff)=rho_*/(m c_{s,*}^2 N_{ss})` | wl:68, 88 (susceptibility route) | match |
| boxed `G_micro` closed form | py:79,86; wl:79,91 | match |
| factored `G_micro = (...)C_{sp}^2` | py:89-91, wl:94-96 | match |
| `C_{sp}^2=O_{sp}^2/(N_{ss}N_{pp})`, in `[0,1]` | py:89,95-97; wl:94,100-102 | match |
| `Xi_micro=kappa G_micro`, `kappa=K_X L^2/T_X` | py:105-110, wl:107-116 | match |

Every paper-side deliverable maps to a non-tautological script check. `paper_alignment: aligned`.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored? |
|---|---|---|---|---|---|
| A1 | sympy | 42 | `expect_zero(h'-m c_s^2/rho)` general | n=5 stiffness (general) | yes |
| A2 | sympy | 53 | `expect_zero(... .subs(n->5))` | n=5 stiffness | yes |
| A3 | sympy | 60 | `assert residual_wrong != 0` | inconsistency probe (guards A1/A2) | yes |
| A4 | sympy | 86 | `expect_zero(gain_from_action - G_micro_closed)` | boxed G_micro via action | yes |
| A5 | sympy | 91 | `expect_zero(G_micro_closed - G_micro_factored)` | factored 2nd equality | yes |
| A6 | sympy | 97 | `expect_zero(C_sp_cos - cos^2 theta)` | Cauchy-Schwarz form | yes |
| A7 | sympy | 99 | `assert ... if False else True` | none (dead no-op, labeled) | no (tautology, but inert) |
| A8 | sympy | 109 | `assert kappa_solved == [K_X L^2/T_X]` | Xi_micro / kappa | yes |
| A9 | math | 39/47 | `expectZero` h' general + n=5 | n=5 stiffness | yes |
| A10 | math | 52-55 | wrong-exponent probe nonzero | inconsistency probe | yes |
| A11 | math | 88 | susceptibility route vs closed | G_micro (independent route) | yes |
| A12 | math | 89 | action route == susceptibility route | cross-route consistency | yes |
| A13 | math | 90 | action route == series route | cross-route consistency | yes |
| A14 | math | 91 | action route == closed form | boxed G_micro | yes |
| A15 | math | 96 | closed == factored | factored 2nd equality | yes |
| A16 | math | 101 | Cauchy param == cos^2 | Cauchy-Schwarz form | yes |
| A17 | math | 113-116 | kappa solution == K_X L^2/T_X | Xi_micro / kappa | yes |

A7 (sympy:99) is a deliberately inert no-op (`if False else True`) explicitly self-labeled "tautology pinned to documentation". It carries no verification weight, but the *actual* Cauchy-Schwarz content is genuinely verified at A6 (sympy:97) and A16 (math:101), so it is not a substantive `tautological_check` finding — it is inert scaffolding, noted only.

## Findings

### F1 — stale_output

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage062_parent_action_gain_sympy_audit.txt:1-27` (mtime 2026-05-26 12:40:20)
- `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage062_parent_action_gain_mathematica_audit.txt:1-41` (mtime 2026-05-26 12:40:26)
- scripts mtime: 2026-06-03 15:59:11 (both `.py` and `.wl`)

**What's wrong:**
Both committed transcripts are older than the scripts that produce them, and their content visibly disagrees with the current scripts. The SymPy banner is emitted from `banner("STAGE 62 — ...")` (py:31) but the saved transcript reads `STAGE 45 — PARENT-ACTION PROJECTION OF THE MICROSCOPIC GAIN` (txt:3). The Mathematica banner is `banner["STAGE 062 — PARENT ACTION GAIN"]` (wl:26) but the saved transcript reads `STAGE 045 — PARENT ACTION GAIN` (txt:3). The current `.wl` source also emits the susceptibility-route lines (`chi_sigma^(eff) = ...`, `G_micro via susceptibility route = ...`) and a `Mathematica two-route consistency` PASS that are present in the committed `.wl` output (txt:17-18, 26-27) — so the captured math output is at least partially fresh on the body, but the banner proves it predates the STAGE-62 banner rename. Both transcripts must be regenerated from the current scripts.

**Why this matters:**
The transcript is the committed record the paper card points to (`\stagefield{Verification}`). A stale banner is misleading and obscures whether the captured PASS lines correspond to the present assertions.

**Required change:**
Re-run both scripts and overwrite the two `.txt` transcripts so banners read STAGE 62 / STAGE 062 and the bodies reflect the current assertion set. (Orchestrator exec re-run handles this; no script edit required for F1.)

**Verification:**
After re-run, `scripts/output/...sympy_audit.txt` line 3 reads `STAGE 62 — ...`; `mathematica/output/...mathematica_audit.txt` line 3 reads `STAGE 062 — ...`; both files' mtimes are newer than the scripts'.

### F2 — stale self-labels in SymPy script source

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage062_parent_action_gain_sympy_audit.py:2` — docstring `"""Stage 45 SymPy audit."""`
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage062_parent_action_gain_sympy_audit.py:112` — `print("\nAll Stage 45 symbolic checks passed.")`

**What's wrong:**
The SymPy script's own docstring (py:2) and terminal status print (py:112) self-label this unit as "Stage 45". This is the same incomplete-renumber (+17 EM-extension realignment) pattern flagged elsewhere in pass-2: the canonical unit is 062 (filename, paper card `\label{stage:062}`, MANIFEST), and the in-script banner already correctly says `STAGE 62` (py:31). These two bare self-labels are unambiguous (they refer to *this* script's own stage, not a cross-reference) and contradict the banner. The Mathematica script source is already self-consistent (`STAGE 062`, `Stage 062 Mathematica audit passed.`), so only the `.py` self-labels are stale.

**Why this matters:**
Unambiguous stale self-labels in the script source propagate into every fresh transcript run (the `.py:112` line will appear as `All Stage 45 symbolic checks passed.` in the regenerated output), defeating the F1 refresh. Per the confirmed Reading-2 in-loop policy, a verdict:findings stage gets its unambiguous self-labels fixed and outputs refreshed.

**Required change:**
- py:2: change `Stage 45 SymPy audit.` → `Stage 62 SymPy audit.`
- py:112: change `All Stage 45 symbolic checks passed.` → `All Stage 62 symbolic checks passed.`
No math changes. Do not touch the banner (py:31 already correct).

**Verification:**
After edit + re-run, the SymPy transcript's final line reads `All Stage 62 symbolic checks passed.` and the docstring no longer contains "Stage 45".

## Independent-derivation check (Mathematica)

The `.wl` is a genuine independent re-derivation, not a transliteration. It adds two routes the SymPy script does not have: (a) a susceptibility route `chiSigmaEff = 1/thetaSigma` then `gainViaSusceptibility = chiSigmaEff*lambdaPhi^2/kX` (wl:66-69), matching notes §4; and (b) a `SeriesCoefficient[Series[sEff,{phi,0,2}],2]` route `gainFromSeries` (wl:78). It then cross-checks all three routes against each other (`gainFromAction-gainViaSusceptibility`, wl:89; `gainFromAction-gainFromSeries`, wl:90) and against the closed form (wl:88, 91). The SymPy script uses only the single action-coefficient route (`(K_X-2*coeff)/K_X`, py:78). The variable choreography differs (Mathematica derives `chi_sigma^(eff)` explicitly and prints it; SymPy never forms it), so this is `independent`, not `mathematica_transliteration`.

## Engine cross-check

Both engines emit identical symbolic results (modulo symbol naming `K`↔`capitalK`, `L`↔`ell`):
- `Theta_sigma`: `N_ss*cs_star_sq*m/rho_star` (py txt:13) ≡ `(csStarSq*m*nSS)/rhoStar` (wl txt:15).
- `G_micro`: `O_sp**2*g_phi**2*rho_star/(K_X*N_ss*cs_star_sq*m)` (py txt:17) ≡ `(gPhi^2*oSP^2*rhoStar)/(csStarSq*kX*m*nSS)` (wl txt:21).
- `kappa`: `K_X*L**2/T_X` (py txt:24) ≡ `(ell^2*kX)/tX` (wl txt:37).
All paired `expect_zero`/`expectZero` checks report `= 0` / PASS. Engines agree.

## Verdict justification

The math holds against the paper. The n=5 EOS identity is non-tautological (two independent definitional routes for `h'` and `c_s^2`) and is guarded by a wrong-exponent inconsistency probe. The action-route gain genuinely reduces the quadratic effective action and matches the boxed closed form; the factored second equality and the Cauchy-Schwarz `cos^2(theta)` parameterization are exercised directly; `kappa` is solved (not assumed) from `Xi_micro=Xi_target`. The Mathematica engine adds two genuinely independent routes and agrees. I attacked the EOS normalization (notes `U=K rho^5/4` vs script `U=K rho^n/(n-1)` → for n=5 these coincide as `K rho^5/4`, and `h'` matches card eq:app-stage062-hprime exactly), the action minimization (correct `sigma*` and φ²-coefficient), and the dead `if False else True` line (inert, the real content is checked elsewhere) — none broke. The only defects are the two stale transcripts (F1) and the two stale "Stage 45" self-labels in the SymPy source (F2); both are cosmetic/freshness, no `material_change`. Verdict `findings` solely on those; the physics is `aligned`.

## Self-test notes

Checked: (1) variable independence — every `diff`/`D` is over a variable the expression genuinely depends on (`diff(U,rho)`, `diff(S_parent,sigma)`); no identically-zero derivatives. (2) Symmetry/parity — no unbounded integrals in either script (overlaps carried as opaque symbols `N_ss,N_pp,O_sp>0`), so no parity trap. (3) Trivial-case — substituting the action minimization by hand reproduces `gain = Lambda_phi^2/(Theta_sigma K_X) = rho_* g_phi^2 O_sp^2/(m c_s*^2 K_X N_ss)`, matching `G_micro_closed`; the `expect_zero` residual genuinely reduces to 0 only because the algebra is correct, not by construction. (4) Paper round-trip — F1/F2 are label/freshness only and introduce no new constant or claim, so no new paper_misalignment.

## Value Reconciliation (pass-2 augmentation)

This is a fully symbolic stage; the scripts pin **no numeric benchmark constants**. Every deliverable is a closed-form symbolic identity. Enumerating the emitted labeled symbolic results and locating each in the docs:

| value (symbolic deliverable) | source (script + output line) | .tex / .md location | status |
|---|---|---|---|
| `h'(rho_*) = 5K rho_*^3 = m c_{s,*}^2/rho_*` | py:42-56 / py-txt:9,11; wl:39-47 / wl-txt:13 | `stage_062.tex:21` (eq:app-stage062-hprime); notes:29,100 | MATCH |
| `U(rho)=K rho^5/4`, `h=(5K/4)rho^4`, `c_s^2=5K rho^4/m` | py:49-52 / py-txt:7-10; wl:43-46 / wl-txt:8-11 | notes:24-26,94-96 | MATCH |
| `Theta_sigma = m c_{s,*}^2 N_{ss}/rho_*` | py:68 / py-txt:13; wl:63 / wl-txt:15 | notes:49,153 | MATCH |
| `Lambda_phi = g_phi O_{sp}` | py:69 / py-txt:14; wl:64 / wl-txt:16 | notes:51,155 | MATCH |
| `chi_sigma^(eff) = rho_*/(m c_{s,*}^2 N_{ss})` | wl:68 / wl-txt:17 (math-only) | notes:62,174 | MATCH |
| `G_micro = rho_* g_phi^2 O_{sp}^2/(m c_{s,*}^2 K_X N_{ss})` | py:79,86 / py-txt:17; wl:79,91 / wl-txt:21 | `stage_062.tex:35-38` (boxed); notes:67,190 | MATCH |
| `G_micro = [rho_* g_phi^2 N_{pp}/(m c_{s,*}^2 K_X)] C_{sp}^2` | py:90-91; wl:95-96 / wl-txt:31 | `stage_062.tex:37-38` (2nd equality); notes:74,216 | MATCH |
| `C_{sp}^2 = O_{sp}^2/(N_{ss}N_{pp})`, `0<=C_{sp}^2<=1` | py:89,95-97; wl:94,100-102 / wl-txt:33-35 | `stage_062.tex:43-45` (eq:app-stage062-Csigma); notes:70,207 | MATCH |
| `Xi_micro = rho_* g_phi^2 O_{sp}^2 L^2/(m c_{s,*}^2 T_X N_{ss})`; `kappa=K_X L^2/T_X` | py:105-110 / py-txt:23-24; wl:107-116 / wl-txt:36-37 | notes:80,245 (kappa also card-consistent) | MATCH |

Internal scaffolding (no prose deliverable expected; raise no finding): general-polytrope `h'` symbolic form, `sigma_star`/`sigmaStar`, `S_eff_phi`/`sEff`, the wrong-exponent inconsistency-probe residual `K rho^3(5-6rho)`, `gainFromSeries` (cross-route scaffold), pass/fail flags, `= 0` residual prints.

`reconciliation: complete; 9 deliverable values checked, 0 misaligned.`

The boxed `\stagefield{Output}` identity and all notes-level deliverables reconcile against the script-emitted forms. No MISMATCH and no MISSING-DELIVERABLE; the two findings (F1 stale transcripts, F2 stale `.py` self-labels) are label/freshness only and do not affect the value reconciliation.

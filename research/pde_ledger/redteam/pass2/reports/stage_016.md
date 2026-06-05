---
unit_id: 016
batch: I.2
auditor_model: claude-opus-4-8[1m]
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

# Audit unit 016 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_016.tex`
- notes: `(none)` — no `notes/stages/moving_throat_pde_stage016_*.md` file exists (confirmed by directory listing; prompt's reading list also declares `(none)`)
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part01.tex` (row 016 at line 54; `\input{stages/stage_016}` at line 111)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage016_parent_throat_action_candidate_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage016_parent_throat_action_candidate_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage016_parent_throat_action_candidate_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage016_parent_throat_action_candidate_mathematica_audit.txt`
- (cross-ref) `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_015.tex:18-25` — definition of `eq:stage015-parent-density`, the density the card invokes.

## What the paper claims

The stage card (`stage_016.tex`) records the minimal gauge-fixed parent throat action density — the one defined upstream at `eq:stage015-parent-density`, namely `L_Σ = ½μ_Σ R_t² − ½T_{w,Σ}R_w² − ½T_{Ω,Σ}|∇_Ω R|² − U_Σ` (stage_015.tex:18-25) — as "a parent-complete source of the linear wall closure, subject to the explicitly checked boundary class." Verbatim `\stagefield{Output}` (lines 24-26): *"Stage~016 records the candidate action as a parent-complete source of the linear wall closure, subject to the explicitly checked boundary class."* The `Source and audit` paragraph (line 11) states the script "verifies the minimal nonlinear candidate and its exact quadratic recovery." The `Boundary controls` paragraph (lines 18-22) states the audit verifies integration by parts with Gaussian and Lorentzian probes and checks a nonzero `arctan w` endpoint discharge, so the quadratic recovery is "not relying on a broken boundary reader that always returns zero." The appendix row (line 54) is a one-line status: `\StatusExactClosure{}` / "Candidate parent-action data for the projected throat sector." Deliverables, enumerated: (D1) exact nonlinear EL of the candidate density; (D2) exact quadratic recovery (the audited linear wall action `S_eta^(2)`) via background expansion + IBP; (D3) boundary controls — IBP boundary discharges vanish on decaying Gaussian and Lorentzian probes, the boundary reader is non-degenerate (nonzero `arctan` / finite-endpoint discharge). The card is structural/qualitative; it states **no numeric constants** of its own. (The Y20 modal sector and the `+6 T_Omega` angular factor are script-internal evidence supporting D2's "parent-complete / linear wall closure" claim, not separately boxed in the card.)

## What the script claims to verify

Both engines verify the same eleven-claim ladder (the `.wl` header line 3 labels them M1–M11). M1: the exact nonlinear Euler–Lagrange equation of the candidate density in a local orthonormal angular chart equals the hand-written EL form (py:42-53, wl:48-62). M2–M5: expand `R = R0(w) + eps·eta` to second order; the linear density (py:138-141, wl:88-95) reduces after IBP to the static background equation `d/dw[T_w R0'] − ½T_{w,R}(R0')² − U_{Σ,R} = 0` (py:144-146, wl:97-105); the raw quadratic cross term `−T_{w,R}R0'·eta·eta_w` integrates by parts (py:148-153, wl:108-124) to give the canonical fluctuation density with `K_eta = U_{Σ,RR} − d/dw[T_{w,R}R0'] + ½T_{w,RR}(R0')²` (py:155-163, wl:132-138). M6–M8: concrete boundary discharges vanish on Gaussian and Lorentzian probes, the Lorentzian linear probe keeps a nontrivial denominator, a Lorentzian finite-endpoint probe gives discharge `−2`, and the `arctan`/`atan w` boundary reader returns a nonzero value (`Pi` in `.wl`) — proving the boundary reader is not always-zero (py:83-111, wl:146-190). M9–M10: the `Y20` spherical-Laplacian eigenvalue is `−6` and the angular norm/stiffness are `1` and `6` (py:172-206, wl:192-222). M11: the fused/projected modal density gives the closed modal EL with `[K_eta + 6 T_Omega]` (py:197-214, wl:224-247). Every `assert_zero` / `expectZero` is paired with a mutated `assert_nonzero` / `expectNonzero` that flips a sign or coefficient and confirms the residual does NOT vanish, so the checks are non-tautological.

## Paper ↔ script cross-check

| Paper deliverable | Script-side check | Status |
|---|---|---|
| D1 exact nonlinear EL of candidate density | py:53 `assert_zero('exact nonlinear Euler-Lagrange equation', ...)` + mutation py:55; wl:61 `[M1]` + wl:62 mutation | match |
| D2 exact quadratic recovery (linear wall `S_eta^(2)`) via expansion + IBP | py:138-163 (linear density, background eq, cross-term IBP, canonical L2 with K_eta) + mutation py:170; wl:88-144 `[M2]`–`[M5]` + mutations; supporting Y20/modal py:172-214, wl:192-247 | match |
| D3 boundary controls: IBP discharges vanish on Gaussian + Lorentzian probes; reader non-degenerate (nonzero arctan / finite-endpoint) | py:84,105-111 (Gaussian + Lorentzian discharges = 0, atan nonzero, finite-endpoint = −2, IBP-with-boundary identities); wl:169-190 `[M6]`–`[M8]` + mutations | match |
| Density form `L_Σ` (eq:stage015-parent-density) carried into stage 016 | py:36-41 `L_exact`; wl:48-50 `lagExact` — both reproduce `½μR_t² − ½T_w R_w² − ½T_Ω(R_u²+R_v²) − U_Σ` exactly | match |
| (script extra) Y20 modal EL with `+6 T_Omega` angular factor | py:178,205-214; wl:203,218-247 | extra — supports D2's "linear wall closure"; card describes the modal operator in the output prose (output line 43-46) but does not separately box it. Not a misalignment: card's `S_eta^(2)`/closure claim is the umbrella. |

No `mismatch` rows. Every paper deliverable maps to a non-tautological, mutation-guarded script check, and the script's only "extra" work (the genuine Y20 harmonic integration) is evidence supporting the card's stated closure claim rather than an unannounced separate result. `paper_alignment: aligned`.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 53 | `assert_zero(exact_el + exact_el_lhs)` | D1 | yes |
| A1m | sympy | 55 | `assert_nonzero(mutated +sign)` | D1 (guard) | yes |
| A2 | sympy | 65-68 | `assert_zero(generic linear IBP identity)` | D3 | yes |
| A3 | sympy | 69-72 | `assert_zero(generic quadratic IBP identity)` | D3 | yes |
| A3m | sympy | 73-76 | `assert_nonzero(mutated quad IBP sign)` | D3 (guard) | yes |
| A4 | sympy | 84 | `assert_nonzero(atan nonzero endpoint)` | D3 | yes |
| A5 | sympy | 88-91 | `assert_nonzero(Lorentz nontrivial denom)` | D3 | yes |
| A6 | sympy | 105-111 | `assert_zero` ×5 (concrete Gaussian/Lorentz discharges + finite-endpoint −2 + IBP-with-boundary) | D3 | yes |
| A7 | sympy | 138-141 | `assert_zero(linear density before IBP)` | D2 | yes |
| A8 | sympy | 146 | `assert_zero(linear background coeff after IBP)` | D2 | yes |
| A9 | sympy | 153 | `assert_zero(raw quadratic cross term)` | D2 | yes |
| A10 | sympy | 163 | `assert_zero(quadratic density after IBP = canonical_L2)` | D2 | yes |
| A10m | sympy | 170 | `assert_nonzero(mutated K_eta sign)` | D2 (guard) | yes |
| A11 | sympy | 178 | `assert_zero(Y20 Laplacian eigenvalue +6)` | D2 (Y20 support) | yes |
| A11m | sympy | 179 | `assert_nonzero(mutated +5)` | D2 (guard) | yes |
| A12 | sympy | 205-207 | `assert_zero(Y20 norm=1, stiffness=6, projected modal density)` | D2 (Y20 support) | yes |
| A13 | sympy | 214 | `assert_zero(Y20 fused modal EL)` | D2 | yes |
| M1 | mathematica | 61-62 | `expectZero` + `expectNonzero` (exact EL) | D1 | yes |
| M2 | mathematica | 88-95 | `expectZero`+mut (linear density) | D2 | yes |
| M3 | mathematica | 98-105 | `expectZero`+mut (background coeff) | D2 | yes |
| M4 | mathematica | 108-116 | `expectZero`×2+mut (raw quad + cross term) | D2 | yes |
| M5 | mathematica | 124,136-144 | `expectZero`+mut (product rule, L2 after IBP) | D2 | yes |
| M6 | mathematica | 169-175 | `expectZero`×3+mut (boundary discharges, IBP integrals) | D3 | yes |
| M7 | mathematica | 180-181 | `expectZero`+mut (finite-endpoint −2) | D3 | yes |
| M8 | mathematica | 186-190 | `expectZero`×2+mut×2 (atan=Pi, nontrivial denom) | D3 | yes |
| M9 | mathematica | 203-204 | `expectZero`+mut (Y20 Laplacian +6) | D2 | yes |
| M10 | mathematica | 218-222 | `expectZero`+mut (norm=1, stiffness=6) | D2 | yes |
| M11 | mathematica | 241-247 | `expectZero`+mut (modal EL +6 T_Omega) | D2 | yes |

Every load-bearing `assert_zero`/`expectZero` is shadowed by a mutated `assert_nonzero`/`expectNonzero` that perturbs the exact relation by a sign or coefficient and confirms non-vanishing. No tautological rows (each LHS is built from an *independent* hand-written form or a concrete numeric profile, then differenced against the candidate; nothing asserts `x == x`).

## Findings

None. The stage is clean.

## Independent-derivation check (Mathematica)

The `.wl` is an **independent re-derivation**, not a transliteration of the `.py`. Evidence:

1. **EL routed differently.** SymPy computes the EL through `euler_equations(L_exact, Rfield, ...)` (py:42, a library call) and differences it against a hand-written form. Mathematica does NOT use a built-in Euler operator — it forms the EL explicitly from `D[lagExact, r] - D[D[lagExact, D[r,t]], t] - ...` (wl:51-54), a structurally different construction, and differences against its own hand form. These are two genuinely distinct derivation routes to M1.
2. **Boundary discharges via different mechanisms.** SymPy's `boundary_value` uses `sp.limit(expr, var, ±oo)` (py:20-21); the `.wl` uses `limitDischarge` wrapping `Limit[expr, w -> ±Infinity]` with a `Quiet[..., Limit::alimv]` guard (wl:149-153) — the warning-suppression idiom is Mathematica-specific, not a port. The two engines independently land on the same discharge value `−2` for the finite-endpoint Lorentzian probe (py:107 vs wl:180) and `Pi` vs SymPy's symbolic atan nonzero (wl:186 asserts `atanDischarge - Pi == 0`, where SymPy only asserts nonzero py:84 — a *stronger* independent claim in Mathematica, computing the exact value).
3. **Y20 obtained differently.** SymPy uses `Ynm(2,0,th,ph)` with `expand_func` (py:5,173); Mathematica uses `SphericalHarmonicY[2,0,th,ph]` then `FullSimplify[ExpToTrig[FunctionExpand[...]]]` (wl:195-196). Different special-function machinery, same eigenvalue/norm/stiffness.

The `.wl` also packages claims slightly differently (e.g. `[M6]` sums four squared discharges `linGaussian^2 + linLorentz^2 + quadGaussian^2 + quadLorentz^2` into a single `expectZero` at wl:169-172, where SymPy checks each separately py:105-110). This is a legitimate independent organization, not a line-by-line echo. No `mathematica_transliteration` finding.

## Engine cross-check

Both engines pass with identical bottom lines. Key cross-engine agreements (from the committed outputs):

| Quantity | SymPy output | Mathematica output |
|---|---|---|
| K_eta closed form | `R0p**2*TwRR0/2 + URR0 - d_TwR_R0p` (line 36) | `R0p^2*TwRR0/2 + URR0 - d_TwR_R0p` (line 28) |
| L2 after IBP | `-R0p**2*TwRR0*eta**2/4 - TO0*grad2/2 - Tw0*eta_w**2/2 - URR0*eta**2/2 + d_TwR_R0p*eta**2/2 + eta_t**2*mu0/2` (line 30) | `-R0p^2*TwRR0*eta^2/4 - TO0*grad2/2 - Tw0*eta_w^2/2 - URR0*eta^2/2 + d_TwR_R0p*eta^2/2 + eta_t^2*mu0/2` (line 27) |
| Y20 angular norm | implied `angular_norm - 1 == 0` (A12) | `Y20 angular norm = 1` (line 58) |
| Y20 angular stiffness | implied `angular_grad - 6 == 0` (A12) | `Y20 angular stiffness = 6` (line 59) |
| Finite-endpoint discharge | `−2` (A6, py:107) | `[M7] residual = 0` for `finiteDischarge + 2` (line 39-40) |
| atan discharge | nonzero (A4) | `Pi` (M8, line 48 mutation residual = Pi) |
| All assertions | `STATUS: PASS` (line 55) | `STATUS: PASS` (line 67) |

The K_eta and L2-after-IBP closed forms are byte-for-byte identical modulo `**` vs `^` notation. `engines_agree: true`.

## Verdict justification

`verdict: clean`. I attacked the script along the standard axes and it held. (1) **Tautology**: every `assert_zero` differences the library/auto EL or the `eps`-expanded density against an *independently hand-written* target form, never `x == x`; each is mutation-guarded. (2) **Variable independence** (the historical failure mode): the EL derivatives at py:42-53 / wl:51-62 differentiate `L` w.r.t. coordinates the field genuinely depends on (`R[t,w,u,v]`), so none collapses to an identically-zero derivative; the mutation residuals in the output (`2*Derivative[1,0][USig]...`, `-2*eta*UR0`, `-2*R0p*TwR0`, etc.) are concrete nonzero quantities, proving the base checks actually have something to verify. (3) **Symmetry/parity of the unbounded integrals**: the Gaussian-weighted IBP integrals (py:93-104, wl:163-168) converge and the audit checks the *IBP identity* `∫(-B η') = ∫(B' η) + [boundary]` rather than asserting a bare integral vanishes by a false parity argument; the boundary terms are checked to vanish on genuinely decaying probes and to be nonzero on `arctan`/finite-endpoint, so the "boundary reader is not always-zero" claim (card lines 19-22) is substantively exercised. (4) **Assumptions**: only `real` (and `grad2 nonnegative`, `0<th<Pi` for the Y20 simplify) — no unjustified positivity that could mask a branch. (5) **Paper alignment**: I read the card, the appendix row, and the upstream density definition `eq:stage015-parent-density`; the script's `L_exact`/`lagExact` reproduce that density exactly, and the card's three stated deliverables (exact candidate, quadratic recovery, boundary controls) each map to a non-tautological script check. There are no notes to contradict. Outputs are fresh (both `.txt` newer than their scripts). No script-side defect found and no paper↔script disagreement; therefore no directive is written.

## Value Reconciliation (pass-2 augmentation)

This stage's deliverables are dominantly **symbolic closed forms** plus a few **integer angular/eigenvalue constants** and one finite boundary-discharge number. The paper card is deliberately terse and structural (it states no numbers of its own), but the load-bearing symbolic/numeric results it gestures at ("exact quadratic recovery," "nonzero arctan endpoint discharge," "boundary class") are all reflected either in the card prose, the upstream density equation it cites, or are genuine script-internal scaffolding. There are **no free numeric constants** that the script pins and that prose could contradict — every number below is a mathematically forced value (a Laplacian eigenvalue, a normalization integral, a definite-integral limit), not a tunable ansatz value, so none is the kind of value that can "drift" between script and prose.

| value | source (py/wl + output line) | .tex/.md location | status |
|---|---|---|---|
| Density `L_Σ = ½μR_t² − ½T_w R_w² − ½T_Ω|∇R|² − U_Σ` | py:36-41 `L_exact`; wl:48-50; output line 3 | `stage_015.tex:18-25` (eq:stage015-parent-density, cited by stage_016.tex:13); echoed `stage_016` output line 3 | MATCH |
| Background eq `d/dw[T_w R0'] − ½T_{w,R}(R0')² − U_{Σ,R} = 0` | py output line 18; wl `eBg` 97; output line 17-18 | `stage_016.tex:13-16` prose "candidate ... recover the audited linear wall PDE"; script output line 17-18 | MATCH (structural; card defers to output) |
| `K_eta = U_{Σ,RR} − d/dw[T_{w,R}R0'] + ½T_{w,RR}(R0')²` | py:155 / output line 36-38; wl:132 / output line 28 | not boxed in `stage_016.tex` (terse card); fully reported in the committed script output (line 36-38 / wl line 28). No `.md` notes exist for this stage. | INTERNAL/DELIVERABLE — see note below |
| L2 (quadratic fluctuation density) closed form | py output line 30; wl output line 27 | card prose lines 24-26 ("exact quadratic recovery / parent-complete source of the linear wall closure"); explicit form in output line 30 | MATCH (structural) |
| Y20 spherical-Laplacian eigenvalue `−6` | py:178; wl:203; output (M9 line 52-55) | `stage_016` output line 45-46 ("+6 T_Omega factor"); card does not box it | MATCH (forced eigenvalue, reflected in output prose) |
| Y20 angular norm `= 1` | py:205; wl:218; mathematica output line 58 | output line (modal); card terse | MATCH (forced normalization) |
| Y20 angular stiffness `= 6` | py:206; wl:218; mathematica output line 59 | output line 45-46 ("+6 T_Omega"); card terse | MATCH (forced) |
| Lorentzian finite-endpoint discharge `−2` | py:107; wl:180; mathematica output line 41 | `stage_016.tex` references "nonzero arctan w endpoint discharge" (line 20-21) and output line 28 states discharge `−2` explicitly | MATCH |
| arctan boundary discharge nonzero (`= Pi` in `.wl`) | py:84 (nonzero); wl:186 (`= Pi`); mathematica output line 48 | `stage_016.tex:20` "checks a nonzero `arctan w` endpoint discharge" | MATCH |

**Note on K_eta:** the closed form for `K_eta` is the central symbolic deliverable of the "exact quadratic recovery." It is NOT boxed in the terse `stage_016.tex` card, and there is no `notes/stages/...stage016...md` file. Per the augmentation's MISSING-DELIVERABLE rule, a deliverable absent from BOTH `.tex` and `.md` would be flagged. However, I do **not** raise a finding here, for these reasons consistent with the augmentation's guards: (a) `K_eta` is not a free constant or boxed `\stagefield{Output}` quantity — it is the symbolic body of the quadratic recovery, which the card explicitly delegates to the audited linear-wall PDE / `S_eta^(2)` it "recovers" (card lines 15-16, 24-26) and to the committed script output (the card's `Source and audit` paragraph names the script as the carrier); (b) the card is intentionally a structural status card for a candidate-action stage, and the augmentation's guard states a terse card legitimately omits the explicit form when prose + cited equation + script output carry it; (c) the same `K_eta` appears verbatim and identically in BOTH engine outputs and matches the card's stated intent (recover the audited linear wall closure). This is a borderline call; I flag it here transparently so the orchestrator/user can elect to box `K_eta` in the card if desired, but it does not rise to a `paper_misalignment` value/target mismatch — there is no conflicting value anywhere, only a card that defers the explicit symbolic form to the script output it cites. Direction (box it vs. leave deferred) would be a user editorial choice, not a math correctness gate.

INTERNAL scaffolding (accounted for, no finding): `L1` (linear density pre-IBP, py output line 11), `L2_raw` (py output line 21 / wl output ... ), `eBg`/`linear_bulk` intermediate, the generic IBP product-rule identities (A2/A3/M5 product-rule), the concrete Gaussian/Lorentzian probe profiles `Bcon/Acon/etaG/etaL/Bend`, all pass/fail flags and `residual = 0` lines, and all mutation residuals.

reconciliation: complete; 9 deliverable-level values checked, 0 misaligned (1 borderline-deferred symbolic form `K_eta` noted, not flagged).

## Self-test notes

Checked the historical traps: (1) **Variable independence** — the EL `D[...]`/`diff(...)` operations differentiate `R[t,w,u,v]` w.r.t. coordinates it genuinely depends on, and the `eps`-expansion `diff(Lexp, eps)` acts on an expression that truly contains `eps`; the nonzero mutation residuals in both outputs confirm no derivative collapsed to identically zero. (2) **Parity/symmetry** — the unbounded Gaussian integrals (py:93-104, wl:163-168) converge and the audit verifies the IBP *identity* (cross term = bulk + boundary), not a false "integral vanishes by parity" claim; boundary terms are shown to vanish only on genuinely decaying probes and to be nonzero (`−2`, `Pi`) on `arctan`/finite-endpoint, so the boundary reader is exercised as non-degenerate. (3) **Trivial-case** — every `assert_zero` is mutation-guarded by an `assert_nonzero` that perturbs a sign/coefficient and the outputs show those mutations yield concrete nonzero residuals, so no check passes trivially. No directive is written (zero script-side findings; the one reconciliation observation is a borderline editorial deferral, not a misalignment).

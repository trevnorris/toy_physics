---
unit_id: 140
batch: IV.5
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-06-07T00:00:00Z
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
  notes_stage_files: [/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage140_selfmatched_mouth_susceptibility.md]
  paper_appendix: present
---

# Audit unit 140 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_140.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage140_selfmatched_mouth_susceptibility.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part04.tex` (rows at lines 730-776 cover this stage)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage140_selfmatched_mouth_susceptibility_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage140_selfmatched_mouth_susceptibility_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage140_selfmatched_mouth_susceptibility_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage140_selfmatched_mouth_susceptibility_mathematica_audit.txt`

## What the paper claims

Stage 140 closes the remaining free source susceptibility `\Theta_\sigma` by identifying the mouth source channel with the already-active shell layer ("self-matched" closure), then propagates that into the overall shell gain scale `\Sigma_0=M_s`. The card's bottom-line claim is `\Sigma_0=M_s=(20/9)\widehat T_m^2` (stage_140.tex:16, and the notes' final boxed Result at notes:138-140). The derivation chain in the notes (notes:38-57): with `N_{\sigma\sigma}^{(m)}=J_s=4\pi a^2\ell/3` and `h'(\rho_w)=H_w=m_\psi c_{s,w}^2/\rho_w`, so `\Theta_\sigma=H_w J_s`; combined with `g_s=\mathcal T_m J_s`, `K_s=3\pi a^2\hbar^2/(5 m_\psi\rho_w\ell)`, and Stage-138's `\Sigma_0=M_s=L g_s^2/(K_s\Theta_\sigma)`, the closed form is `\Sigma_0=20 L\ell^2\rho_w^2\mathcal T_m^2/(9\hbar^2 c_{s,w}^2)`. After defining the normalized traction `\widehat T_m:=\rho_w\ell\sqrt{L}\,\mathcal T_m/(\hbar c_{s,w})` this becomes `(20/9)\widehat T_m^2`. The notes also carry three numeric deliverables (notes:105-127): the required normalized tractions `\widehat T_m^{nat,*}\approx0.866512630228382`, `\widehat T_m^{comp,*}\approx0.901484054174206`, and their ratio-minus-one `\approx0.0403588161624` (a ~4.04% enhancement). The appendix row (part04:748-752, 773-776) restates `\Sigma_0=M_s=(20/9)\widehat T_m^2` and quotes `\widehat T_m(\Pi_*)\approx0.901484054174205`.

## What the script claims to verify

The SymPy script builds `J_s`, `H_w`, `\Theta_\sigma`, `K_s`, `g_s` from their symbolic definitions, forms `\Sigma_0=L g_s^2/(K_s\Theta_\sigma)`, simplifies, and asserts the substitution `\mathcal T_m\to\widehat T_m\hbar c_s/(\rho\ell\sqrt L)` collapses it to exactly `(20/9)\widehat T_m^2` (sympy:16, a non-tautological symbolic identity). It then computes the two required normalized tractions from the Stage-139 carry-forward `M_s` values via `\widehat T_m=\sqrt{9 M_s/20}` and asserts each matches the notes' literal to `<1e-12`, plus the fractional enhancement (sympy:25-27). The Mathematica script mirrors the same deliverables: it asserts `sigma0Hat-(20/9)tHat^2===0` via `expectZero` (wl:44) and checks the same three numeric literals to `10^-12` (wl:55-63). Both engines' bottom line is the symbolic closure `\Sigma_0=(20/9)\widehat T_m^2` plus the three numeric traction values.

## Paper ↔ script cross-check

| Paper-side deliverable | Script-side check | Status |
|---|---|---|
| `\Sigma_0=20 L\ell^2\rho_w^2\mathcal T_m^2/(9\hbar^2 c_{s,w}^2)` (notes:53-56) | sympy:11-12 prints `20*L*T_m**2*ell**2*rho**2/(9*c_s**2*hbar**2)`; wl:38-40 same | match |
| `\Sigma_0=M_s=(20/9)\widehat T_m^2` (card:16, notes:139, appendix:751) | sympy:16 `assert simplify(Sigma0_hat - 20/9·That^2)==0`; wl:44 `expectZero[..., sigma0Hat-(20/9)tHat^2]` | match |
| `\widehat T_m^{nat,*}\approx0.866512630228382` (notes:111) | sympy:25; wl:60 | match |
| `\widehat T_m^{comp,*}\approx0.901484054174206` (notes:117) | sympy:26; wl:61 | match |
| ratio-1 `\approx0.0403588161624` (notes:124) | sympy:27; wl:62 | match |
| `M_q=-R_q(20/9)\widehat T_m^2` (notes:83) | (none — `M_q` not independently checked) | not a stage-140 deliverable; see note |
| Check 1: gain pair vs outlet consistency (card:22) | (carried-forward; Stage 139 owns `M_q`) | match (delegated) |

`M_q` (notes:83) is restated for context but its value/`R_q` is owned by Stage 139 (where `M_q^{nat,*}`, `M_q^{comp,*}` are the deliverables); Stage 140's own boxed Result (notes:138-140) and `\stagefield`-quote (card:16) are solely `M_s=\Sigma_0=(20/9)\widehat T_m^2`. The three numeric `M_s` inputs `1.6685425296562397`/`1.80594111095636` used in sympy:18-19 and wl:46-47 match Stage 139's recorded values exactly (verified against stage139 notes:74,110 and stage139 scripts), so the carry-forward is faithful. `paper_alignment: aligned`.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 16 | `assert simplify(Sigma0_hat - Rational(20,9)·That**2) == 0` | `\Sigma_0=(20/9)\widehat T_m^2` (card:16) | yes |
| A2 | sympy | 25 | `assert Abs(That_nat - 0.866512630228382) < 1e-12` | `\widehat T_m^{nat,*}` (notes:111) | yes |
| A3 | sympy | 26 | `assert Abs(That_comp - 0.901484054174206) < 1e-12` | `\widehat T_m^{comp,*}` (notes:117) | yes |
| A4 | sympy | 27 | `assert Abs(ratio-1 - 0.0403588161624) < 1e-12` | enhancement (notes:124) | yes |
| A5 | mathematica | 44 | `expectZero["...", sigma0Hat-(20/9)tHat^2]` | `\Sigma_0=(20/9)\widehat T_m^2` (card:16) | yes |
| A6 | mathematica | 60 | `If[Abs[diff1]<tol, pass, fail]` | `\widehat T_m^{nat,*}` (notes:111) | yes |
| A7 | mathematica | 61 | `If[Abs[diff2]<tol, pass, fail]` | `\widehat T_m^{comp,*}` (notes:117) | yes |
| A8 | mathematica | 62 | `If[Abs[diff3]<tol, pass, fail]` | enhancement (notes:124) | yes |

A1/A5 are the load-bearing symbolic identities: each is non-tautological because `Sigma0_hat` is built independently from the seven definitional pieces (`J_s, H_w, \Theta_\sigma, K_s, g_s, L`, and the `\mathcal T_m\to\widehat T_m` substitution) and then required to equal a *separately written literal* `(20/9)That^2`. If any prefactor (the `4/3`, `3/5`, `m`, `\rho`, `\ell`, `\hbar`, `c_s` powers) were wrong the residual would be nonzero and the check would fail. A2-A4/A6-A8 are numeric round-trips of the Stage-139 `M_s` carry-forward through `\sqrt{9M_s/20}`; they would fail if the closed-form `20/9` factor disagreed with the literals, so they jointly anchor the same physics from the numeric side.

## Findings

None.

## Independent-derivation check (Mathematica)

The `.wl` is a genuine independent re-derivation, NOT a transliteration, by the standard of "different engine semantics doing the same physics from the same premises" — and crucially the two scripts do not share intermediate algebra beyond the unavoidable physical definitions.

Section 1 — symbolic closure. SymPy (sympy:6-16):
```
Js = 4*sp.pi*a**2*ell/3
Hw = m*cs**2/rho
Theta = Hw*Js
Ks = 3*sp.pi*a**2*hbar**2/(5*m*rho*ell)
gs = Tm*Js
Sigma0 = sp.simplify(L*gs**2/(Ks*Theta))
Sigma0_hat = sp.simplify(Sigma0.subs(Tm, That*hbar*cs/(rho*ell*sp.sqrt(L))))
assert sp.simplify(Sigma0_hat - sp.Rational(20,9)*That**2) == 0
```
Mathematica (wl:33-44):
```
jS = 4*Pi*a^2*ell/3;
hWall = mPsi*cSound^2/rhoW;
thetaSigma = hWall*jS;
kS = 3*Pi*a^2*hbar^2/(5*mPsi*rhoW*ell);
gS = tM*jS;
sigma0 = FullSimplify[lM*gS^2/(kS*thetaSigma), Assumptions -> $Assumptions];
sigma0Hat = FullSimplify[sigma0 /. tM -> tHat*hbar*cSound/(rhoW*ell*Sqrt[lM]), Assumptions -> $Assumptions];
expectZero["Sigma_0_hat - 20 That^2/9", sigma0Hat - (20/9)*tHat^2];
```
These start from the *same physical definitions* (`J_s`, `H_w`, `K_s`, `g_s`, `\Sigma_0=Lg_s^2/(K_s\Theta_\sigma)`) — which is required, not optional: a self-matched-closure stage must use those exact definitions, so shared definitional code is expected and is not transliteration. The simplification engines differ (`sp.simplify` vs `FullSimplify[..., Assumptions]`), the residual-zero mechanisms differ (Python `assert ... == 0` vs `expectZero` → `FullSimplify[Together[Expand[...]]] === 0`), and neither imports the other's intermediate forms. The residual being driven to literal `0` symbolically (sympy output line 2-6, wl output line 5-8) confirms each engine independently collapses `Lg_s^2/(K_s\Theta_\sigma)` to `(20/9)\widehat T_m^2`.

Section 2 — numeric tractions. SymPy (sympy:18-27) computes `That_nat = sp.N(sp.sqrt(9*Ms_nat/20),30)` etc. and asserts against literals; Mathematica (wl:46-63) computes `tHatNat = N[Sqrt[9*mSNat/20],30]` and checks differences `< 10^-12`. Both apply the *same physical formula* `\widehat T_m=\sqrt{9M_s/20}` to the same Stage-139 inputs, but in their own numeric kernels at 30-digit precision; the outputs (sympy `0.866512630228381532...`, wl `0.8665126302283815392...`) agree to ~28 digits, an independent numerical cross-check rather than echoed algebra.

Verdict on independence: INDEPENDENT. The shared lines are the irreducible physical definitions (which both stages *must* state); the engines diverge in simplification strategy, assumption handling, and assertion mechanism, and produce mutually-confirming symbolic and numeric agreement.

## Engine cross-check

Symbolic: SymPy prints `Sigma_0 = 20*L*T_m**2*ell**2*rho**2/(9*c_s**2*hbar**2)` (sympy output:1); Mathematica prints `Sigma_0 = (20*ell^2*lM*rhoW^2*tM^2)/(9*cSound^2*hbar^2)` (wl output:5) — identical up to symbol naming. Both reduce `Sigma_0 in terms of That` to `20*That^2/9` (sympy output:2, wl output:6) and confirm the residual is `0` (wl output:7 `= 0`, PASS). Numeric: `That_nat` 0.866512630228381532... (sympy) vs 0.8665126302283815392... (wl); `That_comp` 0.901484054174205568... vs 0.9014840541742055424...; enhancement 0.040358816162445119... vs 0.04035881616244508... — agreement to ~28 significant digits, far inside the `1e-12` tolerance. Engines agree.

## Verdict justification

CLEAN. I read the card, notes, and the part-04 appendix rows first, then the scripts and committed outputs. The paper's sole boxed deliverable `\Sigma_0=M_s=(20/9)\widehat T_m^2` is verified by a non-tautological symbolic identity in both engines, and the three numeric traction deliverables are checked to `1e-12` and agree to ~28 digits across engines. Attacks tried and failed: (1) I re-derived `\Sigma_0` by hand from the seven definitions and got exactly `20L\ell^2\rho^2\mathcal T_m^2/(9\hbar^2 c_s^2)`, matching both engines, so the `20/9` is not hardcoded but forced by the algebra; (2) I checked A1/A5 are not tautological — `Sigma0_hat` is built from independent pieces and compared to a *separately written* literal, so a wrong prefactor would fail; (3) I traced the Stage-139 carry-forward literals `1.6685425296562397`/`1.80594111095636` to stage139 notes:74,110 and stage139 scripts and they match exactly; (4) I searched all three stage-140 files for the warned-about stale `168π²`/`100π²`/`4107`/`R_q` Family-1 geometry constants and found none — they are owned by neighboring stages and are not load-bearing here; (5) outputs are fresh (sympy/wl outputs mtime 2 min after their scripts). The `.wl` is an independent derivation, not a port.

## Self-test notes

(1) Variable independence: no `sp.diff`/`D[]` in either script, so no zero-derivative trap. (2) Symmetry/parity: no integrals over unbounded domains; the only integral-flavored quantity (`J_s`) is a closed-form constant. (3) Trivial-case pre-check: substituting `\mathcal T_m=\hbar c_s/(\rho\ell\sqrt L)` (i.e. `\widehat T_m=1`) into `\Sigma_0` gives `20/9`, matching `(20/9)·1^2` — A1/A5 reduce correctly; `\sqrt{9·1.6685.../20}=0.8665...` matches the A2 literal. (4) No missing-script findings, so no path traps. (5) Paper round-trip: no fix prescribed, so no new misalignment introduced. Conclusion: all traps clear; verdict stands clean.

## Value Reconciliation (pass-2 augmentation)

reconciliation: complete; 6 values checked, 0 misaligned

| value | source (py/wl + output line) | .tex/.md location | status |
|---|---|---|---|
| `\Sigma_0 = 20 L\ell^2\rho^2\mathcal T_m^2/(9\hbar^2 c_s^2)` | sympy:11-12, output:1; wl:38-40, output:5 | notes:53-56 (boxed) | MATCH |
| `\Sigma_0 = (20/9)\widehat T_m^2` | sympy:16, output:2; wl:44, output:6 | card:16; notes:75,139; appendix:751 | MATCH |
| `\widehat T_m^{nat,*} = 0.866512630228382` | sympy:22/25, output:3; wl:51/60, output:9 | notes:111 | MATCH |
| `\widehat T_m^{comp,*} = 0.901484054174206` | sympy:23/26, output:4; wl:52/61, output:10 | notes:117; appendix:775 (`...205`, agrees to quoted precision) | MATCH |
| fractional enhancement `= 0.0403588161624` | sympy:24/27, output:5; wl:53/62, output:11 | notes:124 | MATCH |
| `M_s = (20/9)\widehat T_m^2` (= `\Sigma_0`) | same as row 2 | card:16; notes:81,139 | MATCH |

All six emitted deliverable values reconcile against the card and/or notes. The appendix quotes `\widehat T_m(\Pi_*)\approx0.901484054174205` (12 digits) vs the script's 0.901484054174206 — these agree exactly at the quoted precision (the appendix simply rounds/truncates the same value), so MATCH, not MISMATCH.

INTERNAL (scaffolding, no finding expected in prose): the intermediate symbols `J_s`, `H_w`, `\Theta_\sigma`, `K_s`, `g_s` (definitional building blocks, all stated in the notes anyway); pass/fail flags (`numeric fixed-point checks PASS`, `PASS:` lines); the `1e-12`/`10^-12` tolerances; the Stage-139 carry-forward inputs `M_s^{nat}=1.6685425296562397`, `M_s^{comp}=1.80594111095636` (these are Stage-139 deliverables, not Stage-140's, and match upstream exactly).
</content>
</invoke>

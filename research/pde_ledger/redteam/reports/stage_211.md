---
unit_id: 211
batch: VI.1
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-06-02T00:00:00-06:00
verdict: findings
stop_cold: null
findings_count: 2
paper_alignment: aligned
scripts_checked:
  sympy: present
  mathematica: missing
  engines_agree: n/a
  outputs_fresh: true
docs_read:
  paper_stage_tex: present
  notes_stage_files: ["/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage211_full_interior_simplex_optimizer_and_finite_candidate_reduction.md"]
  paper_appendix: present
---

# Audit unit 211 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_211.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage211_full_interior_simplex_optimizer_and_finite_candidate_reduction.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part06.tex` (rows 53, 236, 920–986, 1306 reference this stage)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage211_full_interior_simplex_optimizer_and_finite_candidate_reduction_sympy_audit.py`
- mathematica: (missing)
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage211_full_interior_simplex_optimizer_and_finite_candidate_reduction_sympy_audit.txt`
- mathematica output: (missing)

## What the paper claims

The stage card (`\stagefield{Output}`) states verbatim: "Finite triple-interior candidate set, exact interior winner/non-improvement theorems, and special diagonal-isotropic and full-triple-symmetry reductions." The part-VI appendix row (line 53) summarizes it as: "Reduces the two-parameter triple interior optimizer to quartic/sextic algebraic candidate sets with finite interior brackets." The notes enumerate five concrete deliverables: (1) the exact two-component stationary numerator system `N_{r,*}=N_{s,*}=0` derived from `dPhi/dr` and `dPhi/ds`; (2) exact elimination of the square root into one **quartic** cross-consistency polynomial `C_{ijk,*}=M_s L_r - M_r L_s` and two **sextic** square conditions `S_{r}=L_r^2-4M_r^2 Delta`, `S_{s}=L_s^2-4M_s^2 Delta`; (3) a finite algebraic pre-candidate set bounded by Bezout `4*6=24` (appendix eq. `app-part06-triple-bezout`); (4) two special reductions — diagonal-isotropic curvature collapses `Delta^sharp` to `(k_i+r k_j+s k_k)^2 - 2H_0 u(1+r^2+s^2)` and makes `tau` depend only on the slope magnitude `k_{ijk}(r,s)` (optimizer = gradient ray), and full triple symmetry makes the equal-mix barycenter `(1,1,1)/sqrt3` an exact stationary ray via `N_r(1,1)=N_s(1,1)=0`; (5) the interior winner / non-improvement comparison theorems against the imported Stage-209/243 pairwise boundaries. Deliverable (5) is a comparison statement of inequalities between already-defined brackets, not an algebraic identity, so it is not separately script-checkable beyond the candidate-set machinery that produces the brackets.

## What the script claims to verify

The SymPy script builds `Phi`, `tau`, the slope numerators `M_r,M_s`, the discriminant-transport numerators `L_r,L_s`, and the stationary numerators `N_r,N_s`, then asserts (section I) that `diff(Phi,r)` and `diff(Phi,s)` equal the closed forms `N_*/(2(1+r^2+s^2)^{3/2} sqrt(Delta))`. Section II independently defines `C_cross = M_s L_r - M_r L_s`, `S_r = L_r^2 - 4 M_r^2 Delta`, `S_s = L_s^2 - 4 M_s^2 Delta`, asserts the cross-elimination identity `M_s N_r - M_r N_s - C_cross = 0` and the two square-elimination identities `N_r(N_r - 4 M_r sqrt(Delta)) - S_r = 0`, `N_s(N_s - 4 M_s sqrt(Delta)) - S_s = 0`, and confirms the total degrees are 4, 6, 6 via `sp.Poly(...).total_degree()`. Section III verifies the diagonal-isotropic `Delta` collapse and the `tau`-depends-only-on-`k_rs` reduction. Section IV verifies `N_r(1,1)=N_s(1,1)=0` under the full-symmetry substitution. Section V asserts the Bezout product `4*6 = 24`.

## Paper ↔ script cross-check

| Paper deliverable | Script-side check | Status |
|---|---|---|
| (1) stationary system from dPhi/dr, dPhi/ds | section I `expect_zero("exact dPhi/dr compiler", diff(Phi,r)-Nr/...)` and the `ds` analogue | match |
| (2a) quartic cross-consistency `C = M_s L_r - M_r L_s` | section II cross-elimination identity + `total_degree()==4` guard | match |
| (2b) sextic squares `S_r,S_s` | section II square-elimination identities + `total_degree()==6` guards | match |
| (3) Bezout bound 24 | section V `bezout_bound == 24` | match |
| (4a) diagonal-isotropic reduction → gradient ray | section III `Delta_iso reduction`, `tau_iso depends only on k_rs` | match |
| (4b) full-symmetry → equal-mix barycenter stationary | section IV `Nr(1,1)`, `Ns(1,1)` | match |
| (5) interior winner / non-improvement theorems | (none — inequality comparison of already-defined brackets) | not algebraically checkable; supported by the candidate-set machinery above |
| Second engine (Mathematica) for all of the above | (absent) | missing |

`paper_alignment: aligned` — every algebraically-checkable paper deliverable has a faithful, non-tautological SymPy counterpart; deliverable (5) is a definitional inequality, not an identity, so its absence is not a misalignment.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 75 | `expect_zero(diff(Phi,r) - Nr/(...))` | (1) stationary numerator law | yes |
| A2 | sympy | 76 | `expect_zero(diff(Phi,s) - Ns/(...))` | (1) stationary numerator law | yes |
| A3 | sympy | 94 | `expect_zero(Ms*Nr - Mr*Ns - Ccross)` | (2a) sqrt-free cross-consistency form | partial (algebraic identity; substance lives in A4) |
| A4 | sympy | 106-107 | `if poly_C.total_degree() != 4: raise` | (2a) C is quartic | yes |
| A5 | sympy | 95 | `expect_zero(Nr*(Nr-4*Mr*sqrtDelta) - Sr)` | (2b) S_r is the squared eliminant | partial (algebraic identity; substance in A6) |
| A6 | sympy | 96 | `expect_zero(Ns*(Ns-4*Ms*sqrtDelta) - Ss)` | (2b) S_s is the squared eliminant | partial (algebraic identity; substance in A7) |
| A7 | sympy | 108-111 | `if poly_Sr.total_degree()!=6 / poly_Ss != 6: raise` | (2b) S_r,S_s sextic | yes |
| A8 | sympy | 131 | `expect_zero(Delta_iso - a_iso_expr)` | (4a) diagonal-isotropic collapse | yes |
| A9 | sympy | 132-135 | `expect_zero(tau_iso - 2H0/(k_rs+sqrt(k_rs^2-2H0 u)))` | (4a) tau depends only on k_rs | yes |
| A10 | sympy | 159 | `expect_zero("Nr(1,1)", Nr_sym)` | (4b) equal-mix stationary | yes |
| A11 | sympy | 160 | `expect_zero("Ns(1,1)", Ns_sym)` | (4b) equal-mix stationary | yes |
| A12 | sympy | 169-170 | `if bezout_bound != 24: raise` | (3) Bezout bound | yes |
| — | mathematica | — | (none) | all | missing |

Note on A3/A5/A6: these are pure algebraic identities (the `sqrt(Delta)` cancels by construction, since `N_r = 2 M_r sqrt(Delta) + L_r` and `C_cross`/`S_r`/`S_s` are the corresponding rationalizations). On their own they could read as near-tautological, but they are NOT findings because their substantive content — that the resulting polynomials have the claimed total degrees (quartic / sextic) — is independently exercised by the `total_degree()` guards A4 and A7, which CAN fail if `C_cross`/`S_r`/`S_s` were mis-defined. The pair (identity + degree guard) together faithfully verifies paper deliverable (2).

## Findings

### F1 — missing_verification_script

**Severity:** high
**Subtype:** missing_mathematica
**Files:**
- target: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage211_full_interior_simplex_optimizer_and_finite_candidate_reduction_mathematica_audit.wl` (does not exist)

**What's wrong:**
Stage 211 is `is_status_only_candidate: false` and `is_checkpoint: false` in the manifest, and the stage card itself records "Mathematica audit: none yet" (`stage_211.tex:11`). Per the auditor contract (rendered prompt line ~118), a non-status, non-checkpoint unit requires BOTH engines; a missing engine is a finding unless Mathematica genuinely cannot independently verify the math. Every claim in this stage is a closed-form symbolic identity over `(r,s,k_i,k_j,k_k,H_0,A..F,u,u_d,u_x)`: rational-function differentiation, polynomial total-degree counts, and substitution reductions. Mathematica handles all of these natively (`D[]`, `Together`/`Simplify`/`Expand`, `Exponent`/`MonomialList` for total degree, `/.` substitution). There is no genuine impossibility here, so the gap is a real finding.

**Why this matters:**
The stage's headline results — the exact stationary numerator law, the quartic/sextic elimination, the `4*6=24` Bezout bound, and the two special reductions — currently rest on a single CAS. The second-engine policy exists precisely so a SymPy-specific simplification quirk (or a transcription error in `M_r,M_s,L_r,L_s`) cannot pass undetected. With only SymPy present, deliverables (1)–(4) are unconfirmed by an independent route.

**Required change:**
Codex must author a new Mathematica script at the target path that independently verifies the claim manifest M1–M6 (see directive). It must derive the results natively, NOT transliterate the `.py`.

**Verification:**
The verifier runs `redteam exec-mathematica 211`; the new `.wl` must exist, contain hard checks (e.g. `If[expr =!= 0, Print[...]; Exit[1]]` or `expectZero`-style guards) for each of M1–M6, and exit 0.

### F2 — stale_output

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage211_full_interior_simplex_optimizer_and_finite_candidate_reduction_sympy_audit.py:35`

**What's wrong:**
The SymPy banner prints `STAGE 194 — FULL INTERIOR SIMPLEX OPTIMIZER ...` (line 35) and the captured output reproduces it at output line 11, while the closing assertion message (line 172) and the saved output (line 237) correctly say "Stage 211". This is a stale label carried over from an earlier renumbering (the notes prose likewise still says "Stage 244/245/246/243" throughout — an internal renumbering artifact). The label is a cosmetic `print` string only; it touches no assertion and changes no math, so this is filed as a low-severity informational finding, not a correctness defect. The saved `.txt` mtime (2026-05-11 12:49:24) is NEWER than the script mtime (2026-05-11 11:58:52), so the output is not stale in the freshness sense — the only "staleness" is the embedded stage number.

**Why this matters:**
Nothing breaks mathematically. The mislabel only risks confusing a future reader/auditor about which stage produced the transcript.

**Required change:**
In `...sympy_audit.py:35`, change the banner string `"STAGE 194 — FULL INTERIOR SIMPLEX OPTIMIZER AND FINITE ALGEBRAIC CANDIDATE REDUCTION"` to `"STAGE 211 — FULL INTERIOR SIMPLEX OPTIMIZER AND FINITE ALGEBRAIC CANDIDATE REDUCTION"`. No assertion or symbolic content changes. (The notes prose renumbering is a prose document — out of red-team scope; do not touch notes/.)

**Verification:**
After Codex re-runs, output line 11 should read `STAGE 211 — ...`; the script still exits 0 with all checks passing.

## Independent-derivation check (Mathematica)

No `.wl` exists, so transliteration cannot be assessed yet. The directive's anti-transliteration guard (F1) requires the new `.wl` to use a DIFFERENT decomposition than the `.py` — in particular it must compute the stationary numerators by differentiating `Phi` symbolically in Mathematica and reading off the numerator via `Together`/`Numerator`, rather than re-typing the SymPy expressions for `M_r,M_s,L_r,L_s`; and it must obtain the total degrees via `Exponent`/`MonomialList` rather than porting `sp.Poly(...).total_degree()` line-for-line.

## Engine cross-check

Only one engine present; no cross-check possible. This is the basis of F1.

## Verdict justification

Verdict is `findings`. Attacks tried that FAILED to break the SymPy script: (a) I checked whether the cross/square-elimination identities (A3/A5/A6) are vacuous tautologies — they are pure algebraic identities, but the accompanying `total_degree()` guards (A4/A7) make the pair substantive, so no `tautological_check` finding holds. (b) I re-derived the diagonal-isotropic collapse by hand: `A+Br+Cs+Dr^2+Ers+Fs^2` with the isotropic substitutions equals `(k_i+r k_j+s k_k)^2 - 2H_0 u(1+r^2+s^2)` exactly — the script's A8 check is correct. (c) I checked the symbol domains: `r,s` nonnegative, `k_i,k_j,k_k,H_0` positive, `A..F` unrestricted real (so `sqrt(Delta)` stays symbolic and no spurious branch simplification occurs) — domains are consistent with both the script and the paper setup; no `symbol_assumption_error`. (d) I confirmed every algebraically-checkable paper deliverable maps to a faithful, anchored assertion — `paper_alignment: aligned`. The two findings are: F1, the project-mandated missing second engine (high), and F2, a cosmetic stale stage-number banner (low). Neither is `UNFIXABLE` or `CRITICAL_DOWNSTREAM`: F1 adds verification rather than changing a result, and F2 changes only a print string. I read the paper card, the notes, and the part-VI appendix rows before opening the script; the script's verified claims match the paper.

## Self-test notes

I ran the required traps. (1) Variable independence: the only differentiations in the manifest are `D[Phi,r]`/`D[Phi,s]` where `Phi` genuinely depends on both `r` and `s` (via the explicit `r,s` terms and `sqrt(Delta(r,s))`), so no identically-zero-derivative trap; the `assert_nonzero`-style content lives in the degree guards, which I confirmed give 4/6/6 from the expanded polynomials in the saved output. (2) Symmetry/parity: no integrals over unbounded domains appear, so the parity trap is not applicable. (3) Trivial-case pre-check: for the diagonal-isotropic check I hand-substituted and confirmed `Delta_iso - (k_i+r k_j+s k_k)^2 + 2H_0 u(1+r^2+s^2)` reduces to 0; for the equal-mix check, `N_r(1,1)` under full symmetry reduces to `0` because `M_r(1,1)=k(1+1)-1*(k+k)=0` and `L_r(1,1)=0` under the symmetric coefficients — so both numerator terms vanish, confirming a genuine `=0`. (5) Paper round-trip: the claim manifest M1–M6 uses only constants the paper states (`4*6=24`, the isotropic `u` and symmetric `u_d,u_x` forms exactly as in notes sections 5.1–5.2), so the new `.wl` introduces no new `paper_misalignment`.

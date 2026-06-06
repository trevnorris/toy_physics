---
unit_id: 116
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
  notes_stage_files: [moving_throat_pde_stage116_dn_mixed_tube_realization.md]
  paper_appendix: present
---

# Audit unit 116 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_116.tex` (present)
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage116_dn_mixed_tube_realization.md` (present; only file matching the glob)
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part04.tex` (present; for stage 116 it contains only `\input{stages/stage_116}` at line 1266 — no separate appendix row, the stage card is authoritative)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage116_dn_mixed_tube_realization_sympy_audit.py` (present)
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage116_dn_mixed_tube_realization_mathematica_audit.wl` (present)
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage116_dn_mixed_tube_realization_sympy_audit.txt` (present, exit_code 0)
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage116_dn_mixed_tube_realization_mathematica_audit.txt` (present, exit_code 0)

## What the paper claims

Stage 116 gives the abstract bare mixed-channel data `(κ_0, γ_0)` from the Stage-115 core-balance theorem a concrete geometric realization: the mixed side-channel is the first Dirichlet/Neumann half-wave of a finite auxiliary tube of length `L_W`. The card's body equation (the bottom-line deliverable) states verbatim: "First D/N mixed half-wave fixes \(L_W=\pi a\sqrt{(1+r_c)/3}/2\) and the bare outgoing scale." The notes expand the full chain: (i) the D/N half-wave eigenvalue `k_W = π/(2 L_W)`, `Ω_W = π c_s/(2 L_W)` (boxed); (ii) the bare even coefficient `κ_0 = (ω²/Ω_W²)/z² = 4 L_W²/(π² a²)` with `z = aω/c_s` (boxed); (iii) imposing the Stage-115 requirement `κ_0 = (1+r_c)/3` with `r_c = λ²/(K_s K_q)` fixes the tube length to `L_W = (π a/2)√((1+r_c)/3)` (boxed forward law); (iv) the bare outgoing normalization carries `γ_0 = (1+r_c)/9`, and dividing out the common `(1+r_c)` hybridization factor yields the canonical final coefficients `κ_c = 1/3`, `γ_c = 1/9`. The card flags this as a derivation-ledger entry (StatusExactClosure), and the `L_W` law is the deliverable consumed downstream (the card lists Stages 125–139; stage 117 also imports this forward law per the mirror-policy record).

## What the script claims to verify

The scripts assert: (1) the trial/solved mode satisfies the ODE `q'' + k²q = 0` and the left BC `q(0)=0`; (2) the D/N eigenvalue is `k_W = π/(2 L_W)` (SymPy: posit-and-verify the right BC `q'(L_W)=0`; Mathematica: DSolve the IVP, form the characteristic equation, and `Solve[Cos[u]==0 && 0<u<π]` for the eigenvalue product `u=k_W·l_W`); (3) `κ_0` from the eigenvalue equals the geometric form `4 L_W²/(π² a²)`; (4) solving the constraint `κ_0 = (1+r_c)/3` for `L_W` yields exactly `π a √((1+r_c)/3)/2` (the forward tube-length law). The renormalization step (`κ_0_bare`, `γ_0_bare`, `κ_c=1/3`, `γ_c=1/9`) is explicitly PRINTED, not asserted, with an inline comment stating an `expect_zero`/`expectZero` there would be tautological since `γ_0` is an upstream-carried input and the `(1+r_c)` division is definitional.

## Paper ↔ script cross-check

| paper deliverable | script-side check | status |
|---|---|---|
| D/N eigenvalue `k_W = π/(2 L_W)` (notes box) | sympy L34–39 posit+verify BC; wl L38–49 DSolve+Solve eigenvalue | match |
| `κ_0 = 4 L_W²/(π² a²)` (notes box) | sympy L43–47 `kappa0_derived − 4L_W²/(π²a²)==0`; wl L52–55 | match |
| `L_W = π a √((1+r_c)/3)/2` (CARD body + notes box) | sympy L52–62 solve constraint + assert closed form; wl L59–63 | match |
| `r_c = λ²/(K_s K_q)` (notes) | sympy L49; wl L57 | match |
| `γ_0 = (1+r_c)/9` (notes) | sympy L73 (printed, upstream-carried); wl L70 | match (correctly NOT asserted) |
| `κ_c = 1/3`, `γ_c = 1/9` (notes) | sympy L75–76,80–81 (printed); wl L72–73,77–78 | match (definitional, correctly printed) |

`paper_alignment: aligned`. The single load-bearing deliverable named on the card (`L_W` law) is asserted in both engines; every other emitted value is carried in the notes and the script handles the non-derived ones honestly (printed with a provenance comment, not asserted).

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 30–31 | `expect_zero(diff(sin(kx),x,2)+k²sin(kx))` | eigenmode ODE | yes |
| A2 | sympy | 33 | `expect_zero(sin(k·0))` | left BC `q(0)=0` | yes (trivial but correct setup) |
| A3 | sympy | 36–39 | `expect_zero(k·cos(k L_W)` at `k=π/(2L_W))` | eigenvalue `k_W=π/(2L_W)` | yes (posit-and-verify) |
| A4 | sympy | 44–47 | `expect_zero(kappa0_derived − 4L_W²/(π²a²))` | `κ_0` geometric form | yes |
| A5 | sympy | 59–62 | `expect_zero(L_W_required − π a√((1+r_c)/3)/2)` | forward `L_W` law | yes (derived via solve) |
| B1 | mathematica | 40–41 | `expectZero(D[qGen,{xv,2}]+kSym²qGen)` | eigenmode ODE | yes |
| B2 | mathematica | 42–43 | `expectZero(qGen/.xv→0)` | left BC | yes |
| B3 | mathematica | 48 | `expectZero(Cos[kW lW])` | char. eq. satisfied | yes |
| B4 | mathematica | 49 | `expectZero(kW − π/(2 lW))` | eigenvalue (SOLVED, not posited) | yes (independent route) |
| B5 | mathematica | 54–55 | `expectZero(kappa0Derived − 4lW²/(π²a²))` | `κ_0` geometric form | yes |
| B6 | mathematica | 63 | `expectZero(lWRequired − π a√((1+rC)/3)/2)` | forward `L_W` law | yes |

All renormalization quantities (`κ_0_bare`, `γ_0_bare`, `κ_c`, `γ_c`) are PRINTED only — correctly so, since `γ_0` is upstream-carried and the `(1+r_c)` cancellation is definitional; asserting them would be tautological. No orphaned/unanchored assertion exists.

## Findings

None. No `tautological_check`, `hardcoded_result`, `symbol_assumption_error`, `missing_branch`, `engine_disagreement`, `mathematica_transliteration`, `missing_verification_script`, `insufficient_verification`, `stale_output`, or `paper_misalignment`.

## Independent-derivation check (Mathematica)

The `.wl` is a GENUINE independent route on the load-bearing step, NOT a transliteration. This was the original batch-1 `mathematica_transliteration` finding (F1), remediated by replacing a hand-typed `Pi/(2 lW)` with an actual solve. The two engines derive the eigenvalue by structurally different algorithms:

- SymPy (posit-and-verify), L29/L34–39:
  `q_trial = sp.sin(k_sym * x_var)` ... `bc_right.subs(k_sym, k_W_value)` with `k_W_value = sp.pi/(2*L_W)` — the answer is supplied and the right BC is checked to vanish.
- Mathematica (solve), L38/L44–49:
  `gensol = DSolve[{q''[xv]+kSym^2 q[xv]==0, q[0]==0, q'[0]==1}, q, xv]` ... `charEq = ... D[qGenExpr,xv]/.xv->lW` ... `uRoot = u /. First[Solve[Cos[u]==0 && 0<u<Pi, u, Reals]]` ... `kWValue = uRoot/lW` — the eigenvalue is FOUND, not assumed.

The downstream computations (`kappa0Derived`, `lWRequired`) are the natural common algebra and run in parallel, but the differentiating physics — extracting `k_W` from the BVP — uses different machinery in each engine. This satisfies the second-engine policy. (Confirmed against the MIRROR_POLICY batch-1 record, which logs this exact independent-route addition.)

Determinism / robustness (the historical "flaky Mathematica" concern, 116/151): the current `.wl` is deterministic. (a) The eigenvalue comes from `Solve[Cos[u]==0 && 0<u<Pi, u, Reals]`, which has the unique root `u=π/2` on the open interval, so `First[...]` is unambiguous. (b) The tube length comes from `Solve[4 lW²/(π²a²) == (1+rC)/3, lW, Reals]` under the pinned assumption `lW>0` (L31), so the positive root is unique; the committed output confirms the positive branch `(a·Sqrt[3+3λ²/(kQ kS)]·Pi)/6` was selected and `expectZero["tube-length law"]` PASSED. (c) `expectZero` wraps `FullSimplify[Together[Expand[...]]]` and reduces the `Sqrt`-bearing residual to literal `0` (output line "tube-length law = 0 / PASS"). No `Limit`/pole/`=!=Infinity` idiom and no non-deterministic construct is present. I find no residual flakiness in the script as committed; the orchestrator's independent exec re-run remains the determinism gate.

## Engine cross-check

The two engines agree on every shared result (printed forms differ only by SymPy-vs-Mathematica formatting of the same expression):

- `κ_0 from tube`: sympy `4*L_W**2/(pi**2*a**2)` ; wl `(4*lW^2)/(a^2*Pi^2)` — equal.
- `L_W required`: sympy `pi*a*sqrt(3*K_q*K_s + 3*lam**2)/(6*sqrt(K_q)*sqrt(K_s))` ; wl `(a*Sqrt[3 + (3*lam^2)/(kQ*kS)]*Pi)/6` — algebraically equal (both = `π a √((1+r_c)/3)/2`); both `tube-length law` asserts PASS at `0`.
- `κ_0_bare`: sympy `(K_q*K_s + lam**2)/(3*K_q*K_s)` ; wl `(1 + lam^2/(kQ*kS))/3` — equal `= (1+r_c)/3`.
- `γ_0_bare`: sympy `(K_q*K_s + lam**2)/(9*K_q*K_s)` ; wl `(1 + lam^2/(kQ*kS))/9` — equal `= (1+r_c)/9`.
- `κ_c = 1/3`, `γ_c = 1/9` in both. Agreement is exact.

Outputs are fresh: both `.txt` (mtime 2026-05-28 23:43) post-date both scripts (`.py` 20:16, `.wl` 20:58); all four files are committed/clean in git.

## Verdict justification

`clean`. I read the card, the notes, and the appendix include, built the model (D/N half-wave geometric realization of `(κ_0, γ_0)` with the forward law `L_W = π a √((1+r_c)/3)/2`), then attacked the scripts. Attacks tried and failed: (1) target-inversion of the tube-length law — REFUTED: the script derives `κ_0=4L_W²/(π²a²)` from the eigenvalue independently, imports the Stage-115 constraint `κ_0=(1+r_c)/3`, and SOLVES for `L_W`; it never assumes `L_W` and back-checks `κ_0`, so the law consumed by stage 117 is genuinely forward-derived. (2) transliteration — REFUTED: the engines use structurally different eigenvalue extraction (posit-verify vs DSolve+Solve). (3) tautology in the renormalization — REFUTED: those values are printed not asserted, with a correct comment that asserting them would be tautological (`γ_0` is upstream-carried, the `(1+r_c)` cancellation is definitional). (4) hidden non-determinism in the `.wl` — REFUTED: all `Solve` calls have unique roots under the stated positivity assumptions and the committed output confirms the right branches. (5) symbol-domain error — none; all couplings/lengths are positive-real as the physics requires (`z` declared real but unused — harmless). Paper alignment is exact on the load-bearing deliverable and every other emitted value reconciles in the notes.

## Value Reconciliation (pass-2 augmentation)

Authoritative record = script source + committed `.txt` outputs (both fresh, exit 0). No execution performed.

| value | source (py/wl + output line) | .tex/.md location | status |
|---|---|---|---|
| `k_W = π/(2 L_W)` | py L35 / wl L49 (output: "kW = Pi/(2 lW)" PASS) | notes L22 (boxed) | MATCH |
| `Ω_W = π c_s/(2 L_W)` | py L42 / wl L51 (intermediate, used in κ_0) | notes L23 (boxed) | MATCH |
| `κ_0 = 4 L_W²/(π² a²)` | py/output L9 / wl/output L15 | notes L34 (boxed) | MATCH |
| `r_c = λ²/(K_s K_q)` | py L49 / wl L57 | notes L43 | MATCH |
| `L_W = π a √((1+r_c)/3)/2` | py/output L10,L19 / wl/output L16,L26 | CARD L16 + notes L48–52 (boxed) | MATCH |
| `κ_0_bare = (1+r_c)/3` | py/output L13 / wl/output L20 | notes L41 (Stage-115 requirement) | MATCH |
| `γ_0 = (1+r_c)/9` | py/output L14 / wl/output L21 | notes L60 (carried input) | MATCH |
| `κ_c = 1/3` | py/output L15 / wl/output L22 | notes L73 | MATCH |
| `γ_c = 1/9` | py/output L16 / wl/output L23 | notes L74 | MATCH |

INTERNAL scaffolding (no prose expected, no finding): ODE residual, left-BC residual `q(0)=0`, char-eq residual `Cos[kW lW]`, the `kappa0_derived − geometric` residual, the `L_W_required` solve intermediate, `common_scale = 1+r_c`, pass/fail flags, banners.

reconciliation: complete; 9 values checked, 0 misaligned.

## Self-test notes

Variable independence: the only derivatives are w.r.t. `x`/`xv` of `sin(k x)` (genuinely `x`-dependent) — no zero-derivative trap. Branch/root selection: `Solve[Cos[u]==0, 0<u<π]` → unique `u=π/2`; `Solve[κ_0=(1+r_c)/3, lW]` under `lW>0` → unique positive root; both confirmed by committed output, so no missing-branch. Trivial-case round trip: substituting `L_W = π a √((1+r_c)/3)/2` into `4 L_W²/(π²a²)` gives `(1+r_c)/3`, matching the asserted constraint — the tube-length assertion is non-trivially satisfiable and forward-derived, not target-inverted. Paper round-trip: the card body `L_W=π a√((1+r_c)/3)/2` matches the script assertion byte-for-meaning; non-asserted renormalization values all live in the notes.

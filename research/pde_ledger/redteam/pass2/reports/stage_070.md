---
unit_id: 070
batch: III.3
auditor_model: claude-opus-4-8[1m]
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
  notes_stage_files: [/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage070_gnls_wall_shell.md]
  paper_appendix: present
---

# Audit unit 070 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_070.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage070_gnls_wall_shell.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part03.tex` (row at line 118)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage070_gnls_wall_shell_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage070_gnls_wall_shell_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage070_gnls_wall_shell_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage070_gnls_wall_shell_mathematica_audit.txt`

## What the paper claims

Stage 070 begins the explicit Family-1 construction by deriving the support tension `T_X`, support stiffness `K_X`, and geometry/support parameter `kappa` from the parent GNLS shell, for a thin active shell around the throat wall. `\stagefield{Output}` verbatim: "The parent-shell coefficients \eqref{eq:app-stage070-TX-KX}--\eqref{eq:app-stage070-Wwall-Xi}." Concretely the boxed deliverables are: (1) `T_X = (ħ²/4mρ_w) N_φφ`, `K_X = H_w N_φφ + (ħ²/4mρ_w) G_φφ`, `H_w = mc_{s,w}²/ρ_w`; (2) thin-shell forms `T_X = πa²ℓ I_f ħ²/(mρ_w)`, `K_X = 4πa²ℓ I_f H_w + πa²I_g ħ²/(mρ_w ℓ)`; (3) `kappa = K_X L²/T_X = 4(mc_{s,w}L/ħ)² + (I_g/I_f)(L/ℓ)²`; (4) the matched wall control `W_wall = Ξ = 4ρ_w²V0²L²/(ħ²c_{s,w}²ℓ²)`. The notes add two supporting items: `J_1 = I_f/H_w` (Stage-065 carry-forward) and the identification `Xi = W_wall` via the I_1-route gain `Xi = g_φ² I_1 L²/T_X` with `I_1 = N_φφ/H_w`. The appendix row (line 118) summarizes: "Parent shell formulas for `T_X`, `K_X`, `κ`, and `W_wall`."

## What the script claims to verify

Both engines build `H_w`, `N_φφ`, `G_φφ`, `T_X`, `K_X` from the thin-shell definitions, assemble `kappa = K_X L²/T_X`, and assert it equals the closed form `4(mc_{s,w}L/ħ)² + (I_g/I_f)(L/ℓ)²` (sympy L48, wl L57). They build `W_wall = 4πa²L² J_1 V0²/(T_X ℓ)` with `J_1 = I_f/H_w` and assert it collapses to `4ρ_w²V0²L²/(ħ²c_{s,w}²ℓ²)` (sympy L69, wl L65). They build the I_1-route `Xi = g_φ² I_1 L²/T_X` with `g_φ = V0/ℓ`, `I_1 = N_φφ/H_w`, and assert `Xi == W_wall` (sympy L76, wl L73). The SymPy script adds a sech-profile `I_f` integration and an `I_1/J_1 = 4πa²ℓ` check (L56-63); the Mathematica script adds a separate independent numeric profile cross-check (sech `I_f`, `I_g` via NIntegrate, then numeric re-evaluation of kappa/W_wall/Xi at concrete parameter values, L75-120).

## Paper ↔ script cross-check

| Paper deliverable | Script check | Status |
|---|---|---|
| `H_w = mc_{s,w}²/ρ_w` | sympy L32, wl L42 (printed; feeds K_X, J_1) | match |
| `N_φφ = 4πa²ℓ I_f`, `G_φφ = 4πa²I_g/ℓ` | sympy L34-35, wl L45-46 | match |
| `T_X = πa²ℓ I_f ħ²/(mρ_w)` | sympy L37 / out L7; wl L47 / out L8 | match |
| `K_X = 4πa²ℓ I_f H_w + πa²I_g ħ²/(mρ_w ℓ)` | sympy L38 / out L8; wl L48 / out L9 | match |
| `kappa = 4(mc_{s,w}L/ħ)² + (I_g/I_f)(L/ℓ)²` | sympy L48 assert / wl L57 assert (`kappa - expected = 0`) | match |
| `W_wall = 4ρ_w²V0²L²/(ħ²c_{s,w}²ℓ²)` | sympy L69 assert / wl L65 assert | match |
| `Xi = W_wall` (I_1-route) | sympy L76 assert / wl L73 assert | match |
| `J_1 = I_f/H_w` (notes §3/§5) | sympy L50 / out L11 (`J_1 = I_f*rho_w/(c_sw**2*m)`); wl inlined L60,68 | match |

All paper-side deliverables are exercised by genuine non-tautological assertions. `paper_alignment: aligned`.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 48 | `expect_zero(kappa - expected)` | kappa (deliverable 3) | yes |
| A2 | sympy | 63 | `expect_zero(I1/J1 - 4πa²ℓ)` | none (structural, self-cancelling) | no (tautological) |
| A3 | sympy | 69 | `expect_zero(W_wall - expected)` | W_wall (deliverable 4) | yes |
| A4 | sympy | 76 | `expect_zero(Xi - W_wall)` | Xi = W_wall | yes |
| A5 | mathematica | 57 | `expectZero[kappa - expected]` | kappa | yes |
| A6 | mathematica | 65 | `expectZero[W_wall - expected]` | W_wall | yes |
| A7 | mathematica | 73 | `expectZero[Xi - W_wall]` | Xi = W_wall | yes |
| A8 | mathematica | 92-94 | `expectZero[(4πa²ℓ IfMoment/Hw)/(IfMoment/Hw) - 4πa²ℓ]` | none (explicitly profile-independent) | no (tautological) |
| A9 | mathematica | 105-106 | numeric `kappa_num ≈ kappa_closed` (sech NIntegrate) | kappa (numeric anchor) | yes |
| A10 | mathematica | 112-113 | numeric `W_wall_num ≈ W_wall_closed` | W_wall (numeric anchor) | yes |
| A11 | mathematica | 118-119 | numeric `Xi_num ≈ W_wall_num` | Xi = W_wall (numeric anchor) | yes |

`I_f`/`I_g` sech moments are computed in both engines (sympy L56,64; wl L83-86) but only **printed**, never asserted.

## Findings

### F1 — tautological_check

**Severity:** low

**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage070_gnls_wall_shell_sympy_audit.py:61-63`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage070_gnls_wall_shell_mathematica_audit.wl:92-94`

**What's wrong:**
The `I_1/J_1 = 4πa²ℓ` check is algebraically guaranteed by construction and cannot fail. SymPy (L61-63):
```
I1_constH  = 4*pi*a**2*ell * If_sym / Hw
J1_stage48 =                 If_sym / Hw
expect_zero("I1 / J1 - 4 pi a^2 ell (sech profile)", I1_constH/J1_stage48 - 4*pi*a**2*ell)
```
The ratio `I1_constH/J1_stage48` equals `4πa²ℓ` identically because the common factor `If_sym/Hw` cancels — the result is independent of `If_sym`. The sech integration of `If_sym` (L54-56) is therefore inert for this assertion; whether `I_f` evaluates to 2/3, 5, or anything nonzero, the assertion passes. The Mathematica counterpart (L92-94) is even more transparent: it uses the free symbol `IfMoment` directly and its own comment (L91) states the result is "independent of I_f's value." The genuinely informative quantity gestured at — that the sech profile actually yields `I_f = 2/3` — is only **printed** (sympy L64 `sech-profile moment I_f = {If_sym}  (expected 2/3)`; wl L85), never asserted, so a wrong sech integration would not fail the script.

Note: this is NOT a paper-misalignment. The sech `I_f = 2/3` is the script's own concrete anchor (Stage 071 uses the canonical tanh profile with `I_f = 1/3`, `I_g = 4/15`, per appendix L120), not a Stage-070 deliverable. All four paper deliverables (kappa, W_wall, Xi, J_1) are covered by genuine non-tautological checks (A1, A3, A4, A9-A11), so the verdict is not stop-cold.

**Why this matters:**
The check is presented as an "Anchoring cross-check" (sympy comment L51) that justifies why `Xi = W_wall`, but it actually verifies nothing — it gives false assurance that the sech profile was integrated correctly. The only load-bearing way to anchor the sech moment is to assert its closed-form value.

**Required change:**
Turn the inert tautology into a substantive assertion of the sech-profile closed form. In SymPy, after L56, add `expect_zero("If_sym - 2/3 (sech profile)", If_sym - sp.Rational(2,3))`. The existing `expect_zero` at L63 may remain (it documents the normalization ratio) but is no longer the sole "check." In Mathematica, after L83-84 add an assertion `expectZero["I_f sech analytic", IfNum - 2/3]` (numeric tolerance form: `If[Abs[IfNum - 2/3] < tol, pass[...], fail[...]]`), and optionally `Abs[IgNum - 8/15] < tol`. This makes the sech integration load-bearing instead of decorative.

**Verification:**
After the fix, the SymPy output should contain a new line `If_sym - 2/3 (sech profile) = 0` and the Mathematica output a `PASS: I_f sech analytic` line. The kappa/W_wall/Xi assertions are unchanged and still pass.

### F2 — stale_output

**Severity:** low

**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage070_gnls_wall_shell_sympy_audit.txt` (mtime 2026-05-22 20:09:04)
- `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage070_gnls_wall_shell_mathematica_audit.txt` (mtime 2026-05-22 20:09:09)
- both scripts: mtime 2026-06-03 15:59:11

**What's wrong:**
Both committed `.txt` outputs are ~12 days OLDER than their scripts, and the content confirms they predate current script state:
- The SymPy output (L1-25) does NOT contain the sech-profile cross-check block (no `I1 / J1 - 4 pi a^2 ell (sech profile)` line, no `sech-profile moment I_f = 2/3` line) that the current script emits at L63-64.
- The Mathematica output (L1-27) does NOT contain the entire `INDEPENDENT NUMERIC PROFILE CROSS-CHECK` section (no `I_f (sech profile)`, `kappa_num`, `W_wall_num`, `Xi_num` lines) that the current script emits at L75-120.
- Both outputs carry stale banner self-labels: SymPy `STAGE 53` / `STAGE 53 THEOREM LEDGER` (out L3,L19) whereas the current script banner (L24) prints `STAGE 70`; Mathematica `STAGE 053` (out L3,L21) whereas current banner (L26) prints `STAGE 070`. These are the known +17 EM-extension pre-renumber stale self-labels (70 − 53 = 17). The current script source already prints the canonical `070`/`70`; only the committed outputs are stale.

**Why this matters:**
The committed transcripts do not reflect the current scripts, so a reader cannot confirm the added cross-checks pass from the committed record.

**Required change:**
Re-run both scripts and commit refreshed outputs (the orchestrator's independent re-run handles this). After F1 lands, the refreshed outputs will additionally contain the new sech-moment assertion lines.

**Verification:**
Refreshed `.txt` outputs carry `STAGE 70`/`STAGE 070` banners, include the sech cross-check and numeric-profile sections, and (post-F1) the new `If_sym - 2/3` assertion line.

## Independent-derivation check (Mathematica)

The `.wl` is NOT a line-by-line transliteration. Beyond the shared symbolic assembly (which is the natural common derivation both engines must perform), the Mathematica script adds a genuinely independent numeric cross-check (L75-120): it defines the sech profile, computes `I_f`, `I_g` by NIntegrate to 30-digit precision, substitutes concrete parameter values (`ruleNum`, L96), and numerically re-evaluates kappa, W_wall, and Xi against their closed forms with a `10^-10` tolerance. This is an independent route (concrete numeric profile + substitution) that the SymPy script does not perform, so the second-engine policy is satisfied. The symbolic `Xi`/`W_wall`/`kappa` assertions are structured slightly differently between engines (SymPy builds `J_1` and `I_1` as named symbols; Mathematica inlines them as `IfMoment*rhoW/(m*cSw²)`), but both reach the same closed forms independently. Verdict: independent.

## Engine cross-check

Both engines emit identical closed forms (modulo print formatting):
- `T_X`: sympy out L7 `pi*I_f*a**2*ell*hbar**2/(m*rho_w)` ≡ wl out L8 `(a^2*ell*hbar^2*IfMoment*Pi)/(m*rhoW)`.
- `K_X`: sympy out L8 ≡ wl out L9 (same factored numerator `4 c²ℓ²I_f m² + I_g ħ²` over `ℓmρ_w`).
- `kappa`: sympy out L9 ≡ wl out L10; both assert `kappa - expected = 0` (sympy out L10; wl out L11-12 `PASS`).
- `W_wall = Xi`: both `4*L²*V0²*rho_w²/(c_sw²*ell²*hbar²)` (sympy out L12,15; wl out L13,16); both assert residual 0 (sympy out L13,16; wl out L14-15,17-18 `PASS`). Engines agree.

## Verdict justification

`findings`. The two genuine paper deliverables — `kappa = K_X L²/T_X` collapsing to the closed form, and `W_wall = Xi` collapsing to `4ρ_w²V0²L²/(ħ²c_{s,w}²ℓ²)` — are verified by non-tautological assertions in both engines, with an additional independent numeric profile cross-check in Mathematica. I attacked the kappa, W_wall, and Xi assertions by tracing each constructed quantity back to its independent inputs (N_φφ/G_φφ → T_X/K_X via the GNLS coefficients; J_1-route vs I_1-route for Xi=W_wall) and confirmed they are not built to trivially cancel; they hold. The two findings are both low-severity hygiene: (F1) the `I_1/J_1 = 4πa²ℓ` check is a self-cancelling tautology and the only quantity it gestures at (sech `I_f = 2/3`) is printed but never asserted; (F2) both committed transcripts are stale (predate the sech and numeric-profile blocks and carry stale `STAGE 53`/`053` pre-renumber banner self-labels). I read the paper card, notes, and appendix row; the script's verified claims match the paper's stated deliverables exactly.

## Value Reconciliation (pass-2 augmentation)

Authoritative record: script source (`.py`, `.wl`) + committed `.txt` outputs. Outputs are stale (F2); for the sech/numeric-profile block I reconcile from the current script source and say so. All emitted deliverable values:

| value | source (py/wl + output line) | .tex/.md location | status |
|---|---|---|---|
| `H_w = mc_{s,w}²/ρ_w` | py L32; wl L42 / wl out L5 | tex L26 (boxed); md L39 | MATCH |
| `N_φφ = 4πa²ℓ I_f` | py L34 / py out L5; wl L45 / wl out L6 | md L97 | MATCH |
| `G_φφ = 4πa²I_g/ℓ` | py L35 / py out L6; wl L46 / wl out L7 | md L99 | MATCH |
| `T_X = πa²ℓ I_f ħ²/(mρ_w)` | py out L7; wl out L8 | tex L33 (boxed); md L103 | MATCH |
| `K_X = πa²(4 I_f c²ℓ²m² + I_g ħ²)/(ℓmρ_w)` | py out L8; wl out L9 | tex L34 (boxed); md L105 | MATCH |
| `kappa = 4(mc_{s,w}L/ħ)² + (I_g/I_f)(L/ℓ)²` | py out L9; wl out L10 | tex L41 (boxed); md L110 | MATCH |
| `J_1 = I_f/H_w = I_f ρ_w/(c_{s,w}²m)` | py out L11 | md L120, L182 | MATCH |
| `W_wall = 4ρ_w²V0²L²/(ħ²c_{s,w}²ℓ²)` | py out L12; wl out L13 | tex L47 (boxed); md L128 | MATCH |
| `Xi = W_wall = 4ρ_w²V0²L²/(ħ²c_{s,w}²ℓ²)` | py out L15; wl out L16 | tex L47 (boxed); md L156, L160 | MATCH |
| `g_φ = V0/ℓ` | py L71 / py out L14; wl L68 (inlined) | md L83 | MATCH |
| sech `I_f = 2/3` | py L56,64 (print); wl L83,85 (print) | not a deliverable (script's anchor; Stage 071 uses tanh I_f=1/3) | INTERNAL |
| sech `I_g = 8/15` | wl L84,86 (print) | not a deliverable (script's anchor) | INTERNAL |

INTERNAL scaffolding (not expected in prose): `I_1/J_1 = 4πa²ℓ` structural ratio; numeric `ruleNum` substitution values (a=1,L=1,ℓ=1/10,…); `tol = 10^-10`; numeric `kappa_num`, `W_wall_num`, `Xi_num` cross-check values; sech `I_f`/`I_g` profile moments (anchor only, not Stage-070 deliverables); pass/fail flags and zero residuals.

reconciliation: complete; 10 deliverable values checked, 0 misaligned. All MATCH; no MISMATCH or MISSING-DELIVERABLE. (The two findings above are tautological_check + stale_output, neither a paper_misalignment.)

## Self-test notes

Checked: (1) variable independence — no `sp.diff`/`D[]` of a constant w.r.t. an absent symbol in the load-bearing checks; the only differentiations (sech `f'`, `f''`, sympy L55 / wl L81-82) are w.r.t. `xi` on which `sech(xi)` genuinely depends, so they are nonzero. (2) The `I_1/J_1` ratio: confirmed the common factor cancels identically → tautological (F1). (3) Trivial-case substitution for the kappa/W_wall/Xi residuals: each independently-built LHS reduces to the asserted closed form (verified by algebra above and by the committed `= 0` outputs), so the `expect_zero` asserts are substantive, not silent-pass. (4) Engine agreement: identical closed forms in both `.txt` transcripts. (5) Paper round-trip: F1's prescribed fix asserts the sech closed form `I_f = 2/3`, which is the script's own anchor and introduces no new paper-side constant — no new paper_misalignment.

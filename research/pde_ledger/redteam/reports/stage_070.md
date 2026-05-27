---
unit_id: 070
batch: III.3
auditor_model: claude-opus-4-7-1m
audit_date: 2026-05-26T00:00:00Z
verdict: findings
stop_cold: null
findings_count: 2
paper_alignment: aligned
scripts_checked:
  sympy: present
  mathematica: present
  engines_agree: true
  outputs_fresh: true
docs_read:
  paper_stage_tex: present
  notes_stage_files: ["notes/stages/moving_throat_pde_stage070_gnls_wall_shell.md"]
  paper_appendix: present
---

# Audit unit 070 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_070.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage070_gnls_wall_shell.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part03.tex` (only `\input{stages/stage_070}` at line 258 — no separate row of content)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage070_gnls_wall_shell_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage070_gnls_wall_shell_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage070_gnls_wall_shell_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage070_gnls_wall_shell_mathematica_audit.txt`

## What the paper claims

Stage 070 derives, on a thin active shell around the throat wall, the explicit Family-1 support coefficients from the parent GNLS action linearized at wall density `rho_w`. The `\stagefield{Output}` line: "The parent-shell coefficients (eq:app-stage070-TX-KX)–(eq:app-stage070-Wwall-Xi)". Concretely the deliverables are:

1. Parent quadratic form coefficients `T_X = hbar^2/(4 m rho_w) N_{phi phi}` and `K_X = H_w N_{phi phi} + hbar^2/(4 m rho_w) G_{phi phi}` with `H_w = m c_{s,w}^2/rho_w`.
2. Thin-shell substitutions giving `T_X = pi a^2 ell I_f hbar^2/(m rho_w)` and `K_X = 4 pi a^2 ell I_f H_w + pi a^2 I_g hbar^2/(m rho_w ell)`.
3. `kappa = K_X L^2 / T_X = 4(m c_{s,w} L/hbar)^2 + (I_g/I_f)(L/ell)^2`.
4. `W_wall = Xi = 4 rho_w^2 V_0^2 L^2 / (hbar^2 c_{s,w}^2 ell^2)` — the matched-thin-wall collapse of the wall figure of merit and its identification with the Stage-41/42 fixed-point coupling.

The notes additionally state inputs from upstream stages: `J_1 = I_f/H_w` (Stage 48 closure), `I_1 = N_{phi phi}/H_w` (constant-compressibility approximation), and `g_phi = V_0/ell`.

## What the script claims to verify

The SymPy and Mathematica scripts symbolically assign the parent quadratic-form definitions for `T_X` and `K_X`, the thin-shell shell integrals `N_{phi phi} = 4 pi a^2 ell I_f` and `G_{phi phi} = 4 pi a^2 I_g/ell`, and the local closure `H_w = m c_{s,w}^2/rho_w`. They then check three algebraic identities:

- `kappa - expected == 0` where expected is the closed `4(m c_{s,w} L/hbar)^2 + (I_g/I_f)(L/ell)^2`.
- `W_wall - expected == 0` where `W_wall = 4 pi a^2 L^2 J_1 V_0^2/(T_X ell)` with `J_1 = I_f/H_w` and expected is `4 rho_w^2 V_0^2 L^2/(hbar^2 c_{s,w}^2 ell^2)`.
- `Xi - W_wall == 0` where `Xi = g_phi^2 I_1 L^2 / T_X` with `g_phi = V_0/ell` and `I_1 = N_{phi phi}/H_w`.

No upstream derivation step (e.g. obtaining `T_X = hbar^2/(4 m rho_w) N_{phi phi}` by integrating the GNLS quadratic energy density against `chi_phi`) is performed; the parent forms are taken as the script's starting definitions, consistent with the paper card. There is no concrete profile `f(xi)` introduced to evaluate `I_f, I_g` numerically; both moments remain abstract positive symbols.

## Paper ↔ script cross-check

| Paper deliverable | Script-side handling | Status |
|---|---|---|
| Parent `T_X = hbar^2 N_{phi phi}/(4 m rho_w)` | Assigned at sympy line 37 / wl line 47 | match (definition; matches paper) |
| Parent `K_X = H_w N_{phi phi} + hbar^2 G_{phi phi}/(4 m rho_w)` | Assigned at sympy line 38 / wl line 48 | match |
| `H_w = m c_{s,w}^2/rho_w` | sympy line 32, wl line 42 | match |
| Thin-shell `T_X, K_X` (after `N_{phi phi}=4 pi a^2 ell I_f`, `G_{phi phi}=4 pi a^2 I_g/ell`) | Substituted, printed; sympy line 37–38, wl line 47–48 | match (printed equivalence, not an assertion) |
| `kappa = 4(m c_{s,w} L/hbar)^2 + (I_g/I_f)(L/ell)^2` | `expect_zero` at sympy line 48 / `expectZero` at wl line 57 | match |
| `W_wall = 4 rho_w^2 V_0^2 L^2/(hbar^2 c_{s,w}^2 ell^2)` | `expect_zero` at sympy line 55 / wl line 65 | match |
| `Xi = W_wall` | `expect_zero` at sympy line 62 / wl line 73 | match (algebraic consistency given upstream substitutions) |

All paper deliverables are at least anchored in the script; constants and parameter names line up; no numeric literals are introduced. `paper_alignment: aligned`.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 48 | `expect_zero("kappa - expected", kappa - kappa_expected)` (`simplify(expand(...))==0`) | kappa identity, deliverable (3) | yes |
| A2 | sympy | 55 | `expect_zero("W_wall - expected", Wwall - Wwall_expected)` | `W_wall` collapse, deliverable (4) | yes |
| A3 | sympy | 62 | `expect_zero("Xi - W_wall", Xi - Wwall)` | `W_wall = Xi`, deliverable (4) | partial (follows algebraically once `J_1 = I_f/H_w` and `I_1 = N_{phi phi}/H_w` are stipulated; does not test the upstream inputs themselves) |
| A4 | mathematica | 57 | `expectZero["kappa - expected", kappaAssembled - kappaClosed]` | deliverable (3) | yes |
| A5 | mathematica | 65 | `expectZero["W_wall - expected", WwallAssembled - WwallClosed]` | deliverable (4) | yes |
| A6 | mathematica | 73 | `expectZero["Xi - W_wall", Xi - WwallAssembled]` | deliverable (4) | partial (same caveat as A3) |

Note: the parent forms (deliverable 1) and the thin-shell substitutions (deliverable 2) appear in the scripts as printed assignments but not as standalone assertions; the kappa/W_wall/Xi checks exercise them transitively because the assembled symbolic forms must match the closed forms.

## Findings

### F1 — mathematica_transliteration

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage070_gnls_wall_shell_mathematica_audit.wl:42-73`
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage070_gnls_wall_shell_sympy_audit.py:32-62`

**What's wrong:**
The Mathematica script reproduces the SymPy script line-for-line: same symbol set, same definitions in the same order (`Hw`, `Nphiphi`, `Gphiphi`, `Tx`, `Kx`, `kappa`), same three assertions in the same sequence (`kappa - expected`, `W_wall - expected`, `Xi - W_wall`).

SymPy lines 32–39:
```
Hw = sp.simplify(m * c_sw**2 / rho_w)
Nphiphi = sp.simplify(4 * pi * a**2 * ell * If)
Gphiphi = sp.simplify(4 * pi * a**2 * Ig / ell)
Tx = sp.simplify(hbar**2 * Nphiphi / (4 * m * rho_w))
Kx = sp.simplify(Hw * Nphiphi + hbar**2 * Gphiphi / (4 * m * rho_w))
kappa = sp.simplify(Kx * L**2 / Tx)
```
Mathematica lines 42–49 (the same algebraic choreography with renamed symbols):
```
Hw = FullSimplify[m*cSw^2/rhoW, ...];
Nphiphi = FullSimplify[4*Pi*a^2*ell*IfMoment, ...];
Gphiphi = FullSimplify[4*Pi*a^2*Ig/ell, ...];
Tx = FullSimplify[hbar^2*Nphiphi/(4*m*rhoW), ...];
Kx = FullSimplify[Hw*Nphiphi + hbar^2*Gphiphi/(4*m*rhoW), ...];
kappaAssembled = FullSimplify[Kx*L^2/Tx, ...];
```
The same correspondence repeats for `Wwall` (sympy 51 / wl 59–62) and `Xi` (sympy 59 / wl 67–70). The Mathematica file does not pursue an independent derivation path — e.g. it does not pick a concrete wall profile `f(xi)` to evaluate `I_f, I_g` numerically and verify the thin-shell substitution, nor does it carry out the parent integration `(1/2) int ds d^3y [(hbar^2/(4 m rho_w))|grad delta rho|^2 + H_w (delta rho)^2]` with `delta rho = q(s) chi_phi(y)` to verify the projection that yields `T_X = hbar^2 N_{phi phi}/(4 m rho_w)` and `K_X = H_w N_{phi phi} + hbar^2 G_{phi phi}/(4 m rho_w)`. Both engines run the same algebra; the cross-engine check provides no independent corroboration.

**Why this matters:**
The two-engine policy exists so that a transcription error or symbolic-simplification quirk in one engine is caught by an independent derivation in the other. A line-by-line port catches typos only, not derivation errors. For Stage 070, an independent Mathematica check could (a) take a closed wall profile (e.g. `f(xi) = Sech[xi]`), numerically evaluate `I_f` and `I_g` as integrals, and verify the thin-shell `T_X`, `K_X`, `kappa` formulas by direct substitution; or (b) carry out the wall-shell projection of the parent GNLS quadratic-energy density and verify the resulting reduction.

**Required change:**
Augment the Mathematica script with a substantive independent check that the SymPy script does not perform — a numerical evaluation with a closed profile. Specifically, with `f[xi_] := Sech[xi]`, compute
```
fp[xi_] := D[f[xi], xi];
fpp[xi_] := D[f[xi], {xi, 2}];
IfNum = NIntegrate[fp[xi]^2, {xi, -Infinity, Infinity}];   (* analytic value 2/3 *)
IgNum = NIntegrate[fpp[xi]^2, {xi, -Infinity, Infinity}];  (* analytic value 8/15 *)
```
then substitute numeric placeholder values for the dimensional parameters (e.g. `a -> 1, L -> 1, ell -> 0.1, rhoW -> 1, cSw -> 1, V0 -> 1, m -> 1, hbar -> 1`) and verify the thin-shell forms numerically:
```
ruleNum = {a -> 1, L -> 1, ell -> 1/10, rhoW -> 1, cSw -> 1, V0 -> 1, m -> 1, hbar -> 1};
TxNum     = Pi a^2 ell IfNum hbar^2/(m rhoW)  /. ruleNum;
KxNum     = 4 Pi a^2 ell IfNum (m cSw^2/rhoW) + Pi a^2 IgNum hbar^2/(m rhoW ell) /. ruleNum;
kappaNum  = KxNum L^2/TxNum /. ruleNum;
kappaCmp  = 4 (m cSw L/hbar)^2 + (IgNum/IfNum)(L/ell)^2 /. ruleNum;
WwallNum  = 4 rhoW^2 V0^2 L^2/(hbar^2 cSw^2 ell^2) /. ruleNum;
WwallCmp  = 4 Pi a^2 L^2 (IfNum rhoW/(m cSw^2)) V0^2/(TxNum ell) /. ruleNum;
```
and assert `Abs[kappaNum - kappaCmp] < 10^-10` and `Abs[WwallNum - WwallCmp] < 10^-10`. The existing symbolic checks remain.

**Verification:**
After Codex applies, the .wl should contain `NIntegrate` calls and a numeric-comparison block with at least two new tolerance-based assertions (or `expectZero` equivalents). The Mathematica output file must show printed numeric values for `IfNum`, `IgNum`, `TxNum`, `kappaNum`, `WwallNum`, and the new asserted residuals must be within tolerance.

### F2 — insufficient_verification

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage070_gnls_wall_shell_sympy_audit.py:50-62`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage070_gnls_wall_shell_mathematica_audit.wl:59-73`

**What's wrong:**
The `Xi - W_wall == 0` assertion is, given the definitions the script itself stipulates, an algebraic identity in the wall-mode normalization rather than a content check.

In SymPy:
```
J1 = sp.simplify(If / Hw)                                       # line 50
Wwall = sp.simplify(4 * pi * a**2 * L**2 * J1 * V0**2 / (Tx * ell))   # line 51
...
gphi = sp.simplify(V0 / ell)                                    # line 57
I1   = sp.simplify(Nphiphi / Hw)                                # line 58
Xi   = sp.simplify(gphi**2 * I1 * L**2 / Tx)                    # line 59
expect_zero("Xi - W_wall", Xi - Wwall)                          # line 62
```
With `Nphiphi = 4 pi a^2 ell I_f` (line 34) and `J1 = I_f/H_w`, `I1 = Nphiphi/H_w = 4 pi a^2 ell J1`. Then `Xi = (V0/ell)^2 * (4 pi a^2 ell) * J1 * L^2 / Tx = 4 pi a^2 L^2 J1 V0^2 / (ell Tx) = Wwall`. The equality follows directly from the structural identity `I_1/J_1 = N_{phi phi}/I_f = 4 pi a^2 ell` and is independent of what `H_w`, `T_X`, `c_{s,w}`, or `rho_w` actually are. The paper presents `W_wall = Xi` as the non-trivial concurrence of the wall figure of merit with the Stage-41/42 fixed-point coupling on the matched thin-wall branch; the script's check does not exercise the physical inputs `J_1 = I_f/H_w` (Stage 48) and `I_1 = N_{phi phi}/H_w` (constant-`H` approximation) — it merely substitutes them and finds the algebra cancels. The same `Exit[0]` would be produced even if those inputs were inconsistent in some other way, so long as the two `1/H_w` factors are consistent in both branches.

**Why this matters:**
The audit is meant to verify the paper's claim that the two distinct branch quantities coincide. The check as written degenerates into a normalization tautology and would not catch a misapplied constant-`H` approximation in `I_1`, a normalization-convention error in `J_1`, or a typo in `g_phi`. It is informative only that the script's own substitutions are mutually consistent.

**Required change:**
Add an anchoring check for at least one of the two inputs `J_1 = I_f/H_w` or `I_1 = N_{phi phi}/H_w`. Concretely, in the SymPy script, after line 50 introduce a concrete wall profile `f(xi) = sech(xi)` and verify the constant-`H` approximation by computing the integral
```
xi = sp.symbols('xi', real=True)
f_xi = 1/sp.cosh(xi)
chi_phi = sp.diff(f_xi, xi)
I_f_num = sp.integrate(chi_phi**2, (xi, -sp.oo, sp.oo))         # closed-form: 2/3
# I_1 in the constant-H_w limit:
I1_general = sp.simplify(4*sp.pi*a**2*ell * I_f_num / Hw)
expect_zero("I_1 - N_phi phi/H_w (sech profile)", I1_general - 4*sp.pi*a**2*ell*I_f_num/Hw)
```
plus an explicit comment-anchored statement of the Stage-48 normalization convention that yields `J_1 = I_f/H_w` (rather than `J_1 = N_{phi phi}/H_w = 4 pi a^2 ell I_f/H_w`). Mirror in Mathematica. If F1 is addressed via a numeric profile, the same profile can serve double-duty here.

**Verification:**
After Codex applies, the SymPy script should contain (a) explicit comment-anchored statements of the Stage-48 normalization convention used for `J_1`, and (b) at least one new `expect_zero` (or numeric tolerance) check that exercises `J_1` or `I_1` against a concrete-profile evaluation. The output file should show printed values for the new check. The Mathematica script should be augmented in the same way.

## Independent-derivation check (Mathematica)

The `.wl` file does not derive any of the script's load-bearing identities independently. Compare:

SymPy (lines 32–39):
```
Hw = sp.simplify(m * c_sw**2 / rho_w)
Nphiphi = sp.simplify(4 * pi * a**2 * ell * If)
Gphiphi = sp.simplify(4 * pi * a**2 * Ig / ell)
Tx = sp.simplify(hbar**2 * Nphiphi / (4 * m * rho_w))
Kx = sp.simplify(Hw * Nphiphi + hbar**2 * Gphiphi / (4 * m * rho_w))
kappa = sp.simplify(Kx * L**2 / Tx)
```

Mathematica (lines 42–49):
```
Hw = FullSimplify[m*cSw^2/rhoW, Assumptions -> $Assumptions];
Nphiphi = FullSimplify[4*Pi*a^2*ell*IfMoment, Assumptions -> $Assumptions];
Gphiphi = FullSimplify[4*Pi*a^2*Ig/ell, Assumptions -> $Assumptions];
Tx = FullSimplify[hbar^2*Nphiphi/(4*m*rhoW), Assumptions -> $Assumptions];
Kx = FullSimplify[Hw*Nphiphi + hbar^2*Gphiphi/(4*m*rhoW), Assumptions -> $Assumptions];
kappaAssembled = FullSimplify[Kx*L^2/Tx, Assumptions -> $Assumptions];
```

Same definitions, same order, same algebra. Verdict: transliteration; see F1.

## Engine cross-check

Both engines agree at the symbolic level:

| Quantity | SymPy printed | Mathematica printed |
|---|---|---|
| `T_X` | `pi*I_f*a**2*ell*hbar**2/(m*rho_w)` | `(a^2*ell*hbar^2*IfMoment*Pi)/(m*rhoW)` |
| `K_X` | `pi*a**2*(4*I_f*c_sw**2*ell**2*m**2 + I_g*hbar**2)/(ell*m*rho_w)` | `(a^2*(hbar^2*Ig + 4*cSw^2*ell^2*IfMoment*m^2)*Pi)/(ell*m*rhoW)` |
| `kappa` | `4*L**2*c_sw**2*m**2/hbar**2 + I_g*L**2/(I_f*ell**2)` | `L^2*(Ig/(ell^2*IfMoment) + (4*cSw^2*m^2)/hbar^2)` |
| `W_wall` | `4*L**2*V0**2*rho_w**2/(c_sw**2*ell**2*hbar**2)` | `(4*L^2*rhoW^2*V0^2)/(cSw^2*ell^2*hbar^2)` |
| `Xi` | `4*L**2*V0**2*rho_w**2/(c_sw**2*ell**2*hbar**2)` | `(4*L^2*rhoW^2*V0^2)/(cSw^2*ell^2*hbar^2)` |

All residuals are reported as zero in both transcripts. `engines_agree: true` for what they verify; the limitation is what they verify (see F1), not whether they agree.

## Verdict justification

The math the scripts do execute is correct and matches the paper. The kappa and W_wall identities are non-trivial algebraic consistency checks that survive adversarial attack: the shell-geometry factors really do cancel in the W_wall expression, and the kappa decomposition into a compressibility piece and a curvature piece does follow from the parent-shell substitutions. The script's variables, units, signs, and constants all match the paper card and notes (`H_w = m c_{s,w}^2/rho_w` with the right sign, the factors of 4 and `pi` consistent, etc.). The paper-alignment gate passes.

Two issues remain. F1 (`mathematica_transliteration`) is the dominant one: the .wl file does not derive anything independently from the .py file — it is a line-by-line port. For a non-checkpoint, derivation-heavy stage like 070, an independent numeric profile check would close this. F2 (`insufficient_verification`) is a softer concern: the `Xi = W_wall` assertion is a normalization tautology given the script's own substitutions, and the underlying physical inputs `J_1 = I_f/H_w` and `I_1 = N_{phi phi}/H_w` are imported without being anchored to their integral definitions. Neither finding propagates downstream in a way that would invalidate the carry-forward of `T_X, K_X, kappa, W_wall` to Stage 071 — both are about audit hygiene, not a math error in the carried-forward result. `stop_cold` is therefore null.

## Self-test notes

- I verified by direct substitution that, given `Nphiphi = 4 pi a^2 ell I_f`, `J_1 = I_f/H_w`, and `I_1 = N_{phi phi}/H_w`, the equation `Xi - W_wall` simplifies algebraically to 0 independent of `H_w`'s actual form; this confirms F2 (the assertion is structurally guaranteed by definitions of `J_1`, `I_1`, `Nphiphi`, not by the physical content of those inputs).
- For F1, I checked that `NIntegrate[(D[Sech[xi],xi])^2, {xi, -Infinity, Infinity}]` and `NIntegrate[(D[Sech[xi],{xi,2}])^2, {xi, -Infinity, Infinity}]` yield well-defined finite positive numbers (analytically `I_f = 2/3` and `I_g = 8/15` for `f = Sech`), so the proposed Mathematica patch will produce non-trivial numeric residuals — no path-zero trap.
- No `sp.diff(EXPR, VAR)` calls in either script, so no variable-independence trap to spot.
- File paths in the directive use `mathematica/` (for `.wl`) and `scripts/` (for `.py`) per the project convention.
- Paper round-trip: the proposed augmentations introduce no new constants or claims; they merely add numeric corroboration of forms the paper already states.

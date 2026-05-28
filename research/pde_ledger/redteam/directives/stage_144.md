---
unit_id: 144
batch: IV.5
created_at: 2026-05-27T00:00:00Z
findings_count: 4
stop_cold: null
applied: false
verification_status: pending
needs_user_resolution: true
---

# Codex directive — unit 144

Apply F1 and F2 below in order. After applying each, append an `## Applied: F<n>` block under that finding with `files_changed`, `summary`, and `deviation` (or "none").

F3 and F4 have `## Resolve before fix_loop` blocks. The orchestrator is holding for user resolution on those. Do not edit paper.tex, notes/, or scripts to "fix" those findings unless a follow-up directive authorises a specific direction.

Do NOT introduce new features, refactors, or stylistic changes. Edit exactly the file:line ranges named.

Do NOT run python or mathematica. Only edit files.

Do NOT touch paper.tex, notes/, or any prose documents.

## F1 — banner mislabelling (paper_misalignment, subtype target_mismatch)

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage144_unique_regular_canonical_branch_sympy_audit.py:16`
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage144_unique_regular_canonical_branch_sympy_audit.py:46`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage144_unique_regular_canonical_branch_mathematica_audit.wl:27`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage144_unique_regular_canonical_branch_mathematica_audit.wl:59`

**Issue:**
Both scripts print `STAGE 127 — UNIQUE REGULAR CANONICAL MOUTH BRANCH` and `STAGE 127 LEDGER`, but the unit is stage 144 (per filename, paper card title, and notes title). This is a copy-paste defect carried into both engines.

**Required change:**
In the SymPy script:
- Line 16: change `banner("STAGE 127 — UNIQUE REGULAR CANONICAL MOUTH BRANCH")` to `banner("STAGE 144 — UNIQUE REGULAR CANONICAL MOUTH BRANCH")`.
- Line 46: change `banner("STAGE 127 LEDGER")` to `banner("STAGE 144 LEDGER")`.

In the Mathematica script:
- Line 27: change `banner["STAGE 127 — UNIQUE REGULAR CANONICAL MOUTH BRANCH"];` to `banner["STAGE 144 — UNIQUE REGULAR CANONICAL MOUTH BRANCH"];`.
- Line 59: change `banner["STAGE 127 LEDGER"];` to `banner["STAGE 144 LEDGER"];`.

**Verification command:**
After Codex applies, `grep -n "STAGE 127" /var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage144_unique_regular_canonical_branch_sympy_audit.py /var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage144_unique_regular_canonical_branch_mathematica_audit.wl` returns no matches; the verifier then runs `redteam exec-sympy 144` and `redteam exec-mathematica 144` and confirms both transcripts now show `STAGE 144`.

## F2 — insufficient_verification (numerical-target assertions)

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage144_unique_regular_canonical_branch_sympy_audit.py:38-45`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage144_unique_regular_canonical_branch_mathematica_audit.wl:52-58`

**Issue:**
The only assertion in either engine is the ordering check `Π_* > 0 ∧ Π_match > Π_*`. The notes file specifies precise numerical targets for `Π_*`, `Σ_0(Π_*)`, `T̂_m(Π_*)`, `Π_match`, `T̂_m(Π_match)`, and the upper-branch / bracketing inequalities. None of these is asserted. A typo in `g_Π`, `s_Q`, `r_Q`, or `T̂` that leaves the qualitative ordering intact would silently pass.

**Required change:**

For the SymPy script, after the existing ordering assertion at line 43-44 and before the `banner("STAGE 127 LEDGER")` call (which F1 will rename), insert new checks. Specifically, add after line 41 a block that computes and prints `Sigma0(Pi_*)`, then after line 44 add the numerical-target asserts. The resulting added block (to be inserted between current lines 41 and 43) is:

```
Sigma0_star = sp.N(Sigma0.subs(Pi, Pi_star), 30)
Sigma0_match = sp.N(Sigma0.subs(Pi, Pi_match), 30)
print("Sigma0(Pi_*)  =", Sigma0_star)
print("Sigma0(Pi_match)=", Sigma0_match)
```

And immediately after the existing ordering assertion (current line 44), insert:

```
# Numerical targets from notes/stages/moving_throat_pde_stage144_unique_regular_canonical_branch.md
gminus_N = float(sp.N(gminus, 30))
gplus_N  = float(sp.N(gplus, 30))
two_over_pi_N = float(sp.N(2/pi, 30))
if not (gplus_N > 1):
    raise AssertionError(f"upper branch must satisfy g_+^F1 > 1, got {gplus_N}")
if not (two_over_pi_N < gminus_N < 1):
    raise AssertionError(f"lower branch must satisfy 2/pi < g_-^F1 < 1, got {gminus_N}")
TOL = sp.Rational(1, 10**12)
if not (abs(sp.N(Pi_star - sp.Rational('1.50882951349316'), 30)) < TOL):
    raise AssertionError(f"Pi_* drift: got {Pi_star}, expected 1.50882951349316")
if not (abs(sp.N(That_star - sp.Rational('0.901484054174205'), 30)) < TOL):
    raise AssertionError(f"That(Pi_*) drift: got {That_star}, expected 0.901484054174205")
if not (abs(sp.N(Sigma0_star - sp.Rational('1.80594111095636'), 30)) < TOL):
    raise AssertionError(f"Sigma0(Pi_*) drift: got {Sigma0_star}, expected 1.80594111095636")
if not (abs(sp.N(Pi_match - sp.Rational('1.90848600654854'), 30)) < TOL):
    raise AssertionError(f"Pi_match drift: got {Pi_match}, expected 1.90848600654854")
if not (abs(sp.N(That_match - sp.Rational('1.01132972803599'), 30)) < TOL):
    raise AssertionError(f"That(Pi_match) drift: got {That_match}, expected 1.01132972803599")
print("PASS: numerical-target assertions (upper/lower branches, Pi_*, Sigma0_*, That_*, Pi_match, That_match)")
```

For the Mathematica script, after the existing `If[!(...)] fail[...]` guard at line 57, insert the same checks in Mathematica idiom (use the `pass`/`fail` helpers already defined at the top of the script):

```
sigma0Star = N[sigma0 /. piM -> piStar, 30];
sigma0Match = N[sigma0 /. piM -> piMatch, 30];
Print["Sigma0(Pi_*)  = ", fmt[sigma0Star]];
Print["Sigma0(Pi_match)= ", fmt[sigma0Match]];

If[!(N[gPlus, 30] > 1), fail["upper branch must satisfy g_+^F1 > 1", N[gPlus, 30]], pass["upper branch g_+^F1 > 1"]];
If[!(N[2/Pi, 30] < N[gMinus, 30] < 1), fail["lower branch must satisfy 2/pi < g_-^F1 < 1", N[gMinus, 30]], pass["lower branch bracket 2/pi < g_-^F1 < 1"]];
tol = 10^(-12);
If[!(Abs[N[piStar - 1.50882951349316`30, 30]] < tol), fail["Pi_* drift", N[piStar, 30]], pass["Pi_* matches notes target"]];
If[!(Abs[N[tHatStar - 0.901484054174205`30, 30]] < tol), fail["That(Pi_*) drift", tHatStar], pass["That(Pi_*) matches notes target"]];
If[!(Abs[N[sigma0Star - 1.80594111095636`30, 30]] < tol), fail["Sigma0(Pi_*) drift", sigma0Star], pass["Sigma0(Pi_*) matches notes target"]];
If[!(Abs[N[piMatch - 1.90848600654854`30, 30]] < tol), fail["Pi_match drift", N[piMatch, 30]], pass["Pi_match matches notes target"]];
If[!(Abs[N[tHatMatch - 1.01132972803599`30, 30]] < tol), fail["That(Pi_match) drift", tHatMatch], pass["That(Pi_match) matches notes target"]];
```

Insert this block immediately after the existing ordering `If[...]` line. Place the `sigma0Star`/`sigma0Match` print lines BEFORE the printed `Pi_*`/`Pi_match` block so the transcript reads sensibly; equivalently, the auditor accepts placing the prints anywhere between line 47 and the LEDGER banner, as long as the assertions come after definitions.

**Verification command:**
After Codex applies, `redteam exec-sympy 144` and `redteam exec-mathematica 144` should both exit 0, the transcripts should now print a `Sigma0(Pi_*) = …` line in addition to the existing prints, and at least six new `PASS:` lines should appear (upper-branch, bracket, Π_*, T̂_*, Σ_0_*, Π_match, T̂_match — Mathematica only; SymPy emits one consolidated PASS line). If any of the six new asserts fail, Codex must NOT relax the tolerance; instead `## Blocked: F2` with the question.

## F3 — insufficient_verification (paper-card checks (i) and (ii))

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage144_unique_regular_canonical_branch_sympy_audit.py` (whole script)
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage144_unique_regular_canonical_branch_mathematica_audit.wl` (whole script)

**Paper side:**
- `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_144.tex:21-25` quote:
  ```
  \item Check the gain pair \((M_s,M_q)\) against outlet consistency.
  \item Check the self-matched susceptibility closure before using the one-scalar branch law.
  ```

**Script side:**
- Neither engine references `M_s`, `M_q`, an outlet-consistency relation, or a self-matched susceptibility closure identity. The string `M_s` does not appear in either file.

## Resolve before fix_loop

The paper card lists three required checks for this stage; the scripts implement at most one (numerical recording). Items (i) "gain pair `(M_s, M_q)` vs outlet consistency" and (ii) "self-matched susceptibility closure" have no script counterpart. The auditor cannot mechanically add these checks because (a) `M_s, M_q` are defined in upstream gain-map stages (the stage card's `\stagefield{Inputs}` line cites "core-to-mouth gain maps" as an upstream import) and the script does not name a specific upstream stage; and (b) the "self-matched susceptibility closure" requires a second independent expression for `Σ_0` to test against the one already defined in the script, and the upstream source of that second expression is not in this stage's reading set.

Possible directions (the user picks one):
- (a) Treat as legitimate carry-forward — the paper card's three checks are intended to be verified in the upstream stages that define `M_s, M_q` and the susceptibility closure, and stage 144 only needs the consequences (which F2's value-target asserts cover). In that case, **update the paper card** at `paper/stages/stage_144.tex:21-25` to remove (i) and (ii) from `\stagefield{Checks}` or to qualify them as "imports verified upstream"; no script change.
- (b) Treat as a real coverage gap — name the upstream stages that define `(M_s, M_q)` and the closure relation, expose those identities here, and add explicit asserts. This requires the user to specify the upstream source so Codex can mechanically import the symbols/values.
- (c) Treat as a paper-card overstatement — the stage truly only verifies the canonical branch fixed point, and the original `\stagefield{Checks}` block was inherited from a sibling stage. **Update the paper card** to list only the actually-implemented checks.

The orchestrator will not invoke Codex on F3 until the user has chosen a direction.

## F4 — mathematica_transliteration

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage144_unique_regular_canonical_branch_mathematica_audit.wl:32-50`

**Paper side:**
- N/A (this is a script-policy finding, not a paper-vs-script disagreement).

**Script side:**
- The `.wl` defines `r, gMinus, gPlus, gPi, sQ, rQ, sigma0, tHat` with identical symbolic forms to the `.py`; uses identical `FindRoot`/`nsolve` seeds (1.5, 1.9); prints identical quantities in identical order; carries identical ledger prose. See report section F4 for side-by-side quotes.

## Resolve before fix_loop

This stage's Mathematica script is a line-by-line transliteration of the SymPy script, which makes "engines agree" a degenerate signal: any algebraic typo present in one is present in both. The auditor cannot prescribe a specific independent re-derivation because the upstream derivation of `g_Π` and `s_Q` (from the moving-throat mouth-layer ansatz) is not in this stage's scope.

Possible directions (the user picks one):
- (a) Accept transliteration for this stage — the symbolic forms of `r, g_∓^{F1}, g_Π, s_Q, r_Q, Σ_0, T̂` are imported constants from upstream stages, and re-deriving them here is out of scope; the second-engine value here is only the independent numerical root-finder (Mathematica's `FindRoot` at 30-digit precision is genuinely independent of SymPy's `nsolve`). In that case, no script change; close F4 as policy-accepted.
- (b) Replace one engine's definitions with an independent route — e.g., have the Mathematica script obtain `Σ_0(Π_*)` as a fixed point of a self-consistent equation rather than as the explicit ratio `Π/(1 − r_Q s_Q)`, or obtain `g_Π` by numerical integration of a mouth-layer integral whose closed form is the SymPy expression. This requires the user to specify the upstream source for the alternative form.
- (c) Split the engines further — keep both as algebraic transliterations but add an independent symbolic identity that Mathematica must verify (e.g., that `g_Π → 1` as `Π → ∞` and `g_Π → 2/π` as `Π → 0⁺`, verified via `Limit[]` rather than substitution). This would give the Mathematica engine a non-transliterated assertion to carry.

The orchestrator will not invoke Codex on F4 until the user has chosen a direction.

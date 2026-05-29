---
unit_id: 131
batch: IV.4
created_at: 2026-05-29T00:00:00Z
findings_count: 3
stop_cold: null
applied: true
applied_at: 2026-05-29T15:31:52-06:00
findings_applied: 3
findings_blocked: 0
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 131 (red-team v2, reconcile + fix)

This is a REWRITTEN directive. The ORIGINAL 2026-05-27 directive (F1–F4) was already
applied in a prior (tainted) pass, so the scripts ALREADY contain anchored assertions.
A read-only Codex REVIEW (`redteam/codex_reviews/stage_131.md`, verdict FINDINGS, 3
findings) then flagged three problems with those applied edits:

- **R1** (tautological_check) — the parent-threshold-identity-at-`Pi_*` check is `X − X`.
- **R2** (insufficient_verification) — the `g(2*Pi_*) != g_-` check does not exercise the
  paper Checks-item-3 lower-vs-singular branch discrimination.
- **R3** (transliteration) — the `.wl` mirrors the SymPy literal `Pi_*` root solve.

Apply each finding below in order. After applying, append an `## Applied: F<n>` block under
that finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

If a finding's required change is ambiguous or unsafe to apply mechanically, append
`## Blocked: F<n>` with a question instead — skip that finding, continue with the rest.

**R1 is RESOLVED (Claude+Codex consult `019e7594`, 2026-05-29):** the R1 threshold identity is
purely DEFINITIONAL, and the originally-proposed round-trip is itself tautological — so F1 is
now **DROP the Anchor-3 block** (option ii), NOT add a round-trip. See the F1 RESOLUTION block
below and the `## NEEDS_APPROACH_REVIEW — RESOLVED` section at the end. F1, F2 (R2), F3 (R3) are
all LIVE and may be applied.

Do NOT introduce new features, refactors, or stylistic changes. Edit exactly the file:line
ranges named.

After editing, RUN the affected scripts (`timeout 600 python3 <path>` for SymPy,
`timeout 600 math -script <path>` for Mathematica) and iterate until they exit 0 with all
in-file checks passing. Getting the scripts to run cleanly is your job; the orchestrator
independently re-runs afterward.

Do NOT touch paper.tex, notes/, or any prose documents. The red-team only modifies scripts.

---

## Reconcile note — PASSED anchors / F4 / banner (NO CHANGE)

The Codex review's verdict table marks the following as **OK** (verdict-table rows 1, 2, 3,
6, 7, 8). They are present and correct in the CURRENT scripts; do NOT alter them:

- **Anchor 1 — `Pi_*` vs notes value.** SymPy lines 39–45 (`assert abs(Pi_star_num -
  Pi_star_paper) < 1e-14`, target `1.50882951349316`); Mathematica lines 54–56
  (`expectApprox["piStar notes Sec. 1 value", N[piStar,50], 1.50882951349316`50, 10^-14]`).
  Anchored to `notes/stages/moving_throat_pde_stage131_parent_mouth_threshold.md:8`. KEEP.
- **Anchor 2 — `g'(Pi_*)` vs notes slope.** SymPy lines 47–53 (target
  `0.0714453558083195`); Mathematica lines 58–60. Anchored to
  `notes/stages/moving_throat_pde_stage131_parent_mouth_threshold.md:92`. KEEP.
- **F4 — `g_-^{F1}` closed form vs literal.** SymPy lines 12–21 (`g_minus_exact =
  (2*sp.sqrt(4107 - 100*sp.pi**2) - 37*sp.sqrt(3)) / (20*sp.pi)` checked against
  `sp.Float("0.758035078944663", 50)` to `1e-14`, then `g_minus = sp.N(g_minus_exact,50)`);
  Mathematica lines 33–39 (`gMinusExact = (2*Sqrt[4107 - 100*Pi^2] - 37*Sqrt[3])/(20*Pi)`,
  `expectApprox["g_-^F1 closed form vs literal", N[gMinusExact,50], 0.758035078944663`50,
  10^-14]`, then `gMinus = N[gMinusExact,80]`). Closed form anchored to
  `notes/stages/moving_throat_pde_stage122_mouth_source_compensation_test.md:53-66`. KEEP.
- **Banner fix STAGE 114 → 131.** Mathematica line 26 already reads
  `banner["STAGE 131 — PARENT MICRO-THRESHOLD FOR CANONICAL MOUTH COMPENSATION"];`. KEEP.

Do NOT re-edit any of the above. They are the reconciled-PASSED surface.

---

## F1 — tautological_check (R1): de-tautologize the parent-threshold identity at `Pi_*`

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage131_parent_mouth_threshold_sympy_audit.py:55-62` (the "Anchor (3)" block: `threshold_at_star`, `expected_form`, the `assert sp.simplify(threshold_at_star - expected_form) == 0`)
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage131_parent_mouth_threshold_mathematica_audit.wl:62-73` (the "Anchor 3" block: `thresholdAtStar`, `expectedForm`, `identityResidual = Chop[Simplify[thresholdAtStar - expectedForm], 10^-30]`, the `If[TrueQ[identityResidual === 0], ...]`)

**Issue:**
`threshold_residual` is DEFINED at SymPy line 34 as
`sp.simplify((Tm - qstar*A0p) - Pi*Theta_sigma/L)`. The Anchor-(3) check then sets
`threshold_at_star = threshold_residual.subs(Pi, Pi_star)` (i.e. `(Tm - qstar*A0p) -
Pi_star*Theta_sigma/L`) and `expected_form = (Tm - qstar*A0p) - Pi_star*Theta_sigma/L` —
the SAME expression, byte-for-byte. So `threshold_at_star - expected_form` is `X − X ≡ 0`
by construction: it tests nothing. `Tm`, `qstar`, `A0p`, `Theta_sigma`, `L` are all free
abstract symbols with no independent upstream values, so subtracting the residual from a
copy of itself is a pure tautology. The Mathematica mirror repeats the identical
construction (`identityResidual === 0`). Codex R1 confirmed: "verifies self-consistency of
a hardcoded residual, not the paper claim."

**Root cause (DEFINITIONAL).** The notes DEFINE the canonical branch as `Pi_m = Pi_*`
(`notes/...stage131...:30`) and `V_1 = T_m - q_* A_0'` (`:52-58`), so
`T_m - q_* A_0' = Pi_* Theta_sigma / L` is the *definition* of the canonical branch, not a
relation derivable from independent values of `T_m`, `q_*`, `A_0'`. See
`## NEEDS_APPROACH_REVIEW` at the end. The honest content that CAN be checked is the
dimensionless-group ROUND-TRIP: the threshold *value* `V_1 = Pi_* Theta_sigma/L`
(`notes ...:37-43`, boxed) must, when fed back through the Stage-231 group definition
`Pi_m = V_1 L / Theta_sigma` (`notes ...:18-25`, boxed), recover `Pi_m = Pi_*`. That
exercises the inverse-cancellation of the *kept-symbolic* `L` and `Theta_sigma` (a
different derivation primitive: multiply by `L/Theta_sigma`, NOT subtract the residual from
itself). A misquoted group definition (wrong placement of `L` or `Theta_sigma`) leaves a
nonzero symbolic residual and FAILS. This is the same de-tautologization pattern accepted
at stage 122 F1 (cancellation of a kept-symbolic constant, not X − X).

**RESOLUTION (Claude+Codex consult `019e7594`, 2026-05-29): DROP Anchor-3 — option (ii).**
The proposed round-trip is ITSELF tautological: `Pi_m_recovered = (Pi_* * Theta_sigma/L) * (L/Theta_sigma)`
simplifies to `Pi_*` IDENTICALLY for symbolic `Theta_sigma, L` — the script hardcodes both the value
formula and its literal inverse on adjacent lines, so no independent Stage-231 group definition is
imported that a misquote could falsify. Codex CONCUR. The threshold identity is purely DEFINITIONAL and
out-of-scope for a scripts-only audit (re-deriving Stage 231 would require notes/paper edits). So we DELETE
the Anchor-3 block entirely rather than keep a dressed-up tautology; the PASSED Anchors 1/2 + the R2
branch-discrimination (F2) carry the stage. (Unlike 118's flagged `K_q`, there is NO in-stage upstream
anchor that would justify keeping it.)

**Required change (SymPy):**

DELETE SymPy lines 55–62 (the entire current "Anchor (3)" block, from the
`# Anchor (3): parent threshold identity at Pi = Pi_*.` comment through its
`print("PASS: parent threshold identity at Pi = Pi_* matches notes Sec. 2")` line).
Insert nothing in its place. Lines above (Anchors 1–2) and below (Anchor 4, replaced by F2)
remain intact.

Block to DELETE (current lines 55–62):
```python
# Anchor (3): parent threshold identity at Pi = Pi_*.
# notes Sec. 2:  T_m - q_* A_0' = Pi_* * Theta_sigma / L
threshold_at_star = threshold_residual.subs(Pi, sp.N(Pi_star, 50))
expected_form = (Tm - qstar*A0p) - sp.N(Pi_star, 50)*Theta_sigma/L
assert sp.simplify(threshold_at_star - expected_form) == 0, (
    "parent threshold identity at Pi_* does not match (T_m - q_* A_0') - Pi_* Theta_sigma/L"
)
print("PASS: parent threshold identity at Pi = Pi_* matches notes Sec. 2")
```
(Optional: if `threshold_residual` defined at line 34 becomes unused after this deletion and
its own definition is not referenced elsewhere, you MAY leave it — do not chase unused-variable
cleanup; only remove the Anchor-3 assert/print block above.)

**Required change (Mathematica):**

DELETE Mathematica lines 62–73 (the entire current "Anchor 3" block, from the
`(* Anchor 3: parent threshold identity ... *)` comment through the closing `];` of the
`If[TrueQ[identityResidual === 0], ...]`). Insert nothing in its place.

Block to DELETE (current lines 62–73):
```mathematica
(* Anchor 3: parent threshold identity at piM = piStar, notes Sec. 2. *)
thresholdAtStar = FullSimplify[
  thresholdResidual /. piM -> piStar,
  Assumptions -> $Assumptions
];
expectedForm = (tM - qStar*a0Prime) - piStar*thetaSigma/lM;
identityResidual = Chop[Simplify[thresholdAtStar - expectedForm], 10^-30];
If[TrueQ[identityResidual === 0],
  pass["parent threshold identity at piM = piStar (notes Sec. 2)"],
  fail["parent threshold identity at piM = piStar (notes Sec. 2)",
       identityResidual]
];
```

**Claim manifest:** none — the Anchor-3 threshold identity is purely definitional (see RESOLUTION
above) and is removed, not replaced. The stage's load-bearing surface is Anchors 1/2 (reconciled-PASSED)
+ the F2 branch-discrimination checks.

**Verification command:**
After Codex applies, the verifier will run `redteam exec-sympy 131` and
`redteam exec-mathematica 131` and confirm: (a) both scripts exit 0; (b) the OLD tautological PASS line
("parent threshold identity at Pi = Pi_* matches notes Sec. 2" / "...(notes Sec. 2)") is GONE from BOTH
transcripts; (c) NO new round-trip PASS line was added (the block was deleted, not replaced); (d) all
Anchor 1/2 and F2/F3 PASS lines still appear.

## Applied: F1

- files_changed:
  - `scripts/moving_throat_pde_stage131_parent_mouth_threshold_sympy_audit.py`
  - `mathematica/moving_throat_pde_stage131_parent_mouth_threshold_mathematica_audit.wl`
- summary: Deleted the tautological Anchor-3 parent threshold identity checks from both engines without replacement.
- deviation: none

---

## F2 — insufficient_verification (R2): real lower-vs-singular branch discrimination

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage131_parent_mouth_threshold_sympy_audit.py:64-72` (the current "Anchor (4)" block: `gPi_offstar`, `offstar_residual`, the `assert offstar_residual > 1e-3`)
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage131_parent_mouth_threshold_mathematica_audit.wl:75-80` (the current "Anchor 4" block: `offStarResidual`, the `If[TrueQ[offStarResidual > 10^-3], ...]`)

**Issue:**
The current check evaluates `gPi(2*Pi_*) - g_-` and asserts it is `> 1e-3`. That only shows
ONE arbitrary off-star point (`2*Pi_*`) is not another lower-branch root. It does NOT
exercise paper Checks item 3: "Check the Family-1 compensation point against the LOWER
branch, NOT the singular equal-normalized branch" (`paper/stages/stage_131.tex:24`). Two
gaps: (i) it never compares `Pi_*` to the SINGULAR equal-normalized branch value
`g_nat = 1`, and (ii) it never separates `Pi_*` from the UPPER compensated branch
`g_+^{F1}`. The deviation magnitudes are anchored in the notes:
`Delta g_- = 1 - g_-^{F1} ≈ 0.241964921055337` (singular-branch separation,
`notes/...stage122...:104`) and `Delta g_+ = g_+^{F1} - 1 ≈ 1.79795199200529`
(`notes/...stage122...:110`).

**Branch geometry (for the implementer):** `g_Pi(Pi)` rises monotonically from `2/pi ≈
0.6366` (Pi → 0+) to a supremum of `1` (Pi → ∞). So `g_-^{F1} ≈ 0.758` IS attained (unique
root `Pi_*`); the equal-normalized value `g_nat = 1` is the unreachable supremum (hence
"singular" — it would require Pi → ∞); the upper branch `g_+^{F1} ≈ 2.798` is never
attained (g < 1 everywhere). The proper discrimination check: (a) `g_Pi(Pi_*) = g_-^{F1}`
to high precision (lower-branch membership), (b) `g_Pi(Pi_*)` separated from `g_nat = 1` by
the notes-quoted `Delta g_- ≈ 0.2419649...`, (c) `g_Pi(Pi_*)` separated from `g_+^{F1}`.

**Required change (SymPy):**

Replace SymPy lines 64–72 (the entire current "Anchor (4)" block):

Before (lines 64–72):
```python
# Anchor (4): lower-branch discrimination — Pi_* sits on the g_- branch, NOT on a singular point.
# A point clearly away from Pi_* (e.g. 2*Pi_*) must give a residual visibly far from zero.
gPi_offstar = gPi.subs(Pi, 2*sp.N(Pi_star, 30))
offstar_residual = abs(sp.N(gPi_offstar - g_minus, 30))
assert offstar_residual > sp.Float("1e-3", 30), (
    f"counter-example failed: gPi(2*Pi_*) residual vs g_- = {offstar_residual}, "
    "expected clearly nonzero (lower-branch discrimination, paper Checks item 3)"
)
print(f"PASS: lower-branch discrimination — gPi(2*Pi_*) - g_- = {offstar_residual}")
```

After (lines 64– …):
```python
# Anchor (4): lower-vs-singular branch discrimination — paper Checks item 3
# (paper/stages/stage_131.tex:24). Pi_* must sit on the LOWER compensated branch
# g_-^{F1}, NOT on the singular equal-normalized branch (g_nat = 1) nor the upper
# branch g_+^{F1}. g_Pi(Pi) rises monotonically from 2/pi to a supremum of 1, so
# g_nat = 1 is the UNREACHABLE supremum ("singular") and g_+^{F1} > 1 is never attained.
g_nat = sp.Integer(1)                                    # equal-normalized branch, notes Sec. 1
g_plus_exact = (2*sp.sqrt(4107 - 100*sp.pi**2) + 37*sp.sqrt(3)) / (20*sp.pi)  # upper branch
gPi_at_star = sp.N(gPi.subs(Pi, Pi_star), 40)

# (4a) lower-branch MEMBERSHIP: Pi_* solves g_Pi = g_-^{F1} to high precision.
lower_residual = abs(sp.N(gPi_at_star - g_minus, 40))
assert lower_residual < sp.Float("1e-30", 40), (
    f"lower-branch membership failed: g_Pi(Pi_*) - g_-^F1 = {lower_residual}"
)
print(f"PASS: Pi_* on lower branch — |g_Pi(Pi_*) - g_-^F1| = {lower_residual}")

# (4b) SINGULAR equal-normalized branch EXCLUDED: separation matches notes Delta g_-.
sing_sep = sp.N(g_nat - gPi_at_star, 30)
delta_g_minus_notes = sp.Float("0.241964921055337", 30)   # notes stage122 Sec. 4 line 104
assert abs(sp.N(sing_sep - delta_g_minus_notes, 30)) < sp.Float("1e-12", 30), (
    f"singular-branch separation mismatch: g_nat - g_Pi(Pi_*) = {sing_sep}, "
    f"notes Delta g_- = {delta_g_minus_notes}"
)
assert sing_sep > sp.Float("1e-3", 30), (
    f"singular equal-normalized branch NOT excluded: g_nat - g_Pi(Pi_*) = {sing_sep}"
)
print(f"PASS: singular equal-normalized branch excluded — g_nat - g_Pi(Pi_*) = {sing_sep}")

# (4c) UPPER branch EXCLUDED: g_Pi(Pi_*) is far below g_+^{F1}.
upper_sep = abs(sp.N(gPi_at_star - g_plus_exact, 30))
assert upper_sep > sp.Float("1", 30), (
    f"upper branch NOT excluded: |g_Pi(Pi_*) - g_+^F1| = {upper_sep}"
)
print(f"PASS: upper branch excluded — |g_Pi(Pi_*) - g_+^F1| = {upper_sep}")
```

**Required change (Mathematica):**

Replace Mathematica lines 75–80 (the entire current "Anchor 4" block):

Before (lines 75–80):
```mathematica
(* Anchor 4: lower-branch discrimination, gPi at 2*piStar is far from gMinus. *)
offStarResidual = Abs[N[(gPi /. piM -> 2*piStar) - gMinus, 30]];
If[TrueQ[offStarResidual > 10^-3],
  pass["lower-branch discrimination (paper Checks item 3)"],
  fail["lower-branch discrimination (paper Checks item 3)", offStarResidual]
];
```

After (lines 75– …):
```mathematica
(* Anchor 4: lower-vs-singular branch discrimination, paper Checks item 3. *)
(* gPi rises from 2/Pi to a supremum of 1, so gNat = 1 is the unreachable singular *)
(* equal-normalized branch and gPlus > 1 is never attained. Wrap each Rule in parens. *)
gNat = 1;
gPlusExact = (2*Sqrt[4107 - 100*Pi^2] + 37*Sqrt[3])/(20*Pi);
gPiAtStar = N[(gPi /. piM -> piStar), 40];

(* 4a: lower-branch MEMBERSHIP. *)
lowerResidual = Abs[N[gPiAtStar - gMinus, 40]];
If[TrueQ[lowerResidual < 10^-30],
  pass["Pi_* on lower branch (membership)"],
  fail["Pi_* on lower branch (membership)", lowerResidual]
];

(* 4b: SINGULAR equal-normalized branch EXCLUDED, separation matches notes Delta g_-. *)
singSep = N[gNat - gPiAtStar, 30];
deltaGMinusNotes = 0.241964921055337`30;
If[TrueQ[Abs[N[singSep - deltaGMinusNotes, 30]] < 10^-12 && singSep > 10^-3],
  pass["singular equal-normalized branch excluded (notes Delta g_-)"],
  fail["singular equal-normalized branch excluded (notes Delta g_-)", singSep]
];

(* 4c: UPPER branch EXCLUDED. *)
upperSep = Abs[N[gPiAtStar - N[gPlusExact, 40], 30]];
If[TrueQ[upperSep > 1],
  pass["upper branch excluded"],
  fail["upper branch excluded", upperSep]
];
```

**Claim manifest:**
- M4a — `Pi_*` lies on the lower compensated branch `g_-^{F1}` (root of `g_Pi = g_-^{F1}`):
  paper Checks item 3 `paper/stages/stage_131.tex:24`; `g_-^{F1}` value
  `notes/stages/moving_throat_pde_stage122_mouth_source_compensation_test.md:53-66`.
- M4b — the SINGULAR equal-normalized branch `g_nat = 1`
  (`notes/...stage122...:24-28`, boxed) is EXCLUDED, with separation
  `Delta g_- = 1 - g_-^{F1} ≈ 0.241964921055337` (`notes/...stage122...:102-106`, boxed).
- M4c — the upper branch `g_+^{F1} ≈ 2.79795199200529`
  (`notes/...stage122...:61-66`, boxed) is EXCLUDED.

**Verification command:**
After Codex applies, the verifier runs `redteam exec-sympy 131` /
`redteam exec-mathematica 131` and confirms: (a) both exit 0; (b) each transcript contains
the THREE new PASS lines (lower membership, singular excluded, upper excluded); (c) the OLD
`gPi(2*Pi_*)` / `lower-branch discrimination (paper Checks item 3)` single PASS line is
GONE; (d) the singular-branch PASS ties to `0.241964921055337` (notes Delta g_-).

## Applied: F2

- files_changed:
  - `scripts/moving_throat_pde_stage131_parent_mouth_threshold_sympy_audit.py`
  - `mathematica/moving_throat_pde_stage131_parent_mouth_threshold_mathematica_audit.wl`
- summary: Replaced the off-star residual check with lower-branch membership, singular-branch exclusion, and upper-branch exclusion checks.
- deviation: none

---

## F3 — transliteration (R3): give the Mathematica `Pi_*` an independent derivation primitive

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage131_parent_mouth_threshold_mathematica_audit.wl:40-41` (`gPi = FullSimplify[...]; piStar = piM /. FindRoot[gPi == gMinus, {piM, 1.5}, ...]`)

**Issue:**
Codex R3: the `.wl` is a direct transliteration of the SymPy — same `g_-` closed form, same
literals, same root solve (`FindRoot[gPi == gMinus]` mirroring SymPy's
`sp.nsolve(gPi - g_minus, 1.5)`), same threshold structure, same off-star check. The
load-bearing quantity is `Pi_*`, a numerical root of the transcendental equation
`g_Pi(Pi) = g_-^{F1}`. Both engines currently solve the IDENTICAL rational equation with the
IDENTICAL seed via the IDENTICAL primitive (Newton/FindRoot from 1.5), so a wrong target
copied into both could pass.

**Mirror-policy DECISION: STRENGTHEN with a genuinely independent Mathematica route for
`Pi_*` (do NOT accept as a pure policy-mirror).**

Justification (per `notes/MATHEMATICA_MIRROR_POLICY.md`): the default "transliteration is
expected" exemption applies to CONSOLIDATION/STATUS cards whose `.wl` is pure mechanical
algebra (e.g. 117 F1, 169/175 policy-mirrors). Stage 131 is NOT such a card — its claim
status is `\StatusExactClosure` and its load-bearing payload is a NUMERICAL ROOT of a
transcendental equation, exactly the situation where the IV.5 precedent (139, 143, 144 F4)
forced an independent second-engine solver: "the second-engine value here is the
independent numerical root-finder." The fix is low-cost and directly answers R3's "same
root solve" complaint. (`Pi_*` is ALSO defended against a wrong copy by the passed Anchor 1,
which checks both engines against the independent notes literal `1.50882951349316`; the
independent route below makes the Mathematica DERIVATION itself structurally distinct, not
merely the validation.)

The route: derive `Pi_*` in Mathematica from the CLEARED-DENOMINATOR residual — a different
equation FORM (polynomial-in-(piM, Exp[piM]) rather than the rational `gPi == gMinus`) —
solved via `FindRoot` over a BRACKETING seed pair, a structurally different solver path than
SymPy's `nsolve(gPi - g_minus, 1.5)`. The root is unique on (0, ∞) (g_Pi is monotone), so
the bracket robustly isolates it. The transcendental nature precludes a closed-form
`Solve`/`Reduce`, so a different numerical formulation is the correct independent primitive.

**Required change (Mathematica):**

Replace Mathematica line 41:
```mathematica
piStar = piM /. FindRoot[gPi == gMinus, {piM, 1.5}, WorkingPrecision -> 80, AccuracyGoal -> 30, PrecisionGoal -> 30, MaxIterations -> 100];
```
with an INDEPENDENT cleared-denominator route (keep line 40 `gPi = FullSimplify[...]`
unchanged — it is reused by Anchors 2 and 4):
```mathematica
(* INDEPENDENT Pi_* route (not a transliteration of SymPy's nsolve on gPi == gMinus): *)
(* clear denominators so the root equation is polynomial-in-(piM, Exp[piM]) rather than *)
(* the rational gPi == gMinus, and isolate the unique positive root with a bracketing *)
(* seed pair. g_Pi is monotone on (0, Infinity), so the bracket robustly fixes the root. *)
gThresholdResidual[p_] := 40*Pi*p*(2*p*Exp[p] + Pi) - 20*Pi*gMinus*(4*p^2 + Pi^2)*(Exp[p] - 1);
piStar = piM /. FindRoot[gThresholdResidual[piM] == 0, {piM, 1.4, 1.6},
  WorkingPrecision -> 80, AccuracyGoal -> 30, PrecisionGoal -> 30, MaxIterations -> 100];
```

Derivation of the cleared-denominator residual (so Codex can confirm it, not copy blindly):
`gPi - gMinus = 0` with `gPi = 2*piM*(2*piM*Exp[piM] + Pi)/((4*piM^2 + Pi^2)*(Exp[piM] - 1))`.
Multiply both sides by `(4*piM^2 + Pi^2)*(Exp[piM] - 1)`:
`2*piM*(2*piM*Exp[piM] + Pi) - gMinus*(4*piM^2 + Pi^2)*(Exp[piM] - 1) = 0`.
Multiply through by `20*Pi`:
`40*Pi*piM*(2*piM*Exp[piM] + Pi) - 20*Pi*gMinus*(4*piM^2 + Pi^2)*(Exp[piM] - 1) = 0`,
which is exactly the residual coded above (root at piStar; verified numerically ≈ 0 at
piM ≈ 1.50883). The `20*Pi*gMinus` factor keeps `gMinus` symbolic so its value flows in from
F4, not a re-typed literal. (Sign corrected per Claude+Codex consult `019e7594`: an earlier
draft wrote `- (1 - Exp[p])*(...)*(20*Pi*gMinus)`, which flips the sign of the gMinus term —
that residual evaluates to ≈ 6366, NOT 0, at piStar.)

Note (Mathematica hazards): `gThresholdResidual` is a fresh symbol (not a bound name); the
`Pi` inside is Mathematica's built-in π (matches `gPi`'s `Pi`), `piM` is the variable.
Comment bodies contain no `*)` substrings. Do NOT touch line 40.

**SymPy side:** NO CHANGE. SymPy keeps its `sp.nsolve(gPi - g_minus, 1.5, ...)` (line 24).
The two engines now use structurally different `Pi_*` derivations (SymPy: Newton on the
rational `gPi - g_minus` from a single seed; Mathematica: FindRoot on the cleared-denominator
polynomial-in-Exp residual over a bracketing seed pair), and both validate against the
independent notes anchor `1.50882951349316` (Anchor 1, reconciled-PASSED).

**Claim manifest:**
- M5 — `Pi_*` is the unique positive root of `g_Pi(Pi) = g_-^{F1}`
  (`notes/stages/moving_throat_pde_stage131_parent_mouth_threshold.md:6-9`, value
  `1.50882951349316`); the cleared-denominator residual is an algebraic rearrangement of
  the same equation, solved by a structurally distinct numerical primitive.

**Verification command:**
After Codex applies, the verifier runs `redteam exec-mathematica 131` and confirms: (a) the
script exits 0; (b) `Pi_* = 1.50882951349316...` still prints (Anchor 1 still PASSES against
the notes literal); (c) a grep shows the `.wl` `Pi_*` derivation uses `gThresholdResidual`
(cleared-denominator) and a bracketing `{piM, 1.4, 1.6}` seed, NOT `FindRoot[gPi == gMinus,
{piM, 1.5}]`; (d) all Anchor 2/3/4 PASS lines (which reuse `piStar` and `gPi`) still appear.

## Applied: F3

- files_changed:
  - `mathematica/moving_throat_pde_stage131_parent_mouth_threshold_mathematica_audit.wl`
- summary: Changed the Mathematica Pi_* solve to a cleared-denominator residual with a bracketing seed pair.
- deviation: none

---

## NEEDS_APPROACH_REVIEW — RESOLVED (Claude+Codex consult `019e7594`, 2026-05-29): option (ii) DROP

**RESOLUTION:** Codex CONCUR that even the dimensionless-group round-trip is tautological
(`(Pi_* * Theta_sigma/L) * (L/Theta_sigma) ≡ Pi_*` for symbolic `Theta_sigma, L`, since the script
hardcodes both the value formula and its literal inverse — no independent group definition is imported
to falsify). Per option **(ii)**, the Anchor-3 block is DELETED in both engines (see the F1 RESOLUTION
above); the threshold identity is recorded as purely definitional and out-of-scope for a scripts-only
audit. The stage's load-bearing checks are Anchors 1/2 (reconciled-PASSED) + the F2 branch discrimination.
This was a framing/altitude call settled Claude+Codex (NOT escalated to the user; not a conceptual change).
The analysis that led here is preserved below for the record.

---

(Original analysis — superseded by the RESOLUTION above:)

De-tautologizing R1 (above) revealed that the parent-threshold identity
`T_m - q_* A_0' = Pi_* Theta_sigma / L` is **purely DEFINITIONAL**, not an independently
checkable physics relation:

- The notes DEFINE the canonical branch as `Pi_m = Pi_*` (`notes/...stage131...:30`).
- They DEFINE `Pi_m = V_1 L / Theta_sigma` (Stage 231, `:18-25`) and
  `V_1 = T_m - q_* A_0'` (`:52-58`).
- Therefore `T_m - q_* A_0' = Pi_* Theta_sigma / L` is the *definition* of "being on the
  canonical branch," and `T_m`, `q_*`, `A_0'`, `Theta_sigma`, `L` are FREE abstract symbols
  in this stage with NO independent upstream numerical values to substitute. There is
  literally nothing to derive the identity FROM — any non-trivial check would have to
  re-derive the Stage-231 group definition itself, which lives upstream and is out of scope
  for this audit (scripts only, no notes/paper edits).

**Proposed resolution (encoded in F1 above):** check the dimensionless-group ROUND-TRIP —
the boxed threshold VALUE `V_1 = Pi_* Theta_sigma/L` fed through the group definition
`Pi_m = V_1 L/Theta_sigma` must recover `Pi_m = Pi_*`. This is falsifiable against a
MISQUOTED group definition (the kept-symbolic `L`, `Theta_sigma` fail to cancel), exactly
the stage-122-F1-accepted "cancellation of a kept-symbolic constant" pattern. It is the
strongest honest check available without re-deriving Stage 231.

**However** — like stage 118's accepted-but-flagged `K_q closed form` assert (see
`notes/MATHEMATICA_MIRROR_POLICY.md`, "CAVEAT: 118's `K_q closed form` ... effectively X−X
... non-blocking because ... independently anchored upstream") — the round-trip's substance
is modest: it verifies the *dimensional inverse-consistency* of the Sec.1 group definition
with the Sec.1/2 threshold value, NOT a numerical fact about the mouth layer. The
LOAD-BEARING physics of this stage is carried by Anchors 1, 2 (PASSED) and the R2 branch
discrimination (F2 above), not by the threshold identity.

**Decision needed from orchestrator + Codex (Claude+Codex math consult, per memory
`feedback_claude_codex_resolve_math`):** EITHER
- (i) accept the round-trip as the de-tautologized F1 (proposed; consistent with the
  118 `K_q`/122-F1 precedent), recording that the threshold identity is definitional and the
  stage's load-bearing checks are Anchors 1/2 + R2; OR
- (ii) if even the round-trip is judged too thin, DROP the threshold-identity assert
  entirely (delete the Anchor-3 block in both engines) and let Anchors 1/2 + R2 carry the
  stage — rather than keep a definitional check dressed up as a verification.

This is a framing/altitude question, not a conceptual-nature change, so it should be settled
Claude+Codex and NOT escalated to the user. F2 and F3 are unaffected and may be applied
immediately. Apply F1 only after this section is marked resolved with the chosen option.

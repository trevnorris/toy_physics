---
unit_id: 137
batch: IV.4
auditor_model: claude-opus-4-7[1m]
audit_date: 2026-05-27T00:00:00Z
verdict: findings
stop_cold: null
findings_count: 4
paper_alignment: partial
scripts_checked:
  sympy: insufficient
  mathematica: insufficient
  engines_agree: true
  outputs_fresh: true
docs_read:
  paper_stage_tex: present
  notes_stage_files:
    - /var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage137_core_to_mouth_gain_map.md
  paper_appendix: present
---

# Audit unit 137 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_137.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage137_core_to_mouth_gain_map.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part04.tex` (only `\input{stages/stage_137}` and the group summary line on L30 reference this unit)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage137_core_to_mouth_gain_map_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage137_core_to_mouth_gain_map_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage137_core_to_mouth_gain_map_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage137_core_to_mouth_gain_map_mathematica_audit.txt`

## What the paper claims

The card titled "Explicit Core-to-Mouth Gain Map" delivers the explicit form of the coupled mouth-layer gain pair `(M_s, M_q)` as a function of the parent throat-core data. The card's centerpiece block (stage_137.tex L15-17) states verbatim: "Explicit gains are \(M_s=Lg_s^2/(K_s\Theta_\sigma)\), \(M_q=-L(K_sg_q-\lambda g_s)^2/[K_s(K_sK_q+\lambda^2)\Theta_\sigma\]." The supporting notes (Sections 1-2, L23-89) derive that result in two non-trivial steps: (a) varying the explicit electrochemical mouth free energy `F_mouth[sigma, U_s, U_q]` w.r.t. `sigma` yields `Pi = (L/Theta_sigma)[rho_c U_s'(0) - sigma_c U_q'(0)]`, identifying the gain coefficients with `(L/Theta_sigma) rho_c` and `-(L/Theta_sigma) sigma_c`; and (b) the Stage 97 Schur complement supplies the concrete `rho_c = g_s^2/K_s` and `sigma_c = (K_s g_q - lambda g_s)^2 / [K_s(K_s K_q + lambda^2)]`. The card's `Checks` block additionally enumerates three sub-checks (outlet-consistency of `(M_s, M_q)`, self-matched susceptibility closure before the one-scalar branch law, and recording numerically located fixed points as numerical rather than closed-form).

## What the script claims to verify

Both scripts hand-assign `rho_c`, `sigma_c`, `M_s`, `M_q` to the literal boxed expressions from the notes (sympy L9-13; mathematica L33-37). They then `print` those expressions and execute exactly one substantive assertion: that an alternate algebraic packaging of `sigma_c`, namely `(K_s g_q - lambda g_s)^2 / [K_s^2 K_q (1 + lambda^2/(K_s K_q))]`, simplifies to the previously assigned `sigma_c` (sympy L21-23; mathematica L44-46 via `expectZero`). No variational derivation from the free energy is exercised; no link to the Stage 97/99 Schur complement (`delta Lambda_core(z)`) is exercised; no outlet-consistency or self-matched susceptibility check is exercised. The scripts' bottom line is a single algebraic identity between two equal-by-construction packagings of `sigma_c`.

## Paper ↔ script cross-check

| Paper-side deliverable | Script-side check | Status |
|---|---|---|
| Closed form `M_s = L g_s^2 / (K_s Theta_sigma)` (card L16) | `Ms = sp.simplify(L*rho_c/Theta)` with `rho_c = gs**2/Ks` (sympy L9, L12) — assigned, never derived or asserted against an upstream artifact | `partial` (printed, no non-trivial check) |
| Closed form `M_q = -L (K_s g_q - lambda g_s)^2 / [K_s(K_s K_q + lambda^2) Theta_sigma]` (card L16) | `Mq = sp.simplify(-L*sigma_c/Theta)` with `sigma_c` hand-assigned (sympy L10, L13) | `partial` (printed, no non-trivial check) |
| `(rho_c, sigma_c)` come from the Stage 97 Schur complement of `delta Lambda_core(z)` (notes L57-71) | No use of the Schur complement; values hand-assigned (sympy L9-10; .wl L33-34) | `missing` |
| `Pi = (L/Theta_sigma)[rho_c U_s'(0) - sigma_c U_q'(0)]` follows from varying `F_mouth[sigma, U_s, U_q]` (notes L23-49) | No variation; no free energy ever appears in the script | `missing` |
| `M_s = (L/Theta_sigma) rho_c`, `M_q = -(L/Theta_sigma) sigma_c` (notes boxed L75-89) | Direct hand assignment matches form (sympy L12-13) | `partial` (assignment, not an assertion) |
| Card `Checks` item 1: gain pair `(M_s, M_q)` consistent with outlet | No check | `missing` |
| Card `Checks` item 2: self-matched susceptibility closure | No check | `missing` |
| Card `Checks` item 3: numerically located fixed points recorded as numerical | No check (no numerics computed) | `missing` |

`paper_alignment` is set to `partial`: the script's printed expressions textually match the paper's closed form, but no script-side assertion exercises the derivation chain (free-energy variation + Schur complement) the paper card and notes invoke. The card's three enumerated `Checks` are entirely absent on the script side.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 23 | `assert sp.simplify(sigma_c - expr_rc) == 0` where `expr_rc = (K_s g_q - lam g_s)^2 / [K_s^2 K_q (1 + lam^2/(K_s K_q))]` | None — both sides are equal-by-construction algebraic packagings of `sigma_c`; this never touches `M_s`, `M_q`, the free energy, the Schur complement, or any of the card's `Checks` items | no (tautological algebraic rearrangement) |
| A2 | mathematica | 46 | `expectZero["sigma_c equivalence with r_c form", sigmaC - sigmaCRc]` (same as A1 in Mathematica syntax) | None — same tautological rearrangement | no |

Both scripts contain exactly one assertion each. Neither anchors any of the paper's eight enumerated deliverables. There is no second engine in any meaningful sense (see F2).

## Findings

### F1 — insufficient_verification

**Severity:** high
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage137_core_to_mouth_gain_map_sympy_audit.py:9-23`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage137_core_to_mouth_gain_map_mathematica_audit.wl:33-46`

**What's wrong:**
The script's only assertion is `assert sp.simplify(sigma_c - expr_rc) == 0` (sympy L23), where `sigma_c` and `expr_rc` are both hand-built from the same numerator `(K_s g_q - lam g_s)^2` with denominators `K_s(K_s K_q + lam^2)` and `K_s^2 K_q (1 + lam^2/(K_s K_q))`. These two denominators are algebraically identical by distributivity (`K_s^2 K_q + K_s lam^2 = K_s(K_s K_q + lam^2)`), so the assertion cannot fail and exercises no physics. The card's actual deliverable (`M_s = L g_s^2/(K_s Theta_sigma)`, `M_q = -L(K_s g_q - lam g_s)^2/[K_s(K_s K_q + lam^2) Theta_sigma]`, stage_137.tex L16) is only `print`-ed, never asserted against anything independent. The card's enumerated `Checks` (stage_137.tex L21-25: outlet consistency of `(M_s, M_q)`, self-matched susceptibility closure, numerical recording of fixed points) are completely absent.

**Why this matters:**
A regression that broke the `M_s` or `M_q` formula on the script side would still PASS the only assertion in the script, because `sigma_c - expr_rc == 0` survives any change to `M_s`, `M_q`, `rho_c`, `Theta`, or `L`. The script gives a PASS that does not certify the card's closed form. The card's `Checks` block explicitly states three sub-checks the audit should exercise; none of them appear.

**Required change:**
Add substantive assertions in both engines that exercise the paper's gain map without trivializing into a rearrangement of the same hand-assigned expression. Concretely:
- (S1) Derive `M_s` and `M_q` from a separate symbolic computation `(L/Theta) * rho_c` and `-(L/Theta) * sigma_c` (do not pre-multiply into the answer), then assert `simplify(M_s - L*g_s**2/(K_s*Theta)) == 0` and `simplify(M_q + L*(K_s*g_q - lam*g_s)**2/(K_s*(K_s*K_q + lam**2)*Theta)) == 0` against the literal forms quoted in stage_137.tex L16.
- (S2) Anchor `rho_c, sigma_c` to the Schur complement structure described in notes L57-71: build `delta_Lambda_core(z) = rho_c - sigma_c/(1 - kappa_c z^2 - I*gamma_c z^5)` as a sympy expression, then assert that the residue / leading static term at `z=0` reproduces `rho_c - sigma_c` and that the `O(z^0)` Laurent extraction yields the hand-assigned `rho_c, sigma_c`. (At minimum, assert the static limit `z -> 0` returns `rho_c - sigma_c`.)
- (S3) Optionally add the free-energy variation check from notes L23-49: define `F_mouth` symbolically, take `diff(F_mouth, sigma)`, and confirm the stationarity equation gives `Pi = (L/Theta)[rho_c U_s'(0) - sigma_c U_q'(0)]` so that `M_s = (L/Theta) rho_c` and `M_q = -(L/Theta) sigma_c` fall out by inspection.

Mirror in the `.wl` script using `FullSimplify`/`expectZero`.

**Verification:**
After the fix, the SymPy script must contain at least three independent assertions (one each for `M_s`, `M_q`, and the Schur or free-energy anchor) where the LHS of each `simplify(... - ...) == 0` is built from a different starting expression than the RHS. The Mathematica script must mirror those assertions via `expectZero`. Output transcripts must show three (or more) PASS lines, not just a single equivalence line.

### F2 — mathematica_transliteration

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage137_core_to_mouth_gain_map_mathematica_audit.wl:33-46`
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage137_core_to_mouth_gain_map_sympy_audit.py:9-23`

**What's wrong:**
The `.wl` file is a line-by-line port of the `.py` file:

| sympy L9-10 | mathematica L33-34 |
|---|---|
| `rho_c = sp.simplify(gs**2 / Ks)` | `rhoC = FullSimplify[gS^2/kS, ...]` |
| `sigma_c = sp.simplify((Ks*gq - lam*gs)**2 / (Ks*(Ks*Kq + lam**2)))` | `sigmaC = FullSimplify[(kS*gQ - lam*gS)^2/(kS*(kS*kQ + lam^2)), ...]` |

| sympy L12-13 | mathematica L36-37 |
|---|---|
| `Ms = sp.simplify(L * rho_c / Theta)` | `mS = FullSimplify[lM*rhoC/thetaSigma, ...]` |
| `Mq = sp.simplify(-L * sigma_c / Theta)` | `mQ = FullSimplify[-lM*sigmaC/thetaSigma, ...]` |

| sympy L21-23 | mathematica L44-46 |
|---|---|
| `expr_rc = sp.simplify((Ks*gq - lam*gs)**2 / (Ks**2*Kq*(1 + lam**2/(Ks*Kq))))` | `sigmaCRc = FullSimplify[(kS*gQ - lam*gS)^2/(kS^2*kQ*(1 + lam^2/(kS*kQ))), ...]` |
| `assert sp.simplify(sigma_c - expr_rc) == 0` | `expectZero["sigma_c equivalence with r_c form", sigmaC - sigmaCRc]` |

Every variable choreography, every intermediate, every numerator/denominator packaging is identical. The `.wl` is not an independent derivation — it is the same algebra rewritten in Mathematica syntax. This violates the second-engine policy that both engines must derive the result independently from the physical premises.

**Why this matters:**
A bug in the algebraic packaging (sign flip, missing factor) would propagate to both engines identically, so the "engine cross-check" provides no independent confirmation. The whole point of dual-engine audits is to catch authorial transcription errors; a transliterated mirror catches nothing.

**Required change:**
Once F1 lands, ensure the Mathematica re-derivation uses different algebraic packaging than the SymPy one. Concretely, in the Schur-complement anchor (S2 above), derive `sigma_c` in the `.wl` by computing the residue of `delta_Lambda_core(z)` at the appropriate pole via `Residue[...]` or by series-expanding via `Series[...]`, rather than rebuilding the same `(kS*gQ - lam*gS)^2/(kS*(kS*kQ + lam^2))` symbol product. Mathematica should arrive at the same closed form by a route the SymPy script does not use.

**Verification:**
The verifier compares the two scripts side-by-side after F1+F2 land. At least one of the three new assertions in the `.wl` must construct its LHS by an algebraic route different from the `.py` script's corresponding LHS (e.g., `Residue`, `SeriesCoefficient`, or `Solve` of the variational equation, rather than direct simplification of the same product).

### F3 — hardcoded_result

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage137_core_to_mouth_gain_map_sympy_audit.py:9-13`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage137_core_to_mouth_gain_map_mathematica_audit.wl:33-37`

**What's wrong:**
`rho_c`, `sigma_c`, `M_s`, `M_q` are each written down as the answer at script L9-13 (sympy) / L33-37 (.wl). There is no in-script derivation from the upstream Schur complement of Stages 97-99 (notes L57-71) or from the free-energy variation (notes L23-49). The boxed final result of the notes is the input to the script. The single assertion (A1/A2) only rearranges `sigma_c` algebraically, so the "result" never has to survive a non-trivial computation.

**Why this matters:**
A hardcoded result with only a tautological rearrangement check provides no defense against an authorial typo when the values were copied in. If the upstream Schur complement of Stage 97 yields a different `sigma_c` (e.g., due to a sign convention change), this script would still PASS because nothing connects it to the upstream artifact.

**Required change:**
Replace the hand assignment of `rho_c, sigma_c` with a derivation step. Two acceptable routes:
- (a) Define `delta_Lambda_core(z) = rho_c_sym - sigma_c_sym/(1 - kappa_c*z**2 - sp.I*gamma_c*z**5)` symbolically using independent `rho_c_sym, sigma_c_sym` symbols, then match the leading static / static-limit Laurent extraction against the closed forms `g_s**2/K_s` and `(K_s*g_q - lam*g_s)**2/(K_s*(K_s*K_q + lam**2))` via `simplify(... - ...) == 0`.
- (b) Build the explicit two-channel core stiffness matrix `[[K_s, lam], [lam, K_q]]` with mouth coupling vector `[g_s, g_q]`, compute its Schur complement / inverse-coupling-quadratic-form, and assert that the resulting `sigma_c` equals `(K_s*g_q - lam*g_s)**2/(K_s*(K_s*K_q + lam**2))`. Similarly for `rho_c`.

Either route makes the assertion non-trivial: it can fail if `rho_c` or `sigma_c` are misquoted.

**Verification:**
The verifier inspects the SymPy script and confirms that `rho_c` and `sigma_c` no longer appear as direct expression literals on a single RHS; instead they are extracted via `Series`, `limit`, `solve`, or matrix inversion of an independently posed object. The Mathematica script must do the same via an independent route (per F2).

### F4 — paper_misalignment

**Subtype:** notes_contradicts_script

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage137_core_to_mouth_gain_map_mathematica_audit.wl:26`

**What's wrong:**
The Mathematica script's banner reads `banner["STAGE 120 — EXPLICIT CORE-TO-MOUTH GAIN MAP"]` (L26), but this is unit 137 (paper card title: "Stage 137: Explicit Core-to-Mouth Gain Map"). The stage number in the banner is stale (120 vs 137).

The Mathematica transcript also records this banner verbatim: `STAGE 120 — EXPLICIT CORE-TO-MOUTH GAIN MAP` (mathematica_audit.txt L11).

This is a low-severity script-side label bug, not a paper-side dispute — the paper card is correct and the script's substantive content matches the paper's topic. It is filed under `paper_misalignment` only because the script's self-labeling contradicts the unit number the paper assigns. No user resolution is needed; the fix direction is unambiguous (correct the script banner to "STAGE 137").

**Why this matters:**
Mislabeled banners produce stale output transcripts that the orchestrator's downstream pattern-matching (e.g., "find the stage 137 transcript") may use to attribute results to the wrong unit. Low impact, but trivially fixable.

**Required change:**
At `mathematica/moving_throat_pde_stage137_core_to_mouth_gain_map_mathematica_audit.wl:26`, change `banner["STAGE 120 — EXPLICIT CORE-TO-MOUTH GAIN MAP"];` to `banner["STAGE 137 — EXPLICIT CORE-TO-MOUTH GAIN MAP"];`.

**Verification:**
After re-running the `.wl`, the transcript at `mathematica/output/moving_throat_pde_stage137_core_to_mouth_gain_map_mathematica_audit.txt:11` reads `STAGE 137 — EXPLICIT CORE-TO-MOUTH GAIN MAP`.

## Independent-derivation check (Mathematica)

The Mathematica script is a transliteration of the SymPy script (see F2 for the line-by-line correspondence). The variable choreography is identical: declare ρ_c, σ_c with identical numerators/denominators; build M_s, M_q by the same multiplication; assert the same algebraic rearrangement of σ_c. Mathematica idioms (`FullSimplify`, `expectZero`, `Assumptions`) are used, but the algebra never deviates from the Python path. No independent derivation (e.g., `Residue`, `Series`, `Solve` of a stationarity condition) appears.

## Engine cross-check

Both transcripts produce identical closed forms (modulo variable naming `Ks` vs `kS`, etc.):

- sympy_audit.txt L9-13:
  - `rho_c = g_s**2/K_s`
  - `sigma_c = (K_s*g_q - g_s*lam)**2/(K_s*(K_q*K_s + lam**2))`
  - `M_s = L*g_s**2/(K_s*Theta)`
  - `M_q = -L*(K_s*g_q - g_s*lam)**2/(K_s*Theta*(K_q*K_s + lam**2))`
  - `sigma_c (r_c form) = (K_s*g_q - g_s*lam)**2/(K_s*(K_q*K_s + lam**2))`
- mathematica_audit.txt L13-18:
  - `rho_c = gS^2/kS`
  - `sigma_c = (gQ*kS - gS*lam)^2/(kS*(kQ*kS + lam^2))`
  - `M_s = (gS^2*lM)/(kS*thetaSigma)`
  - `M_q = -(((gQ*kS - gS*lam)^2*lM)/(kS*(kQ*kS + lam^2)*thetaSigma))`
  - `sigma_c (r_c form) = (gQ*kS - gS*lam)^2/(kS*(kQ*kS + lam^2))`
  - `sigma_c equivalence with r_c form = 0` -> `PASS`

The outputs agree — but the agreement is trivial because the two engines run the same algebra. This is not a genuine cross-check (see F2).

## Verdict justification

The scripts trivially pass on a hand-assigned identity that exercises none of the paper card's substantive content. The printed closed forms textually match the card, but no script-side assertion would fail if `M_s`, `M_q`, `rho_c`, `sigma_c`, `L`, or `Theta_sigma` were mistyped — the only assertion rearranges `sigma_c` with itself. The Mathematica script is a transliteration of the Python script, so the engine cross-check is no cross-check. The paper card additionally enumerates three explicit `Checks` (outlet consistency, self-matched susceptibility closure, numerical fixed-point status); none of those appear in either script. The banner in the Mathematica file is mislabeled (`STAGE 120`). No `stop_cold` is warranted: the math is not internally inconsistent, and the boxed result the paper claims is in fact algebraically correct (it factors `K_s K_q + lam^2` consistently); the failure is in verification depth, not in the underlying claim. Verdict: `findings`.

## Self-test notes

- Variable independence: no `sp.diff` / `D[]` calls in either script, so the "differentiation against an inert variable" trap does not apply here.
- Symmetry/parity: no integrals; trap not applicable.
- Trivial-case pre-check: the proposed S1 fix would assert `simplify(M_s - L*g_s**2/(K_s*Theta)) == 0` — mentally, with `M_s` defined as `L*rho_c/Theta` and `rho_c = g_s**2/K_s`, the difference is `L*g_s**2/(K_s*Theta) - L*g_s**2/(K_s*Theta) = 0`; non-trivially failable only if the LHS is changed (e.g., `M_q` swapped). For S2 (Schur anchor), at `z -> 0` the Lorentzian factor `1/(1 - kappa_c z^2 - I gamma_c z^5)` reduces to `1`, so `delta_Lambda_core(0) = rho_c - sigma_c`, matching the hand-assigned `rho_c - sigma_c` — assertion can fail if the upstream Schur form is misquoted.
- Path specifications: F1 and F3 edit existing files at named line ranges; F2 modifies the existing `.wl` at L33-46; F4 edits the existing `.wl` at L26. No missing-script directive shape needed.
- Paper round-trip: the proposed S1 assertions quote the exact literals from stage_137.tex L16; no new constant or symbol is introduced; no new paper_misalignment risk.

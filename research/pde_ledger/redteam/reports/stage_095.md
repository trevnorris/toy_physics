---
unit_id: 095
batch: IV.1
auditor_model: claude-opus-4-7[1m]
audit_date: 2026-05-27T00:00:00Z
verdict: findings
stop_cold: null
findings_count: 3
paper_alignment: partial
scripts_checked:
  sympy: present
  mathematica: present
  engines_agree: true
  outputs_fresh: true
docs_read:
  paper_stage_tex: present
  notes_stage_files: [moving_throat_pde_stage095_second_order_geometry_contamination.md]
  paper_appendix: present
---

# Audit unit 095 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_095.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage095_second_order_geometry_contamination.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part04.tex` (only the `\input{stages/stage_095}` line at 1224 references this unit; the file is otherwise generic part-04 framing)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage095_second_order_geometry_contamination_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage095_second_order_geometry_contamination_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage095_second_order_geometry_contamination_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage095_second_order_geometry_contamination_mathematica_audit.txt`

## What the paper claims

The stage card (`stage_095.tex:15-17`) quotes the bottom-line claim verbatim:

> Weak \(l=0\leftrightarrow l=2\) mixing contaminates the grouped module only at \(O(\chi^2)\).

The notes (`notes/stages/moving_throat_pde_stage095_second_order_geometry_contamination.md:23-44`) make the claim quantitative by starting from the bilinear action `L = (1/2) q D_q q + (1/2) g D_g g + chi M_0 q g + J q` and *deriving* the Schur-complement effective kernel `D_eff(omega) = D_q(omega) - chi^2 M_0^2 / D_g(omega)`. From there the low-frequency expansion yields `K_(g,2)^eff = chi^2 M_0^2 G_2 / G_0^2` and `K_(g,4)^eff = chi^2 M_0^2 (G_0 G_4 - G_2^2) / G_0^3`, hence `eps_2, eps_4 = O(chi^2)` and `c_pole = 1/4 [1 + eps_4 - 2 eps_2 + O(chi^4)]`. The card's `\stagefield{Checks}` (lines 21-25) lists three deliverables: (i) static limit `eps_2 = eps_4 = 0` returns `c_pole = 1/4`; (ii) `l=0`/`l=2` orthogonality before applying the firewall; (iii) the support/source statement still carries the minimal-module hypothesis (a provenance tag, not a math check).

## What the script claims to verify

The SymPy script (`scripts/moving_throat_pde_stage095_second_order_geometry_contamination_sympy_audit.py`) starts directly from `Dg = G0 + G2*w**2 + G4*w**4` and writes the Schur-complement correction as `corr = series(-chi**2 * M0**2 / Dg, w, 0, 6)`. It then extracts `K0corr, K2corr, K4corr` (the w^0, w^2, w^4 coefficients), asserts each carries a factor `chi**2`, builds `eps2, eps4, cpole`, and asserts (line 44) that `delta = cpole - 1/4` vanishes at `chi -> 0`, plus (line 51) that `d c_pole / d chi |_{chi=0} == 0`. The closing print (line 53) declares "contamination begins at O(chi^2)." The Mathematica script (`mathematica/.../audit.wl`) mirrors this structure but adds three closed-form residual checks (lines 44-46) confirming the series coefficients equal the notes' closed forms (`K2corr/chi^2 == G2*M0^2/G0^2`, etc.). Both scripts therefore presume the Schur form of `D_eff` rather than derive it from the L action in the notes.

## Paper ↔ script cross-check

| Paper-side deliverable | Script-side check | Status |
|---|---|---|
| Bilinear action `L = (1/2) q D_q q + (1/2) g D_g g + chi M_0 q g` → Schur `D_eff = D_q - chi^2 M_0^2 / D_g` (notes §1) | Neither script derives `D_eff` from the L action; both write `-chi^2 M_0^2 / D_g` ab initio | `missing` (derivation skipped) |
| Closed forms `K_(g,2)^eff = chi^2 M_0^2 G_2 / G_0^2`, `K_(g,4)^eff = chi^2 M_0^2 (G_0 G_4 - G_2^2)/G_0^3` (notes §2) | sympy prints series coefficients but never asserts equality with these closed forms; mathematica `expectZero` lines 44-46 do | `partial` (Mathematica only) |
| `eps_2 = eps_4 = O(chi^2)` (notes §2; card quote line 16) | Built into starting expression by writing `-chi^2 M_0^2 / D_g`; the `has(chi**2)` asserts (sympy lines 27-29) and `d c_pole / dchi |0 == 0` (sympy line 51; mathematica line 59) follow tautologically from that construction | `mismatch`/tautological — see F1 |
| Static-limit check `eps_2=eps_4=0 ⇒ c_pole=1/4` (card line 22) | sympy line 44 asserts `delta.subs(chi,0) == 0`; mathematica only prints `cPoleSeries` but does not assert the static limit explicitly | `partial` (Mathematica missing explicit assert) |
| `l=0`/`l=2` orthogonality before firewall (card line 23) | Neither script tests orthogonality of the two modes; orthogonality is implicit in the diagonal kinetic/dispersive sector but never exercised | `missing` |
| Minimal-module hypothesis carried with status tag (card line 24) | Not a math claim; provenance/policy item | n/a |

Dominant pattern: `partial`. Paper-alignment field set to `partial`.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 27 | `assert sp.factor(K0corr).has(chi**2)` | static-renorm part of contamination | no — tautological (chi^2 is a factor of the input by construction) |
| A2 | sympy | 28 | `assert sp.factor(K2corr).has(chi**2)` | `eps_2 = O(chi^2)` | no — tautological |
| A3 | sympy | 29 | `assert sp.factor(K4corr).has(chi**2)` | `eps_4 = O(chi^2)` | no — tautological |
| A4 | sympy | 44 | `assert sp.simplify(delta.subs(chi, 0)) == 0` | static limit `c_pole = 1/4` | partial — passes by construction (eps2,eps4 each carry `chi^2`) but does check the static limit |
| A5 | sympy | 51 | `assert chi1 == 0` (i.e. `d c_pole/dchi |0 = 0`) | first nonzero term is even in chi | no — tautological (cpole is even in chi by construction) |
| M1 | mathematica | 44 | `expectZero["K0corr / chi^2 static factor", k0Corr/chi^2 + m0^2/g0]` | closed form of `K_(g,0)^eff` | yes — substantive residual against the closed form |
| M2 | mathematica | 45 | `expectZero["K2corr / chi^2 dynamic factor", k2Corr/chi^2 - g2*m0^2/g0^2]` | closed form of `K_(g,2)^eff` | yes — substantive |
| M3 | mathematica | 46 | `expectZero["K4corr / chi^2 dynamic factor", k4Corr/chi^2 - m0^2*(g0*g4 - g2^2)/g0^3]` | closed form of `K_(g,4)^eff` | yes — substantive |
| M4 | mathematica | 59 | `expectZero["d c_pole / dchi |0", D[cPole,chi]/.chi->0]` | first nonzero term is even in chi | no — tautological |

The Mathematica script carries three substantive residual checks (M1-M3) that the SymPy script does not mirror. The SymPy script has zero substantive assertions: every assert is true by construction.

## Findings

### F1 — insufficient_verification

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage095_second_order_geometry_contamination_sympy_audit.py:27-29,44,51`

**What's wrong:**
Every SymPy assertion is tautological given how the script *constructs* its starting expression. The script never derives the Schur form `D_eff = D_q - chi^2 M_0^2 / D_g` from the L action of the notes; it simply writes `corr = series(-chi**2 * M0**2 / Dg, w, 0, 6)`. From that starting point:
- Lines 27-29 `assert sp.factor(K{0,2,4}corr).has(chi**2)`: every term in `corr` was multiplied by `chi**2` in the very first line, so every coefficient inherits the factor. The asserts cannot fail.
- Line 44 `assert sp.simplify(delta.subs(chi, 0)) == 0`: `delta = cpole - 1/4` where `cpole = (1 + eps4)/(4*(1+eps2)^2)` and `eps2, eps4` each carry a literal `chi**2` factor, so `chi=0` makes them vanish and reduces `cpole` to `1/4` algebraically. The assert cannot fail.
- Line 51 `assert chi1 == 0`: `cpole` depends on `chi` only via `chi**2`, so it is an even function of `chi` and its odd-order derivatives at `chi=0` are zero by parity. The assert cannot fail.

The Mathematica mirror (lines 44-46) adds three non-tautological residual checks against the notes' closed forms (`K2corr/chi^2 - g2*m0^2/g0^2`, etc.). These are substantive — if SymPy's `series(...).coeff(w, 2)` had a bug, these checks would catch it. The SymPy script has no such check.

**Why this matters:**
The stage's stated claim (paper quote: "Weak l=0↔l=2 mixing contaminates the grouped module only at O(chi^2)") is built into the starting expression `-chi^2 M_0^2 / D_g`. The scripts thus verify that the construction is consistent with itself, not that the construction is correct. The Mathematica residual checks (M1-M3) anchor the coefficients to the explicit closed forms from notes §2; the SymPy script does not do this, so a regression in SymPy's series expansion would slip through.

**Required change:**
Mirror the Mathematica closed-form residual checks in SymPy. After lines 20-25 (where `K0corr, K2corr, K4corr` are computed), add three asserts:
```python
assert sp.simplify(K0corr - (-M0**2*chi**2/G0)) == 0
assert sp.simplify(K2corr - (G2*M0**2*chi**2/G0**2)) == 0
assert sp.simplify(K4corr - (M0**2*chi**2*(G0*G4 - G2**2)/G0**3)) == 0
```
These are the closed forms quoted verbatim in the notes (`notes/.../stage095_*.md:74-76` and the static renorm line 65). They are not tautological because they pin the series-extracted coefficients against an independent algebraic answer.

**Verification:**
After the patch, `redteam exec-sympy 095` should still exit 0 (the closed forms match the series), and the saved output should now print at least three new lines confirming the coefficient equalities. If SymPy's series expansion were to drift, these asserts would fail.

### F2 — mathematica_transliteration

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage095_second_order_geometry_contamination_sympy_audit.py:12-53`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage095_second_order_geometry_contamination_mathematica_audit.wl:28-61`

**What's wrong:**
The two scripts follow the same algebraic choreography step-for-step:

SymPy lines 16-18:
```
Dg = G0 + G2*w**2 + G4*w**4
corr = sp.expand(sp.series(-chi**2 * M0**2 / Dg, w, 0, 6).removeO())
```
Mathematica lines 33-34:
```
dG = g0 + g2*w^2 + g4*w^4;
corr = Expand[Normal[Series[-chi^2*m0^2/dG, {w, 0, 5}]]];
```

SymPy lines 20-22 use `coeff(w, 0/2/4)`; Mathematica lines 35-37 use `Coefficient[corr, w, 0/2/4]`. SymPy lines 32-33 build `eps2 = OQ**2 * K2corr / Kpole`; Mathematica lines 48-49 build `eps2 = oQ^2*k2Corr/kPole`. SymPy line 38 builds `cpole = (1 + eps4)/(4*(1+eps2)**2)`; Mathematica line 53 the identical form. There is no independent derivation in Mathematica from the L action; both engines compute the same `D_eff` via the same Schur formula and the same series.

Mitigating note: the Mathematica script does add closed-form residual checks (lines 44-46) that the SymPy script lacks, so the .wl is *not* purely a transliteration — it does extra substantive work the .py does not. But the load-bearing derivation path is identical.

**Why this matters:**
A genuine second-engine cross-check would derive `D_eff` from the L action by integrating out `g` (e.g., using Mathematica's `Eliminate` or explicit Lagrangian completion), not by writing the Schur complement by hand and running the same series. Both engines could share a common conceptual bug in the Schur step and neither would catch it.

**Required change:**
Add to the Mathematica script (after line 31, before line 33) an independent derivation of the Schur complement from the bilinear action. One acceptable shape: define the quadratic form `Q = (1/2)*q^2*dQ + (1/2)*g^2*dG + chi*m0*q*g` (treating `dQ, dG` as scalars for fixed `omega`), apply `Solve[D[Q, g] == 0, g]` to integrate out `g`, substitute back, collect the `q^2` coefficient, and identify it with `(1/2)*(dQ + corrDerived)`. Assert (via `expectZero`) that `corrDerived` equals `-chi^2*m0^2/dG`. Only then proceed to the existing series expansion of `corrDerived`.

**Verification:**
The Mathematica saved output should gain a `PASS: D_eff Schur derivation` line that did not previously exist. The SymPy script remains the algebraic check; the Mathematica script becomes the independent derivation.

### F3 — paper_misalignment

**Subtype:** script_missing_paper_claim
**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_095.tex:23` quote: "Check `l=0` and `l=2` orthogonality before applying the geometry firewall."
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage095_second_order_geometry_contamination_sympy_audit.py` (no corresponding check)
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage095_second_order_geometry_contamination_mathematica_audit.wl` (no corresponding check)

**What's wrong:**
The stage card's `\stagefield{Checks}` lists three deliverables; the second one — orthogonality of the `l=0` and `l=2` modes before applying the firewall — has no counterpart in either script. The L action (notes §1) assumes diagonal kinetic/dispersive sectors and writes the cross-coupling as `chi M_0 q g`; orthogonality is implicit in that diagonality but never verified. The script writes a 2x2 dispersion with diagonal `D_q, D_g` and an off-diagonal `chi M_0`, which assumes (rather than checks) orthogonality.

Whether this should be resolved by adding a check or by removing the bullet from the card is the user's call. The notes do not describe a concrete orthogonality test for the reduced-model setup; this bullet may be a carry-over from a prior, larger-mode-basis version of the stage.

**Why this matters:**
A paper-side deliverable with no script-side check is a documentation gap either way (script under-verifies or card over-promises). Since the notes do not state a concrete orthogonality test that maps onto the present reduced 2-mode model, both directions of resolution are plausible.

**Required change:**
See `## Resolve before fix_loop` in the directive.

**Verification:**
After user resolution: either (a) a new orthogonality assertion appears in both scripts and the output shows it passing, or (b) the bullet at `stage_095.tex:23` is removed/clarified.

## Independent-derivation check (Mathematica)

The Mathematica script is structurally a line-by-line port of the SymPy choreography for the central derivation: identical `Dg = g0 + g2 w^2 + g4 w^4`, identical `Series[-chi^2*m0^2/dG, ...]`, identical coefficient extraction, identical `eps2/eps4` construction, identical `cPole` formula. Sample correspondences quoted above in F2. The Mathematica script does, however, add three closed-form residual checks (`expectZero` lines 44-46) that the SymPy script lacks; these are independently authored, not transliterated. Net: the load-bearing derivation is a transliteration (F2), but the script is not *purely* a port.

## Engine cross-check

Both engines agree on the final symbolic forms:

SymPy output:
```
K0corr = -M0**2*chi**2/G0
K2corr = G2*M0**2*chi**2/G0**2
K4corr = M0**2*chi**2*(G0*G4 - G2**2)/G0**3
d^2 c_pole / dchi^2 |0 = M0**2*OQ**2*(-G0*G2 + OQ**2*(G0*G4 - G2**2)/2)/(G0**3*Kpole)
```

Mathematica output:
```
K0corr = -((chi^2*m0^2)/g0)
K2corr = (chi^2*g2*m0^2)/g0^2
K4corr = (chi^2*(-g2^2 + g0*g4)*m0^2)/g0^3
d^2 c_pole / dchi^2 |0 = -1/2*(m0^2*(2*g0*g2*oQ^2 + (g2^2 - g0*g4)*oQ^4))/(g0^3*kPole)
```

The two `d^2 c_pole / dchi^2 |0` expressions are equal after common-denominator collection:
- SymPy: `M0^2*OQ^2*(-G0*G2 + OQ^2*(G0*G4 - G2^2)/2)/(G0^3*Kpole)`
- Multiply numerator: `M0^2*OQ^2*(-G0*G2 + (OQ^2*G0*G4 - OQ^2*G2^2)/2)/(G0^3*Kpole)`
- Mathematica: `-(m0^2*(2*g0*g2*oQ^2 + (g2^2 - g0*g4)*oQ^4))/(2 g0^3 kPole)` = `m0^2*(-2 g0 g2 oQ^2 + (g0 g4 - g2^2) oQ^4)/(2 g0^3 kPole)`

These match. Engines agree.

## Verdict justification

The paper's bottom-line claim (contamination first appears at `O(chi^2)`) is *built into* the starting expression `-chi^2 M_0^2 / D_g` rather than derived from the L action stated in notes §1. The SymPy script's five asserts are all tautological under that construction; the Mathematica script adds three substantive residual checks (M1-M3) but still does not derive `D_eff` from the action. The stage is therefore under-verified — not wrong, but the verifications do not actually exercise the claim, only its consequences after the claim has been assumed. F1 patches the SymPy gap by mirroring the Mathematica closed-form residual checks. F2 raises the bar on the Mathematica script to derive `D_eff` independently. F3 routes the `l=0`/`l=2` orthogonality bullet to the user for resolution. Engines agree on all final forms; outputs are fresh; no stop-cold.

## Self-test notes

- Variable-independence trap (sympy d/dchi at chi=0): `cpole` depends on chi only through `chi**2`, so `d/dchi` at chi=0 vanishes by parity — confirmed; the assert at line 51 cannot fail. Flagged as tautological in F1.
- Symmetry/parity (series in chi): `corr` is even in chi (only `chi**2`), so `cpole` is even in chi; this fact is what makes A5/M4 tautological, not a math error.
- Trivial-case substitution: For F1's proposed asserts, set `chi=1, M0=1, G0=1, G2=0, G4=0`: then `Dg=1`, `corr = -1`, so `K0corr = -1`, `K2corr = 0`, `K4corr = 0`. Closed forms give `-M0^2 chi^2/G0 = -1` ✓, `G2 M0^2 chi^2/G0^2 = 0` ✓, `M0^2 chi^2 (G0 G4 - G2^2)/G0^3 = 0` ✓. The proposed asserts pass on this concrete profile and reduce to nontrivial identities for general symbolic G2, G4. Non-tautological.
- Paper round-trip: F1's proposed closed forms are quoted verbatim from notes lines 65, 74, 76; F2's proposed derivation step uses only quantities already named in the script and in notes §1; F3 is a user-resolution finding so no script edit is prescribed. No new paper_misalignment is introduced.

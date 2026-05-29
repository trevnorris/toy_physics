---
unit_id: 116
batch: IV.3
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-05-28T00:00:00Z
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
  notes_stage_files: [moving_throat_pde_stage116_dn_mixed_tube_realization.md]
  paper_appendix: present
---

# Audit unit 116 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_116.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage116_dn_mixed_tube_realization.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part04.tex` (subsec "Core-balance theorem and D/N tube", lines 511-573; eq. `app-part04-LW-compensation` line 544-549)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage116_dn_mixed_tube_realization_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage116_dn_mixed_tube_realization_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage116_dn_mixed_tube_realization_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage116_dn_mixed_tube_realization_mathematica_audit.txt`

## What the paper claims

The stage gives a concrete geometric realization of the abstract mixed-channel data `(kappa_0, gamma_0)`. The mixed side-channel lives on a finite auxiliary tube of length `L_W` carrying the first Dirichlet/Neumann half-wave of `q'' + k^2 q = 0` with `q(0)=0`, `q'(L_W)=0`, giving `k_W = pi/(2 L_W)` and `Omega_W = pi c_s/(2 L_W)`. The bare even coefficient is `kappa_0 = (omega^2/Omega_W^2)/z^2 = 4 L_W^2/(pi^2 a^2)` with `z = a*omega/c_s`. The Stage-98 compensation requirement `kappa_0 = (1+r_c)/3`, `r_c = lambda^2/(K_s K_q)`, then fixes the tube length (the card's boxed Output: `L_W = pi a sqrt((1+r_c)/3)/2`; appendix eq. line 547). The compensation also requires `gamma_0 = (1+r_c)/9`, realized by the pure-scale bare outgoing form `D_W^bare(z) = (1+r_c)(1 - z^2/3 - i z^5/9) + O(z^6)`; removing the `(1+r_c)` hybridization factor (Stage-97 renormalization) leaves the canonical coefficients `kappa_c = 1/3`, `gamma_c = 1/9`. The notes flag the `D_W^bare` form explicitly as "a simple concrete realization", i.e. a posited ansatz, not an in-stage derivation.

## What the script claims to verify

Both scripts assert: (i) `sin(k x)` solves the ODE and `q(0)=0`; (ii) the Neumann BC `q'(L_W)=0` holds at `k = pi/(2 L_W)` (genuine eigenvalue check); (iii) the dimensionless combination `(omega/Omega_W)^2/z^2` collapses to `4 L_W^2/(pi^2 a^2)`; (iv) solving `4 L_W^2/(pi^2 a^2) = (1+r_c)/3` for `L_W` gives the boxed closed form `pi a sqrt((1+r_c)/3)/2`; (v) a battery of follow-on checks: round-tripping the geometric `kappa_0` back through the solved `L_W`, dividing `gamma_0=(1+r_c)/9` by `(1+r_c)` to get `1/9`, dividing `kappa_0=(1+r_c)/3` by `(1+r_c)` to get `1/3`, and reading the `z^2`/`z^5` coefficients back off the hand-built `D_bare`. The genuine content is (ii)-(iv); checks in (v) are constructed to be unfalsifiable.

## Paper ↔ script cross-check

| Paper deliverable | Script check | Status |
|---|---|---|
| D/N eigenvalue `k_W = pi/(2 L_W)` | sympy 34-39 / wl 43-46 (Neumann BC at k_W) | match (genuine) |
| `kappa_0 = (omega/Omega_W)^2/z^2 = 4 L_W^2/(pi^2 a^2)` | sympy 41-47 / wl 48-52 | match (genuine) |
| Tube-length law `L_W = pi a sqrt((1+r_c)/3)/2` (card Output, appendix 547) | sympy 52-62 / wl 56-60 (solve + closed-form check) | match (genuine) |
| `gamma_0 = (1+r_c)/9` | sympy 70,79-83 / wl 65,73-75 | partial (hardcoded ansatz + coefficient read-back, no independent derivation) |
| Renormalization -> `kappa_c = 1/3`, `gamma_c = 1/9` | sympy 71-74,84-93 / wl 66-70,76-80 | partial (divide-by-own-scale tautologies) |

`paper_alignment` = aligned: every paper-side target has a corresponding script-side target with the correct constants (`1/3`, `1/9`, `(1+r_c)/3`, `(1+r_c)/9`, tube-length closed form all match the card/notes/appendix). The deficiency is verification strength on the renormalization block, not a target/value mismatch -- hence `insufficient_verification`, not `paper_misalignment`.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 30-31 | `simplify(q''+k^2 q)==0` | ODE (claim i) | yes |
| A2 | sympy | 32-33 | `q(0)==0` | left BC (claim i) | yes (trivial but genuine) |
| A3 | sympy | 34-39 | `cos(k_W L_W)==0` at k_W | eigenvalue (claim i) | yes |
| A4 | sympy | 43-47 | `kappa0_derived - 4L_W^2/(pi^2 a^2)==0` | claim ii | yes |
| A5 | sympy | 52-62 | `solve(...)[0] - pi a sqrt((1+r_c)/3)/2==0` | tube-length law (claim iii / Output) | yes |
| A6 | sympy | 64-69 | `4 L_W_required^2/(pi^2 a^2) - (1+r_c)/3==0` | claim iii | no (solve-then-recheck round-trip) |
| A7 | sympy | 73 | `kappa_c - 1/3==0` | renorm (claim v) | no (divide by own scale) |
| A8 | sympy | 74 | `gamma_c - 1/9==0` | renorm (claim v) | no (divide hardcoded value by own scale) |
| A9 | sympy | 79-83 | `i*D_bare.coeff(z,5) - (1+r_c)/9==0` | gamma_0 (claim iv) | no (read back constructed coefficient) |
| A10 | sympy | 84-88 | `D_bare/(1+r_c) - (1-z^2/3-i z^5/9)==0` | renorm (claim v) | no (divide by own scale) |
| A11 | sympy | 89-93 | `-D_bare.coeff(z,2) - (1+r_c)/3==0` | kappa_0 (claim iii) | no (read back constructed coefficient) |
| B1-B11 | mathematica | 39-80 | one-to-one mirror of A1-A11 | same | same (see transliteration finding) |

A1-A5 (genuine) carry the stage's load-bearing physics -- the eigenvalue and the tube-length law that the card names as its Output. A6-A11 are unfalsifiable.

## Findings

### F1 — mathematica_transliteration

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage116_dn_mixed_tube_realization_mathematica_audit.wl:38-81`

**What's wrong:**
The `.wl` is a line-by-line port of the `.py`, not an independent re-derivation. Corresponding sections:

SymPy 29-39:
```
q_trial = sp.sin(k_sym * x_var)
ode_residual = sp.simplify(sp.diff(q_trial, x_var, 2) + k_sym**2 * q_trial)
...
bc_right = sp.simplify(sp.diff(q_trial, x_var).subs(x_var, L_W))
k_W_value = sp.pi / (2 * L_W)
```
Mathematica 38-46:
```
qTrial[xVar_, kArg_] := Sin[kArg*xVar];
odeRes = FullSimplify[D[qTrial[x, kSym], {x, 2}] + kSym^2 * qTrial[x, kSym], ...];
...
kWValue = Pi/(2*lW);
bcRightAtKW = FullSimplify[(D[qTrial[x, kSym], x] /. x -> lW) /. kSym -> kWValue, ...];
```
Both *posit* the same trial function `sin(k x)` rather than solving the ODE+BCs independently (e.g. via `DSolve`), and verify it identically.

SymPy 52-62 vs Mathematica 56-60:
```
L_W_required = sp.solve(sp.Eq(kappa0_from_tube, (1 + r_c) / 3), L_W)[0]
expect_zero("tube-length law...", sp.simplify(L_W_required - sp.pi * a * sp.sqrt((1 + r_c) / 3) / 2))
```
```
lWRequired = FullSimplify[lW /. First[Solve[kappa0FromTube == (1 + rC)/3, lW, Reals]], ...];
expectZero["tube-length law", lWRequired - (Pi*a*Sqrt[(1 + rC)/3])/2];
```
`sp.solve(...)[0]` <-> `First[Solve[...]]`, identical variable choreography (`kappa0_derived`<->`kappa0Derived`, `r_c`<->`rC`, `L_W_required`<->`lWRequired`), identical step ordering through the entire file, and identical coefficient read-backs (`D_bare.coeff(z,5)` <-> `Coefficient[dBare, z, 5]`). This is a transliteration, not a second independent engine.

**Why this matters:**
The two-engine policy exists so a sign/branch error baked into one derivation is caught by a structurally different one. A transliteration cannot catch such an error because it reproduces the same algebra. The stage's load-bearing claim (the tube-length law) is therefore confirmed by only one genuine derivation, mirrored.

**Required change:**
In the `.wl`, derive the eigenvalue independently rather than positing `Sin[k x]`: use `DSolve[{q''[x] + k^2 q[x] == 0, q[0] == 0, q'[lW] == 0}, q, x]`, or solve the characteristic condition `Cos[k lW] == 0` for the smallest positive root via `Reduce[Cos[k lW] == 0 && 0 < k lW < Pi, k]`, to obtain `kW = Pi/(2 lW)` from first principles, then feed that into the `kappa0`/tube-length checks. Keep the asserted targets identical (they already match the paper). See directive F1.

**Verification:**
The new `.wl` must obtain `kW = Pi/(2 lW)` from an ODE/BC solve (not from a posited `Sin`) and still print all PASS lines and exit 0; the `kappa0`/tube-length/renorm assertions are unchanged.

### F2 — insufficient_verification

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage116_dn_mixed_tube_realization_sympy_audit.py:64-93`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage116_dn_mixed_tube_realization_mathematica_audit.wl:62-80`

**What's wrong:**
A block of checks cannot fail regardless of the physics:

- Solve-then-recheck round-trip (sympy 64-69 / wl 62-64): `L_W_required` was obtained at line 52-55 by solving `4 L_W^2/(pi^2 a^2) = (1+r_c)/3`. Lines 65-69 then compute `kappa0_bare_geom = 4 L_W_required^2/(pi^2 a^2)` and assert it equals `(1+r_c)/3`. By construction this is identically true; the assertion adds nothing over the genuine tube-length check at 52-62.
- Divide-by-own-scale (sympy 70-74, 84-88 / wl 65-70, 76-77): `gamma0_bare` is *defined* as `(1+r_c)/9` at line 70 and then divided by `(1+r_c)` at line 72 to assert `gamma_c == 1/9` at line 74. Likewise `D_final = D_bare/(1+r_c)` at line 84 is asserted to equal the canonical form at 85-88, but `D_bare` was *built* as `(1+r_c)*(canonical)` at line 78. Dividing a quantity by the literal factor it was constructed with is unfalsifiable.
- Read-back of a constructed coefficient (sympy 79-83, 89-93 / wl 73-75, 78-80): `D_bare` is the hand-written ansatz `(1+r_c)(1 - z^2/3 - i z^5/9)`; lines 79-83 read its `z^5` coefficient and assert it gives `(1+r_c)/9`, and lines 89-93 read its `z^2` coefficient and assert `(1+r_c)/3`. These read back exactly the literals typed into line 78; they confirm Python's `coeff` extractor, not the physics.

The genuine content of the stage -- the D/N eigenvalue and the tube-length law (A1-A5) -- is correctly and non-tautologically verified, so this is a verification-strength gap, not a wrong target or paper misalignment.

**Why this matters:**
The renormalization claim (that removing `(1+r_c)` yields canonical `1/3`, `1/9`) and the `gamma_0`/`kappa_0` coefficient claims are presented in the transcript as passing checks, giving false assurance that something physical was tested. As written they would pass even if the canonical coefficients, the bare scale factor, or the `D` ansatz were wrong, because each quantity is divided by / read back from the very expression that defines it.

**Required change:**
Tie the renormalization checks to the *genuinely derived* `kappa_0` rather than to round-tripped or hardcoded quantities, and add a single falsifiable joint-uniformity check, so a wrong canonical value would actually surface. See directive F2 for the exact edits. Codex should apply what is mechanically safe and `## Blocked` any sub-step it deems unsafe.

**Verification:**
After the fix: sympy line range 64-93 no longer contains a `kappa0_bare_geom - (1+r_c)/3` round-trip, no `D_bare/(1+r_c)`-equals-canonical read-back, and no `D_bare.coeff(...)`-equals-literal read-back; instead a `3*kappa_0 - 9*gamma_0` (uniform-scale) check appears and the renormalization assertion references `(1+r_c)/3` as the required `kappa_0`. Both scripts still print all PASS lines and exit 0. The `.wl` shows the same structural change.

## Independent-derivation check (Mathematica)

Not independent. The `.wl` posits the identical trial function `Sin[k x]` (line 38) the SymPy script posits (line 29), performs the identical `First[Solve[...]]` the SymPy `sp.solve(...)[0]` performs (line 56 vs 52-55), and reads the identical `Coefficient[dBare, z, 5]` / `Coefficient[dBare, z, 2]` the SymPy `D_bare.coeff(z, 5)` / `D_bare.coeff(z, 2)` reads (lines 73,78 vs 79,89). Every intermediate variable maps one-to-one. This is `mathematica_transliteration` (F1).

## Engine cross-check

Both engines pass and agree. SymPy output line 10 / Mathematica output line 14 give the tube length in algebraically equivalent forms: `pi*a*sqrt(3*K_q*K_s + 3*lam**2)/(6*sqrt(K_q)*sqrt(K_s))` vs `(a*Sqrt[3 + (3*lam^2)/(kQ*kS)]*Pi)/6`, both equal to `pi a sqrt((1+r_c)/3)/2`. Final coefficients `kappa_c = 1/3`, `gamma_c = 1/9` agree (sympy out line 22, wl out line 33). No `engine_disagreement`.

## Verdict justification

Verdict is `findings`. The stage's central, card-named deliverable -- the tube-length law `L_W = pi a sqrt((1+r_c)/3)/2` derived from the D/N half-wave eigenvalue -- is genuinely and non-tautologically verified, with constants matching the card, notes, and appendix (so `paper_alignment` = aligned, no `paper_misalignment`). Two real defects remain: (F1) the Mathematica script is a line-by-line transliteration of the SymPy script and so is not an independent second engine; (F2) the renormalization / coefficient block (sympy 64-93, wl 62-80) is built from solve-then-recheck round-trips, divide-by-own-scale operations, and read-backs of a hand-constructed ansatz -- all unfalsifiable. Attacks that failed: the eigenvalue check is not trivial (Neumann BC genuinely forces `cos(k_W L_W)=0`); the `kappa0` dimensionless collapse genuinely cancels `omega` and `c_s`; the `solve(...)[0]` index happens to return the positive root and the closed-form assertion would catch a wrong branch; outputs are fresh (both `.txt` newer than their scripts); the symbol positivity assumptions (`a, L_W, K_s, K_q, lam > 0`) are justified by the geometric/physical setup. No stop-cold condition: neither finding propagates a sign/constant change downstream (the verified targets are unchanged).

## Self-test notes

Checked the variable-independence trap: the proposed F2 uniform-scale check `3*kappa_0 - 9*gamma_0` uses no `diff`, so no identically-zero-derivative trap applies; the existing genuine `diff` calls (ODE residual, BC) depend on `x` and are non-trivial. Trivial-case substitution: with `r_c -> 0`, `3*kappa_0 = 9*gamma_0 = 1` so the proposed check gives `1 - 1 = 0` correctly, and against a hypothetically-wrong `gamma_0 = (1+r_c)/8` it gives `(1+r_c) - 9(1+r_c)/8 = -(1+r_c)/8 != 0`, confirming the new check is falsifiable. Paper round-trip: the F2 replacement keeps `kappa_0 = (1+r_c)/3` and `gamma_0 = (1+r_c)/9` exactly as the notes/appendix state, introducing no new constant, so no new `paper_misalignment`; F1 keeps every asserted target identical and only changes how the eigenvalue is obtained.

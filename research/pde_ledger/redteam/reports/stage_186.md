---
unit_id: 186
batch: V.2
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-05-30T00:00:00Z
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
  notes_stage_files: [/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage186_similarity_orbit_closure.md]
  paper_appendix: present
---

# Audit unit 186 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_186.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage186_similarity_orbit_closure.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part05.tex` (Table row line 103; "Similarity-orbit theorem" section lines 1049-1085; M_* matrix lines 1025-1045)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage186_similarity_orbit_closure_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage186_similarity_orbit_closure_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage186_similarity_orbit_closure_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage186_similarity_orbit_closure_mathematica_audit.txt`

## What the paper claims

Stage 186 (`\stagefield{Output}`): "Shows the zero-defect equations are the tangent equations of a five-parameter monomial-preserving similarity orbit \(\mathcal G_*\)." The notes and appendix (lines 1049-1085, eqs `app-part05-Keta/TU/muW-orbit-law`) make this concrete: (1) define a finite five-parameter multiplicative family `G_*` with free log-scalings `(Λ,C,Γ,U,W)` on `(λ_W, c_ηU, γ, K_U, K_W^eff)` and three *dependent* scalings forced by monomial preservation — `K_η^eff → e^{2C−U}`, `T_U → e^{U − ((1+δ_*)/(1+χ_*))(Γ+C−U)}`, `μ_W → e^{M(Λ,C,Γ,U,W)}` with the explicit M-exponent (lines 1077-1083); (2) `G_*` preserves the three direct monomials `(C_tr,*, C_nt,*, ε_η)` exactly; (3) the monomial-drift map `M_*` (3×8, lines 1028-1032) has rank 3 via the dependent minor `det M_*^{(τ,κ_η,μ)} = 1+χ_0,* > 0`, hence `dim ker M_* = 5`; (4) the linearization of `G_*` at the identity reproduces the Stage 185/253 compatibility ledger exactly, so `ker M_* = T_id G_*`. Every load-bearing constant (the matrix entries, the dependent-scaling exponents, the minor determinant) is stated in the paper.

## What the script claims to verify

The SymPy script (docstring lines 2-16) claims to: build the exact 3×8 matrix `M_*`; verify the convenient 3×3 minor determinant equals `1+χ`; solve the linear compatibility ledger for `(τ_1, κ_η, μ_1)` and match the Stage-185 closed forms; construct the finite five-parameter orbit and verify exact preservation of the three monomials; and verify the orbit's linearization reproduces the compatibility formulas. The Mathematica script (`.wl`) makes the identical claim with the same checks. Both exit 0 with all checks passing.

## Paper ↔ script cross-check

| Paper deliverable | Script-side check | Status |
|---|---|---|
| M_* matrix (lines 1028-1032) | hardcoded `M` (py 52-56) / `mMat` (wl 40-44) — entry-by-entry identical to paper | match |
| rank 3 / 5-dim kernel via `det M_*^{(τ,κ_η,μ)}=1+χ_*` | `minor.det()==1+chi` (py 60-64) / `mMat[[All,{8,5,7}]]` (wl 48-52) | match |
| dependent scaling `K_η→e^{2C−U}` | `Eta_exp=2*C-U` + `Eta_orbit` check (py 101,112,116) | partial (check is tautological — see F1) |
| dependent scaling `T_U` | `Tau_exp` + `Ctr_orbit` check (py 102,110,114) | match (substantive cancellation) |
| dependent scaling `μ_W` | `Mu_exp` + `Cnt_orbit` check (py 103,111,115) | match (substantive cancellation) |
| orbit preserves 3 monomials | `Ctr_orbit/Cnt_orbit/Eta_orbit==0` (py 114-116) | match for C_tr, C_nt; tautological for ε_η (F1) |
| `ker M_* = T_id G_*` (compatibility match) | `solve()` + linearization checks (py 67-122) | match |
| independent second-engine derivation | `.wl` is a line-by-line transliteration of `.py` (F2) | mismatch (policy) |

`paper_alignment: aligned` — every paper-side constant and identity is faithfully encoded; the two findings are about check quality (one tautological assertion) and second-engine independence, not about the script verifying a wrong claim.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 63-64 | `det_minor != 1+chi` → raise | rank-3 minor (claim 3) | yes |
| A2 | sympy | 81-84 | `expect_zero(tau1 - StageForm)` | compatibility τ_1 (claim 4) | yes |
| A3 | sympy | 85 | `expect_zero(kEta - (2c1-kU))` | compatibility κ_η (claim 4) | yes |
| A4 | sympy | 86-96 | `expect_zero(mu1 - StageForm)` | compatibility μ_1 (claim 4) | yes |
| A5 | sympy | 114 | `expect_zero(Ctr_orbit)` | orbit preserves C_tr (claim 2) | yes (Tau_exp nontrivial) |
| A6 | sympy | 115 | `expect_zero(Cnt_orbit)` | orbit preserves C_nt (claim 2) | yes (Mu_exp nontrivial) |
| A7 | sympy | 116 | `expect_zero(Eta_orbit)` | orbit preserves ε_η (claim 2) | **no — tautological (F1)** |
| A8 | sympy | 120 | `expect_zero(Tau_exp.subs - tau_formula)` | ker=T_id G_* (claim 4) | yes |
| A9 | sympy | 121 | `expect_zero(Eta_exp.subs - kEta_formula)` | ker=T_id G_* (claim 4) | yes |
| A10 | sympy | 122 | `expect_zero(Mu_exp.subs - mu_formula)` | ker=T_id G_* (claim 4) | yes |
| B1-B10 | mathematica | 49-103 | identical structure to A1-A10 | same | same (A7↔wl 97 tautological) |

## Findings

### F1 — tautological_check

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage186_similarity_orbit_closure_sympy_audit.py:101,112,116`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage186_similarity_orbit_closure_mathematica_audit.wl:83,93,97`

**What's wrong:**
The "finite orbit preserves ε_η" check is structurally `X − X` and cannot fail. The script sets the dependent K_η scaling and the ε_η drift form from the *same literal* `2*C - U`:
```
Eta_exp  = 2*C - U          # py line 101  (the chosen K_η^eff scaling)
Eta_orbit = 2*C - U - Eta_exp   # py line 112  = (2*C - U) - (2*C - U) ≡ 0
expect_zero("finite orbit preserves eps_eta", Eta_orbit)   # py 116
```
Because `Eta_orbit` substitutes `Eta_exp` (= `2*C - U`) into a drift form whose only other term is the identical `2*C - U`, the residual is identically zero by construction — independent of whether `2C−U` is actually the ε_η-preserving scaling. The companion C_tr and C_nt checks (A5/A6) are NOT tautological: their dependent scalings `Tau_exp` and `Mu_exp` are non-trivial, so their cancellations genuinely certify the chosen scalings. Only the ε_η check degenerates. The Mathematica check at wl:97 has the identical defect (`etaExp = 2*cSym - uSym`, `etaOrbit = 2*cSym - uSym - etaExp`).

The paper-side claim this is supposed to exercise is that `G_*` preserves `ε_η = c_ηU² / (K_U K_η^eff)` exactly (notes eq lines 282-287; appendix `eq:app-part05-epsilon-eta-monomial` line 989; orbit law `eq:app-part05-Keta-orbit-law` line 1062). The script never grounds the `2` (from `c_ηU²`) and the two `−1`s (from `K_U`, `K_η^eff`) in the monomial's actual exponent structure; it asserts the cancellation of two textually-identical copies.

**Why this matters:**
A transcription error in the ε_η monomial structure (e.g. if a future edit mistyped the exponent on `c_ηU` or the K_η scaling) would be invisible to this check — it passes no matter what the shared literal is. The check provides zero verification coverage for the one of the three monomial-preservation deliverables that it nominally covers.

**Required change:**
Replace the trivial `Eta_orbit = 2*C - U - Eta_exp` with a check that grounds the exponents in the actual ε_η monomial. Build the ε_η log-drift from independent per-variable scalings of `c_ηU`, `K_U`, `K_η^eff` (with the K_η scaling = `Eta_exp`), so the `2`, `−1`, `−1` exponents are load-bearing:
```python
# ε_η = c_ηU^2 / (K_U * K_η^eff); free scalings C on c_ηU, U on K_U, Eta_exp on K_η^eff
eps_eta_drift = 2*C - U - Eta_exp
```
is still the same expression, so instead make the per-variable exponents explicit and assert the chosen `Eta_exp` is the unique scaling that zeroes the *monomial* drift — i.e. solve for the preserving K_η scaling and confirm it equals the paper value `2C − U`:
```python
eta_scaling = sp.symbols("eta_scaling", real=True)
# ε_η = c_ηU^2 * K_U^{-1} * K_η^{-1}; log-drift under (C, U, eta_scaling):
eps_eta_logdrift = 2*C - U - eta_scaling
solved_eta = sp.solve(sp.Eq(eps_eta_logdrift, 0), eta_scaling)[0]
expect_zero("K_eta preserving scaling matches paper 2C-U", solved_eta - (2*C - U))
expect_zero("chosen Eta_exp solves eps_eta preservation",
            eps_eta_logdrift.subs(eta_scaling, Eta_exp))
```
The first new check is non-tautological: `solve()` derives the preserving scaling from the monomial's `(2, −1, −1)` exponent vector independently of the hardcoded `Eta_exp`, then compares against the paper-stated `2C − U`. Apply the analogous change at wl:83,93,97.

**Verification:**
After the fix the sympy output should show a new line `K_eta preserving scaling matches paper 2C-U = 0` (and the Mathematica `PASS:` analogue), and the script still exits 0. The new `solve`-derived scaling must equal `2*C - U`.

### F2 — mathematica_transliteration

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage186_similarity_orbit_closure_mathematica_audit.wl` (whole file)
- compared against `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage186_similarity_orbit_closure_sympy_audit.py`

**What's wrong:**
The `.wl` is a line-by-line port of the `.py`, not an independent re-derivation. The variable choreography, intermediate steps, and check order are identical; only the syntax is rewritten. Corresponding sections:

1. Drift expressions — `.py:44-46`
```python
dlog_Ctr = (1 + delta) * (gam1 + c1 - kU) + (1 + chi) * (tau1 - kU)
dlog_Cnt = (2 + 2*E) * lam1 + 2*E * gam1 + mu1 + (F - E) * kU - kEta - (2 + E) * kW - F * tau1
```
`.wl:32-33` — byte-for-byte identical algebra (`E`→`eParam`, `F`→`fParam`):
```wolfram
dlogCtr = (1 + delta)*(gam1 + c1 - kU) + (1 + chi)*(tau1 - kU);
dlogCnt = (2 + 2*eParam)*lam1 + 2*eParam*gam1 + mu1 + (fParam - eParam)*kU - kEta - (2 + eParam)*kW - fParam*tau1;
```

2. Hardcoded matrix — `.py:52-56` `M = sp.Matrix([[0, 1+delta, ...], ...])` ↔ `.wl:40-44` `mMat = {{0, 1+delta, ...}, ...}` — identical entries; same minor column choice `M[:, [7,4,6]]` ↔ `mMat[[All, {8,5,7}]]`.

3. Orbit exponents — `.py:101-103`
```python
Eta_exp = 2*C - U
Tau_exp = U - (1 + delta) * (Gam + C - U) / (1 + chi)
Mu_exp = Eta_exp + 2*W - 2*Lam - E * (2*Gam + 2*Lam - U - W) + F * (Tau_exp - U)
```
`.wl:83-85` — identical line-for-line transliteration.

Both engines hardcode the *same* `M_*` and the *same* dependent-scaling exponents, then run `solve`/`Solve` on the *same* equations. They cannot disagree because they encode identical algebra; the Mathematica run is not an independent witness. Both even share the copy-paste banner artifact `"STAGE 169 — EXACT MICROSCOPIC SIMILARITY ORBIT"` (`.py:33` / `.wl:26`) on what is Stage 186.

**Why this matters:**
The second-engine policy requires the Mathematica audit to derive the result independently from the physical premises so that an algebra error in one engine is caught by the other. A transliteration provides no such cross-check: a wrong matrix entry or wrong orbit exponent typed once and copied to both engines passes in both. The engine "agreement" in the outputs is therefore not evidence of correctness.

**Required change:**
Re-derive at least the load-bearing objects independently in the `.wl` rather than copying the `.py` algebra:
- Build `M_*` rows as the logarithmic exponent vectors of the three monomials from their factored monomial forms (`C_tr,*`, `C_nt,*`, `ε_η` as given in notes eqs lines 108-123 / appendix lines 971-990), i.e. read each matrix entry off the monomial exponents rather than transcribing the `.py` literal matrix. Then assert this independently-built matrix equals the expected entries.
- Obtain the dependent scalings by `Solve`-ing the monomial-preservation equations `M_* · (orbit log-vector) == 0` for `(K_η, T_U, μ_W)` scalings, rather than transcribing `Eta_exp/Tau_exp/Mu_exp` from the `.py`. Compare the solved scalings against the paper-stated closed forms (notes lines 229-259 / appendix lines 1060-1083).
- Fix the banner string to read Stage 186 (`.wl:26`, and `.py:33`).

This keeps the same physical premises but makes the Mathematica path a genuine re-derivation. Do not introduce new claims; the verified identities stay exactly those in the paper.

**Verification:**
After the fix, the `.wl` should construct `M_*` from monomial exponents (not a transcribed literal) and `Solve` for the dependent scalings; the verifier confirms the independently-built matrix matches the paper entries, the solved scalings match the paper closed forms, and the script exits 0. Banner line reads "STAGE 186".

## Independent-derivation check (Mathematica)

Transliteration confirmed — see F2. The `.wl` mirrors the `.py` step for step (drift expressions, hardcoded matrix, minor column selection, orbit exponents, basis substitution, and even the stale "STAGE 169" banner). It is not an independent derivation of the Stage 186 claim.

## Engine cross-check

Both engines produce identical results and all checks pass:
- `minor det(tau1,kEta,mu1)` = `chi + 1` (sympy line 24) ↔ `1 + chi` (mathematica line 19) — agree.
- `kEta = 2*c1 - kU` in both (sympy 30, mathematica 25) — agree.
- `tau1`, `mu1` solved forms are algebraically equal (sympy 29,31; mathematica 24,26) — agree (different surface form, same content).
- All six orbit-preservation / linearization checks read `= 0` / `PASS:` in both.

The agreement is real at the surface but, per F2, is not independent evidence because the two engines encode the same hardcoded algebra. No `engine_disagreement` finding (both agree); the concern is that they are too identical to constitute a cross-check.

## Verdict justification

The script faithfully encodes every paper-side constant and identity for Stage 186: the M_* matrix entries match line-for-line, the rank-3 minor determinant `1+χ` is checked, the compatibility solve for `(τ_1, κ_η, μ_1)` independently reproduces the Stage-185 closed forms (genuine, non-tautological), and the C_tr and C_nt orbit-preservation checks substantively certify the non-trivial dependent scalings `Tau_exp` and `Mu_exp`. Paper alignment is exact — no `paper_misalignment`. Two real script-side defects remain: (F1) the ε_η orbit-preservation check is tautological (`X−X` by construction, present in both engines), giving zero coverage for that one monomial-preservation deliverable; and (F2) the Mathematica script is a line-by-line transliteration of the SymPy script, so the second-engine "agreement" is not independent. Attacks that failed: I tried to break the M·δx ↔ dlog consistency (rows match exactly), the minor column selection (correct in both 0- and 1-indexed forms), the symbol domains (chi,delta>0 matches the constructive branch; division by 1+chi is safe), and the rank argument (a nonzero 3×3 minor of a 3-row matrix rigorously certifies rank 3 ⇒ kernel dim 5). Outputs are fresh (both .txt mtimes postdate their scripts). Verdict: findings.

## Self-test notes

Variable-independence: the proposed F1 fix uses `sp.solve(eps_eta_logdrift==0, eta_scaling)`; `eps_eta_logdrift = 2*C − U − eta_scaling` depends on `eta_scaling`, so the solve is well-posed and yields `2*C − U` (a nonzero literal), not a degenerate identity. Symmetry/parity: no integrals in this unit — N/A. Trivial-case pre-check: substituting `Eta_exp = 2*C − U` into the new `eps_eta_logdrift` gives `2C−U−(2C−U)=0` (assert_zero holds), and the `solve`-derived scaling `2C−U` minus paper `2C−U` is `0` while a hypothetical wrong `Eta_exp` (e.g. `C−U`) would yield a nonzero residual `C` — so the rewritten check genuinely can fail. Paper round-trip: both fixes keep the verified objects (M_* entries, dependent scalings, monomials) exactly equal to the paper/notes values (notes lines 229-259; appendix lines 989, 1060-1083), introducing no new constant or `paper_misalignment`.

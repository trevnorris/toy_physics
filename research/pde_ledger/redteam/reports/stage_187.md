---
unit_id: 187
batch: V.2
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-05-30T00:00:00-06:00
verdict: findings
stop_cold: null
findings_count: 2
paper_alignment: partial
scripts_checked:
  sympy: present
  mathematica: present
  engines_agree: true
  outputs_fresh: true
docs_read:
  paper_stage_tex: present
  notes_stage_files: [moving_throat_pde_stage187_orbit_quotient_closure.md]
  paper_appendix: present
---

# Audit unit 187 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_187.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage187_orbit_quotient_closure.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part05.tex` (rows at lines 105, 1087-1136 reference this unit)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage187_orbit_quotient_closure_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage187_orbit_quotient_closure_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage187_orbit_quotient_closure_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage187_orbit_quotient_closure_mathematica_audit.txt`

## What the paper claims

Stage 187 proves the **converse** of the Stage 186/237 orbit construction: that the three branch monomials \((\mathfrak C_{{\rm tr},*},\mathfrak C_{{\rm nt},*},\epsilon_\eta)\) are *complete* finite orbit invariants. `\stagefield{Output}` states verbatim: "Proves finite level sets of \((\mathfrak C_{{\rm tr},*},\mathfrak C_{{\rm nt},*},\epsilon_\eta)\) are exactly \(\mathcal G_*\)-orbits." The distinctive new content over Stage 186 (notes §2 line 149, §3 line 211) is the **finite-level** statement: because the three invariants are *monomials* in the eight microscopic variables (notes §1, lines 83-98), equality of the invariant *values* between two arbitrary positive states \(x,\widetilde x\) is **exactly linear** in the finite log-ratios \(\Delta\mathbf x=\ln(\widetilde x/x)\), giving \(M_*\Delta\mathbf x=0\) with the rank-3 Stage 237 matrix (notes §3, lines 198-207). Deliverables: (D1) the three finite monomial-equality conditions equal the three linear rows \(M_*\Delta\mathbf x=0\) **derived from the §1 monomials**; (D2) the selected dependent minor is nonsingular, \(\det = 1+\chi_{0,*}>0\) (notes §4 line 228); (D3) solving for \((\Delta_\eta,\Delta_T,\Delta_\mu)\) reproduces the three boxed orbit laws (notes §4 lines 235-256); (D4) hence fibres = \(\mathcal G_*\)-orbits and \(\mathcal M_+/\mathcal G_*\cong(\mathbb R_{>0})^3\) (notes §5-6).

## What the script claims to verify

The SymPy docstring (lines 7-15) claims: (1) the exact finite log-ratio equalities of the three monomial invariants are exactly \(M_*\Delta\mathbf x=0\); (2) solving for \((\Delta_\eta,\Delta_T,\Delta_\mu)\) reproduces the Stage 186 orbit laws; (3) back-substitution annihilates the equations; (4) the monomials are therefore complete finite orbit invariants. What the assertions actually test: the hand-typed matrix `M` reproduces the hand-typed rows (`Mx[i] - row_i == 0`); the linear solve of the three posited rows reproduces the three closed-form orbit laws; back-substitution of the solve zeroes the rows. The Mathematica script adds tautological `Log[Exp[row]] - row == 0` ratio checks and `ratio-after-solve == 1` checks, plus an explicit assertion that the minor det = \(1+\chi\).

## Paper ↔ script cross-check

| Paper deliverable | Script-side coverage | Status |
|---|---|---|
| D1: three finite monomial-equality conditions equal \(M_*\Delta\mathbf x=0\), **derived from the §1 monomials** | rows `row_tr/nt/eta` are **posited as given linear expressions** (sympy 46-48, wl 34-36); never built from \(\mathfrak C_{tr,*}=(\gamma c_{\eta U}/K_U)^{1+\delta}(\pi^2T_U/L^2K_U)^{1+\chi}\) etc.; the only "monomial" touch is `ctrRatio=Exp[row]` then `Log[ctrRatio]-row` (wl 38-55), which is tautological | partial / mismatch |
| D2: \(\det\)(dependent minor) = \(1+\chi>0\) | sympy only **prints** det (line 73, no assert); wl asserts `Det[minor]-(1+chiStar)==0` (line 75); `minor` is hand-typed, not extracted from `M` | partial |
| D3: solve reproduces the three orbit laws | sympy `expect_zero` on `sol[DEta/DT/DM] - *_expected` (lines 89-91); wl analogues (lines 94-96) | match |
| D4: fibres = orbits, quotient \(\cong(\mathbb R_{>0})^3\) | asserted only in print statements (sympy 99-113, wl 105-121); no quotient/projector check (those are Stage 191-192 per appendix line 1138) | partial (acceptable — D4 is the prose conclusion of D1-D3) |

Set `paper_alignment: partial`: the solve (D3) is faithfully exercised, but the *distinctive finite-linearity content* (D1) — the whole reason Stage 187 exists separately from Stage 186 — is asserted, not derived, in either engine.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 64-66 | `Mx[i] - row_i == 0` | D1 (matrix↔rows) | no — tautological; `M` hand-typed to match the rows |
| A2 | sympy | 73 | `print(minor.det())` (no assert) | D2 | no — print only, not asserted |
| A3 | sympy | 89-91 | `sol[DEta/DT/DM] - *_expected == 0` | D3 | yes — genuine linear solve vs closed forms |
| A4 | sympy | 94-96 | `row_i.subs(sol) == 0` | D3 (consistency) | partial — follows from A3 |
| A5 | wl | 53-55 | `Log[Exp[row]] - row == 0` | D1 (finite ratio) | no — tautological (exp∘log of same row) |
| A6 | wl | 65-67 | `mx[i] - row_i == 0` | D1 (matrix↔rows) | no — tautological, mirrors A1 |
| A7 | wl | 75 | `Det[minor] - (1+chiStar) == 0` | D2 | partial — `minor` hand-typed, not extracted from `m` |
| A8 | wl | 94-96 | `(deta/dt/dm /. sol) - *Expected == 0` | D3 | yes — genuine solve vs closed forms |
| A9 | wl | 98-103 | `row /. sol == 0`, `(ratio /. sol) - 1 == 0` | D3 | partial — follows from A8; ratio checks tautological w/ A5 |

The substantive, non-tautological content of both engines reduces to D3 (the linear solve reproduces the orbit laws). D3 is real and matches the paper exactly. D1 — the finite monomial-to-linear-row construction that is the stage's reason for existing — is not exercised by any assertion in either engine.

## Findings

### F1 — insufficient_verification

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage187_orbit_quotient_closure_sympy_audit.py:46-48`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage187_orbit_quotient_closure_mathematica_audit.wl:38-55`

**What's wrong:**
The stage's distinctive deliverable (D1) is the *finite-level* claim that equality of the actual monomial **values** between two states is exactly the linear log-ratio condition. Notes §1 give the monomials explicitly, e.g.
> `\mathfrak C_{{\rm tr},*} = (\gamma c_{\eta U}/K_U)^{1+\delta_{U,*}} (\pi^2 T_U / L^2 K_U)^{1+\chi_{0,*}}` (notes lines 83-87)

and notes §2 line 149 / §3 line 211 assert this equality is "exactly linear in these finite log-ratios", yielding `row_tr/nt/eta`. The SymPy script never builds these monomials; it posits the rows directly:
> `row_tr = (1 + deltaU) * (DG + DC - DU) + (1 + chi) * (DT - DU)` (sympy line 46)

The Mathematica script's only gesture toward the monomials is
> `ctrRatio = FullSimplify[Exp[(1 + deltaStar)*(dg + dc - du) + (1 + chiStar)*(dt - du)]]` then `expectZero["log C_tr ratio - row_tr", Log[ctrRatio] - rowTr]` (wl 38-41, 53),

which is `Log[Exp[rowTr]] - rowTr == 0` — tautological. It defines the ratio *as* `Exp[row]` rather than from the §1 monomial structure, so it cannot detect a transcription error in any exponent. Consequently a wrong coefficient in any row (e.g. an exponent `1+\delta` mistyped as `1+\chi`) would survive: A1/A6 are tautological against the same rows, and the A3/A8 solve simply propagates whatever rows were typed. **No script in either engine checks that the monomials of notes §1 actually produce these rows at finite level.** (The upstream Stage 185 script verifies only the *infinitesimal* first-order drift via `first_ratio_drift = diff(ratio,eps)|_0`, not the finite log-ratio, so the finite claim is not covered upstream either.)

**Why this matters:**
Stage 187's entire reason for existing separate from Stage 186 is the finite-level upgrade. As written, both scripts only re-run Stage 186's linear-algebra solve on hand-copied rows; the load-bearing monomial→linear-row step is unguarded. A checkpoint-grade stage (`is_checkpoint: false` here, but the appendix flags this as `\StatusExactClosure`) should have its central claim exercised, not asserted.

**Required change:**
In **both** scripts, add a block that constructs the three invariant ratios *from the §1 monomial forms* in terms of positive primitive ratios, then verifies `log(ratio) - row == 0`. Concretely (SymPy): declare positive symbols `rL, rC, rG, rU, rEta, rW, rM, rT` for the eight microscopic ratios \(\widetilde x_i/x_i\); set `DL=log(rL), …, DT=log(rT)`; build
`Ctr_ratio = (rG*rC/rU)**(1+deltaU) * (rT/rU)**(1+chi)`,
`Cnt_ratio = (rL**2*rM/(rEta*rW**2)) * (rG**2*rL**2/(rU*rW))**E * (rT/rU)**(-F)`,
`eps_ratio = rC**2/(rU*rEta)`,
then `expect_zero("log C_tr ratio - row_tr", sp.expand_log(sp.log(Ctr_ratio), force=True) - row_tr)` and likewise for the other two. Mirror the structure in Mathematica independently (build the ratios from the primitive ratios, not from `Exp[row]`). The constant factors \(\pi^2/L^2\) and \(\sigma\) drop out of the ratios, so no new constants are introduced. (See `## Resolve`/directive for the exact monomial forms reconstructed from notes §1.)

**Verification:**
After the fix, both transcripts must show three new non-tautological PASS lines (`log C_tr ratio - row_tr = 0`, `log C_nt ratio - row_nt = 0`, `log epsilon_eta ratio - row_eta = 0`) where the ratio is built from primitive ratios, and both scripts still exit 0.

### F2 — mathematica_transliteration

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage187_orbit_quotient_closure_mathematica_audit.wl:34-103`
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage187_orbit_quotient_closure_sympy_audit.py:46-96`

**What's wrong:**
The `.wl` is a line-by-line port of the `.py` rather than an independent re-derivation. Corresponding sections:
- Rows: wl 34-36 (`rowTr = (1 + deltaStar)*(dg + dc - du) + (1 + chiStar)*(dt - du)` …) vs sympy 46-48 (`row_tr = (1 + deltaU) * (DG + DC - DU) + (1 + chi) * (DT - DU)` …) — identical posited expressions, same term order.
- Matrix: wl 57-61 vs sympy 56-60 — identical 3×8 literal entries.
- Solve + expected forms: wl 77-96 vs sympy 76-91 — same `Solve`/`solve` over `{deta,dt,dm}`, same closed-form `*Expected`.
The only `.wl` additions (the `Exp/Log` ratio checks, A5) are tautological and do not constitute independent derivation. The second-engine policy requires each engine to derive the rows from the physical premises (the §1 monomials); neither does, and the `.wl` simply echoes the `.py`'s linear algebra.

**Why this matters:**
A transliterated second engine provides no real cross-check: a coefficient error copied into both files passes both. Independent derivation is the entire point of the two-engine policy.

**Required change:**
The F1 fix resolves this if applied **independently in each engine**: have the Mathematica script build the invariant ratios from the §1 monomials and take their `Log` to *derive* `rowTr/rowNt/rowEta` (rather than positing them and rederiving via `Exp[row]`), and assert the derived rows equal the matrix action `m.dx`. Once each engine derives the rows from the monomials by its own algebra, the cross-check becomes genuine. No separate edit beyond F1 is required if F1's Mathematica side is implemented as an independent derivation rather than a copy of the SymPy block.

**Verification:**
After the fix, the `.wl` must contain a monomial-based derivation of the rows (not `Exp[row]`), and its PASS lines for the new ratio checks must reference ratios built from primitive ratios `rL,rC,…`. The two engines then independently arrive at the same rows.

## Independent-derivation check (Mathematica)

The `.wl` is a transliteration (see F2). Beyond the structural copy, its distinctive lines (38-55, `ctrRatio = Exp[row]; Log[ctrRatio] - row`) are tautological and add no independent content. The solve, matrix, minor, and expected forms (lines 57-96) mirror the SymPy script entry-for-entry.

## Engine cross-check

Both engines agree at the level they claim. SymPy output (lines 23-25) and Mathematica output (lines 34-36) report the same `Delta_eta = 2*dc - du`, the same `Delta_T = (-((1+deltaStar)(dc+dg)) + (2+chiStar+deltaStar)du)/(1+chiStar)` (algebraically equal to the SymPy form, verified by hand), and the same `Delta_mu`. All `expect_zero`/`expectZero` checks report 0 and both exit 0. The agreement is genuine for the solve (D3) but, per F2, is partly an artifact of shared posited rows rather than independent derivation.

## Verdict justification

The substantive content that both engines exercise — the linear solve of the three rows reproduces the three orbit laws (D3) — holds up and matches the notes §4 boxed formulas and the appendix orbit-fibre theorem exactly. I verified the rows against notes §3 (exact match, including the `2(1+E)`, `(F-E)`, `-(2+E)`, `-F` coefficients of row_nt), verified the `Delta_T` SymPy/Mathematica forms are algebraically identical, and verified the hand-typed minor determinant equals `1+chi` (magnitude; the script's column ordering η,μ,T differs from the paper's T,η,μ by one swap, a sign convention, not a defect). The verdict is `findings` because the stage's *distinctive* deliverable — the finite monomial-equality ⇔ linear-row claim (D1) — is asserted but never derived from the §1 monomials in either engine (the Mathematica "ratio" check is tautological), and the `.wl` is a transliteration of the `.py`. Both are script-side, fixable by Codex (build the finite ratios from the monomials, independently per engine); there is no paper↔script disagreement (rows and solve match the paper), so no `paper_misalignment` and no user routing. Not stop-cold: the underlying math is correct and the fix only adds the missing derivation.

## Self-test notes

I checked: (1) Variable-independence/derivative traps — my proposed fix uses no `diff`; for monomials `log(∏(r_i)^{a_i}) = Σ a_i log r_i = Σ a_i Δ_i = row` exactly, with positive ratio symbols (`expand_log(..., force=True)` collapses it). (3) Trivial-case pre-check: all `r_i=1` ⇒ Δ=0, row=0, residual 0; `r_γ=e` only ⇒ `log(C_tr ratio)=1+δ=row_tr`, residual 0 — both pass. (5) Paper round-trip: the monomials are reconstructed verbatim from notes §1 (the `\pi^2/L^2` and `\sigma` constants cancel in the ratios, introducing no new constant), so the fix introduces no new paper_misalignment. No integrals, so the parity trap is N/A; no new files, so the path-spec trap is N/A.

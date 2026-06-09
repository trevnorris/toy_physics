---
unit_id: 165
batch: V.1
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-06-08T00:00:00Z
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
  notes_stage_files: [/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage165_exact_branch_drifts.md]
  paper_appendix: present
---

# Audit unit 165 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_165.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage165_exact_branch_drifts.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part05.tex` (rows: line 61 status row; lines 334-338 narrative; eq:app-part05-LW-co-transport)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage165_exact_branch_drifts_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage165_exact_branch_drifts_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage165_exact_branch_drifts_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage165_exact_branch_drifts_mathematica_audit.txt`

## What the paper claims

The card `\stagefield{Output}` states verbatim: *"Shows \(\delta\ln L_W=\delta\ln a\) and reduces the lower compensated branch to four irreducible drifts \((\delta\ln\mathcal Z_q,\delta\ln\rho_w,\delta\ln c_s,\delta\ln a)\)."* The notes go far beyond the terse card. They establish: (1) the D/N co-transport law `δln L_W = δln a` (notes §2); (2) the exact `δln v_w0` and `δln T_m` log-drift laws solved from the fixed-`r` and fixed-`g` conditions inherited from Stage 164 (notes §3-4); (3) the product/ratio split `δln(v_w0/T_m)=2δln c_s−δln a` and `δln(v_w0 T_m)=δln Z_q+3δln c_{s,w}−δln ρ_w−4δln a` (notes §6); (4) the `n=5` wall-EOS identity `δln c_{s,w}=2δln ρ_w` and the collapse to four coordinates (notes §7); and (5) the identical vanishing of the Stage-164 off-family coordinate `δ_⊥=0` (notes §8). Notes §5-6 also report closed-form *amplitude* laws for `T_m` and `v_w0` together with five numeric branch constants (`r_*≈1.77799…`, `g_*≈0.758035…`, `T_m≈1.27158…`, `v_w0≈−1.14288…`, `v_w0/T_m≈−0.89878…`, `v_w0 T_m≈−1.45328…`). These amplitude/numeric deliverables are NOT mentioned on the card; the card's load-bearing claim is the log-drift collapse.

## What the script claims to verify

Both engines verify the same five log-drift identities (docstring items 1-5). Concretely: (A) at fixed `r`, `a·d/da[ln L_W] = 1`, i.e. `δln L_W = δln a`; (B) the 2×2 linear system `eqR`/`eqG` (the fixed-`r` and fixed-`g` conditions) solves to the stated `dv`/`dT` log-drift laws; (C) after `dLW→da`, the difference `dvDN−dTDN` equals `2dcs−da` and the sum equals `dZ+3dcsw−drho−4da`; (D) after `dcsw→2drho`, the `n=5` product law equals `dZ+5drho−4da`; and (E) the Stage-164 channels evaluated at the solved `(dv,dT)` print as `0` — explicitly labeled in-comment as a solver-consistency print, *not* asserted. The scripts operate entirely on the symbolic log-drift coordinates; they do NOT verify the §5-6 closed-form amplitudes or numeric branch constants.

## Paper ↔ script cross-check

| Paper deliverable | Script-side check | Status |
|---|---|---|
| `δln L_W = δln a` (card Output; notes §2; appendix eq:app-part05-LW-co-transport) | py:44 / wl:33 `expect_zero("d ln L_W - d ln a at fixed r_*", dlog_LW - 1)` | match |
| `δln v_w0`, `δln T_m` log-drift laws (notes §3-4) | py:49 / wl:37 solve of eqR/eqG; printed py:61-62 / wl:49-50 | match (printed, exercised downstream by C/D) |
| ratio law `2dcs−da` (notes §6) | py:67 / wl:54 expectZero | match |
| product law `dZ+3dcsw−drho−4da` (notes §6) | py:68 / wl:55 expectZero | match |
| `n=5` EOS `δln c_{s,w}=2δln ρ_w` + collapse (notes §7) | py:71-79 / wl:58-64 subst + expectZero | match |
| `δ_⊥=0` Stage-164 collapse (notes §8) | py:86-89 / wl:71-74 — solver-consistency print only, NOT asserted (so labeled) | partial (intentional; correctly disclaimed) |
| Closed-form `T_m`, `v_w0` amplitudes + 5 numeric branch constants (notes §5-6) | none | missing (not a card claim; see F2) |

`paper_alignment: partial` — all card-level (log-drift) deliverables are faithfully exercised; the notes-only amplitude/numeric deliverables have no script check, and the second engine is a transliteration (F1).

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored? |
|---|---|---|---|---|---|
| A1 | sympy | 44 | `expect_zero(dlog_LW - 1)` | `δln L_W=δln a` | yes |
| A2 | sympy | 67 | `expect_zero(ratio_law - (2dcs - da))` | ratio law §6 | yes |
| A3 | sympy | 68 | `expect_zero(prod_law - (dZ+3dcsw-drho-4da))` | product law §6 | yes |
| A4 | sympy | 76-79 | `expect_zero(n=5 product - (dZ+5drho-4da))` | n=5 collapse §7 | yes |
| A5 | sympy | 88-89 | `print(channel_g/channel_r)` — no assert | δ_⊥=0 §8 | no (print only, disclaimed) |
| A6 | math | 33 | `expectZero[a*D[Log[lwLaw],a] - 1]` | `δln L_W=δln a` | yes |
| A7 | math | 54 | `expectZero[ratioLaw - (2dcs - da)]` | ratio law §6 | yes |
| A8 | math | 55 | `expectZero[prodLaw - (dZ+3dcsw-drho-4da)]` | product law §6 | yes |
| A9 | math | 62-65 | `expectZero[n=5 product - (dZ+5drho-4da)]` | n=5 collapse §7 | yes |
| A10 | math | 71-74 | `Print[channelG/channelR]` — no assert | δ_⊥=0 §8 | no (print only, disclaimed) |

A1-A4 / A6-A9 are non-tautological: each subtracts the *independently transcribed* notes target from a *solver-derived* quantity. The `dv`/`dT` laws come out of `solve`/`Solve` on eqR/eqG; the ratio/product/n=5 targets are written separately and must match. If the notes' boxed RHS were wrong by any factor, these would fail. A5/A10 are correctly identified by the authors as tautological (substituting a linear solution back into its own system) and are deliberately demoted to non-asserting prints with an explanatory comment — this is honest, not a defect.

## Findings

### F1 — mathematica_transliteration

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage165_exact_branch_drifts_mathematica_audit.wl:35-74`
- (mirror) `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage165_exact_branch_drifts_sympy_audit.py:47-89`

**What's wrong:**
The `.wl` is a line-by-line port of the `.py`, not an independent re-derivation. The same variable choreography, same intermediate steps, same ordering, and even the same comment text appear in both. Corresponding sections:

1. The linear system + solve are identical step-for-step:
   - py:47-49 `eq_r = sp.Eq(dZ + 2*dcs + 3*dcsw - drho - 2*dv - 2*da - 3*dLW, 0)` … `sol = sp.solve([eq_r, eq_g], [dv, dT], dict=True)[0]`
   - wl:35-37 `eqR = dZ + 2*dcs + 3*dcsw - drho - 2*dv - 2*da - 3*dLW == 0;` … `sol = Solve[{eqR, eqG}, {dv, dT}, Reals][[1]];`

2. The DN substitution and product/ratio construction match exactly:
   - py:57-66 `dv_DN = ...subs(dLW, da)` … `ratio_law = sp.simplify(dv_DN - dT_DN)` / `prod_law = sp.simplify(dv_DN + dT_DN)`
   - wl:45-53 `dvDN = FullSimplify[dvSol /. dLW -> da, ...]` … `ratioLaw = ...[dvDN - dTDN, ...]` / `prodLaw = ...[dvDN + dTDN, ...]`

3. The Stage-164 channel block is identical down to the comment wording:
   - py:82-85 comment `they vanish by construction (substituting a linear solution into its own system), so this is a consistency print of the solver, not an independent verification of delta_perp = 0. Reported, not asserted.`
   - wl:68-70 comment `they vanish by construction, so this is a solver-consistency print, not an independent verification of delta_perp = 0. Reported, not asserted.`

**Why this matters:**
The second-engine policy requires both engines to confirm the result by an independent route, not echo the same algebra. Here the `.wl` adds essentially no independent confirmation: it solves the same 2×2 system from the same two transcribed equations. NOTE on severity: the stage's mathematical content is pure 8-variable linear algebra with one canonical solution path, so the *scope* for genuine independence is small — the only independent thing a second engine could do is re-derive eqR/eqG from the §3/§4 structural log-derivatives of `K_s K_q/λ²` and `g_q K_s/(g_s λ)` rather than copying the pre-assembled linear coefficients. That is why this is rated **low**, not medium: the result is robustly true and engine-agreement holds; the finding is a policy-fidelity flag, not a correctness flag.

**Required change:**
Make the `.wl` reach `eqR`/`eqG` by its own route rather than copying the coefficient vectors. Concretely, define the two Stage-164 imbalance channels symbolically as log-derivative sums and let Mathematica assemble the coefficients, e.g. build `chanR = δln(K_s K_q/λ²)` and `chanG = δln(g_q K_s/(g_s λ))` from the §3/§4 substitutions for `K_s, K_q, λ, J_s` in terms of `(Z_q, ρ_w, c_s, c_{s,w}, v_w0, a, L_W)`, then set each to zero to *obtain* eqR/eqG instead of writing them as literal coefficient sums at wl:35-36. Then solve as before. This gives the second engine a genuinely different derivation of the same linear system.

**Verification:**
The verifier confirms wl:35-36 no longer hardcode the coefficient vectors verbatim from py:47-48, that the assembled `eqR`/`eqG` still match the SymPy ones (printed for comparison), and that the script still exits 0 with all four `expectZero` checks PASS.

### F2 — insufficient_verification

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage165_exact_branch_drifts_sympy_audit.py:1-98`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage165_exact_branch_drifts_mathematica_audit.wl:1-89`

**What's wrong:**
Notes §5-6 box closed-form amplitude laws and report five numeric branch constants as results of the stage:
- `T_m ≈ 1.2715890393387603 · …` (notes §5, line 232)
- `v_w0 ≈ -1.1428896163056477 · …` (notes §5, line 241)
- `v_w0/T_m ≈ -0.8987885086678338/q_* · c_s²/a` (notes §6, line 310)
- `v_w0 T_m ≈ -1.4532859092683434/q_* · …` (notes §6, line 318)
built from `r_* ≈ 1.77799353547498` and `g_* ≈ 0.758035078944663` (notes §5, lines 198-200).
Neither engine evaluates any of these. The scripts verify only the *logarithmic* drift relations; the closed-form prefactor algebra (e.g. `3√10·3^{3/4}/(5π g_* (1+r_*²)^{1/4})` → `1.27158…`) and the product/ratio numeric constants are entirely unchecked. A transcription slip in any of those notes constants would not be caught.

**Why this matters:**
This is NOT a card-Output deliverable (the card's claim is the log-drift collapse, which IS fully verified), so it does not lower the card-level verdict. But the notes present these amplitudes as boxed stage results, and there is no engine that confirms the closed-form prefactor evaluates to the quoted decimal. Severity is **low** because (a) the load-bearing card claim is covered, and (b) these are downstream-numeric niceties rather than the irreducible-drift theorem.

**Required change:**
Add to BOTH engines a small numeric block that (1) pins `r_* = 1.77799353547498`, `g_* = 0.758035078944663`; (2) evaluates the §5/§6 closed-form prefactors `T_m_pref = 3*sqrt(10)*3^(3/4)/(5*pi*g_**(1+r_*^2)^(1/4))`, `v_pref = 9*sqrt(10)*3^(1/4)*r_*/(20*(1+r_*^2)^(3/4))`, `ratio_pref = sqrt(3)*pi*g_*r_*/(4*sqrt(1+r_*^2))`, `prod_pref = 81*r_*/(10*pi*g_**(1+r_*^2))`; and (3) `expect_close`/numeric-`expectZero` each against the notes decimals `1.2715890393387603`, `1.1428896163056477`, `0.8987885086678338`, `1.4532859092683434` to a tight tolerance (e.g. `1e-12`). This makes the §5-6 numeric deliverables script-verified.

**Verification:**
New numeric checks appear in both scripts; the SymPy and Mathematica outputs each print the four prefactor values and PASS the four numeric comparisons against the notes decimals; both scripts exit 0.

## Independent-derivation check (Mathematica)

The `.wl` is a transliteration of the `.py` (see F1 for three quoted corresponding section pairs, including verbatim-shared comment text). Confidence: **high** that it is a port. Mitigating context: the stage is single-path linear algebra, so the achievable independence is limited; this drives the F1 severity to low rather than medium. Engine outputs agree exactly (see below), so there is no correctness concern — only a second-engine-policy concern.

## Engine cross-check

The two saved outputs agree on every emitted quantity (modulo SymPy vs Mathematica term ordering / normalization, which is cosmetic):

| Quantity | SymPy output | Mathematica output |
|---|---|---|
| `L_W` law | `sqrt(3)*pi*a*sqrt(r**2 + 1)/6` | `(a*Pi*Sqrt[1 + r^2])/(2*Sqrt[3])` (equal forms) |
| `δln v_w0` (pre-DN) | `-3*dLW/2 + dZ/2 - da + dcs + 3*dcsw/2 - drho/2` | `(-2*da + 2*dcs + 3*dcsw - 3*dLW - drho + dZ)/2` (equal) |
| `δln v_w0` (post-DN) | `dZ/2 - 5*da/2 + dcs + 3*dcsw/2 - drho/2` | `(-5*da + 2*dcs + 3*dcsw - drho + dZ)/2` (equal) |
| `δln v_w0` (n=5) | `dZ/2 - 5*da/2 + dcs + 5*drho/2` | `(-5*da + 2*dcs + 5*drho + dZ)/2` (equal) |
| ratio/product/n=5 checks | all `= 0` | all `= 0`, PASS |
| Stage-164 channels | `0`, `0` (print) | `0`, `0` (print) |

`engines_agree: true`.

## Verdict justification

`findings`. Attacks attempted that FAILED to break the scripts: (1) tested A1-A4/A6-A9 for tautology — each subtracts a separately-transcribed notes target from a solver-derived quantity, so a wrong notes coefficient would fail; non-tautological. (2) Re-derived the 2×2 solve by hand from eqR/eqG and confirmed `dv = dZ/2 + dcs + 3dcsw/2 − drho/2 − da − 3dLW/2` matches the notes §3 box and output line 10. (3) Checked symbol domains: `r,a` positive-real (correct for a tube length/radius), drift variables real-unrestricted (correct for linearized log-drifts); the `a*D[Log[lwLaw],a]` derivative genuinely depends on `a` (not identically zero), so A1/A6 are live. (4) Checked the §8 `δ_⊥=0` block is correctly disclaimed as a non-asserting solver-consistency print rather than masquerading as a proof — honest. The card-level Output claim (`δln L_W=δln a` + four-drift collapse) is fully and faithfully verified by both engines. The two findings are: F1 (low) the `.wl` is a transliteration adding no independent route, and F2 (low) the notes-only §5-6 closed-form amplitudes and five numeric branch constants are not verified by either engine. I confirm I read the card, the notes, and the appendix rows before the scripts; the script's verified claim matches the card's claim (the gap is at the notes-detail level, not the card level). No `paper_misalignment` (no value differs between script and docs; the unverified §5-6 values live correctly in the notes and so reconcile as MATCH — see below).

## Value Reconciliation (pass-2 augmentation)

The scripts emit symbolic log-drift relations and the `L_W` law form; they emit NO numeric result constants. Every emitted symbolic deliverable is reflected in the notes (and, for the load-bearing ones, the card). The numeric branch constants in notes §5-6 are NOT emitted by the scripts (that is the F2 gap), so they are not part of the script→doc reconciliation — they appear here only to confirm the notes carry them.

| value (script-emitted) | source (py/wl + output line) | .tex/.md location | status |
|---|---|---|---|
| `δln L_W = δln a` (`L_W = √3 π a √(r²+1)/6`) | py:40-44 / wl:31-33; sympy.txt:5-7, math.txt:5-7 | notes §2 lines 66-74 (`L_W=π a/2 √((1+r²)/3)`, box `δln L_W=δln a`); card Output line 15; appendix eq:app-part05-LW-co-transport (line 337) | MATCH |
| `δln v_w0 = dZ/2 + dcs + 3dcsw/2 − drho/2 − da − 3dLW/2` | py:49,54 / wl:37,42; sympy.txt:10, math.txt:12 | notes §3 box lines 116-123 | MATCH |
| `δln T_m = dZ/2 − dcs + 3dcsw/2 − drho/2 − da − dLW/2` | py:49,55 / wl:37,43; sympy.txt:11, math.txt:13 | notes §4 box lines 168-176 | MATCH |
| `δln(v_w0/T_m) = 2 dcs − da` | py:67 / wl:54; sympy.txt:16, math.txt:20 | notes §6 box lines 269-276; §7 lines 359-364 | MATCH |
| `δln(v_w0 T_m) = dZ + 3dcsw − drho − 4da` | py:68 / wl:55; sympy.txt:17, math.txt:22 | notes §6 box lines 291-303 | MATCH |
| `n=5`: `δln c_{s,w} = 2 δln ρ_w` → product `dZ + 5drho − 4da` | py:71,76-79 / wl:58,62-65; sympy.txt:22, math.txt:30 | notes §7 boxes lines 336-339, 366-375 | MATCH |
| `δ_⊥ = 0` (channels) | py:86-89 / wl:71-74 (print, not asserted); sympy.txt:23-24, math.txt:36-37 | notes §8 box lines 398-402 | MATCH (printed; consistency only) |

Notes-only numeric constants NOT script-emitted (relevant to F2 only, not reconciled here): `r_*≈1.77799353547498`, `g_*≈0.758035078944663`, `T_m≈1.2715890393387603`, `v_w0≈−1.1428896163056477`, `v_w0/T_m≈−0.8987885086678338`, `v_w0 T_m≈−1.4532859092683434` (notes §5-6). No `168π²/100π²/4107` Family-1-radius constants appear anywhere in this stage's scripts or notes; that stale-constant class does not apply here.

INTERNAL scaffolding (no finding, not expected in prose): `expect_zero`/`expectZero`/`banner`/`pass`/`fail`/`fmt` helpers, the `expectZero` residual prints (`= 0`), `eqR`/`eqG` intermediate linear system, `dvDN`/`dTDN`/`dvN5`/`dTN5` intermediate substitutions, PASS flags.

reconciliation: complete; 7 script-emitted deliverable values checked, 0 misaligned.

## Self-test notes

Checked the self-test traps: (1) Variable independence — `a*D[Log[lwLaw],a]` (wl:33) / `a*diff(log(LW_law),a)` (py:41): `LW_law` genuinely depends on `a` (linear factor), so the derivative is `1`, not identically zero; A1/A6 are live, not silent-pass. (2) No unbounded integrals in this stage, so the parity trap is N/A. (3) Trivial-case pre-check — set all drift symbols to 0 except one (e.g. `da=1`): ratio law gives `−1` matching `2·0−1`, product law gives `−4` matching `0+0−0−4`; residuals reduce to 0 as asserted. (4/5) F2's proposed numeric block uses exactly the notes constants/prefactors (no new constant introduced), so it cannot create a fresh `paper_misalignment`; F1 only reroutes the `.wl` derivation of eqR/eqG without changing any verified value. Both findings are low-severity and script-side; a directive is written.

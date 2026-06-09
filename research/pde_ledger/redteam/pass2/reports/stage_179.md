---
unit_id: 179
batch: V.2
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-06-08T00:00:00Z
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
  notes_stage_files: [moving_throat_pde_stage179_transfer_shape_theorem.md]
  paper_appendix: present
---

# Audit unit 179 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_179.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage179_transfer_shape_theorem.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part05.tex` (rows 89, 671-687, 1464 reference this unit / its anchor)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage179_transfer_shape_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage179_transfer_shape_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage179_transfer_shape_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage179_transfer_shape_mathematica_audit.txt`

## What the paper claims

Stage 179 (`MTDC-T9.2`) proves the **wall-normalized transfer-shape theorem**. The card's `\stagefield{Output}` states verbatim: "Factors each outgoing port as \(N_{A,0}^{(r)}=K_A\mathcal T_{A,r}^2\) and gives \(\Xi_1=2\sum_r\rho_r^{(N)}\tau_r\)." The notes expand this into four deliverables: (1) the exact wall-normalized factorization \(N_0^{(r)}/K=[(\widehat G_W+\widehat R\widehat G_U)/(1-\widehat R^2)]^2 = \mathcal T_r^2\); (2) the central slope identity \(\nu_r=\kappa_1+2\tau_r\) with the explicit closed form \(\tau_r=\widehat\alpha_r\mathfrak w_r+\widehat\beta_r(\mathfrak u_r+\mathfrak c_r)+\frac{2\widehat R_r^2}{1-\widehat R_r^2}\mathfrak c_r\); (3) the collapse \(\Xi_1=\sum_r\rho_r^{(N)}(\nu_r-\kappa_1)=2\sum_r\rho_r^{(N)}\tau_r\) given \(\sum_r\rho_r^{(N)}=1\); (4) equivalence to the Stage 176/177/178 slippage language, \(\tau_r=\mathfrak m_r+\frac{\mathcal I_r}{1+\mathcal I_r}\mathfrak i_r+\frac{\mathcal H_r}{1-\mathcal H_r}\mathfrak h_r\). The appendix (eqs. `app-part05-port-transfer-shape`, `app-part05-Xi1-transfer-shape-sum`) restates (1) and (3) identically. This is an exact-closure (non-checkpoint) algebraic-identity stage.

## What the script claims to verify

The SymPy docstring (lines 4-11) lists four checks that map one-to-one onto the four notes deliverables. The assertions are: (A1) `N0/K - T^2 == 0` after building `N0 = P^2/Delta^2` from `P = OU2*GW + R*GU`, `Delta = OU2*OW2 - R^2` and the wall-normalized `T` (line 53); (A2) `nu_direct - (kappa1 + 2*tau) == 0`, where `nu_direct` is independently obtained by `sp.series(sp.log(N0A), eps, 0, 2)` of the perturbed coefficient and `tau` is the notes' closed form (line 91); (A3) `tau - slippage form == 0` and `(nu-kappa1) - 2*tau_slippage == 0` (lines 102-103); (A4) `Xi_1 - 2 weighted tau == 0` for three weights with `rho3 = 1 - rho1 - rho2` (line 112). The Mathematica `.wl` runs the same five `expectZero` checks in the same order with the same constructions.

## Paper ↔ script cross-check

| paper deliverable | script-side check | status |
|---|---|---|
| (1) factorization `N0^(r)=K T_r^2` | A1 `N0/K - T^2 == 0` (py:53) | match |
| (2) slope identity `nu_r = kappa1 + 2 tau_r` with closed-form tau | A2 `nu_direct - (kappa1+2 tau) == 0` (py:91); tau closed form (py:84) | match |
| (4) slippage equivalence `tau = m + I/(1+I) i + H/(1-H) h` | A3 (py:102-103) | match |
| (3) defect collapse `Xi_1 = 2 sum_r rho_r tau_r`, `sum rho=1` | A4 `Xi_1 - 2 weighted tau` (py:112) | partial (see F2) |

`paper_alignment: aligned` — every paper Output quantity (`N_{A,0}^{(r)}=K_A\mathcal T_{A,r}^2`, `\Xi_1=2\sum_r\rho_r^{(N)}\tau_r`) has a faithful script-side counterpart with the exact same form. The only weakness is that A4 is a pure-algebra restatement (F2), not a verification anchored to the physics object; that is a verification-quality issue, not a paper/script disagreement.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 53 | `simplify(N0/K - T^2) == 0` | claim 1 (factorization) | yes |
| A2 | sympy | 91 | `nu_direct - (kappa1+2*tau) == 0` | claim 2 (slope identity) | yes |
| A3a | sympy | 102 | `tau - tau_slippage == 0` | claim 4 (slippage equiv) | yes |
| A3b | sympy | 103 | `(nu-kappa1) - 2*tau_slippage == 0` | claim 4 | yes |
| A4 | sympy | 112 | `Xi - Xi_expected == 0` | claim 3 (defect collapse) | partial (tautological-leaning) |
| M1 | mathematica | 44 | `expectZero[N0/K - T^2]` | claim 1 | yes |
| M2 | mathematica | 75 | `expectZero[nuDirect-(kappa1+2 tau)]` | claim 2 | yes |
| M3a | mathematica | 84 | `expectZero[tau - tauSlippage]` | claim 4 | yes |
| M3b | mathematica | 85 | `expectZero[(nu-kappa1)-2 tauSlippage]` | claim 4 | yes |
| M4 | mathematica | 93 | `expectZero[xi - xiExpected]` | claim 3 | partial |

A1/A2/A3 are genuine: A1 builds `T` from scratch via the wall-normalized substitutions and confronts it with `P^2/Delta^2`; A2 contrasts an **independently obtained** logarithmic slope (`series`/`D` of `log N0A`) against the notes' closed-form `kappa1 + 2*tau` — these are constructed from different routes, so the cancellation is informative. A4 is the soft spot: both `Xi` and `Xi_expected` are written directly in terms of `kappa1 + 2*tau_i`, so the check reduces to `sum_i rho_i (kappa1 + 2 tau_i) - kappa1 == 2 sum_i rho_i tau_i` given `sum rho_i = 1` — an algebraic tautology in `tau1,tau2,tau3` that does not exercise the physical `nu_r = kappa1 + 2 tau_r` or the weights `rho_r^(N)` deriving from the actual port data.

## Findings

### F1 — stale_output

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage179_transfer_shape_sympy_audit.txt:3` (and `:68`)
- `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage179_transfer_shape_mathematica_audit.txt:3` (and `:14,:22,:41`)

**What's wrong:**
The committed outputs are older than the scripts: outputs mtime `May 30 01:39`, both scripts mtime `Jun 3 15:59`. The content confirms staleness — the saved banners read `STAGE 162 — WALL-NORMALIZED TRANSFER-SHAPE THEOREM` (both `.txt:3`) and `Stage 176/160/161 slippage formulas` (sympy `.txt:14` header / mathematica `.txt:22`), whereas the current scripts emit `STAGE 179 — ...` (py:32, wl:26) and `Stage 176/160/161` is only in the banner string (py:93, wl:77). The mathematica `.txt:14` also still prints the old `nu_direct` form. Because the scripts were edited after these outputs were captured, the transcripts no longer reflect the current source.

**Why this matters:**
A reviewer relying on the committed transcript sees stage label `162`, not `179`, and an out-of-date slope expression. The captured PASS lines may not correspond to the current assertion text.

**Required change:**
Re-run both scripts and recommit the transcripts so the banner reads `STAGE 179` and the slope output reflects the current source. (Informational; the orchestrator's independent re-run regenerates these.)

**Verification:**
After re-run, `scripts/output/..._sympy_audit.txt:3` and `mathematica/output/..._mathematica_audit.txt:3` should both read `STAGE 179 — WALL-NORMALIZED TRANSFER-SHAPE THEOREM`.

### F2 — insufficient_verification

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage179_transfer_shape_sympy_audit.py:105-112`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage179_transfer_shape_mathematica_audit.wl:87-93`

**What's wrong:**
The "Weighted defect identity" block constructs both sides from the SAME symbolic template:
```
Xi = rho1*(kappa1 + 2*tau1) + rho2*(kappa1 + 2*tau2) + rho3*(kappa1 + 2*tau3) - kappa1
Xi_expected = 2*(rho1*tau1 + rho2*tau2 + rho3*tau3)
expect_zero("Xi_1 - 2 weighted tau", Xi - Xi_expected)
```
with `rho3 = 1 - rho1 - rho2` and `tau1,tau2,tau3` free symbols. The check therefore reduces to the trivial algebraic identity `sum_r rho_r (kappa1 + 2 tau_r) - kappa1 = 2 sum_r rho_r tau_r` whenever `sum_r rho_r = 1`. It cannot fail for any physics, and it does not connect `tau_r` here to the closed-form `tau` verified by A2, nor `rho_r` to the actual outgoing weights `rho_r^(N)`. The genuine content of deliverable (3) — that the per-port `nu_r - kappa1` equals `2 tau_r` for the ACTUAL port slope — is already carried by A2; A4 is only the weighted-sum bookkeeping. This is `insufficient_verification`, not a paper disagreement.

**Why this matters:**
A4 reads as a substantive verification of `Xi_1 = 2 sum rho tau` but is in fact a `sum rho = 1` algebra restatement. If A2 ever regressed, A4 would still pass, giving false confidence in the collapse step.

**Required change:**
Strengthen A4 so the `tau_i` it sums are the SAME closed-form `tau` object verified in A2 (i.e. substitute distinct port instances of the A2-validated `tau`, or at minimum feed `nu_i = kappa1 + 2*tau_i` from the per-port `nu_direct` construction rather than re-declaring fresh free `tau_i` symbols). Then `Xi_1 = sum_r rho_r (nu_r - kappa1)` collapsing to `2 sum_r rho_r tau_r` would exercise the physical slope, not just the linearity of weighted averages.

**Verification:**
The new A4 should reference the `tau` defined at py:84 / wl:70 (or `nu_direct`) rather than independent `tau1,tau2,tau3` symbols introduced at py:108 / wl:88; the residual must still simplify to 0.

## Independent-derivation check (Mathematica)

The `.wl` is a **transliteration / port** of the `.py`, not an independent derivation. Evidence:

1. **Identical object choreography.** Both build `P/p = ou2*gw + r*gu`, `Delta/delta = ou2*ow2 - r^2`, `N0/n0 = P^2/Delta^2`, then the SAME three wall-normalized substitutions `gWh = gw/(Sqrt[k]*ow2)`, `gUh = gu/(Sqrt[k]*Sqrt[ou2]*Sqrt[ow2])`, `rHat = r/(Sqrt[ou2]*Sqrt[ow2])`, and `t = (gWh + rHat*gUh)/(1 - rHat^2)`. Compare py:40-48 with wl:32-39 line-for-line.

2. **Same five checks, same order, same banners.** `N0/K - T^2`, then `nu_direct - (kappa1+2 tau)`, then `tau - slippage form` + `(nu-kappa1) - 2 tau_slippage`, then `Xi_1 - 2 weighted tau`. Even the banner strings are byte-identical, including the stale `"Stage 176/160/161 slippage formulas"` (py:93 ≡ wl:77) and the entire carry-forward `Print` block (py:114-119 ≡ wl:96-101).

3. **Same slope object.** The one place the engines differ mechanically is the slope: SymPy uses `sp.series(sp.log(N0A), eps, 0, 2).removeO().coeff(eps,1)/lam` (py:73); Mathematica uses `(D[Log[n0A], eps] /. eps -> 0)/lam` (wl:61). But this is the SAME constructed object `Log[N0A]` built from the SAME perturbed `n0A = (ou2A*gwA + rA*guA)^2/(ou2A*ow2A - rA^2)^2` with the SAME `(1 + eps*lam*·)` choreography (py:61-70 ≡ wl:51-60). Series-coefficient vs. first-derivative-at-0 of the same log is the same extraction, not an independent route.

This is a `mathematica_transliteration` candidate. However, the underlying algebra is a closed-form polynomial/rational identity (no special functions, no integration, no branch choices), so a second engine has very little room to "derive it differently" — the realistic independent value of the `.wl` here is engine-cross-validation of SymPy's `simplify`, which it does provide. I record it as a transliteration call below and as a low-severity note, but do not raise it as a blocking finding given the identity-only nature of the stage; see Verdict.

## Engine cross-check

Both engines pass all five checks with residual exactly 0 (sympy `.txt`: `N0/K - T^2 = 0`, `nu_direct - (kappa1 + 2 tau) = 0`, `tau - slippage form = 0`, `(nu-kappa1) - 2*tau_slippage = 0`, `Xi_1 - 2 weighted tau = 0`; mathematica `.txt`: same five with `PASS:` and `Stage 179 Mathematica audit passed.`). The printed `T`, `P`, `Delta`, `nu_direct`, `tau`, `nu_expected` forms agree symbolically between the two transcripts (e.g. `T = (GU*R + GW*OU2)/(sqrt(K)*(OU2*OW2 - R**2))` vs `(gw*ou2 + gu*r)/(Sqrt[k]*(ou2*ow2 - r^2))`). No `engine_disagreement`.

## Verdict justification

`findings` (2, both low severity). The core mathematics is sound and the paper alignment is exact: A1 verifies the factorization `N0^(r)=K T_r^2`, A2 verifies the central identity `nu_r = kappa1 + 2 tau_r` by confronting an independently extracted log-slope against the notes' closed-form `tau`, and A3 verifies the slippage-language equivalence — these are non-tautological and trace cleanly to the Output and notes. I attacked the symbol domains (all positivity assumptions `K,OU2,OW2,GW,GU>0` are physically justified — they are squared frequencies and positive port amplitudes; `R` correctly left unrestricted real so `1 - Rhat^2` and `Delta` can change sign), the slope extraction (genuinely independent of `tau`), and the factorization (built from primitives, not hardcoded). Those held. The two findings are quality issues: F1 stale outputs (banner `162`, regenerated by the orchestrator re-run) and F2 the weighted-defect check being a `sum rho = 1` algebra tautology that does not feed the A2-validated `tau`. The `.wl` is a port of the `.py` (see independent-derivation section), but because the stage is a pure closed-form rational identity, the second engine still provides legitimate cross-validation of SymPy's simplifier; I did not escalate transliteration to a blocking finding, only noted it.

## Self-test notes

Checked: (1) Variable independence — `nu_direct` differentiates `Log[N0A]` w.r.t. `eps`, and `N0A` genuinely depends on `eps` through all six `(1+eps*lam*·)` factors, so the derivative is non-trivially nonzero (the printed `nu_direct` is a large rational, confirming). (2) No unbounded integrals — pure algebra, parity trap n/a. (3) Trivial-case: setting all slopes `gW=gU=oU=oW=rr=0` gives `nu_direct = -kappa1*(...)/(...) `? No — `kappa1` enters `KA` only and cancels in `log(N0A)` since `N0A` has no `KA`; with all port slopes 0 `nu_direct→0` and `kappa1+2*tau→0` consistently, so A2's zero is not accidental. (4) Path specs n/a (no missing-script directive). (5) Paper round-trip — F2's proposed strengthening reuses the existing `tau` (py:84), introduces no new constant, so no new `paper_misalignment`.

## Value Reconciliation (pass-2 augmentation)

The scripts emit only closed-form **symbolic** deliverables (no numeric constants, no figure-of-merit numbers). Each is the load-bearing result of the stage and is checked against the card `\stagefield{Output}` (stage_179.tex:15), the notes boxed equations, and the appendix eqs.

| value | source (py/wl + output line) | .tex/.md location | status |
|---|---|---|---|
| `N0^(r)/K = T^2`, i.e. `N_{A,0}^{(r)}=K_A T_{A,r}^2` | py:53 / wl:44; out `N0/K - T^2 = 0` (sympy .txt:8) | tex:15 Output; md:115-117 (boxed `N0^(r)=K T_r^2`); appendix:674 eq. `app-part05-port-transfer-shape` | MATCH |
| `T_r = (Ghat_W + Rhat Ghat_U)/(1 - Rhat^2)` | py:48 / wl:39; printed `T = ...` (sympy .txt:7) | tex:15 (implicit in `\mathcal T_{A,r}`); md:107-109 (boxed); appendix:676-678 | MATCH |
| `nu_r = kappa1 + 2 tau_r` | py:85,91 / wl:71,75; out `nu_direct - (kappa1+2 tau)=0` | md:222-228 (boxed, "central identity"); (terse on card — lives in notes) | MATCH |
| `tau_r = alpha w + beta(u+c) + 2 Rhat^2/(1-Rhat^2) c` | py:84 / wl:70 | md:195-203 (boxed closed form) | MATCH |
| `tau_r = m + I/(1+I) i + H/(1-H) h` (slippage form) | py:101 / wl:83; out `tau - slippage form=0` | md:298-306 (boxed) | MATCH |
| `Xi_1 = 2 sum_r rho_r^(N) tau_r` | py:111 / wl:92; out `Xi_1 - 2 weighted tau=0` | tex:15 Output; md:248-252 (boxed); appendix:683 eq. `app-part05-Xi1-transfer-shape-sum` | MATCH |

Carry-forward print block (py:114-119 / wl:96-101) restates exactly these six and matches the notes/appendix. Both Output-field quantities (`N_{A,0}^{(r)}=K_A\mathcal T_{A,r}^2` and `\Xi_1=2\sum_r\rho_r^{(N)}\tau_r`) appear verbatim on the card and in the appendix.

INTERNAL (scaffolding, no finding): `P`, `Delta`, `N0` primitives; perturbed `PA/DeltaA/N0A`; `alpha`, `beta`, `w`, `u`, `c`, `Ih`, `Hh`, `m`, `i`, `h`, `nu_expected`, `nu_direct` intermediates; pass/fail flags; the three `rho`/`tau` placeholder symbols in the weighted block.

reconciliation: complete; 6 values checked, 0 misaligned

Note (numbering, NOTE-only — do not fix here): the `\stagefield{Purpose}` carries no `Stage NNN` digit drift, but the script banner string `"Exact equivalence to the Stage 176/160/161 slippage formulas"` (py:93, wl:77) and the docstring "Stage 176/160/161" (py:10) use `160/161`, which appear to be pre-renumber (+? ) labels — the notes call these "Stage 176/177/178" (md:153,286,308). The committed `.txt` banners also still read `STAGE 162` (both `.txt:3`), a stale self-label captured before the May-30→Jun-3 banner fix to `179`. Flagged for the dedicated script/output numbering-band pass, not corrected here.

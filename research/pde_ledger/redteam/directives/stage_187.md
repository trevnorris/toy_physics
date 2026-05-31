---
unit_id: 187
batch: V.2
created_at: 2026-05-30T00:00:00-06:00
findings_count: 2
stop_cold: null
applied: true
applied_at: 2026-05-30T01:34:33-06:00
findings_applied: 2
findings_blocked: 0
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 187

Apply each finding below in order. After applying, append an `## Applied: F<n>` block under that finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

If a finding's required change is ambiguous or unsafe to apply mechanically, append `## Blocked: F<n>` with a question instead — skip that finding, continue with the rest.

Do NOT introduce new features, refactors, or stylistic changes. Edit exactly the file:line ranges named. Do NOT touch paper.tex, notes/, or any prose documents — the red-team only modifies scripts.

After editing, RUN the affected scripts (`python3 <path>` for SymPy, `math -script <path>` for Mathematica) and iterate until they exit 0 with all in-file checks passing.

## F1 — insufficient_verification

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage187_orbit_quotient_closure_sympy_audit.py` — add a derivation block after the row definitions (after line 53, before the matrix `M`).
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage187_orbit_quotient_closure_mathematica_audit.wl` — add an *independently coded* equivalent (after the row definitions, around line 36, before the matrix `m`).

**Issue:**
The stage's distinctive finite-level deliverable — that equality of the three actual monomial *values* between two positive states is exactly the linear log-ratio condition `M_* Delta x = 0` — is never derived from the monomials. The rows `row_tr/nt/eta` are posited directly (sympy 46-48, wl 34-36). The Mathematica `Log[Exp[row]] - row == 0` checks (wl 38-55) are tautological. No assertion in either engine guards the monomial→row exponents.

**Required change:**
Add, in each engine, a block that builds the three invariant *ratios* `Ctilde/C = (xtilde/x)` from the Stage 187 §1 monomial forms in terms of positive primitive ratios, then asserts `log(ratio) - row == 0`. Build the ratio from the primitive ratios (NOT from `Exp[row]`).

Monomial forms (reconstructed from the notes §1; constants `pi^2/L^2` and `sigma` cancel in the ratios and must NOT appear):
- M1: `C_tr,*  = (gamma * c_etaU / K_U)^(1+deltaU) * (T_U / K_U)^(1+chi)`
  ⇒ `Ctr_ratio  = (rG*rC/rU)^(1+deltaU) * (rT/rU)^(1+chi)`,  and `log(Ctr_ratio) == row_tr`.
- M2: `C_nt,*  = (lambda_W^2 * mu_W / (K_eta * K_W^2)) * (gamma^2 * lambda_W^2 / (K_U * K_W))^E * (T_U/K_U)^(-F)`
  ⇒ `Cnt_ratio  = (rL^2*rM/(rEta*rW^2)) * (rG^2*rL^2/(rU*rW))^E * (rT/rU)^(-F)`,  and `log(Cnt_ratio) == row_nt`.
- M3: `epsilon_eta = c_etaU^2 / (K_U * K_eta)`
  ⇒ `eps_ratio  = rC^2/(rU*rEta)`,  and `log(eps_ratio) == row_eta`.

SymPy (insert after line 53):
```python
# Derive the rows from the actual Stage-187 monomials at finite level.
# Positive primitive ratios (xtilde/x) for the eight microscopic variables.
rL, rC, rG, rU, rEta, rW, rM, rT = sp.symbols(
    "r_lambda r_c r_gamma r_U r_eta r_W r_mu r_T", positive=True
)
log_subs = {DL: sp.log(rL), DC: sp.log(rC), DG: sp.log(rG), DU: sp.log(rU),
            DEta: sp.log(rEta), DW: sp.log(rW), DM: sp.log(rM), DT: sp.log(rT)}
Ctr_ratio = (rG * rC / rU) ** (1 + deltaU) * (rT / rU) ** (1 + chi)
Cnt_ratio = (rL ** 2 * rM / (rEta * rW ** 2)) * (rG ** 2 * rL ** 2 / (rU * rW)) ** E * (rT / rU) ** (-F)
eps_ratio = rC ** 2 / (rU * rEta)
expect_zero("log C_tr ratio - row_tr", sp.expand_log(sp.log(Ctr_ratio), force=True) - row_tr.subs(log_subs))
expect_zero("log C_nt ratio - row_nt", sp.expand_log(sp.log(Cnt_ratio), force=True) - row_nt.subs(log_subs))
expect_zero("log eps_eta ratio - row_eta", sp.expand_log(sp.log(eps_ratio), force=True) - row_eta.subs(log_subs))
```
(If `expect_zero`'s `sp.simplify(sp.expand(...))` does not fully collapse the `log`, wrap the argument in `sp.expand_log(..., force=True)` first, as shown.)

Mathematica (insert after line 36, before the existing `ctrRatio`/`cntRatio` block; REPLACE the tautological `Exp[row]` definitions at lines 38-55 with monomial-built ratios so the engine derives independently):
```wolfram
(* Positive primitive ratios (xtilde/x); declare positivity for Log expansion. *)
$Assumptions = $Assumptions && (rL | rC | rG | rU | rEta | rW | rM | rT) > 0;
logSubs = {dl -> Log[rL], dc -> Log[rC], dg -> Log[rG], du -> Log[rU],
           deta -> Log[rEta], dw -> Log[rW], dm -> Log[rM], dt -> Log[rT]};
ctrRatio = (rG*rC/rU)^(1 + deltaStar) * (rT/rU)^(1 + chiStar);
cntRatio = (rL^2*rM/(rEta*rW^2)) * (rG^2*rL^2/(rU*rW))^eStar * (rT/rU)^(-fStar);
epsEtaRatio = rC^2/(rU*rEta);
expectZero["log C_tr ratio - row_tr", PowerExpand[Log[ctrRatio]] - (rowTr /. logSubs)];
expectZero["log C_nt ratio - row_nt", PowerExpand[Log[cntRatio]] - (rowNt /. logSubs)];
expectZero["log epsilon_eta ratio - row_eta", PowerExpand[Log[epsEtaRatio]] - (rowEta /. logSubs)];
```
Then DELETE the old lines 38-55 (`ctrRatio = FullSimplify[Exp[...]]` … through the three `expectZero["log … ratio - row…", Log[…Ratio] - row…]`) so the only ratio checks are the monomial-built ones above. Leave the `ctrRatio /. sol` etc. solve-consistency checks (lines 101-103) intact but rebind them against the new monomial-built ratios with the `logSubs` substitution applied (i.e. `expectZero["C_tr ratio after solve - 1", (ctrRatio /. (logSubs applied to sol)) - 1]`); if that rebinding is awkward, drop the three `… ratio after solve - 1` checks (lines 101-103) entirely, since the `row /. sol == 0` checks at lines 98-100 already cover solve consistency. Note `PowerExpand`/`Log` of negative-power monomials with positive symbols is exact; the `(2+E)`-style exponents are symbolic but the primitive symbols are positive, so `PowerExpand` is valid here.

**Claim manifest:**
- M1: `log(C_tr,* ratio) = (1+deltaU)(Delta_gamma+Delta_c-Delta_U) + (1+chi)(Delta_T-Delta_U) = row_tr`.
- M2: `log(C_nt,* ratio) = 2(1+E)Delta_lambda + 2E*Delta_gamma + (F-E)Delta_U - Delta_eta - (2+E)Delta_W + Delta_mu - F*Delta_T = row_nt`.
- M3: `log(epsilon_eta ratio) = 2*Delta_c - Delta_U - Delta_eta = row_eta`.

**Verification command:**
The verifier will run `redteam exec-sympy 187` and `redteam exec-mathematica 187` and confirm three new PASS lines `log C_tr ratio - row_tr = 0`, `log C_nt ratio - row_nt = 0`, `log epsilon_eta ratio - row_eta = 0` (built from primitive ratios, not `Exp[row]`) appear AND both scripts exit 0.

## Applied: F1

- files_changed:
  - `scripts/moving_throat_pde_stage187_orbit_quotient_closure_sympy_audit.py`
  - `mathematica/moving_throat_pde_stage187_orbit_quotient_closure_mathematica_audit.wl`
- summary: Added primitive-ratio monomial derivations for `C_tr`, `C_nt`, and `epsilon_eta` and checked their logarithms against the finite rows.
- deviation: Dropped the later Mathematica ratio-after-solve checks as allowed because the new ratios are expressed in primitive ratio variables while the retained row-after-solve checks cover solve consistency.

## F2 — mathematica_transliteration

**Target:** `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage187_orbit_quotient_closure_mathematica_audit.wl:38-55`

**Issue:**
The `.wl` is a line-by-line port of the `.py`; its only distinctive content (the `Exp[row]`/`Log` ratio checks at lines 38-55) is tautological, so the second engine provides no real cross-check of the monomial→row step.

**Required change:**
F2 is resolved by implementing the F1 Mathematica change as an **independent** monomial-based derivation (the `PowerExpand[Log[ctrRatio]]` block above), replacing the tautological `Exp[row]` definitions. No additional edit beyond F1 is required: once the Mathematica engine builds the ratios from the §1 monomials by its own algebra (not `Exp[row]`, not copied from the SymPy text), the two engines derive the rows independently. Do NOT merely rename the SymPy variables — the Mathematica monomial expressions must be written natively against the `rL,rC,rG,rU,rEta,rW,rM,rT` primitive ratios.

**Verification command:**
The verifier confirms the `.wl` no longer contains `Exp[(1 + deltaStar)*...]`-style row-echoing ratio definitions and instead builds `ctrRatio/cntRatio/epsEtaRatio` from primitive ratios, and that `redteam exec-mathematica 187` exits 0 with the new PASS lines present.

## Applied: F2

- files_changed:
  - `mathematica/moving_throat_pde_stage187_orbit_quotient_closure_mathematica_audit.wl`
- summary: Replaced the tautological `Exp[row]` Mathematica ratio definitions with native primitive-ratio monomial expressions.
- deviation: none

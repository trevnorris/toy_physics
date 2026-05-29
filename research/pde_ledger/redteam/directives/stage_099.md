---
unit_id: 099
batch: IV.1
created_at: 2026-05-27T16:58:52Z
findings_count: 4
stop_cold: null
applied: false
verification_status: pending
needs_user_resolution: true
---

# Codex directive — unit 099

Apply each non-paper_misalignment finding below in order. After applying, append an `## Applied: F<n>` block under that finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

For `paper_misalignment` findings (F2, F3), do nothing — the orchestrator is holding for user resolution. Do not edit paper.tex, notes/, or scripts to "fix" a paper_misalignment unless the user has explicitly chosen a direction in a follow-up directive.

If a non-paper_misalignment finding's required change is ambiguous or unsafe to apply mechanically, append `## Blocked: F<n>` with a question instead — skip that finding, continue with the rest.

Do NOT introduce new features, refactors, or stylistic changes. Edit exactly the file:line ranges named.

After editing, RUN the affected scripts (`python3 <path>` for SymPy, `math -script <path>` for Mathematica) and iterate until they exit 0 with all in-file checks passing. Getting the scripts to run cleanly is your job; the orchestrator independently re-runs afterward.

Do NOT touch paper.tex, notes/, or any prose documents. The red-team only modifies scripts (except when a follow-up directive explicitly authorizes a paper-side edit after user resolution).

## F2 — paper_misalignment

**Subtype:** script_missing_paper_claim

**Paper side:**
- `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_099.tex:21-25` quote:
  > \stagefield{Checks}{\begin{verificationchecklist}
  > \item Check the static limit \(\epsilon_2=\epsilon_4=0\) returns \(c_{\rm pole}=1/4\).
  > \item Check \(l=0\) and \(l=2\) orthogonality before applying the geometry firewall.
  > \item Check that any support/source success statement still carries the minimal-module hypothesis.
  > \end{verificationchecklist}}

**Script side:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage099_reduced_finish_line_sympy_audit.py` and `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage099_reduced_finish_line_mathematica_audit.wl`: no assertion anywhere exercises (i) the static-limit pole coefficient `c_pole = 1/4`, (ii) `l=0`/`l=2` orthogonality, or (iii) the minimal-module hypothesis carry. The conservative module `Yhat_Q^cons` is bound and printed but not used in any partial-fraction or residue assertion.

## Resolve before fix_loop

The stage card lists three Checks that neither engine exercises. Are those Checks intended to live in stage 099's audit scripts, or did they migrate (during the recent reorder commit `0d09ef6 fully reorder the pde ledger`) from upstream stages 091-098 where the geometry-lane content actually originates?

Possible directions (the user picks one):
- (a) Card is correct → add three assertions to both scripts: (1) residue/pole-coefficient check confirming `c_pole = 1/4` (via `sp.residue(Yhat_Q_cons, omega, Omega_Q)` in sympy and `Residue[yhatCons, {omega, omegaQ}]` in Mathematica, asserting the residue corresponds to `c_pole = 1/4` in the chosen sign convention), (2) explicit `l=0` and `l=2` Legendre orthogonality integrals over `mu in [-1, 1]` (i.e., `Integrate[LegendreP[0,mu]*LegendreP[2,mu], {mu,-1,1}] == 0`, and `Integrate[LegendreP[2,mu]^2, {mu,-1,1}] == 2/5`), (3) an explicit assertion or comment marker noting which upstream stage holds the minimal-module hypothesis being carried in (e.g., reference to stages 091-098). Re-run sympy + Mathematica audits.
- (b) Card lists the wrong checks → strip the three `verificationchecklist` items from `stage_099.tex:22-24` and migrate them to whichever upstream stage card actually tests them (likely one of stages 091-098). No script change.
- (c) The Checks are conceptually correct but already verified upstream → keep them on the card as *carry-forward references*, add explicit upstream citations to the card (e.g., "verified in stage 094 / 097 / etc."), and leave the audit scripts as is.

The orchestrator will not invoke Codex on this unit until the user has chosen a direction.

## F3 — paper_misalignment

**Subtype:** notes_contradicts_script

**Paper side:**
- `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_099.tex:1` quote: `\section[Stage~116]{Stage~116: The Reduced Finish Line After the Geometry-Lane Check}` with `\label{stage:099}` on line 2.

**Script side:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage099_reduced_finish_line_sympy_audit.py:3,25,61` quote: `"""Stage 82 SymPy audit..."""`, `banner("STAGE 82 — REDUCED FINISH LINE")`, `print("\nSTAGE 82 AUDIT PASSED")`.
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage099_reduced_finish_line_mathematica_audit.wl:26,73` quote: `banner["STAGE 082 — REDUCED FINISH LINE"]`, `Print["STAGE 082 AUDIT PASSED"]`.

## Resolve before fix_loop

The stage's three identifiers disagree: the card title says "Stage~116", the audit-unit filename uses `099`, and the script banners say "STAGE 82 / 082" (a stale carry-over from the pre-reorder numbering). Which number is canonical for the script-side banners?

Possible directions (the user picks one):
- (a) Audit-unit number is canonical → rename banners to "STAGE 099" in both engines (sympy lines 3, 25, 61; Mathematica lines 26, 73). Re-run audits to refresh transcripts.
- (b) Card title number is canonical → rename banners to "STAGE 116" in both engines at the same lines.
- (c) Banners should carry both for traceability → use e.g. "STAGE 099 (paper Stage 116, formerly Stage 82)" at the banner lines.

The orchestrator will not invoke Codex on this unit until the user has chosen a direction.

## F1 — tautological_check

**Target:** `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage099_reduced_finish_line_sympy_audit.py:38-55` and `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage099_reduced_finish_line_mathematica_audit.wl:40-66`.

**Issue:**
All four `R_n - (N_Q - 1) == 0` assertions in both engines are guaranteed by the script's own construction. In sympy, `K_n = N_Q * K_n_target` makes `R_n = K_n/K_n_target - 1 = N_Q - 1` identically; for `Gamma_5 = 9 K_2^{5/2}/K_0^{3/2}`, the same multiplicative `N_Q` cancellation gives `R_5 = N_Q^{5/2}/N_Q^{3/2} - 1 = N_Q - 1` identically. In Mathematica, `kbar = nQ*k0Target*yhatCons` followed by series extraction reproduces the same per-construction reductions (`Coefficient[nQ*X*yhatCons, omega, 2n] = nQ * Coefficient[X*yhatCons, omega, 2n] = nQ * k_n_target`). No assertion can fail.

**Required change:**
Replace the four `R_n - (N_Q - 1)` checks with a non-tautological identity actually used in the appendix's factorization. Concretely, in **both** engines, after line 41 (sympy) / line 46 (Mathematica) — but **without** the upstream `K0 = N_Q * K0_target` substitution — keep `K0` as a free positive symbol and assert the canonical odd-coefficient/static-target factorization:

For sympy, edit lines 38-55 to:

```python
# Conservative quadrupole structural relations (forced by Yhat_Q^cons):
K0_sym = sp.symbols("K0_sym", positive=True, real=True)
K2_struct = K0_sym / (4 * Omega_Q**2)
K4_struct = K0_sym / (4 * Omega_Q**4)
Gamma5_struct = 9 * K0_sym / (32 * Omega_Q**5)

# Branch identity input (card stagefield Inputs at stage_099.tex:9):
expect_zero("branch identity K0 K4 = 4 K2^2", K0_sym * K4_struct - 4 * K2_struct**2)

# Equivalent Gamma_5 forms: derived from K_n series, and from canonical odd coeff
expect_zero("Gamma_5 sqrt form matches canonical odd-coeff form",
            9 * K2_struct**sp.Rational(5, 2) / K0_sym**sp.Rational(3, 2) - Gamma5_struct)

# Appendix factorization eq:app-part04-factorized-defect-again (chi_Q = 1 branch):
# Gammabar_5 / (2 G / (5 c^5)) = N_Q  (with K0 = N_Q * K0_target)
expect_zero("Gamma_5 normalization equals N_Q on chi_Q=1 branch",
            Gamma5_struct.subs(K0_sym, N_Q * K0_target) / (2 * G / (5 * c**5)) - N_Q)

# Print the structural moments for the transcript:
print("K2_struct =", sp.factor(K2_struct))
print("K4_struct =", sp.factor(K4_struct))
print("Gamma5_struct =", sp.factor(Gamma5_struct))
print("\nSTAGE 099 AUDIT PASSED")
```

The key change: `K0_sym` is a *free* symbol, so the assertion `K0 K4 - 4 K2^2 = 0` and the `Gamma_5` form-equivalence are forced by the *structural* relations between K_2, K_4, Gamma_5 and K_0 — not by the multiplicative `N_Q` recipe. The final `Gamma_5 / (2 G/(5 c^5)) = N_Q` assertion then exercises the appendix's factorization `Gammabar_5/(2G/(5c^5)) = chi_Q N_Q` at `stage_appendix_part04.tex:286` on the `chi_Q = 1` canonical branch.

For Mathematica, edit lines 40-66 to mirror the same structure: introduce a free symbol `k0Sym` (Element positive Real), define `k2Struct = k0Sym/(4*omegaQ^2)`, `k4Struct = k0Sym/(4*omegaQ^4)`, `gamma5Struct = 9*k0Sym/(32*omegaQ^5)`, then:

```mathematica
expectZero["branch identity K0 K4 = 4 K2^2", k0Sym*k4Struct - 4*k2Struct^2];
expectZero["Gamma_5 sqrt form matches canonical odd-coeff form",
           9*Sqrt[k2Struct^5/k0Sym^3] - gamma5Struct];
expectZero["Gamma_5 normalization equals N_Q on chi_Q=1 branch",
           (gamma5Struct /. k0Sym -> nQ*k0Target) / (2*G/(5*c^5)) - nQ];
```

Also keep (and do not break) the existing geometric-form check at sympy line 36 / Mathematica line 38 (`K0_target` two-form identity).

**Verification command:**
After Codex applies, the verifier will run `redteam exec-sympy 099` and `redteam exec-mathematica 099` and confirm:
- the new check `branch identity K0 K4 = 4 K2^2` appears and exits 0;
- the new check `Gamma_5 sqrt form matches canonical odd-coeff form` appears and exits 0;
- the new check `Gamma_5 normalization equals N_Q on chi_Q=1 branch` appears and exits 0;
- the old `R_n - (N_Q - 1)` checks are removed.

## F4 — insufficient_verification

**Target:** `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage099_reduced_finish_line_sympy_audit.py:31-32` and `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage099_reduced_finish_line_mathematica_audit.wl:33-34`.

**Issue:**
The conservative module expression `Yhat_Q^cons = 3/4 + (1/4)/(1 - omega^2/Omega_Q^2)` is the load-bearing object of the entire reduced finish line (notes §1.2; appendix `eq:app-part04-intro-3quarter-module`), but it is only printed, never asserted against. The sympy script binds `Yhat_Q_cons` at line 31 and never uses it again (the structural `K_n` are written directly from `N_Q * K_0_target` at lines 38-40). The Mathematica script uses it in `kbar` (line 40) but does not assert its static-slot value or its partial-fraction decomposition.

**Required change:**
In **both** engines, immediately after the `Yhat_Q^cons` definition (sympy line 31; Mathematica line 33), add two assertions:

For sympy at line 31 onward — add after line 32 (`print("Yhat_Q^cons(omega) =", Yhat_Q_cons)`):

```python
# Static slot: Yhat_Q^cons(omega=0) = 1
expect_zero("Yhat_Q^cons static slot equals 1",
            Yhat_Q_cons.subs(omega, 0) - 1)

# Pole structure: residue of Yhat_Q^cons at omega = Omega_Q is -Omega_Q/8
# (i.e., the (1/4)/(1 - omega^2/Omega_Q^2) part has residue -Omega_Q/8 at omega = Omega_Q,
#  confirming c_pole = 1/4 in the partial-fraction sense of the card's check (a) at
#  stage_099.tex:22).
expect_zero("Yhat_Q^cons pole residue at omega=Omega_Q is -Omega_Q/8",
            sp.residue(Yhat_Q_cons, omega, Omega_Q) - (-Omega_Q / 8))
```

For Mathematica at line 33 onward — add after line 34:

```mathematica
expectZero["Yhat_Q^cons static slot equals 1", (yhatCons /. omega -> 0) - 1];
expectZero["Yhat_Q^cons pole residue at omega=omegaQ is -omegaQ/8",
           Residue[yhatCons, {omega, omegaQ}] - (-omegaQ/8)];
```

**Verification command:**
After Codex applies, the verifier will run `redteam exec-sympy 099` and `redteam exec-mathematica 099` and confirm:
- the new check `Yhat_Q^cons static slot equals 1` appears and exits 0;
- the new check `Yhat_Q^cons pole residue at omega=Omega_Q is -Omega_Q/8` (Mathematica: `omegaQ`) appears and exits 0.

If perturbing `yhatCons` to `1/2 + 1/(2*(1 - omega^2/Omega_Q^2))` causes the residue assertion to fail (residue becomes `-Omega_Q/4`, not `-Omega_Q/8`), the check is doing real work.

---
## Applied: F2 (orchestrator-direct, post-user-resolution per batch-IV1-paper-alignment Cluster A direction (a))

- files_changed: scripts/moving_throat_pde_stage099_reduced_finish_line_sympy_audit.py
- summary: SymPy docstring annotates the three "Checks" items as upstream carry-ins from Part III + IV.1 (same anchors as 097: orthogonality at 094; static-limit at 091/092/094/096; minimal-module at 088/089/090). The local script also exercises Yhat_Q^cons static-slot and pole-residue as substantive anchors (F4 implementation).
- deviation: none

## Applied: F3 (orchestrator-direct, post-user-resolution per batch-IV1-paper-alignment Cluster C direction (a))

- files_changed: scripts/moving_throat_pde_stage099_reduced_finish_line_sympy_audit.py, mathematica/moving_throat_pde_stage099_reduced_finish_line_mathematica_audit.wl
- summary: Script-side banners and AUDIT PASSED prints relabeled to STAGE 099. Paper card title "Stage~116" deferred to PAPER_CLEANUP_TRACKER per Cluster C (a).
- deviation: paper-side edits explicitly deferred.

## Applied: F1
- files_changed: scripts/moving_throat_pde_stage099_reduced_finish_line_sympy_audit.py, mathematica/moving_throat_pde_stage099_reduced_finish_line_mathematica_audit.wl
- summary: Replaced 4 tautological R_n - (N_Q - 1) checks with 3 non-tautological structural identities using a FREE positive symbol K0_sym: branch identity K0 K4 = 4 K2^2, Gamma_5 sqrt vs canonical odd-coeff form equivalence, and Gamma_5/(2G/5c^5) = N_Q on chi_Q = 1 branch.
- deviation: none

## Applied: F4
- files_changed: scripts/moving_throat_pde_stage099_reduced_finish_line_sympy_audit.py, mathematica/moving_throat_pde_stage099_reduced_finish_line_mathematica_audit.wl
- summary: Added Yhat_Q^cons static-slot (=1) and pole-residue (-Omega_Q/8 at omega = Omega_Q) asserts in both engines.
- deviation: none

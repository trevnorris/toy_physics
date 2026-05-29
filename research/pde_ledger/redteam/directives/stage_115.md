---
unit_id: 115
batch: IV.3
created_at: 2026-05-27T00:00:00-06:00
findings_count: 1
stop_cold: null
applied: false
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 115

Apply each non-paper_misalignment finding below in order. After applying, append an `## Applied: F<n>` block under that finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

For `paper_misalignment` findings, do nothing — the orchestrator is holding for user resolution. Do not edit paper.tex, notes/, or scripts to "fix" a paper_misalignment unless the user has explicitly chosen a direction in a follow-up directive.

If a non-paper_misalignment finding's required change is ambiguous or unsafe to apply mechanically, append `## Blocked: F<n>` with a question instead — skip that finding, continue with the rest.

Do NOT introduce new features, refactors, or stylistic changes. Edit exactly the file:line ranges named.

After editing, RUN the affected scripts (`python3 <path>` for SymPy, `math -script <path>` for Mathematica) and iterate until they exit 0 with all in-file checks passing. Getting the scripts to run cleanly is your job; the orchestrator independently re-runs afterward.

Do NOT touch paper.tex, notes/, or any prose documents. The red-team only modifies scripts (except when a follow-up directive explicitly authorizes a paper-side edit after user resolution).

## F1 — mathematica_transliteration

**Target:** `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage115_core_balance_compensation_mathematica_audit.wl:46-47` (insertion point between the existing `expectZero["sigma_c on balance surface", …]` at line 46 and the `kappa0Can = …` definition at line 48)

**Issue:**
The `.wl` script is a line-by-line rewrite of `scripts/moving_throat_pde_stage115_core_balance_compensation_sympy_audit.py:24-70`. Every algebraic step, the branch pick (`First[gQSolutions]` vs `gq_solutions[0]`), and the substitution sequence are identical. There is no independent derivation of the balance surface or of `sigma_c = sigma_*` in the Mathematica script — it merely echoes the SymPy algebra in Mathematica notation. The second-engine policy requires an independent algebraic route. The Part-IV appendix (`paper/appendices/stage_appendix_part04.tex` lines 530-543) supplies one: the parent-overlap reparametrization `mathfrak r = lambda/sqrt(K_s K_q)`, `mathfrak g = g_q sqrt(K_s)/(g_s sqrt(K_q))` with the equivalent compensation family `1 + mathfrak r^2 = 4 (mathfrak g - mathfrak r)^2` and roots `mathfrak g_+- = mathfrak r +- (1/2) sqrt(1 + mathfrak r^2)`.

**Required change:**
Insert (alongside the existing derivation, not replacing it) an independent verification block in the Mathematica script that derives `sigma_c = sigma_*` on the compensated branch via the parent reparametrization rather than via `Solve[balanceEq == 0, gQ]`.

Edit the file `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage115_core_balance_compensation_mathematica_audit.wl` by inserting the following block as new lines immediately after the existing line 46 (which currently reads `expectZero["sigma_c on balance surface", (sigmaC /. gQ -> gQBranch) - sigmaStar];`) and before the existing line 48 (`kappa0Can = FullSimplify[(1 + rC)/3, Assumptions -> $Assumptions];`):

```mathematica

(* Independent derivation: parent-overlap reparametrization (Part-IV appendix \
   eqs eq:app-part04-r-g-parent-ratios and eq:app-part04-parent-compensation-family). *)
(* Define the dimensionless parent ratios. *)
frakR = lam/Sqrt[kS*kQ];
frakG = gQ*Sqrt[kS]/(gS*Sqrt[kQ]);
(* The parent compensation family condition: 1 + frakR^2 - 4 (frakG - frakR)^2 = 0. *)
parentFamilyResidual = 1 + frakR^2 - 4*(frakG - frakR)^2;
(* Equivalence to the original balance equation: parentFamilyResidual vanishes \
   iff the original balanceEq vanishes (up to a positive multiplicative factor). *)
expectZero[
  "independent: parent family ≡ balance equation",
  FullSimplify[
    parentFamilyResidual - balanceEq*(kS*kQ)/(gS^2),
    Assumptions -> $Assumptions
  ]
];
(* Solve the parent family for the dimensionless mixed coupling frakG. *)
frakGRoots = frakG /. Solve[parentFamilyResidual == 0, frakG];
If[Length[frakGRoots] =!= 2,
  fail["expected two parent-family roots", frakGRoots]];
Print["Parent-family roots for frakG = ", fmt[frakGRoots]];
(* Pick the lower branch frakR - (1/2) Sqrt[1 + frakR^2]. *)
frakGMinus = FullSimplify[frakR - Sqrt[1 + frakR^2]/2,
  Assumptions -> $Assumptions];
expectZero[
  "independent: frakG_- root matches Solve output",
  FullSimplify[(frakGRoots[[1]]) - frakGMinus,
    Assumptions -> $Assumptions] *
  FullSimplify[(frakGRoots[[2]]) - frakGMinus,
    Assumptions -> $Assumptions]
];
(* Translate back to g_q via g_q = frakG * gS * Sqrt[kQ]/Sqrt[kS]. *)
gQFromFrakMinus = FullSimplify[
  frakGMinus*gS*Sqrt[kQ]/Sqrt[kS],
  Assumptions -> $Assumptions
];
(* Verify that sigma_c = sigma_* under this independent substitution. *)
expectZero[
  "independent: sigma_c = sigma_* via parent reparametrization",
  FullSimplify[
    (sigmaC /. gQ -> gQFromFrakMinus) - sigmaStar,
    Assumptions -> $Assumptions
  ]
];

```

Notes for Codex on the inserted block:

1. The `expectZero["independent: parent family ≡ balance equation", …]` check shows that `parentFamilyResidual` and the original `balanceEq` differ only by the positive factor `gS^2/(kS*kQ)` (so their zero sets coincide). The exact factor was derived as follows: `balanceEq = (gS^2 kQ kS - 4 gQ^2 kS^2 + 8 gQ gS kS lam - 3 gS^2 lam^2) / (kS (kQ kS + lam^2))` (per the saved Mathematica output line 13). Multiplying by `(kS*kQ)/gS^2` and clearing the denominator yields `(kQ*kS + lam^2 - 4 (gQ*kS/gS)^2 (kQ/kS) + 8 (gQ*kS/gS)(lam/sqrt(kS*kQ))*sqrt(kQ/kS) − 3 (lam^2/(kS*kQ))*(kS*kQ))/(kS*kQ + lam^2)`, which after substituting `frakG = gQ*sqrt(kS)/(gS*sqrt(kQ))` and `frakR = lam/sqrt(kS*kQ)` and simplifying matches `1 + frakR^2 - 4 (frakG - frakR)^2`. Codex should rely on `FullSimplify` to verify this; if the residual after `FullSimplify` is not literally `0`, Codex must NOT silently change the multiplicative factor — instead append a `## Blocked: F1` block asking the auditor to recompute the equivalence factor.
2. The `expectZero["independent: frakG_- root matches Solve output", …]` check is the product of the two differences `(roots[[1]] - frakGMinus)` and `(roots[[2]] - frakGMinus)`. This product is zero iff `frakGMinus` matches one of the two roots (Solve order is implementation-dependent, so we cannot index a single root reliably). This is non-tautological: it would fail if `Sqrt[1 + frakR^2]/2` were the wrong half-discriminant.
3. The final `expectZero["independent: sigma_c = sigma_* via parent reparametrization", …]` is the load-bearing check: it substitutes `gQ` from the *parent-family* root (not from `gQBranch = First[gQSolutions]`) into `sigmaC` and confirms the result equals `sigmaStar`. The two derivations are now algebraically independent paths to the same conclusion.
4. The block uses only symbols already in scope (`kS, kQ, lam, gS, gQ, sigmaC, sigmaStar, balanceEq, fmt, expectZero, fail`) and the existing `$Assumptions` declaration on line 29-31. Do NOT add new global symbols beyond `frakR`, `frakG`, `parentFamilyResidual`, `frakGRoots`, `frakGMinus`, `gQFromFrakMinus`.
5. Leave the rest of the file (lines 48-74) unchanged. In particular, do NOT modify the `kappa0Can`, `gamma0Can`, `deltaCore`, `targetDelta`, `lambdaOut`, `lambdaEff`, `yEff`, or `yTarget` blocks. Do NOT modify line 41 (`If[Length[gQSolutions] =!= 2, fail[...]];`). Do NOT remove the existing `expectZero` at line 46.

**Verification command:**
After Codex applies, the verifier will run `redteam exec-mathematica 115` and confirm:

- the new lines `independent: parent family ≡ balance equation`, `independent: frakG_- root matches Solve output`, and `independent: sigma_c = sigma_* via parent reparametrization` all appear in the transcript;
- each of those three new `expectZero` calls reports residual `0` (PASS);
- the existing three `expectZero` calls (`sigma_c on balance surface`, `exact collapse to compensated branch`, `normalized outgoing fingerprint preserved`) still report residual `0`;
- the script exits 0;
- `Stage 115 Mathematica audit passed.` still prints at the end.

---
unit_id: 178
batch: V.2
created_at: 2026-05-30T00:00:00-06:00
findings_count: 2
stop_cold: null
applied: true
applied_at: 2026-05-30T01:15:34-06:00
findings_applied: 2
findings_blocked: 0
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 178

Apply each finding below in order. After applying, append an `## Applied: F<n>` block under that finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

If a finding's required change is ambiguous or unsafe to apply mechanically, append `## Blocked: F<n>` with a question instead — skip that finding, continue with the rest.

Do NOT introduce new features, refactors, or stylistic changes beyond the named edits. Do NOT touch paper.tex, notes/, or any prose documents.

After editing, RUN the affected scripts (`python3 <path>` for SymPy, `math -script <path>` for Mathematica) and iterate until they exit 0 with all in-file checks passing. The orchestrator independently re-runs afterward.

## F1 — mathematica_transliteration

**Target:** `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage178_outgoing_port_coloading_mathematica_audit.wl` (add inside Section 4, after the existing `expectZero["nu_r - [kappa1 + sigma_r]", ...]` at line 101)

**Issue:** The `.wl` is a line-by-line port of the `.py`: it builds `pExpected`/`dExpected` exactly as SymPy does and then assembles `nuExpected` from the same slippage formula, so the Mathematica run does not independently confirm the central identity — a shared construction error would be reproduced, not caught. The algebra is correct; the missing piece is an independent path from the port data to `nu_r` that does not retrace the `2*(pExpected - dExpected)` factoring.

**Required change:**
Add an independent cross-derivation of `nu_r` in Section 4 of the `.wl`. Use the component-level expansions `pA` and `dA` already defined at lines 62–63 (these are \(P_{A,r}\) and \(\Delta_{A,r}\) built from the five component drifts `oU, gW, rr, gU, oW`), and obtain `nu_r` via the log-slope of \(N_{A,0}^{(r)}=P_{A,r}^2/\Delta_{A,r}^2\) rather than via `pExpected`/`dExpected`. Insert immediately after line 101:

```
(* Independent route: nu_r as the eps*lam log-slope of N = P^2/Delta^2,
   built straight from the component drifts, not from pExpected/dExpected. *)
nuFromData = FullSimplify[
  Coefficient[Normal[Series[Log[pA^2/dA^2], {eps, 0, 1}]], eps*lam],
  Assumptions -> $Assumptions
];
expectZero["nu via log-data vs slippage", nuFromData - nuExpected];
Print["nuFromData = ", fmt[nuFromData]];
```

Leave the existing §1–§4 checks untouched. Do not delete or rewrite `pExpected`/`dExpected`; just add the independent confirmation above.

**Self-test (auditor-verified):**
- `pA` and `dA` (lines 62–63) both depend on `eps` (through `eps*lam*(...)`), so `Series[..., {eps,0,1}]` is non-trivial and the `eps*lam` coefficient is the genuine first-order log-slope \(\nu_r=2(\mathfrak p_r-\mathfrak d_r)\). No identically-zero-derivative trap.
- `Coefficient[..., eps*lam]` is well-defined because `eps` and `lam` always occur as the product `eps*lam` in `pA`/`dA`.
- `nuFromData - nuExpected` reduces to `0` for the same reason the existing `nuDirect - nuExpected` (line 101) does — both equal \(2(\mathfrak p_r-\mathfrak d_r)\) — so the new `expectZero` passes, but now via a path that never references `pExpected`/`dExpected`.

**Verification command:**
After Codex applies, the verifier runs `redteam exec-mathematica 178` and confirms the new `expectZero["nu via log-data vs slippage", ...]` appears, prints residual `0`, and the script exits 0.

**Note for orchestrator:** this is a policy-driven finding. If an established MATHEMATICA_MIRROR_POLICY for this paper explicitly sanctions a mirrored second engine, the user/orchestrator may waive F1; the auditor cannot read that tracker.

## Applied: F1

- files_changed:
  - `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage178_outgoing_port_coloading_mathematica_audit.wl`
- summary: Added the independent `nu_r` log-slope check from `pA` and `dA` immediately after the existing slippage residual.
- deviation: none

## F2 — stale_output (banner label defect)

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage178_outgoing_port_coloading_sympy_audit.py:34`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage178_outgoing_port_coloading_mathematica_audit.wl:26`

**Issue:** Both scripts print `banner("STAGE 161 — OUTGOING-PORT CO-LOADING THEOREM")`, which propagates a wrong stage number ("STAGE 161") into both saved transcripts. This unit is stage 178; the docstring header and the Mathematica final line already say 178.

**Required change:**
- SymPy line 34: change `banner("STAGE 161 — OUTGOING-PORT CO-LOADING THEOREM")` to `banner("STAGE 178 — OUTGOING-PORT CO-LOADING THEOREM")`.
- Mathematica line 26: change `banner["STAGE 161 — OUTGOING-PORT CO-LOADING THEOREM"];` to `banner["STAGE 178 — OUTGOING-PORT CO-LOADING THEOREM"];`.

No other change. No assertion depends on this string; running the scripts after the edit must still exit 0.

**Verification command:**
After Codex applies, the verifier reruns both engines and confirms the first banner of each refreshed transcript reads `STAGE 178 — OUTGOING-PORT CO-LOADING THEOREM`, both exit 0.

## Applied: F2

- files_changed:
  - `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage178_outgoing_port_coloading_sympy_audit.py`
  - `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage178_outgoing_port_coloading_mathematica_audit.wl`
- summary: Corrected the first audit banner in both engines from stage 161 to stage 178.
- deviation: none

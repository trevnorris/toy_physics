---
unit_id: 239
batch: VII.2
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-06-02T00:00:00-06:00
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
  notes_stage_files: ["/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage239_rigid_mouth_physical_normal_form_exact_physical_to_microscopic_correction_compiler_and_cartesian_orbit_lock_theorem_sympy_audit.md"]
  paper_appendix: present
---

# Audit unit 239 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_239.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage239_rigid_mouth_physical_normal_form_exact_physical_to_microscopic_correction_compiler_and_cartesian_orbit_lock_theorem_sympy_audit.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part07.tex` (row at line 90; narrative at lines 1101-1172)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage239_rigid_mouth_physical_normal_form_exact_physical_to_microscopic_correction_compiler_and_cartesian_orbit_lock_theorem_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage239_rigid_mouth_physical_normal_form_exact_physical_to_microscopic_correction_compiler_and_cartesian_orbit_lock_theorem_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage239_rigid_mouth_physical_normal_form_exact_physical_to_microscopic_correction_compiler_and_cartesian_orbit_lock_theorem_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage239_rigid_mouth_physical_normal_form_exact_physical_to_microscopic_correction_compiler_and_cartesian_orbit_lock_theorem_mathematica_audit.txt`

## What the paper claims

The `\stagefield{Output}` is: "Cartesian orbit-lock theorem: rigid-mouth orbit lock is exactly `U=V=0`, equivalently `T^2=T_ref^2` and `eps_eta=eps_eta_ref`." The stage uses the physical logarithmic variables `U=ln(T^2/T_ref^2)`, `V=ln(eps_eta/eps_eta_ref)` as the rigid-mouth chart and (per the derivation ledger and notes §1-§8) proves: (i) the packet is diagonal in `(U,V)` with `q_nt=U`, `q_eta=V` and compiler matrix `M_phys=I_2`; (ii) the exact target-ratio formula `R_target/R_target_ref = ((1-eps_ref e^V)/(1-eps_ref)) e^{-U}`; (iii) complementary projectors `P_T`, `P_eta` and the commuting transfer/dressing legs; (iv) the exact physical-to-microscopic dependent compiler `(Delta_T, Delta_Keta, Delta_mu)=(0, -V, U-V)` with left inverse `U=Delta_mu-Delta_Keta`, `V=-Delta_Keta`; (v) the two axis images `(U,0)->(0,0,U)` and `(0,V)->-V(0,1,1)`; (vi) the static-only / post-static / full orbit correction split; (vii) support-blindness in `zeta` and `M_mix`; and (viii) the orbit-lock equivalence in finite and first-order form. Constants and matrices are stated identically in the stage card, notes §9, and appendix eqs. `app-part07-uv-def` ... `app-part07-uv-first-order-dependent`. All values are exact-symbolic; there are no free numeric ansatz constants.

## What the script claims to verify

The SymPy docstring (lines 4-18) enumerates eight deliverables that line up one-to-one with the paper's. The 51 assertions test: the log-coordinate definitions, the diagonal `M_phys=I_2`, the target-ratio formula by two independent substitutions, projector idempotence/orthogonality/completeness, the packet decomposition, the transfer/dressing legs and their factorization, the Stage-238 branch identity and the `q_nt=U` / `q_eta=V` identifications carried through an explicit Stage-238 branch formula, the dependent compiler `(0,-V,U-V)`, a genuine `solve`-based left inverse plus `L.C=I`, the two axis images, the three correction compilers and their additive split, support-blindness via abstract functions of `(zeta,M_mix)` with derivative-zero substitution, a `solve`-based orbit-lock `{U:0,V:0}`, and the first-order `D[...,h]/.h->0` compiler. The Mathematica script asserts the same set with the same labels.

## Paper <-> script cross-check

| Paper deliverable | Script check | Status |
|---|---|---|
| Diagonal chart `q_nt=U`, `q_eta=V`, `M_phys=I_2` | `diagonal rigid-mouth packet compiler`; `Stage 238 ... q_nt=U` / `q_eta=V` | match |
| Target ratio `((1-eps_ref e^V)/(1-eps_ref)) e^{-U}` | `target-ratio formula from U,V` + `target-ratio reconstruction in physical chart` | match |
| Projectors + commuting legs | `P_T/P_eta` block, `pure transfer/dressing leg`, `exact commutativity` | match |
| Compiler `(0,-V,U-V)` | `physical-to-microscopic dependent compiler` | match |
| Left inverse `U=Dmu-DK`, `V=-DK` | `left inverse ...`, `left inverse of physical compiler`, `recovery of U,V` | match |
| Axis images `(U,0)->(0,0,U)`, `(0,V)->-V(0,1,1)` | `pure transfer/dressing microscopic image` | match |
| Correction split (static/post-static/full) | section 5 block | match |
| Support-blindness in `zeta`, `M_mix` | section 6 block | match |
| Orbit-lock `U=V=0` (finite + first order) | section 7 block | match |

`paper_alignment: aligned` — every paper-side deliverable has a faithful, non-tautological script-side check, the constants/matrices match exactly, and there are no extra script claims unanchored to the paper.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored? |
|---|---|---|---|---|---|
| A1 | sympy | 71-72 | `exp(U_def)-T2/T2_ref==0` | chart def | yes |
| A2 | sympy | 76 | `M_phys*[U,V]-[U,V]==0` | diagonal chart (i) | partial (I_2 only) |
| A3 | sympy | 84-87 | ratio_identity - ratio_UV ==0 | target ratio (ii) | yes |
| A4 | sympy | 90-98 | ratio_UV(U_def,V_def)-ratio_identity==0 | target ratio (ii) | yes |
| A5 | sympy | 106-115 | projector algebra + decomposition | projectors (iii) | yes |
| A6 | sympy | 118-128 | legs + factorization | commuting legs (iii) | yes |
| A7 | sympy | 138-171 | Stage-238 branch identity, q_nt=U, q_eta=V | chart identification (i) | yes |
| A8 | sympy | 174-215 | compiler, propagation, solve-based left inverse | compiler+inverse (iv) | yes |
| A9 | sympy | 220-241 | axis images + correction split | (v)(vi) | yes |
| A10 | sympy | 259-305 | support-blindness chain rule | (vii) | yes |
| A11 | sympy | 310-352 | orbit-lock solve + first order | (viii) | yes |
| B1-B... | mathematica | 87-370 | same set, same labels | same | yes (but transliterated — see F1) |

A2 is marked partial because it only confirms `M_phys=I_2` acts as identity; that is exactly the paper's claim (the chart is diagonal with `M_phys=I_2`), so it is anchored, just narrow by nature. No assertion is tautological in the `x=expr; assert x==expr` sense; the `solve`/`Solve` calls (left inverse, orbit lock, equilibria) are genuine inversions, not assumptions.

## Findings

### F1 — mathematica_transliteration

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage239_rigid_mouth_physical_normal_form_exact_physical_to_microscopic_correction_compiler_and_cartesian_orbit_lock_theorem_mathematica_audit.wl:84-370`

**What's wrong:**
The `.wl` is a line-by-line port of the `.py`, not an independent second-engine derivation. Concretely:

- Same variable choreography, just camelCased. SymPy `U_def = sp.log(T2/T2_ref)` / `V_def = sp.log(eps_eta/eps_eta_ref)` (py:68-69) become `Udef = Log[T2/T2ref]` / `Vdef = Log[epsEta/epsEtaRef]` (wl:84-85). SymPy `M_phys = sp.eye(2)` (py:74) -> `Mphys = IdentityMatrix[2]` (wl:90). SymPy `S_rm_dep = sp.Matrix([[0,0],[0,-1],[1,-1]])`, `C_phys_dep = S_rm_dep*M_phys` (py:174-175) -> `SrmDep = {{0,0},{0,-1},{1,-1}}`, `CphysDep = SrmDep . Mphys` (wl:185-186), identical hardcoded decomposition.
- Same intermediate-step ordering and identical check label strings, verbatim, e.g. "exact commutativity / factorization of transfer and dressing legs" (py:127 vs wl:135), "Stage 236/221 dependent compiler propagated into the physical chart" (py:190 vs wl:201), "recovery of U,V from dependent correction" (py:215 vs wl:229). The 51 assertions appear in the same order with the same names.
- Same `solve` choreography: SymPy `sp.solve([Eq(yK_var,-q_eta_var), Eq(yMu_var,q_nt_var-q_eta_var)],[q_nt_var,q_eta_var])` with `expected_solution=[{q_nt_var: yMu_var-yK_var, q_eta_var:-yK_var}]` (py:197-204) is transcribed verbatim to `Solve[{yK==-qEtaVar, yMu==qNtVar-qEtaVar},{qNtVar,qEtaVar},Reals]` with `expectedInverse={{qNtVar->yMu-yK, qEtaVar->-yK}}` (wl:208-221).
- Same first-order trick: SymPy `T2_pert = T2*sp.exp(h*dlnT2)`, `diff(log(T2_pert/T2),h).subs(h,0)` (py:335,338) -> `T2pert = T2 Exp[h dlnT2]`, `D[Log[T2pert/T2],h] /. h->0` (wl:354,357).

Per the dual-engine policy, both CAS must derive the result independently from the physical premises, not echo each other's algebra. As written, the `.wl` provides no independent confirmation — it would reproduce any algebra error the `.py` made, because it makes the identical choices (same `S_rm_dep`, same projectors, same propagation subs).

**Why this matters:**
The second engine exists to catch errors the first engine's specific decomposition/assumptions could hide. A transliteration cannot do that; the engine cross-check is illusory.

**Required change:**
Re-author the `.wl` from the paper's premises via a genuinely independent route (different decomposition, native Mathematica primitives). See the directive for an independent-route prescription and the claim manifest M1-M9. Keep the `_mathematica_audit.wl` filename suffix (it is already correct).

**Verification:**
The re-authored `.wl` must not share the `S_rm_dep . M_phys` construction, the `expected_solution` literal-comparison choreography, or the verbatim check-label strings of the `.py`; it must still establish the same nine physical claims and exit 0.

### F2 — stale stage label (script-side, cosmetic)

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage239_..._mathematica_audit.wl:59` (banner `"STAGE 222 — RIGID-MOUTH PHYSICAL NORMAL FORM"`)
- same file: `T2Stage221`, `RtargetStage221`, `qNtRigidStage221`, `UStage221`, `VStage221` (wl:141-143,153-164,168-203), `Stage221` variable suffixes

**What's wrong:**
The Mathematica script's banner says "STAGE 222" and its carried Stage-238 branch variables are suffixed `Stage221`, while the canonical unit is 239 (paper card `stage_239.tex`, script filename `stage239`, appendix row 90). This is the known project-wide incomplete-renumber drift from the EM-extension realignment, surfacing on the script side. The `Print` banner propagates "STAGE 222" into the saved output header (mathematica output line 11). Canonical numbering is ground truth; the labels are simply stale. (The SymPy script's section comments reference "Stage 236/221" in the same drift pattern, e.g. py:11, py:131,173,193 — same root cause, no separate finding.)

NOTE: This is a *script-side* label, not a paper/notes disagreement, so it is safe for Codex to correct. This is NOT a basis for any batch-wide renumber sweep; fix only the labels in this one `.wl` (and, if Codex re-authors per F1, simply use the correct "Stage 239" labels in the new file).

**Why this matters:**
The output header advertises the wrong stage number, which is confusing for downstream trust audits, but the math is unaffected.

**Required change:**
In the re-authored (or, if F1 not yet applied, the existing) `.wl`, change the banner to "STAGE 239 — RIGID-MOUTH PHYSICAL NORMAL FORM" and rename the `Stage221` variable suffixes to `Stage238` (the carried branch is from Stage 238 per the paper card `\stagefield{Inputs}` and notes §0). Do not touch paper/ or notes/.

**Verification:**
`grep -i "STAGE 222\|Stage221" <wl>` returns nothing; output header reads "STAGE 239".

## Independent-derivation check (Mathematica)

The `.wl` is a transliteration, not an independent derivation. The three sections quoted in F1 (chart definitions wl:84-93 vs py:68-76; compiler decomposition wl:185-196 vs py:174-181; solve-based left inverse wl:208-227 vs py:197-211) are mechanical syntax translations preserving every variable name, every intermediate object, and every check label. See F1.

## Engine cross-check

Both engines exit 0 with all 51/50 checks PASS (sympy output 50 `[ok]` lines + final pass; mathematica output identical set). Their final residuals all reduce to `0` / zero matrices, so they agree numerically and symbolically. The agreement is expected but, per F1, not independent — it is the agreement of a script with its own translation. Minor coverage note: the SymPy has two distinct target-ratio checks (`target-ratio formula from U,V` py:84-87 and `target-ratio reconstruction in physical chart` py:90-98); the Mathematica carries only the second (wl:103-109). Not a defect on its own (the second subsumes the first) but the re-authored `.wl` should still cover the target-ratio identity.

## Verdict justification

`verdict: findings`. The mathematics holds up: I attacked the orbit-lock `solve` (genuine inversion of `-V=0`, `U-V=0`, not assumed), the support-blindness checks (real chain-rule propagation through abstract `T2_sb(zeta,M_mix)`, plus a legitimate direct check that the Stage-238 branch formula contains no support variables), the left-inverse `L.C=I` (non-trivial), and the first-order `D[...,h]/.h->0` derivation (correct) — all survive. Paper alignment is exact: constants `(0,-V,U-V)`, matrices `M_phys=I_2`, `C_phys_dep`, `L_phys_dep`, and the target ratio match the card, notes §9, and appendix verbatim. The two findings are: (F1) the Mathematica script is a transliteration of the SymPy and therefore does not provide the independent second-engine confirmation a checkpoint requires; and (F2) cosmetic stale "STAGE 222"/"Stage221" labels in the `.wl` from the known renumber drift. No `paper_misalignment`, no stop-cold: the corrections are script-only and do not propagate downstream (the verified identities are unchanged; only the route and labels change).

## Self-test notes

I checked: (1) variable-independence in the support-blind `D[...]` checks — `T2Stage221` genuinely has no `zeta`/`M_mix`, so those derivatives are correctly identically zero (the claim is exactly that support-blindness); the abstract-function checks `D[Log[T2sb/T2ref],zeta]/.rules` are real chain-rule, non-vacuous. (2) parity/integrals — none in this stage (pure algebra/linear maps). (3) trivial-case: substituting `U=V=0` into `y_dep` and `ratio_from_UV` gives `(0,0,0)` and `1` respectively, as asserted. (4) the F1/F2 fixes change route and labels only, introducing no new constant, so the paper round-trip is preserved. The directive for F1 is a re-authoring requirement (claim manifest only), not a pre-written script.

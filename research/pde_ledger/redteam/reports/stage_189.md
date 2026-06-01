---
unit_id: 189
batch: V.3
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-06-01T00:00:00-06:00
verdict: findings
stop_cold: null
findings_count: 3
paper_alignment: aligned
scripts_checked:
  sympy: present
  mathematica: missing
  engines_agree: n/a
  outputs_fresh: true
docs_read:
  paper_stage_tex: present
  notes_stage_files: [moving_throat_pde_stage189_transfer_shape_prefactor_compiler.md]
  paper_appendix: present
---

# Audit unit 189 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_189.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage189_transfer_shape_prefactor_compiler.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part05.tex` (row at line 109; supporting context lines 47, 644-699, 1466)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage189_transfer_shape_prefactor_compiler_sympy_audit.py`
- mathematica: (missing)
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage189_transfer_shape_prefactor_compiler_sympy_audit.txt`
- mathematica output: (missing)

## What the paper claims

The stage card `\stagefield{Output}` states: "Equates \(\Xi_1\) with the logarithmic transfer-shape and outgoing-prefactor slope, including corrected nontracking form." (verbatim, stage_189.tex:15; the part-05 appendix row at line 109 repeats this verbatim.) The notes expand this into five concrete deliverables: (1) the exact first-order transfer-shape identities `δln T² = Ξ₁`, `δln(1-ε_η) = R₁+Ξ₁`, `δln R_target = R₁`, plus the selected-branch compatibility identity `R_target·T² = Λ₀(1-ε_η)` with `Λ₀ := 27π²Gc_s⁵/(20a⁵c⁵)`; (2) the exact isotropic grouped compiler `(D₀,D₂,D₄,N₀,N₂,N₄) ↦ (u₂,u₄,P₀,P₂,P₄)` with `u₂=-D₂/D₀`, `u₄=(D₂²-D₀D₄)/D₀²`, `P₀=N₀/D₀`, `P₂=(D₀N₂-2D₂N₀)/D₀²`, `P₄=(D₀²N₄-2D₀(D₂N₂+D₄N₀)+3D₂²N₀)/D₀³`; (3) the outgoing-branch coefficients `K₀=P₀`, `K₂=P₂+A·P₀`, `K₄=P₄+A·P₂+B·P₀`, `Γ₅=G₅P₀=a⁵P₀/(27c_s⁵)`, with fingerprint constants `A=a²/(9c_s²)`, `B=4a⁴/(81c_s⁴)`, `G₅=a⁵/(27c_s⁵)`; (4) the normalization equivalence `m̂₀²Γ₅=2G/(5c⁵) ⟺ m̂₀²P₀=54Gc_s⁵/(5a⁵c⁵)`; (5) the carry-forward bridge `Ξ₁ = δln T_A²/(ελ_A) = P₁/P₀` with `P₁=(N₁D₀-N₀D₁)/D₀²`. The notes also assert the rank-2 compatibility relation `δln R_target + δln T² - δln(1-ε_η) = 0` as a derived (not definitional) condition (notes:217-224).

## What the script claims to verify

The script (banner "STAGE 172", ledger "STAGE 172") is organized in six sections mirroring the notes. Section I builds the observable→transfer compiler matrix `C_obs→trf` and checks it reproduces the transfer packet, plus the defect-notation identities for `Ξ`, `R₁`. Section II checks the selected-branch transfer-shape identity `Λ₀(1-ε_η)/R_target = T_A²`. Section III independently series-expands `Y(ω)=D₀/D_cons(ω)` and `Pref(ω)=D₀N(ω)/D_cons(ω)²` and compares to the hand-written `u₂,u₄,P₀,P₂,P₄` compiler formulas, plus the static `P₀=(K_bl/D₀)T_eff²` bridge. Section IV series-expands the weak-axisymmetric `P_A(ε)` and checks `P₁` and `P₁/P₀=N₁/N₀-D₁/D₀`. Section V forms the product `Pref·Ŷ₂^(out)` and checks the `K₀,K₂,K₄,Γ₅` coefficients. Section VI checks the normalization equivalence at the claimed `P₀_target` value and the constant-prefactor branch (`P₂=0,P₄=0` conditions, `K₂=AP₀`, `K₄=BP₀`).

## Paper ↔ script cross-check

| Paper-side deliverable | Script-side check | Status |
|---|---|---|
| Transfer-shape identities `δln T²=Ξ₁`, `δln R_target=R₁`, `δln(1-ε_η)=R₁+Ξ₁` | I, lines 70-72 | mismatch (tautological — see F2) |
| Rank-2 compatibility `δln R_target+δln T²-δln(1-ε_η)=0` (derived) | I, line 63 | mismatch (tautological — baked into definition; F2) |
| Compiler matrix `C_obs→trf` reproduces transfer packet | I, line 62 | match (matrix written independently of component defs) |
| Selected-branch identity `R_target·T²=Λ₀(1-ε_η)`, `Λ₀=27π²Gc_s⁵/(20a⁵c⁵)` | II, line 89 | mismatch (tautological round-trip — F1; and Λ₀'s explicit value never instantiated) |
| Coherent D/N specialization `T²=Z_W(1+χ₀)²/(Ω_W²(1-ε)²)` | II, lines 91-95 | partial (printed, not asserted — print-only; acceptable as display) |
| Grouped compiler `(D...,N...)↦(u₂,u₄,P₀,P₂,P₄)` | III, lines 106-119 | match (independent series vs hand formulas) |
| Static prefactor bridge `P₀=(K_bl/D₀)T_eff²` | III, line 123 | match |
| Weak-axisym slope `P₁`, `P₁/P₀=N₁/N₀-D₁/D₀`, `Ξ₁=P₁/P₀` bridge | IV, lines 132-135 | match (Ξ₁=P₁/P₀ bridge implicit: both equal slope; see note below) |
| Outgoing coefficients `K₀,K₂,K₄,Γ₅` and fingerprint `A,B,G₅` | V, lines 145-158 | match (independent product expansion) |
| Normalization equivalence `m̂₀²Γ₅=2G/(5c⁵) ⟺ m̂₀²P₀=54Gc_s⁵/(5a⁵c⁵)` | VI, lines 166-173 | match (genuine value/target consistency) |
| Constant-prefactor branch `P₂=0,P₄=0`, `K₂=AP₀`, `K₄=BP₀` | VI, lines 177-180 | match |

`paper_alignment` set to `aligned`: every deliverable has a corresponding script-side section and no verified identity contradicts the paper. The findings are about checks that PASS FOR THE WRONG REASON (tautology), not about verifying a wrong identity, so this is not a `paper_misalignment`. Note: the central `Output` claim `Ξ₁=δln T_A²/(ελ_A)=P₁/P₀` is only verified on the `P₁/P₀` half (Section IV); the `Ξ₁=δln T_A²/(ελ_A)` half is asserted in the notes (line 294) but never independently exercised in the script — it is carried as an upstream theorem. This is acceptable carry-forward, not a finding, since the script does verify the `P₁/P₀` algebra that the bridge equates to.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 62 | `expect_zero(C·obs - trf)` | compiler matrix vs component defs | yes (matrix independent) |
| A2 | sympy | 63 | `expect_zero(trf[2]+trf[0]-trf[1])` | rank-2 compatibility | no (definitional — tautology, F2) |
| A3 | sympy | 70 | `expect_zero(Xi - (dn - Bstar*dr))` | `δln T²=Ξ₁` | no (Xi≡dlnT2≡dn-Bstar·dr; tautology, F2) |
| A4 | sympy | 71 | `expect_zero(Rcal - dlnRtarget)` | `δln R_target=R₁` | no (Rcal≡dlnRtarget; tautology, F2) |
| A5 | sympy | 72 | `expect_zero((Rcal+Xi) - dlnOneMinus)` | `δln(1-ε_η)=R₁+Ξ₁` | no (definitional; tautology, F2) |
| A6 | sympy | 89 | `expect_zero(T2_selected - T2_direct)` | selected-branch identity | no (R_target back-defined as Λ₀(1-εη)/T2_direct; tautology, F1) |
| A7 | sympy | 111 | `expect_zero(Y_series - (1+u2ω²+u4ω⁴))` | `u₂,u₄` compiler | yes (independent series) |
| A8 | sympy | 119 | `expect_zero(Pref_series - (P0+P2ω²+P4ω⁴))` | `P₀,P₂,P₄` compiler | yes (independent series; P4 nontrivial) |
| A9 | sympy | 123 | `expect_zero(P0 - (K_bl/D0)T_eff²)` | static prefactor bridge | yes |
| A10 | sympy | 132 | `expect_zero(P_A_series - (P0+ελP1))` | `P₁` slope | yes (independent ε-series) |
| A11 | sympy | 134 | `expect_zero(ln_ratio - ελ(P1/P0))` | log slope `Ξ₁` | yes |
| A12 | sympy | 135 | `expect_zero(P1/P0 - (N1/N0-D1/D0))` | `P₁/P₀` identity | yes |
| A13 | sympy | 155 | `expect_zero(out_trunc - (K0+K2ω²+K4ω⁴+iΓ5ω⁵))` | outgoing coefficients | yes (independent product) |
| A14 | sympy | 166 | `expect_zero((m̂₀²Γ5-2G/5c⁵).subs(P0,P0_target))` | normalization equivalence | yes (value/target consistency) |
| A15 | sympy | 170 | `expect_zero(Γ5 - a⁵P0/(27c_s⁵))` | `Γ₅=G₅P₀` | no (Gamma5≡G5out·P0≡a⁵P0/(27c_s⁵); mild tautology) |
| A16 | sympy | 177 | `expect_zero(P2.subs(N2,N2_const))` | `P₂=0` condition | yes |
| A17 | sympy | 178 | `expect_zero(P4.subs(...))` | `P₄=0` condition | yes |
| A18 | sympy | 179 | `expect_zero(K2.subs(...) - A·P0)` | `K₂=AP₀` collapse | yes |
| A19 | sympy | 180 | `expect_zero(K4.subs(...) - B·P0)` | `K₄=BP₀` collapse | yes |

Section III–VI (A7–A14, A16–A19) are substantive: SymPy independently series-expands / multiplies and compares against hand-written compiler formulas. These genuinely settle the rational identities the stage's deliverables (2)–(5) require. Section I (A2–A5) and Section II (A6) are tautological. A15 is a mild self-restatement (Gamma5 defined as exactly that product).

## Findings

### F1 — tautological_check

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage189_transfer_shape_prefactor_compiler_sympy_audit.py:81-89`

**What's wrong:**
The selected-branch identity check is a self-inverted round-trip. The script defines (line 83)
```
Rtarget_formula = Lambda0 * OmW2 * (1 - epseta) * (1 - epsW)**2 / (Zw * (1 + rho)**2)
```
which is exactly `Lambda0*(1-epseta) / T2_direct` (since `T2_direct = Zw*(1+rho)**2/(OmW2*(1-epsW)**2)`, line 82). Then line 84:
```
T2_selected = sp.simplify(Lambda0 * (1 - epseta) / Rtarget_formula)
```
and line 89:
```
expect_zero("Lambda_0 (1-epseta) / R_target - T_A^2", T2_selected - T2_direct)
```
By construction `Lambda0*(1-epseta)/Rtarget_formula = Lambda0*(1-epseta) * T2_direct/(Lambda0*(1-epseta)) = T2_direct`, so the residual is identically zero no matter what any of `Zw, rho, OmW2, epsW, epseta, Lambda0` are. The paper's boxed claim (notes:139-145) is the *selected-branch identity* `R_target·T² = Λ₀(1-ε_η)` — a physical relation that should be checked by independently constructing `R_target` from its own one-port definition (notes:242 gives `R_{target,A}` as the quantity such that `T_A² = Λ₀(1-ε_{η,A})/R_{target,A}`) and confirming the product equals `Λ₀(1-ε_η)`. As written, `R_target` is defined *by inverting the very identity being tested*, so the check verifies nothing about the physics. Additionally, the explicit value `Λ₀ = 27π²Gc_s⁵/(20a⁵c⁵)` (notes:144) is never instantiated — `Lambda0` stays a free symbol — so the script never confirms the front-factor constant.

**Why this matters:**
A reader of the output sees `Lambda_0 (1-epseta) / R_target - T_A^2 = 0` and concludes the selected-branch identity is verified. It is not — it is an algebraic tautology. If the true one-port form of `R_target` differed from `Λ₀(1-ε_η)/T²` (e.g., a sign or factor error in the continuum-kernel derivation), this check would still pass. The load-bearing `Λ₀` constant is also unverified.

**Required change:**
Replace the self-inverted definition with an independent construction. Define `R_target` from its one-port continuum meaning so the identity is a genuine consequence, OR (simpler and faithful to the notes) instantiate `Lambda0` to its stated value and verify the product `R_target*T2 == Lambda0*(1-epseta)` where `R_target` is an *independent* symbol/expression, not `Λ₀(1-εη)/T²`. Concretely: set `Lambda0_val = sp.Rational(27,20)*sp.pi**2*G*cs**5/(a**5*c**5)` and add a check `expect_zero("Lambda_0 value", Lambda0 - Lambda0_val)` only after binding it; and rewrite the selected-branch check so `R_target` is not algebraically forced to invert `T2_direct`. See directive F1 for the exact edit.

**Verification:**
After the fix, line ~89 region should contain a check whose residual is NOT identically zero by substitution of the `R_target` definition — i.e., the check should reference an independently-defined `R_target` (carried symbol or its own continuum formula) and confirm `R_target*T2_direct - Lambda0*(1-epseta) == 0`. The output line should still print `= 0`, but now the zero is a real cancellation rather than a `x/x` round-trip. The `Λ₀` value check should print `Lambda_0 value = 0`.

### F2 — tautological_check

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage189_transfer_shape_prefactor_compiler_sympy_audit.py:45-72`

**What's wrong:**
The Section I defect-notation and compatibility checks re-assert definitions and cannot fail. The script defines (lines 45-47):
```
dlnT2       = sp.simplify(dn - Bstar * dr)
dlnOneMinus = sp.simplify(-epsetas / (1 - epsetas) * de)
dlnRtarget  = sp.simplify(dlnOneMinus - dlnT2)
```
and (lines 67-69) `Theta = dr`, `Xi = dlnT2`, `Rcal = dlnRtarget`. Then:
- Line 70: `expect_zero("Xi - (dln N_* - B_* dln R_tr)", Xi - (dn - Bstar * dr))` — but `Xi := dlnT2 := dn - Bstar*dr`, so this is `(dn - Bstar*dr) - (dn - Bstar*dr) == 0`, tautological.
- Line 71: `expect_zero("R_1 - dln R_target", Rcal - dlnRtarget)` — `Rcal := dlnRtarget`, so `dlnRtarget - dlnRtarget == 0`, tautological.
- Line 72: `expect_zero("(R_1 + Xi_1) - dln(1-epseta)", (Rcal + Xi) - dlnOneMinus)` — `Rcal + Xi = (dlnOneMinus - dlnT2) + dlnT2 = dlnOneMinus`, so `== 0` by construction.
- Line 63: `expect_zero("selected-branch compatibility row", trf[2] + trf[0] - trf[1])` — `trf[2]=dlnRtarget=dlnOneMinus-dlnT2`, `trf[0]=dlnT2`, `trf[1]=dlnOneMinus`, so `(dlnOneMinus-dlnT2)+dlnT2-dlnOneMinus == 0` identically. The notes present this rank-2 compatibility relation (notes:217-224) as a *derived* consequence of the three identities, but the script bakes `dlnRtarget := dlnOneMinus - dlnT2` into the definition, so the compatibility "check" is guaranteed.

**Why this matters:**
The stage card `Output` is precisely the equating of `Ξ₁` with the transfer-shape/prefactor slope and the corrected-form identities; the appendix checklist (stage_189.tex:19-24) demands "Check that direct transfer-shape drift `Ξ₁` is kept separate from selected-branch dressing residuals." The Section I checks are the script's representation of those identities, yet they are vacuous — they confirm that a symbol equals its own definition. The genuine content (that the *Stage-239 observable* `δln N_*` decomposes as `Ξ₁+B_*Θ₁`, notes:122) is never the thing being tested here; the script just renames `dn - Bstar*dr` to `Xi` and checks the rename. The one non-tautological Section I check is the matrix reproduction (line 62, A1), which is fine.

**Why this is not UNFIXABLE / not paper_misalignment:**
The identities are correct; the defect is only that the checks don't exercise them. The fix is to introduce the independent upstream definitions (`δln N_* = Ξ₁ + B_*Θ₁`, `δln R_tr = Θ₁`, `δln ε_η = Σ_η` from the notes) as separate symbols and verify the transfer identities follow, rather than defining `Xi`/`Rcal` to be the answer.

**Required change:**
Introduce the Stage-239 observable relations as independent inputs and derive the transfer identities from them, so each `expect_zero` is a real cancellation. Specifically: declare `Theta1, Xi1, Sigma_eta` as independent defect symbols; impose the observable definitions `dn == Xi1 + Bstar*Theta1` and `dr == Theta1` and `de == Sigma_eta` (as substitutions or solved relations); then check `dlnT2 == Xi1`, `dlnRtarget == Rcal_expr`, and the compatibility row, where the transfer-packet components are computed from `obs` via the matrix and the defect substitutions — NOT defined to equal the target. See directive F2.

**Verification:**
After the fix, lines 70-72 (or their replacement) should reference independent defect symbols `Xi1, Theta1, Sigma_eta` on one side and the matrix-compiled packet on the other, so the residual is a genuine cancellation using the substitution `dn -> Xi1 + Bstar*Theta1`. The compatibility-row check (line 63) should derive `dlnRtarget` from the two physical identities `δln R_target = R₁` and `δln(1-ε_η) = R₁+Ξ₁` rather than the definitional subtraction, so it is non-vacuous.

### F3 — paper_misalignment

**Severity:** low
**Subtype:** notes_contradicts_script
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage189_transfer_shape_prefactor_compiler_sympy_audit.py:35,182` (banner/ledger say "STAGE 172")
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage189_transfer_shape_prefactor_compiler_sympy_audit.py:38,66,183` (comments/prose say "Stage 188" / "Stage 188 observable packet")
- vs paper `stage_189.tex` (this is Stage 189) and notes (body says "Stage 240" / source "Stage 239")

**What's wrong:**
The script's own labels are internally inconsistent and disagree with every paper-side label. The banner (line 35) and ledger (line 182) read "STAGE 172 — TRANSFER-SHAPE / OUTGOING-PREFACTOR COMPILER". The comments and ledger prose reference "Stage 188 observable packet" (lines 38, 66, 183). The filename and paper card are stage 189. The notes file body calls this "Stage 240" with upstream source "Stage 239". So four different stage numbers (172 / 188 / 189 / 240) appear across the paper-script-notes triple for what is one stage. This is the documented original-vs-restored renumbering convention (card: "restores original stage order"), so the *identity verified* is correct and aligned — but the banner literal "172" is plainly wrong for both the restored (189) and original (240) numbering and is not a renumbering artifact; it is a stale copy-paste. Because this touches only display strings (banner/ledger/comments), not any load-bearing assertion, it does not change the verdict on the math, but it obscures traceability and is exactly the kind of stale-revision signal the prompt flags (tactical reminder: "the script was written against a stale paper revision").

**Why this matters:**
A reviewer cross-referencing the output file (which prints "STAGE 172") against paper stage 189 cannot tell whether this is the right script for the stage. The "Stage 188" upstream reference vs the notes' "Stage 239" upstream reference compounds the confusion about which observable packet is actually being compiled. This is a low-severity hygiene/traceability finding, routed to the user because it spans the script labels and the paper/notes numbering convention (Codex must not unilaterally pick a number).

**Required change (routed to user — see directive Resolve block):**
Determine the canonical display number for this stage's script banner/ledger/comments. Likely the restored number 189 (to match filename + card) or the original 240 (to match notes body). Then the script's banner (line 35), ledger banner (line 182), and the "Stage 188"/"172" prose (lines 38, 66, 183) should be made consistent with whichever convention the pipeline uses for script labels. No assertion changes.

**Verification:**
After resolution, the script banner/ledger and the regenerated output `.txt` should show a single stage number consistent with the chosen convention; no `expect_zero` line changes.

## Independent-derivation check (Mathematica)

No `.wl` script exists for this unit. The `mathematica_transliteration` category does not apply. See "Engine cross-check" / missing-engine judgment below.

## Engine cross-check

Only the SymPy engine is present; `engines_agree` is `n/a`.

On the missing-Mathematica question (prompt line 110-118, this unit is `is_status_only_candidate: False`, `is_checkpoint: False`): I did NOT raise a `missing_mathematica` finding. The substantive deliverables of this stage — the grouped response/prefactor compiler `(u₂,u₄,P₀,P₂,P₄)`, the outgoing coefficients `(K₀,K₂,K₄,Γ₅)`, the weak-axisymmetric slope `P₁/P₀`, the normalization equivalence `m̂₀²P₀=54Gc_s⁵/(5a⁵c⁵)`, and the constant-prefactor branch — are all exact rational/series identities. SymPy settles them genuinely and independently: it computes `series(D₀/D_cons)`, `series(D₀N/D_cons²)`, `series((N₀+ελN₁)/(D₀+ελD₁))`, and `expand(Pref·Ŷ₂)` from the raw definitions and compares to the hand-written compiler formulas (A7–A14, A16–A19). A series expansion of a rational function is a closed, deterministic computation; a second engine would not add confidence to these identities — it would only re-confirm SymPy's exact-arithmetic series. This matches the prompt's line-114 judgment and established pipeline precedent (stages 121/122/123 verified SymPy-only as non-status-only). The genuine defects (F1, F2) are tautologies that would be tautological in either engine; adding a `.wl` would not fix them. There is no specific claimed result here that SymPy fails to verify *in principle* — the failures are in how the script wires the checks, not in SymPy's capability.

## Verdict justification

`verdict: findings`. The substantive core of the stage (Sections III–VI: the prefactor compiler, outgoing coefficients, normalization equivalence, constant-prefactor branch) holds up under attack — I tried to break the `P₄` and `K₄` formulas and the normalization-target consistency, and SymPy's independent series/product expansions genuinely confirm them against the notes' boxed compiler formulas, so those identities are aligned with the paper and real. What does NOT hold up: the Section II selected-branch identity check (F1) is a self-inverted round-trip that verifies nothing about `R_target`'s physics and never instantiates the `Λ₀` constant; and the Section I defect-notation/compatibility checks (F2) re-assert their own definitions and cannot fail. These are PASS-FOR-THE-WRONG-REASON tautologies, which is the single most important failure class for this audit, so the verdict is `findings` (medium severity, fixable in-script). F3 is a low-severity traceability/label inconsistency routed to the user. No `stop_cold`: the identities are correct, nothing is internally inconsistent, and no downstream-propagating constant changes. I confirm I read the stage card, the notes file, and the part-05 appendix rows (line 109 row matches the card `Output` verbatim), and the script's *verified* identities match the paper's claims — the defects are in check construction, not in claimed content, so `paper_alignment: aligned`.

## Self-test notes

I checked: (1) Variable independence / no zero-derivative traps — Sections III–IV use `sp.series` in `omega`/`eps`, and the expanded series genuinely depend on those variables (the printed output shows nonzero ω⁴ and ε terms), so no vacuous-derivative trap. (2) No unbounded-domain integrals appear, so the parity trap is N/A. (3) Trivial-case pre-check: I mentally substituted the F1 round-trip (`Lambda0*(1-epseta)/Rtarget_formula` reduces to `T2_direct` regardless of inputs → residual identically 0, confirming the tautology) and the F2 definitions (`Xi - (dn-Bstar·dr)` with `Xi:=dn-Bstar·dr` → 0 identically), confirming both are tautologies and that the prescribed fixes (independent `R_target`, independent defect symbols `Θ₁,Ξ₁,Σ_η`) make the residuals genuine cancellations. (4) Path: no missing-script finding, so no target-path reconstruction needed. (5) Paper round-trip: the F1 fix instantiates `Λ₀=27π²Gc_s⁵/(20a⁵c⁵)` exactly as the notes state (notes:144) and the F2 fix uses the notes' observable relations (`δln N_*=Ξ₁+B_*Θ₁`, notes:122), so neither fix introduces a new constant or a new paper_misalignment.

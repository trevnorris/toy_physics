---
unit_id: 004
batch: I.1
auditor_model: claude-opus-4-8
audit_date: 2026-06-04T00:00:00Z
verdict: clean
stop_cold: null
findings_count: 0
paper_alignment: aligned
scripts_checked:
  sympy: present
  mathematica: present
  engines_agree: true
  outputs_fresh: true
docs_read:
  paper_stage_tex: present
  notes_stage_files: []
  paper_appendix: present
---

# Audit unit 004 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_004.tex`
- notes: `(none)` — no files matching `notes/stages/moving_throat_pde_stage004_*.md` exist (confirmed by directory listing)
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part01.tex` (row 30 + narrative line 9)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage004_projected_maxwell_bundle_index_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage004_projected_maxwell_bundle_index_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage004_projected_maxwell_bundle_index_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage004_projected_maxwell_bundle_index_mathematica_audit.txt`

## What the paper claims

Stage 004 is an **index / bundle-ordering stage** for the projection-first electromagnetic sector. The card's `\stagefield{Output}` reads verbatim: *"Stage~004 fixes the source ordering for the projection-first EM bundle and anchors the projection-by-parts identity used in the next stage."* The one displayed equation (`eq:stage004-projection-ibp`) is the pointwise product-rule identity \(\partial_w(WQ)-W\partial_wQ-(\partial_wW)Q=0\); the card states that **after integration over \(w\)** this is "the bookkeeping identity that turns the parent \(w\)-flux into a boundary term plus a kernel-gradient leakage term." The "Bundle role" paragraph fixes three ordering rules (projection-by-parts creates a real transverse leakage/source term; homogeneous Maxwell projects in measured fields while inhomogeneous equations carry source-coupled flux; reduction is a matched channel, not the parent derivation). The "Source and audit" paragraph explicitly says the script "imports the bundle-level source index … and checks that the covariant, vector, and projection/reduction comparison audits are present in the ledger script order." The appendix row (line 30) summarizes it as "Bundle ordering and projection-by-parts identity," status `\StatusExact{}`. Distinct deliverables for stage 004 itself: (D1) the projection-by-parts/IBP identity, and (D2) the presence/ordering of the three downstream bundle scripts (005 covariant, 006 vector, 007 projection/reduction comparison). The Gaussian normalization, matched-kernel, and delta-source-ratio results that the script also runs are the substance of the downstream 005–007 stages, previewed here as part of the index role.

## What the script claims to verify

The SymPy docstring states it is a "Bundle-level audit for step_01_projected_maxwell_readme.md … an index note for the first projected-Maxwell bundle," verifying "the core formulas summarized there" and checking "that the three underlying derivation scripts are present." Concretely the assertions test: (1) existence of the three bundle scripts stage005/006/007 (`FileNotFoundError` if missing); (2) the integrated IBP identity \(\int W Q_w\,dw + \int W_w Q\,dw = 0\) with concrete decaying Gaussians \(W=e^{-w^2/\lambda^2}\), \(Q=w e^{-w^2/\lambda^2}\); (3) cyclic Bianchi for \(F=dA\) on the three space-time triples; (4) the three Maxwell-Faraday components built from \(E_i=-F_{0i}\), \(B=(F_{23},F_{31},F_{12})\); (5) Gaussian integrals \(\int Z=\sqrt\pi\,\lambda\), \(\int Z^2=\sqrt{2\pi}\,\lambda/2\), matched overlap \(\int W_{\rm match}Z=\sqrt2/2\); (6) the delta-source projection/reduction coupling ratio \(=\sqrt2\). The Mathematica `.wl` mirrors checks (2)–(6) as M1–M6 (it does not replicate the filesystem inventory check, which is not a math claim).

## Paper ↔ script cross-check

| Paper-side deliverable | Script-side check | Status |
|---|---|---|
| D1: projection-by-parts / IBP identity | py L41–44 integrated IBP; wl M1 L16–26 | match (script verifies the integrated, boundary-vanishing form that the card's prose explicitly invokes; the card's *displayed* equation is the trivial pointwise product rule, of which the integrated form is the load-bearing consequence) |
| D2: presence/ordering of bundle scripts 005/006/007 | py L23–30 file-inventory check | match (all three files confirmed present; check is non-vacuous) |
| (no stage-004 deliverable) Bianchi / Maxwell-Faraday | py L63–84; wl M2 L55–92 | extra — consistent with the EM-sector index role; no paper disagreement |
| (no stage-004 deliverable) Gaussian \(\sqrt\pi\lambda\), \(\sqrt{2\pi}\lambda/2\), \(\sqrt2/2\), ratio \(\sqrt2\) | py L96–99; wl M3–M6 | extra — these are the previewed deliverables of downstream stages 005–007 (the matched-reduction / coupling-mismatch results), run here in the index; not stage-004 outputs |

`paper_alignment: aligned`. Every stage-004 deliverable (D1, D2) has a faithful, non-tautological script check. The additional checks are previews appropriate to the declared "bundle index" role; none disagrees with any paper value, so no `paper_missing_script_claim` finding is warranted (the index stage legitimately runs bundle-wide preview checks, and the augmentation guards exclude flagging values that are deliverables of *other* stages).

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 28–30 | `FileNotFoundError` if any of 005/006/007 missing | D2 (bundle ordering/presence) | yes |
| A2 | sympy | 41–44 | `assert_zero(ibp_lhs - ibp_rhs)` | D1 (IBP, integrated form) | yes |
| A3 | sympy | 63–68 | `assert_zero(cyc)` ×3 triples | none (EM-sector preview) | yes (to its own sign-error claim) |
| A4 | sympy | 76–84 | `assert_zero(mf1/mf2/mf3)` | none (preview) | yes (to E,B↔F sign claim) |
| A5 | sympy | 96 | `assert_zero(Z_int - sqrt(pi)*lam)` | none (005–007 preview) | yes |
| A6 | sympy | 97 | `assert_zero(Z2_int - sqrt(2pi)*lam/2)` | none (preview) | yes |
| A7 | sympy | 98 | `assert_zero(I_WZ - sqrt(2)/2)` | none (preview) | yes |
| A8 | sympy | 99 | `assert_zero(mu_proj_delta/mu_red - sqrt(2))` | none (preview) | yes |
| A9 | math | 23–25 | `If[FullSimplify[m1Left]=!=0, Exit[1]]` | D1 (IBP) | yes |
| A10 | math | 60–70 | `If[residual=!=0, Exit[1]]` ×3 Bianchi | none (preview) | yes |
| A11 | math | 82–92 | `If[residual=!=0, Exit[1]]` ×3 Maxwell-Faraday | none (preview) | yes |
| A12 | math | 102–104 | `If[...=!=0]` Z_int = sqrt(Pi)*lam | none (preview) | yes |
| A13 | math | 114–116 | `If[...=!=0]` Z2_int = sqrt(2Pi)*lam/2 | none (preview) | yes |
| A14 | math | 127–129 | `If[...=!=0]` overlap = sqrt(2)/2 | none (preview) | yes |
| A15 | math | 140–143 | `If[...=!=0]` ratio = sqrt(2) | none (preview) | yes |

Note: A1 (file inventory) is Python-only and non-vacuous (all three named files verified present on disk). It is a filesystem check, not a math identity, so its absence from the `.wl` is not `engine_disagreement`.

## Findings

None. The audit was adversarial across IBP boundary-term vanishing, Gaussian integral closed forms, the Bianchi/Maxwell-Faraday sign maps, the inventory check, and the script↔paper alignment; all held up. See the value reconciliation below — every emitted deliverable value reconciles.

## Independent-derivation check (Mathematica)

The `.wl` is **not** a line-by-line transliteration of the `.py`. It re-derives each result from the physical premises with idiomatic Mathematica structure:
- M1 (IBP, wl L16–26) combines `W*f' + W'*f` into one integrand and lets `Integrate` collapse it to the boundary term, whereas the `.py` (L41–42) computes the two integrals **separately** (`ibp_lhs`, `ibp_rhs`) and subtracts. Different evaluation route, same identity.
- M2 (wl L40–92) precomputes the six field-strength components as integer-literal `fStr[i,j]` expressions and substitutes them literally into each loop body — an explicit workaround for Mathematica's `Module`/`Part` pre-evaluation trap; the `.py` (L55–68) uses a Python closure `F(mu,nu)`. Independent choreography.
- M3–M6 (wl L94–144) recompute the Gaussian normalization, squared norm, matched overlap, and delta-source ratio directly; the `.py` (L86–99) computes them through intermediate symbols `W_match`, `I_WS_delta`, `mu_proj_delta`, `mu_red`. Same targets, independently assembled.
Conclusion: genuine second-engine derivation, no `mathematica_transliteration` finding.

## Engine cross-check

Both saved outputs report `STATUS: PASS`. The `.wl` output prints every residual as literal `0` (M1–M6, including the three Bianchi triples and three Maxwell-Faraday components). The `.py` output prints the summary PASS (all `assert_zero` calls survived). The shared targets agree exactly:

| Quantity | sympy assertion | math residual | agree? |
|---|---|---|---|
| IBP integrated | `ibp_lhs - ibp_rhs == 0` | M1 = 0 | yes |
| Bianchi ×3 | `cyc == 0` | M2 Bianchi = 0,0,0 | yes |
| Maxwell-Faraday ×3 | `mf1/2/3 == 0` | M2 MF = 0,0,0 | yes |
| \(\int Z = \sqrt\pi\lambda\) | `Z_int - sqrt(pi)*lam == 0` | M3 = 0 | yes |
| \(\int Z^2 = \sqrt{2\pi}\lambda/2\) | `Z2_int - sqrt(2pi)*lam/2 == 0` | M4 = 0 | yes |
| matched overlap \(\sqrt2/2\) | `I_WZ - sqrt(2)/2 == 0` | M5 = 0 | yes |
| delta ratio \(\sqrt2\) | `mu_proj_delta/mu_red - sqrt(2) == 0` | M6 = 0 | yes |

No `engine_disagreement`.

## Verdict justification

`clean`. I read the paper card, the absent-notes status, and the appendix row before opening the scripts, and built the deliverable model (D1 IBP identity, D2 bundle-script presence). Attacks attempted and survived: (a) the IBP check is non-tautological — it computes two separate integrals and subtracts, and depends on the Gaussian boundary term genuinely vanishing (a non-decaying probe would break it); (b) the Gaussian closed forms were hand-verified: \(\int e^{-w^2/\lambda^2}=\sqrt\pi\lambda\), \(\int e^{-2w^2/\lambda^2}=\sqrt{\pi/2}\lambda=\sqrt{2\pi}\lambda/2\), matched overlap \(=Z2\_int/Z\_int=\sqrt2/2\), delta ratio \(=(1/(\sqrt\pi\lambda)\cdot 2/\sqrt2)/(1/(\sqrt\pi\lambda))=\sqrt2\) — all match the asserted targets exactly; (c) the file-inventory check is non-vacuous (the three named neighbor scripts exist on disk); (d) the Bianchi/Maxwell-Faraday checks are mathematical identities for any smooth \(A\), so they only discriminate sign/index errors in the \(E,B\leftrightarrow F\) map — which is precisely the claim their comments assert, so they are adequately anchored. The paper↔script alignment is exact for the stage-004 deliverables; the extra preview checks are consistent with the declared "bundle index" role and disagree with no paper value. Positivity/scale assumptions (`lam_ibp>0`, `lambda>0`, `mu0>0`) are physically justified (length scale, coupling) and required for the Gaussian integrals to converge to the asserted forms. Outputs are fresh (both `.txt` mtimes post-date their scripts).

## Self-test notes

Checked the v2 traps: (1) Variable independence — every `sp.diff`/`D[]` is taken w.r.t. a coordinate the integrand/field genuinely depends on (the Gaussians depend on `w`; the `F` and Maxwell-Faraday derivatives are in `t,x,y,z` on functions `A0..A3(t,x,y,z)`), so no derivative is identically zero and no `assert_zero` passes vacuously. (2) Parity/symmetry — the IBP integrand `d/dw(WQ)` with `WQ = w e^{-2w^2/λ^2}` is an odd, decaying function, so its symmetric integral is genuinely 0 by both parity and the vanishing boundary term, matching the assertion. (3) Trivial-case — substituting the concrete Gaussians reproduces the asserted closed forms exactly (verified by hand above). No trap fired; verdict stands `clean`.

## Value Reconciliation (pass-2 augmentation)

Enumerating every RESULT/deliverable value the scripts emit (source `.py`/`.wl` + saved `.txt`), excluding scaffolding. The stage-004 card is, by design, an **index card**: its sole displayed result is the projection-by-parts identity (D1). The Gaussian/coupling constants are previews of the downstream 005–007 bundle deliverables (the script docstring and the card's own "Source and audit" paragraph both frame this stage as a bundle index that imports/checks neighbor results). Per the augmentation guards, a value that is a **deliverable of another stage** and is correctly absent from the terse stage-004 card is not a stage-004 MISSING-DELIVERABLE.

| value | source (py/wl + output line) | .tex/.md location | status |
|---|---|---|---|
| IBP identity \(\partial_w(WQ)-W\partial_wQ-(\partial_wW)Q=0\) (integrated form) | py L41–44; wl M1; sympy.txt L2, math.txt L2–3 | `stage_004.tex:27–33` (`eq:stage004-projection-ibp` + prose) | MATCH |
| presence/ordering of bundle scripts 005/006/007 | py L23–30 | `stage_004.tex:9–12` ("checks that the covariant, vector, and projection/reduction comparison audits are present in the ledger script order") | MATCH |
| \(\int e^{-w^2/\lambda^2}\,dw = \sqrt\pi\,\lambda\) | py L96; wl M3; math.txt L16 | not in stage_004 card/notes — deliverable of downstream 005-bundle (matched reduction); previewed here | INTERNAL (preview; not a stage-004 deliverable) |
| \(\int e^{-2w^2/\lambda^2}\,dw = \sqrt{2\pi}\,\lambda/2\) | py L97; wl M4; math.txt L18 | not in stage_004 card/notes — downstream preview | INTERNAL (preview) |
| matched overlap \(\int W_{\rm match}Z\,dw = \sqrt2/2\) | py L98; wl M5; math.txt L20 | not in stage_004 card/notes — downstream preview | INTERNAL (preview) |
| delta-source projection/reduction ratio \(=\sqrt2\) | py L99; wl M6; math.txt L22 | not in stage_004 card/notes — downstream (coupling-mismatch) preview | INTERNAL (preview) |
| cyclic Bianchi residual \(=0\) (×3) | py L63–68; wl M2; math.txt L4–9 | covered by D1/EM-sector framing; identity, no numeric carrier | INTERNAL (identity check) |
| Maxwell-Faraday residual \(=0\) (×3) | py L76–84; wl M2; math.txt L10–15 | EM-sector framing; identity, no numeric carrier | INTERNAL (identity check) |

INTERNAL scaffolding (accounted for, no finding): pass/fail flags and `STATUS: PASS`, all `residual = 0` print lines, intermediate symbols `W_match`/`I_WS_delta`/`mu_proj_delta`/`mu_red` (py) and `matchedWeight`/`pointSourceCoupling`/`volumeReducedCoupling` (wl), and the `missing` list driving the `FileNotFoundError`.

Reconciliation result: both stage-004 deliverables (D1, D2) reconcile to the card; the remaining emitted values are correctly classified as downstream-bundle previews (deliverables of stages 005–007), legitimately omitted from the terse stage-004 index card under the augmentation guards. No MISMATCH and no stage-004 MISSING-DELIVERABLE.

`reconciliation: complete; 8 deliverable-level values checked, 0 misaligned` (2 stage-004 deliverables MATCH; 6 emitted constants/identities classified INTERNAL as downstream previews).

---
batch_id: IV.1
range: 91-102
created_at: 2026-05-27
status: pending_user
auditor_pass: v2 (first-pass paper-grounded)
clusters: 3
affected_stages: [091, 095, 097, 098, 099, 100, 101, 102]
---

# Batch IV.1 paper-alignment resolution

The IV.1 audit surfaced **three clusters** of `paper_misalignment` findings across 8 stages. The remaining 4 stages (092, 093, 094, 096) had no paper-alignment findings (096's `insufficient_verification` F2 is script-side).

Total `paper_misalignment` findings: 10 (across 8 stages).

The clusters are independent — the user picks a direction for each.

---

## Cluster A — Paper card `\stagefield{Checks}` list items not exercised in scripts

**Stages affected:** 091 (F1), 095 (F3), 097 (F1), 099 (F2)

**Pattern:** Each stage card's `\stagefield{Checks}` block lists items the audit scripts do not exercise. The recurring items are:

1. **"Check `l=0` and `l=2` orthogonality before applying the geometry firewall"** — appears in 091, 095, 097, 099. The full orthogonality is verified at **stage 094** ("Isotropic Geometry-Decoupling Theorem"), which already executes 15 angular integrals plus the Laplace eigenvalue check in both engines.
2. **"Check the static limit `eps_2=eps_4=0` returns `c_pole=1/4`"** — appears in 097, 099. The static-limit reduction is verified at **stages 091, 092, 094, 096** (different facets).
3. **"Check that any support/source success statement still carries the minimal-module hypothesis"** — appears in 097, 099. Minimal-module is the Part III chain, anchored at stage 088 + carried through stages 089/090 (the III.5 checkpoints).

The Checks read as if each stage is supposed to re-run them, but each is a Part-level deliverable upstream of these stages.

### Cluster A options

- **(a) Script-side docstring carry-forward** (Recommended). Each of 091/095/097/099 adds a docstring comment naming the upstream verifying stage(s): "orthogonality carries from stage 094; static-limit carries from stages 091/096; minimal-module hypothesis carries from Part III chain stages 088–090." No new assertions added. Mirrors the III.5 F1 087 consolidation pattern: orchestrator-direct, no paper edits, math content untouched.
- (b) Script-side re-run. Each of 091/095/097/099 adds full orthogonality / static-limit / minimal-module checks. Duplicates work already in 094/091/088-090; substantial new script content (esp. for orthogonality which is 15 integrals).
- (c) Paper-side edit. Reword the Checks blocks in the four cards to read as carry-forward references (e.g., "Inputs: orthogonality from §094"). Defer to a follow-up paper-side directive after user confirms direction.

### Cluster A destination verification

```bash
# Stage 094 orthogonality coverage (the upstream anchor)
grep -n "Orthog\|orthog\|l=0\|l=2\|Y_2\|LegendreP\|spherical harm" \
  scripts/moving_throat_pde_stage094_isotropic_geometry_decoupling_sympy_audit.py \
  mathematica/moving_throat_pde_stage094_isotropic_geometry_decoupling_mathematica_audit.wl | head -30
# Confirm the orthogonality block is substantive (per audit report 094: "15 angular integrals + Laplace eigenvalue")
```

---

## Cluster B — Stage 100/101 DtN-fingerprint + higher-odd-term Checks

**Stages affected:** 100 (F1+F2+F3 coupled), 101 (F4)

**Pattern:** Stage cards 100 and 101 each list three `\stagefield{Checks}`:
- Check 1: keep source / conservative / outgoing factors separate (`mhat_0^2 chi_Q N_Q = 1`)
- Check 2: higher odd terms begin beyond the point-particle 2.5PN coefficient (`O(omega^7)`)
- Check 3: outgoing `l=2` DtN fingerprint against `z = omega a/c_s` expansion (`Yhat_2^out(z) = 1 + z^2/9 + 4 z^4/81 + i z^5/27 + O(z^6)`)

Per notes/stages/moving_throat_pde_stage101_natural_source_map_reduction.md line 41–51: "the canonical compact branch's `chi_Q = 1` identification is from Stage 80" — i.e., **stage 097** owns Check 3. Per appendix part04 line 295: the higher-odd statement is presented in a single sentence without scripted backing; it lives in the **stage 102** "higher odd irrelevance" audit (which does verify omega^5 → omega^7 series structure). The Check 1 closure `mhat_0^2 chi_Q N_Q = 1` is the stage-100 Output: stage 100 must own it.

Stage 100's current A3/A8 verify `Gamma_5/Gamma_5^target - chi_Q*N_Q = 0` — but this is tautological under the script's construction `K_n = N_Q * K_n^target` plus the literal definition `Gamma_5 = 9*K_2^(5/2)/K_0^(3/2)`. The closure isn't *derived* from the observable condition; it's pre-baked by symbol choice.

### Cluster B options

- **(c) Hybrid: script-side docstring carry-forwards + targeted closure assertion in 100** (Recommended).
  - **Stage 100**: Add docstring naming stage 097 as the Check 3 anchor and stage 102 as the Check 2 anchor. Strengthen the Output closure by imposing `mhat_0^2 * Gamma_5 = Gamma_5_target` as the observable condition (not by construction), then deriving `mhat_0^2 * chi_Q * N_Q = 1` non-tautologically. The script already has all the symbolic infrastructure (`Gamma_5`, `Gamma_5_target`, `N_Q`, `chi_Q`, `mhat_0`); the change is replacing `solve(...) - 1/(mhat_0^2 * chi_Q)` with a substantive derivation. F2 (tautological A5/A10) and F3 (algebraic-shadow A4/A9) close automatically.
  - **Stage 101**: Add docstring naming stage 097 (Check 3) and stage 102 (Check 2) as upstream anchors. No new assertions for those Checks.
- (a) Strengthen both scripts fully. Stage 100 adds full DtN-side derivation (`Lambda_2^out(z) = z d/dz ln h_2^(1)(z)`) and a `omega^7` series extension. Stage 101 adds the same DtN and higher-odd-term checks. Large content addition; duplicates stage 097/102 work.
- (b) Paper-side trim. Reword the Checks blocks in 100/101 cards as carry-forward references to 097/102. Defer to a paper-side directive.

### Cluster B destination verification

```bash
# Stage 097 chi_Q = 1 / DtN-fingerprint coverage
grep -n "chi_Q\|chiQ\|Hankel\|DtN\|fingerprint\|z = " \
  scripts/moving_throat_pde_stage097_single_normalization_defect_sympy_audit.py \
  mathematica/moving_throat_pde_stage097_single_normalization_defect_mathematica_audit.wl | head -20
# Stage 102 omega^7 series coverage
grep -n "omega\|tauQ\|chiQ\|Series\|Coefficient" \
  scripts/moving_throat_pde_stage102_higher_odd_irrelevance_sympy_audit.py \
  mathematica/moving_throat_pde_stage102_higher_odd_irrelevance_mathematica_audit.wl | head -20
```

---

## Cluster C — Stage-number inconsistency (banner sweep + paper card titles)

**Stages affected:** all 12 (091–102)

**Pattern:**
- Every paper card section title in IV.1 carries a stale number from a previous renumber pass: 091→108, 092→109, 093→110, 094→111, 095→112, 096→113, 097→114, 098→115, 099→116, 100→117, 101→118, 102→119. Labels (`\label{stage:NNN}`) are correct.
- Every script banner is also stale, ranging from "STAGE 74" through "STAGE 085" plus some bare-banner cases (`.py` without explicit banner exist for 093, 100, 101, 102 — but 093 has no `.py` at all per the status-only carve-out).
- Auditor reports flagged the title mismatch as `paper_misalignment` in 098 (F1), 099 (F3), and 102 (cosmetic note, not filed). The other 9 auditors did not file it but the pattern is identical.

### Cluster C options

- **(a) Script-side banner sweep across all 12; paper-side card titles deferred** (Recommended). Mirrors III.5 banner sweep (12 scripts). Paper card titles are pure presentation typos; they don't affect math content, labels, file paths, or cross-references. Flag for paper-cleanup tracker as a deferred batch item. Same approach as III.5.
- (b) Script + paper sweep. Banner sweep plus 12 paper section title corrections. Larger surface area; requires paper-side authorization.
- (c) Status quo. Skip both. Banners remain stale; cards remain stale.

### Cluster C destination verification

```bash
# Confirm all 12 stale banners
for n in 091 092 093 094 095 096 097 098 099 100 101 102; do
  echo "--- $n ---"
  grep -nE 'STAGE [0-9]+' "scripts/moving_throat_pde_stage${n}"_*sympy*.py 2>/dev/null | head -1
  grep -nE 'STAGE [0-9]+' "mathematica/moving_throat_pde_stage${n}"_*.wl 2>/dev/null | head -1
done
```

---

## Other script-side findings (no user gate needed)

These ride on the orchestrator's normal apply pass after the three clusters resolve:

- **092 F1** (mathematica_transliteration): restructure Mathematica audit in `(eps_2, eps_4)` dimensionless vars from the outset.
- **092 F2** (insufficient_verification): assert the three first-order series coefficients (`1/4`, `-1/2`, `1/4`).
- **094 F1** (engine_disagreement): align the Mathematica cross-term `cCross` definition to the SymPy `Ccross` definition (`tw*overlap` not `tw*gradCross`).
- **094 F2** (insufficient_verification): add static-limit `eps_2=eps_4=0 → c_pole=1/4` assertion block in both engines.
- **095 F1** (insufficient_verification): replace 5 tautological SymPy asserts with three closed-form residual checks mirroring the Mathematica side.
- **095 F2** (mathematica_transliteration): add Schur-derivation block via `Eliminate`/`Solve` before the existing series machinery.
- **096 F1** (tautological_check, checkpoint): delete tautological `eps_2 == 0`/`eps_4 == 0`/`zeta_req - c_pole/c_geom == 0`; keep printed values.
- **096 F2** (insufficient_verification, checkpoint): append "HYPOTHESIS CARRIED" print block in both engines.
- **097 F2** (mathematica_transliteration): rework Mathematica via `SeriesCoefficient[kbarCons[w], {w, 0, n}]` instead of the same direct-definition route as SymPy.
- **098 F2** (hardcoded_result): add provenance comment naming the upstream source of `zeta_max^(F1) = 2.46752922945601` (per [[feedback-mathematica-script-idioms]] pitfall pattern; mirrors III.5 089 F4 resolution via provenance comment).
- **098 F3** (insufficient_verification): mirror Mathematica's `expectApprox` numeric pins on the SymPy side.
- **099 F1** (tautological_check): replace 4 `R_n - (N_Q - 1)` asserts with non-tautological structural checks (branch identity, `Gamma_5` form, `2 G/(5 c^5)` factorization).
- **099 F4** (insufficient_verification): assert `Yhat_Q^cons` static slot (=1) and pole residue (= -Omega_Q/8).
- **100 F4** (mathematica_transliteration): **BLOCKED** in directive — design-level rewrite, surfaced separately if Cluster B direction (c) is approved.
- **100 F5** (symbol_assumption_error): `chiQ` should be declared real, not positive (since it's pinned to 1 by DtN comparison upstream).
- **101 F1** (tautological_check): replace 3 `expectZero` calls with anchors to input equation `mhat0^2 * chi_Q * n_Q - 1`.
- **101 F2** (script_doesnt_cover_claim, high): SymPy script has zero `assert` statements — only `print` calls. Add 4 substantive asserts.
- **101 F3** (insufficient_verification): assert linearization `N_Q - 1 = -Delta_Q + O(Delta_Q^2)` against its paper-stated form.
- **102 F1** (insufficient_verification): SymPy script lacks `assert` statements paralleling Mathematica's `expectZero`. Add three asserts for `D1`, `D2`, `D3`.

Total script-side findings (excluding cluster A/B): 18 across 10 stages.

---

## Process notes

- Per `[[feedback-codex-iterates-until-clean]]` + III.4 stall lesson: Codex is bypassed; orchestrator applies directly.
- Per `[[feedback-no-parallel-exec-sympy]]`: outputs refreshed serially via `python3 …` and `math -script …`, never `$RT exec-*` in parallel.
- Per `[[feedback-mathematica-single-seat]]`: only one `math -script` invocation at a time.
- Per `[[feedback-no-fake-scripts]]`: applies happen by reading + editing, not commentary scripts.
- Per `[[feedback-squash-followup-fixes]]`: small post-batch fixes squash into the prior commit rather than separate commits.

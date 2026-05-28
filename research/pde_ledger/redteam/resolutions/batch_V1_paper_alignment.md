---
batch: V.1
created_at: 2026-05-28
clusters_resolved: 1
mechanical_groups: 6
status: directions_set
---

# Red-team batch V.1 — paper-alignment resolution

12 audits landed (range 164-175). 22 findings raised across all 12 (dirty) stages; 0 clean. One user-gated cluster surfaced (Cluster A, script banners); the standard `(Recommended)` auto-apply direction was confirmed up front and applied. The single paper_misalignment (170 F1) was resolved by explicit user choice. Six mechanical groups applied without a gate.

## Cluster A — script-banner renumbering (mechanical, −17 offset)

**Pattern:** every stage's `.py`/`.wl` opening `banner(...)` printed `STAGE N−17` (the legacy offset):
164→147, 165→148, 166→149, 167→150, 168→151, 169→152, 170→153, 171→154, 172→155, 173→156, 174→157, 175→158.

Unlike IV.x, **no notes-H1 offset and no paper-card-banner offset surfaced** — V.1's Cluster A is entirely script-side print banners.

**Direction:** fix all in place (banner string only; no other content).

- 11 stages had their banners fixed inside the per-finding directives (transliteration / stale_output / mislabel findings that named the banner).
- **4 residual banners were not flagged as findings** (the auditors only named a banner where they raised a finding, and 164/174 transliteration directives were scoped `.wl`-only): `164` `.py`, `165` `.py` + `.wl`, `174` `.py`. These were mass-fixed in place under the pre-authorized Cluster A "fix all" direction (sed, banner string only), and their outputs refreshed.

After the fix, a full sweep confirms zero `STAGE 1[3-5][0-9]` banners remain in any 164-175 script.

## Cluster B — body-text forward-stage citation re-attribution

**None this batch.** No notes prose in stages 164-175 carried legacy upstream stage numbers needing re-attribution. (Documented per the post-batch tracker discipline: absence is recorded, not skipped.)

## Cluster C — paper-card `Checks` downgrade

**None this batch.** Stage 170's `\stagefield{Checks}` item — "Check the weak-axisymmetric signature (1,1/2,-1) before reducing grouped defects to a scalar" — was content-disambiguated against the downstream stages that mention the same signature. Stages 171 (`microscopic_grouped_obstructions`) and 173 (`axisymmetric_loading_mismatch`) verify the `(1,½,−1)` lane *pattern* only on their own quantities (microscopic obstructions K_micro/G_micro; load observable Δ_Q), **not** on stage 170's outlet maps δκ_W/δγ_W. Since the signature inheritance into 170's outlet maps (with amplitudes κ1/γ1) is stage 170's *own* Sec. 5 deliverable, a downstream carry-forward citation (the Cluster C shape, cf. IV.6 stage 158) does **not** fit. Resolved as a script-side addition instead — see the paper_misalignment resolution below.

## Paper_misalignment — stage 170 F1 (user resolution)

**Finding:** `script_missing_paper_claim`. The card's Checks item and notes Sec. 5 enumerate the weak-axisymmetric `(λ_20,λ_21,λ_22)=(1,½,−1)` signature inheritance into δκ_W/δγ_W and the scalar amplitudes κ1, γ1, but neither stage-170 script verified it.

**User direction chosen (2026-05-28): (a) add the missing script check** (no paper-side edit). Rationale: the deliverable is stage 170's own (notes Sec. 5), the maps δκ_W/δγ_W are already verified linear, and the signature inheritance is a low-risk corollary that makes 170 self-verifying against its card.

**Applied (orchestrator-direct, both engines):** a "Weak-axisymmetric signature (1, 1/2, -1)" section that feeds lane-scaled bundle defects `δD_(A,n)=ε·λ_A·D_n^(1)`, `δN_(A,0)=ε·λ_A·N_0^(1)` (λ=(1,½,−1)) through the SAME verified linear outlet maps and asserts:
- `δκ_W^(2A) = ε·λ_A·κ1` with `κ1 = 3(1−σ_*)(D2^(1)+D0^(1)/9)/(σ_* D0)`;
- `δγ_W^(2A) = ε·λ_A·γ1` with `γ1 = −(1−σ_*)(N0^(1)−P0·D0^(1))/(9σ_* N0)`;
- the grouped ratios `(21)=½·(20)`, `(22)=−(20)` for both κ and γ.

+10 checks per engine (sympy 21 `= 0`, mathematica 21 PASS). Non-tautological: a wrong κ1/γ1 literal leaves a nonzero residual. No paper-side file touched; no carried result changed.

## Mechanical groups M1-M6 (applied without gate)

### M1 — `.wl` banner mass-rename
Fixed in 9 stages' Mathematica scripts (166-175 where flagged) plus the 165/164/174 residuals. All `STAGE N−17` → `STAGE N`.

### M2 — `.py` banner mass-rename
Fixed in the SymPy scripts where flagged (166, 167, 168, 169, 170, 171, 172, 173, 175) plus residuals 164, 165, 174. All `STAGE N−17` → `STAGE N`.

### M3 — Notes H1 mass-rename
**None for V.1** (no notes-H1 offsets surfaced). Empty group.

### M4 — Directive substance per stage
Each non-paper_misalignment finding's "Required change" applied orchestrator-direct (Codex bypassed). Substantive script-side edits: 164 (`Series` route), 165 F1 (real log-derivative assert) + F2 (tautology demoted to prints), 166 F1 (4 general-inversion asserts/engine) + F2 (`Inverse[Mmat]` route), 169 F1 (3 numeric coefficient asserts/engine), 170 F2 (`D[…,eps]`+direct-`Solve` route), 171 F1 (collected-literal + `Series` route, after rework), 172 F1 (implicit-diff route), 173 F1 (`Coefficient[Series[…]]` route), 174 F1 (`D[…,eps]` perturbation), 175 F1 (minimal resolution) + F2 (`Ξ_load` Σρ=1 aggregate). Two transliteration findings accepted as policy mirrors (169 F2, 175 F3-step3).

### M5 — Rework loop reruns
Three orchestrator catches required re-runs: 166 (round-trip vector-residual → `Total[(…)^2]`), 175 (F1 cross-check near-trivial → minimal resolution), 171 (directive route tautological → collected-literal + `Series`; this one surfaced from a `needs_rework` verification and was re-verified clean by a fresh agent). Final sweep: 24/24 scripts exit 0, 0 FAIL.

### M6 — Stale output refresh
All 24 output `.txt` files regenerated by direct `python3` / `math -script` (not via `$RT exec-*`, per the no-parallel-exec rule). Freshness audit confirms every output mtime is newer than its script.

## Closing

22 findings resolved (0 blocked_legitimate). 12 stages verified. Material change in 1 stage (170 — additive Sec. 5 deliverable coverage; zero downstream propagation). Thirteenth consecutive batch clear of stop-cold; first non-cluster paper_misalignment (170 F1), resolved additively by user direction (a).

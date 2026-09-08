# S11c-b carrier ablation HARNESSES — review record (two-leg, review-until-clear)

**Artifacts (astra-written, Codex-authored ⇒ legs = a fresh Claude agent + Grok):**
- `scripts/S11c_b_carrier_ablation_harness_sympy.py` (421 ln) — SymPy harness
- `mathematica/S11c_b_carrier_ablation_harness.wl` (233 ln) + driver `scripts/S11c_b_carrier_ablation_harness_wl.py` (96 ln) — WL harness
- run transcripts (evidence): `_measurements/S11c_b_carrier_ablation_harness_{sympy,wl}.md`

**Role.** These harnesses certify that the S11c-b brane operator's **pressure-slot carrier** —
`∂(operator row)/∂(native pressure atom) |_{atoms→0}` over every scalar row × pressure slot — is computed on the
LIVE emitted operator and that each knife (K_A/K_T/K_W) is a genuine one-site × one-FORM structural ablation, not
a fabricated "bite". The carrier is the object c2's N6 binds as its supplied premise. Knife-list cleared at
`21988ca7`; design + expected cones in `_measurements/S11c_b_carrier_ablation_harness_{reauthor_review,adjudication_key}.md`.

## Legs

**Leg 1 — fresh Claude agent (Sonnet), full pair (SymPy + WL): SOUND.**
Ran 6 bounded CAS builds ablating BOTH engines. Findings: a no-op knife mutation drives the carrier diff
**identically 0 in both engines** (harness cannot fabricate a bite); a coefficient rescale is **not** misreported
as a structural bite; the harness's three self-tests are genuine; wrap-engine / print-not-PASS / complete-carrier
compliant; the leg modified no working-tree file (ablated /tmp copies). Raw report: this workstream's session
transcript, 2026-09-07/08 (Agent leg; not separately filed).

**Leg 2 — Grok, independent WL FORM ablation: SOUND.**
Report: `_measurements/S11c_b_carrier_ablation_harness_wl_leg_grok.md`. Comparison instrument (leg-authored):
`/tmp/gwl/compare_carriers.wl`. Six WL kernels, one seat at a time, `timeout 600`, /tmp copies only, working tree
untouched. Grounded results (literal stdout in the report):
- **Wrap-engine (read):** carrier computed on `evaluatedModel["EULERIAN","MATERIAL_ADVECTED","RHO4_CONSTANT"]["OPERATOR"]`
  (harness L183-186) via `extractCarrier` `D[row,atom]/.(atoms→0)` (L121-132); no reimplementation, no hand-typed
  carrier; `modelRecord` after the marker not loaded.
- **One site × one FORM (read):** `patches` L66-80 + `StringReplace spec[[2]]→spec[[3]]` L165-166; `validateSites`
  L91-112 pins each before-string to a single occurrence in-source AND in its named-function window. Tower flux at
  `engine:2862` is a different string after the marker ⇒ not loaded.
- **No-op ⇒ diff identically 0:** with each knife's patch made identity (`spec[[3]]:=spec[[2]]`), K_A/K_T/K_W
  carrier diffs are all identically 0; BASELINE `OPERATOR_SHA256 8cba37be…` unchanged under all three no-ops.
- **Rescale ⇒ coefficient, not structure:** `RESCALE_K_A` adds exactly +1× the baseline `MASS_EVOLUTION_ROW`
  affinity carrier `−λ_A0/(ρ_M(1−iωτ_A))` (both faces) — the same term astra's live K_A bite removes at coeff −1;
  no vanished/new coupling on momentum or thickness.
- **DEAD_PATH (pressure-free) ⇒ carrier unmoved:** deleting the pressure-free `kineticEwLive` from `THICKNESS_ROW`
  moves the operator digest (`8cba37be…→6f101c16…`) but leaves the carrier diff identically 0 — the extractor
  responds only to pressure-bearing structure.
- **Extractor order:** default is `∂` then `→0` (physically correct for `∂/∂p|_{p→0}`); `→0`-first kills the linear
  response (all 0); the self-test prints the disagreement rather than silently swapping.
- **PRINT-not-PASS (read):** no PASS/FAIL/verdict/"bite" payload; `emitTriple` prints `{baseline,corrupted,diff}`
  then guards shape + `$Aborted|$Failed|Indeterminate`; the only `SameQ` stop is the canonical-copy drift guard.

## Orchestrator verification (G4) — mechanical fact-lookups

Grok's independent WL run reproduced astra's committed operator digests exactly (confirms it ran the real harness
against the real engine, agreeing with the committed baseline):
- BASELINE digest `8cba37be5fe30437e6bf09c69a6e51767a635d809d4f803995cebd1808f4334a` — present **3×** in
  `_measurements/S11c_b_carrier_ablation_harness_wl.md` (BASELINE, UNABLATED_COPY, report baseline).
- `RESCALE_K_A` digest `4519e9b9…` and `DEAD_PATH` digest `6f101c16…` — each present **1×**.
- K_A affinity-carrier term `lambdaAZero/(rhoM*(1 - I*frequency*tauA))` — present in the transcript.
- Harness SHA grok cited `9cdce889…52873c6470` — matches the working-tree file byte-for-byte (grok reviewed the
  exact committed harness).
- Harness line references (extractCarrier L121-132, operator L183-184, patches L66-80, validateSites L91-112,
  emitTriple L144-152, drift guard L212-213) — match the artifact on read.

Command + literal stdout for the above: this record's companion Bash lookups (digest counts 3/1/1; term present;
SHA match), run 2026-09-08.

## Verdict

SymPy harness: two-leg-closed (Leg 1 + Grok's prior static+SymPy-ablation pass, cleared earlier). WL harness:
**now two-leg-closed** (Leg 1 full-pair WL ablation + Leg 2 independent WL FORM ablation). Both legs SOUND, no
outstanding finding changes what is computed or may be claimed ⇒ **review-until-clear satisfied. S11c-b carrier
ablation harnesses ACCEPTED.** Committed as the reviewed baseline in the accompanying commit.

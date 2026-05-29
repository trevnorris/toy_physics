---
unit_id: 164
batch: V.1
created_at: 2026-05-28T00:00:00Z
findings_count: 1
stop_cold: null
applied: false
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 164

Apply each finding below in order. After applying, append an `## Applied: F<n>` block under that finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

If a finding's required change is ambiguous or unsafe to apply mechanically, append `## Blocked: F<n>` with a question instead — skip that finding, continue with the rest.

Do NOT introduce new features, refactors, or stylistic changes outside the file:line ranges named.

After editing, RUN the affected scripts (`python3 <path>` for SymPy, `math -script <path>` for Mathematica) and iterate until they exit 0 with all in-file checks passing. Getting the scripts to run cleanly is your job; the orchestrator independently re-runs afterward.

Do NOT touch paper.tex, notes/, or any prose documents. The red-team only modifies scripts.

## F1 — mathematica_transliteration

**Target:** `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage164_microscopic_log_channels_mathematica_audit.wl:32,116-122`

**Issue:**

The Mathematica script mirrors the SymPy script line-for-line: identical variable order, identical `ellLock`/`ksLock` healing-lock substitution choreography, identical hand-typed `firstHeal`/`secondHeal`/`deltaPerpExpected`, and the same copy-paste banner typo "STAGE 147" for a Stage-164 unit. Cross-engine agreement on the load-bearing healing-locked monomials and the `delta_perp` compression is therefore structurally guaranteed rather than independently confirmed. Add an independent route that obtains the two log-channel coefficient vectors directly from the explicit healing-locked product monomials (already asserted at L115-116) by a `Series`/`Coefficient` perturbation — machinery with no counterpart in the existing SymPy choreography — and reconcile it against the hand-written forms and the compressed `delta_perp`. Keep all existing assertions in place.

**Required change:**

Edit only the Mathematica script `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage164_microscopic_log_channels_mathematica_audit.wl`.

### Change A: fix the banner stage number

Current line 32:
```
banner["STAGE 147 — MICROSCOPIC LOG-IMBALANCE CHANNELS"];
```
Change to:
```
banner["STAGE 164 — MICROSCOPIC LOG-IMBALANCE CHANNELS"];
```
(This is a one-token label correction, not a stylistic refactor; the wrong number is the transliteration tell.)

### Change B: independent series-route derivation of the two log-channel coefficient vectors and delta_perp

Insert the following block immediately AFTER the existing healing-product asserts (after current line 116, `expectZero["healing second product exact formula", secondProdHeal - secondHealExpected];`) and BEFORE the `banner["delta_perp on the healing-locked branch"];` line (current line 118).

```mathematica
banner["Independent series route: log channels from explicit monomials"];

(* Take the already-asserted explicit healing-locked product monomials and
   extract each log-channel coefficient vector by a first-order multiplicative
   perturbation p -> p (1 + eps dlnP), reading off the O(eps) coefficient of
   the perturbed/unperturbed product RATIO. Using the ratio (not Log) avoids the
   negative-sign branch of the first product. This route uses Series + Coefficient,
   which the SymPy script does not use, so it is an independent derivation. *)
Clear[eps];
pertRule = {
  zq  -> zq  (1 + eps dlnZ),
  csw -> csw (1 + eps dlncsw),
  rhoW -> rhoW (1 + eps dlnrho),
  tm  -> tm  (1 + eps dlnTm),
  vw0 -> vw0 (1 + eps dlnv),
  a   -> a   (1 + eps dlna),
  lw  -> lw  (1 + eps dlnLw),
  cs  -> cs  (1 + eps dlncs)
};

firstRatio  = (firstHealExpected /. pertRule) / firstHealExpected;
secondRatio = (secondHealExpected /. pertRule) / secondHealExpected;

firstHealSeries  = Coefficient[Normal[Series[firstRatio,  {eps, 0, 1}]], eps];
secondHealSeries = Coefficient[Normal[Series[secondRatio, {eps, 0, 1}]], eps];

firstHealHand  = dlnZ + 3*dlncsw - dlnrho - dlnTm - dlnv - 2*dlna - 2*dlnLw;
secondHealHand = dlnZ + 2*dlncs + 3*dlncsw - dlnrho - 2*dlnv - 2*dlna - 3*dlnLw;

expectZero["first channel via series route",  firstHealSeries  - firstHealHand];
expectZero["second channel via series route", secondHealSeries - secondHealHand];

(* Build delta_perp from the series-derived coefficient vectors and reconcile
   against the compressed A_*/B_*/C_* target. *)
bWeight = 1/(4*Sqrt[1 + rstar^2]);
deltaPerpSeries = Expand[gstar*firstHealSeries + bWeight*secondHealSeries];
deltaPerpSeriesExpected = Expand[
  (gstar + bWeight)*(dlnZ - dlnrho)
  + 3*(gstar + bWeight)*dlncsw
  + 2*bWeight*dlncs
  - gstar*dlnTm
  - (gstar + 2*bWeight)*dlnv
  - 2*(gstar + bWeight)*dlna
  - (2*gstar + 3*bWeight)*dlnLw
];
expectZero["delta_perp via series route", deltaPerpSeries - deltaPerpSeriesExpected];

Clear[eps, pertRule, firstRatio, secondRatio, firstHealSeries, secondHealSeries,
      firstHealHand, secondHealHand, bWeight, deltaPerpSeries, deltaPerpSeriesExpected];
```

Notes for Codex:
- `firstHealExpected` and `secondHealExpected` are already defined at current lines 110-111 and asserted at 115-116, so this block depends only on prior in-scope definitions.
- The trailing `Clear[...]` keeps the later `delta_perp on the healing-locked branch` section (which redefines `b`, `firstHeal`, `secondHeal`, `deltaPerp`, `deltaPerpExpected`) in a clean scope; do not delete or alter that later section.
- Do NOT modify the existing `b = 1/(4*Sqrt[1 + rstar^2]);` at current line 120 — the new block uses a separately named `bWeight` to avoid clobbering it.
- Leave all other lines unchanged.

**Verification command:**

After Codex applies, the verifier will run `redteam exec-mathematica 164` and confirm:

1. The script still exits 0.
2. Three new `PASS` lines appear: `first channel via series route`, `second channel via series route`, `delta_perp via series route`.
3. The banner now reads `STAGE 164 — MICROSCOPIC LOG-IMBALANCE CHANNELS`.
4. All pre-existing `PASS` lines (parent-ratio identities, uniform products, healing products, `delta_perp compressed form`, and the five `expectApprox` numeric checks) still appear.

The captured `.txt` output file is regenerated by the same exec step.

## Applied: F1

files_changed:
  - /var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage164_microscopic_log_channels_mathematica_audit.wl

summary: Corrected the banner stage number from 147 to 164 and inserted an independent Series/Coefficient series-route derivation of the two log-channel coefficient vectors and delta_perp (with reconciliation asserts and a trailing Clear) between the healing-product asserts and the existing delta_perp section.

deviation: none

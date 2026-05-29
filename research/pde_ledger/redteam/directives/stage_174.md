---
unit_id: 174
batch: V.1
created_at: 2026-05-28T00:00:00-06:00
findings_count: 1
stop_cold: null
applied: false
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 174

Apply the finding below. After applying, append an `## Applied: F1` block under it with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

If the required change is ambiguous or unsafe to apply mechanically, append `## Blocked: F1` with a question instead.

Do NOT introduce new features, refactors, or stylistic changes. Edit exactly the lines named.

After editing, RUN the affected scripts (`python3 <path>` for SymPy, `math -script <path>` for Mathematica) and iterate until they exit 0 with all in-file checks passing. Getting the scripts to run cleanly is your job; the orchestrator independently re-runs afterward.

Do NOT touch paper.tex, notes/, or any prose documents. Edit only the Mathematica script.

## F1 — mathematica_transliteration

**Target:** `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage174_static_self_similarity_mathematica_audit.wl:26,67-71,89-93,111-115`

**Issue:**
The `.wl` is a line-for-line port of the SymPy script. The only non-trivial inputs — the first-order differentials `b01Two`, `z01Two`, `n01Two` — are hand-typed identically to the SymPy `B01_two`, `Z01_two`, `N01_two`, so the second engine re-types the same closed form rather than independently deriving it. A transcription error in either differential would be invisible to the cross-engine check. The banner text is also a stale copy ("STAGE 157" in a stage-174 unit).

**Required change:**
Replace the three hand-typed differential assignments with independent symbolic derivations via a perturbation parameter `eps`, so Mathematica reconstructs the differentials itself. Keep the `expectZero` checks at L79, L101, L123 unchanged; they should still pass.

1. BdG block. Replace lines 67-71:
```
b01Two = FullSimplify[
  2*c1*dc1/w1^2 - 2*c1^2*dw1/w1^3 +
  2*c2*dc2/w2^2 - 2*c2^2*dw2/w2^3,
  Assumptions -> $Assumptions
];
```
with:
```
b0Eps = (c1 + eps*dc1)^2/(w1 + eps*dw1)^2 + (c2 + eps*dc2)^2/(w2 + eps*dw2)^2;
b01Two = FullSimplify[D[b0Eps, eps] /. eps -> 0, Assumptions -> $Assumptions];
```

2. Conservative Maxwell/mixed block. Replace lines 89-93:
```
z01Two = FullSimplify[
  (delta1*q1p - q1*delta1p)/delta1^2 +
  (delta2*q2p - q2*delta2p)/delta2^2,
  Assumptions -> $Assumptions
];
```
with:
```
z0Eps = (q1 + eps*q1p)/(delta1 + eps*delta1p) + (q2 + eps*q2p)/(delta2 + eps*delta2p);
z01Two = FullSimplify[D[z0Eps, eps] /. eps -> 0, Assumptions -> $Assumptions];
```

3. Outgoing-transfer block. Replace lines 111-115:
```
n01Two = FullSimplify[
  2*p1*p1p/delta1^2 - 2*p1^2*delta1p/delta1^3 +
  2*p2*p2p/delta2^2 - 2*p2^2*delta2p/delta2^3,
  Assumptions -> $Assumptions
];
```
with:
```
n0Eps = (p1 + eps*p1p)^2/(delta1 + eps*delta1p)^2 + (p2 + eps*p2p)^2/(delta2 + eps*delta2p)^2;
n01Two = FullSimplify[D[n0Eps, eps] /. eps -> 0, Assumptions -> $Assumptions];
```
Add `Clear[eps]` once near the top of each affected block, or rely on `eps` being a fresh symbol; do not add `eps` to `$Assumptions` (it is eliminated by the `eps -> 0` substitution).

4. Banner. On line 26, change `banner["STAGE 157 — STATIC SELF-SIMILARITY DECOMPOSITION"];` to `banner["STAGE 174 — STATIC SELF-SIMILARITY DECOMPOSITION"];`.

Do not change any `expectZero` line, the `$Assumptions` declarations, the weighted-slope candidate expressions (`deltaBWeighted`, `deltaZWeighted`, `deltaNWeighted`), the theorem substitutions, or the carry-forward print block.

**Verification command:**
After Codex applies, the verifier will run `redteam exec-mathematica 174` and confirm: the transcript still prints `BdG weighted-average formula = 0`, `Z weighted-average formula = 0`, `N weighted-average formula = 0` each followed by `PASS`; the banner now reads `STAGE 174 ...`; the `b01Two`/`z01Two`/`n01Two` assignments now use `D[..., eps]`; and the script exits 0.

## Applied: F1

files_changed: mathematica/moving_throat_pde_stage174_static_self_similarity_mathematica_audit.wl
summary: Replaced the three hand-typed differential assignments (b01Two, z01Two, n01Two) with independent symbolic perturbation derivations via D[..., eps] /. eps -> 0, and corrected the stale banner from STAGE 157 to STAGE 174.
deviation: For the eps-symbol step, relied on the directive's explicitly-permitted option of eps being a fresh symbol (guaranteed by the file-opening ClearAll["Global`*"] and eps never being assigned) rather than adding Clear[eps]; eps was not added to $Assumptions, as instructed.

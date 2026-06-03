---
unit_id: 227
batch: VII.1
created_at: 2026-06-02T19:13:50Z
findings_count: 2
stop_cold: null
applied: true
applied_at: 2026-06-02T23:17:54Z
findings_applied: 2
findings_blocked: 0
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 227

Apply each non-paper_misalignment finding below in order. After applying, append an `## Applied: F<n>` block under that finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

The `paper_misalignment` finding(s) (F1) have been RESOLVED by the user (2026-06-02) — see the `## RESOLVED` block at the END of this directive and apply the authorized notes-only edit there as part of this fix loop. (Codex applies notes/*.md; Claude reviews.)

If F2's required change is ambiguous or unsafe to apply mechanically, append `## Blocked: F2` with a question instead.

Do NOT introduce new features, refactors, or stylistic changes beyond authoring the requested `.wl`. Do NOT touch paper.tex or the appendix. The ONLY authorized prose/notes edit is the notes-only change specified in the `## RESOLVED` block at the END of this directive (user-authorized 2026-06-02); apply exactly that and make no other notes/prose edits.

After editing, RUN the affected scripts (`math -script <path>` for Mathematica) and iterate until they exit 0 with all in-file checks passing.

---

## F1 — paper_misalignment

**Subtype:** notes_contradicts_script

**Paper side (notes — authoritative detail source; the stage card and appendix do NOT state the determinant's exact form, only "nonzero"):**
- `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage227_pure_transfer_load_factor_outgoing_rigidity_sieve_and_first_co_loading_no_go_sympy_audit.md:294` quote: `-\frac{19(-25+98\pi^2)(251+215\pi^2)(441\pi^2+4400)}{6(8670000+14894275\pi^2+2117682\pi^4)}`

**Script side:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage227_pure_transfer_load_factor_outgoing_rigidity_sieve_and_first_co_loading_no_go_sympy_audit.py:215` quote: `expected_det = -sp.Integer(19) * (-25 + 98 * sp.pi**2) * (200 + 147 * sp.pi**2) * (441 * sp.pi**2 + 4400) / (6 * (8670000 + 14894275 * sp.pi**2 + 2117682 * sp.pi**4))`

The two closed forms agree on every factor except the middle one: notes `(251+215π²)` vs script `(200+147π²)`. The difference is exactly `51 + 68π²` — the same notes-side drift signature seen in sibling stages 222/223 (notes = script + 68 in the π² coefficient). The script value is what SymPy computes (confirmed in the committed output `...stage227_..._sympy_audit.txt:24`) and is derived in-script from the reduced rows and the corridor nullspace, so the script is internally self-consistent.

## Resolve before fix_loop

The notes give the combined `i=h` rigidity determinant middle factor as `(251+215π²)`; the SymPy script asserts `(200+147π²)`, and the committed run reproduces the script value. Which is correct?

Possible directions (the user picks one):
- (a) Script is correct (expected — the committed output confirms SymPy emits `(200+147π²)`) → Claude corrects the notes line 294 middle factor to `(200+147π²)` (notes are prose, Claude-owned per file-ownership policy); no script change; script continues to exit 0. This matches the 222/223 notes-drift remediation pattern.
- (b) Notes are correct → the script's `expected_det` literal at line 215 must change to `(251+215π²)`; but this would only pass if SymPy actually produces it, which the current committed output shows it does not — so this direction also implies a genuine algebra error somewhere upstream and needs deeper review.
- (c) Both are wrong / third source → flag for deeper review.

F1 is RESOLVED by the user (see the ## RESOLVED block at the end) — apply the authorized notes edit. The co-loading no-go conclusion itself is unaffected: both `200+147π²` and `251+215π²` are strictly positive in π², so the determinant is nonzero either way.

## Applied: F1

- files_changed:
  - `notes/stages/moving_throat_pde_stage227_pure_transfer_load_factor_outgoing_rigidity_sieve_and_first_co_loading_no_go_sympy_audit.md`
- summary: Corrected the notes determinant middle factor from `251+215\pi^2` to `200+147\pi^2`.
- deviation: none

---

## F2 — missing_verification_script (missing_mathematica)

**Target:** `mathematica/moving_throat_pde_stage227_pure_transfer_load_factor_outgoing_rigidity_sieve_and_first_co_loading_no_go_mathematica_audit.wl`

**Issue:** No Mathematica audit exists for unit 227 (the stage card records "Mathematica audit: none yet"). Stage 227 is entirely closed-form symbolic linear algebra over exact rationals and `Pi` — scalar-parameter differentiation of rational functions, a 3×5 mixed matrix and its 2D null space, a 2×2 reduced determinant, one-dimensional null spaces of three reduced rows, and a Gram-projector operator norm. Independent Mathematica verification is clearly POSSIBLE, so the dual-engine rule requires a `.wl`. This stage carries the first exact same-charge co-loading no-go (anchor MTDC-T11.3) cited downstream, so a single-engine result is insufficient — and F1 shows precisely the kind of disagreement a second engine catches.

**Required change:**
Author a NEW file at the exact target path above. It must independently verify the claims below using native Mathematica primitives. Write the output capture to `mathematica/output/moving_throat_pde_stage227_pure_transfer_load_factor_outgoing_rigidity_sieve_and_first_co_loading_no_go_mathematica_audit.txt` consistent with the repo's existing `.wl` output convention.

**Anti-transliteration guard (MANDATORY — do NOT port the .py line-by-line):**
- Do NOT reproduce the SymPy choreography (`exp(eps*x)` dressing then `D[..., eps] /. eps -> 0`). Instead, derive the first-order microscopic slopes directly as `δln` linear forms in the primitive log-variables `{xLU, xLW, xLR, xOU, xOW}` from the closed-form definitions of `P`, `Delta`, `I`, `H` — i.e. take `D[Log[expr], {param}]` symbolically (or apply the log-derivative rule analytically) rather than dressing-and-differentiating an exponential ansatz.
- Use Mathematica-native exact arithmetic with `Pi` (not a float for `pi`), and use `FullSimplify`/`Together`/`Cancel` rather than mirroring SymPy's `simplify` call sequence.
- For the corridor basis and survivors, use `NullSpace` (exact) on the rational matrices directly; do not echo the `.py`'s `nullspace()[0]` index choreography — compute the null spaces independently and compare the resulting *directions* (up to sign/scale) against the expected numeric vectors.
- For the operator norm, build the Gram projector with `Transpose`/`Inverse` natively and use `Norm`/quadratic form on the `Xi_1` coefficient vector; compare the numeric value, not the symbolic intermediate steps.
- Pick a different decomposition where natural (e.g. compute the determinant via `Det` of the 2×2 reduced matrix obtained by restricting `i`,`h` rows to the corridor, and independently cross-check `det != 0` by confirming the 2×2 has full rank via `MatrixRank`).

**Claim manifest** (each must be an independently-derived, asserting check — use `FullSimplify[lhs-rhs]===0` or an `Abs[...] < tol` numeric check followed by `Print`+`Exit[1]` on failure; do NOT only `Print`):

- **M1 — pure-transfer theorem.** On the sample branch + Stage-243 pure-transfer corridor (`xK=xM=xLB=xV=0`), the same-charge scalar satisfies
  `Ξ_1 = N_01/N_0 = 2(P_01/P) − 2(Δ_01/Δ)`,
  where `N_0 = P^2/Δ^2`, `P = Ω_U^2 G_W + R G_U`, `Δ = Ω_U^2 Ω_W^2 − R^2`, and the first-order coefficients `N_01,P_01,Δ_01` are the `δ` (log-derivative) responses to the mixed-sector drifts. Derive `N_01` independently from `N_0`'s definition (not from `P_01/Δ_01`).

- **M2 — one-port factorization (symbolic identity).**
  `Λ = P/Δ = (G_W/Ω_W^2)·(1+I)/(1−H)` with `I = R G_U/(Ω_U^2 G_W)` and `H = R^2/(Ω_U^2 Ω_W^2)`. Verify as an exact symbolic identity in `{G_U,G_W,R,Ω_U,Ω_W}`.

- **M3 — microscopic slope law + sample coefficients.**
  `Ξ_1 = 2[ m + (I/(1+I)) i + (H/(1−H)) h ]` with the log-forms `m = xLW − 2 xOW`, `i = xLR + xLU − xLW − 2 xOU`, `h = 2 xLR − 2 xOU − 2 xOW`; and on the sample branch `I = 3/16`, `H = 25/(98 π^2)`, giving the specialization `Ξ_1 = 2 m + (6/19) i + 50/(98 π^2 − 25) h`. Derive `m,i,h` analytically from the `δln` of the definitions; assert `I===3/16` and `H===25/(98 Pi^2)`.

- **M4 — combined i=h rigidity determinant nonzero (the co-loading no-go).**
  The 2×2 determinant of the `(i,h)` rows reduced onto the pure-transfer corridor basis is nonzero. Assert `det != 0` (e.g. `FullSimplify[det] =!= 0` AND `MatrixRank[reduced2x2] === 2`). NOTE: do NOT hardcode the disputed middle factor `(200+147π²)` vs `(251+215π²)` from F1 — verify only that the determinant is nonzero (which is the paper's actual MTDC-T11.3 claim and is robust to F1). Optionally `Print` the symbolic determinant for the verifier to compare, but the asserting check is `det != 0`, not equality to a literal.

- **M5 — one-dimensional rigid survivors.** Each of `i=0`, `h=0`, `m=0` reduced onto the corridor leaves a one-dimensional null space; the unit survivor directions (up to sign) match
  `v_i ≈ (0.45280825, −0.29424612, −0.82815170, −0.04054866, 0.14458380)`,
  `v_h ≈ (0.66561963, −0.38941932, 0.46712837, 0.03609301, 0.43103536)`,
  `v_m ≈ (0.13386239, −0.10586713, −0.98242900, −0.05389175, −0.05293356)` (tol 1e-8).

- **M6 — same-charge gain scales.** `|Ξ_1(v_i)| ≈ 1.26576248`, `|Ξ_1(v_h)| ≈ 2.04509123`, `|Ξ_1(v_m)| ≈ 0.29342952` (tol 1e-8).

- **M7 — corridor norm and transported 10%-loss ceilings.** `σ_transfer ≈ 2.31561904386057` (tol 1e-11); with `budget_both = 0.367930328492646` and `budget_nonempty = 0.737619063660757`, the ceilings `(budget/σ)` reproduce: transfer `(0.15889070, 0.31854077)`, i `(0.29067881, 0.58274682)`, h `(0.17990900, 0.36067783)`, m `(1.25389678, 2.51378617)` (tol 1e-8).

Sample-branch constants (for reference, derive don't guess): `κ = 2√2/π`, `λ_B=1/2`, `λ_U=3/10`, `λ_W=2/5`, `λ_R=1/4`, `Ω_U=1`, `Ω_W=7/5`, `ϖ=2`, `M=1`, with `G_U=λ_U`, `G_W=κ λ_W`, `R=κ λ_R`. The Stage-243 corridor is the 2D null space of the mixed-sector `{D01,D21,D41}` first-order matrix; the `K` compatibility value is fixed by the Stage-242 shape-preservation relation `K = B0+Z0+3(1+B2+Z2)^2/(B4+Z4)` on the sample branch (see `...sympy_audit.py:120-121` for the relation only — re-derive, do not transliterate). The expected corridor basis (tol 1e-11) is
`t_1 ≈ (−4.359222794718, 3.107402039105, 18.703510605854, 1.0, 0.0)`,
`t_2 ≈ (1.909256655687, −1.163651238154, −0.482414494705, 0.0, 1.0)`.

**Verification command:**
After Codex applies, the verifier runs `redteam exec-mathematica 227` and confirms the new `.wl` exists at the exact target path, all checks M1–M7 appear and assert, and the script exits 0. The verifier also confirms the `.wl` is an independent derivation (different decomposition per the anti-transliteration guard), not a line-by-line port of the SymPy script.

## Applied: F2

- files_changed:
  - `mathematica/moving_throat_pde_stage227_pure_transfer_load_factor_outgoing_rigidity_sieve_and_first_co_loading_no_go_mathematica_audit.wl`
  - `mathematica/output/moving_throat_pde_stage227_pure_transfer_load_factor_outgoing_rigidity_sieve_and_first_co_loading_no_go_mathematica_audit.txt`
- summary: Added an independent Mathematica audit that verifies M1–M7 using exact log-variation, NullSpace, determinant/rank, survivor, gain, and ceiling checks.
- deviation: none

## RESOLVED — F1 paper_misalignment (USER-AUTHORIZED 2026-06-02)

Direction: correct the notes to match the verified SymPy script (authoritative; SymPy output `scripts/output/...stage227..._sympy_audit.txt`). Notes-only; Codex applies, Claude reviews. Do NOT change the script, paper.tex, or appendix.
- In `notes/stages/moving_throat_pde_stage227_..._sympy_audit.md`, the i=h combined-rigidity determinant factor reads `(251 + 215 π²)` (TeX: `251+215\pi^2`); correct it to the script/output value `(200 + 147 π²)` (`200+147\pi^2`) — a `+51 + 68π²` slip.
- Acceptance: the stale `251`/`215` determinant factor no longer appears; `200`/`147` present in that factor. Append `## Applied: F1`.

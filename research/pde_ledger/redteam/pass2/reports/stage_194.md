---
unit_id: 194
batch: V.3
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-06-09T00:00:00Z
verdict: findings
stop_cold: null
findings_count: 2
paper_alignment: aligned
scripts_checked:
  sympy: present
  mathematica: present
  engines_agree: true
  outputs_fresh: false
docs_read:
  paper_stage_tex: present
  notes_stage_files: [moving_throat_pde_stage194_outgoing_l2_fingerprint_and_deformation_algebra.md]
  paper_appendix: present
---

# Audit unit 194 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_194.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage194_outgoing_l2_fingerprint_and_deformation_algebra.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part05.tex` (rows at lines 119, 1263–1304, 1467)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage194_outgoing_l2_fingerprint_and_deformation_algebra_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage194_outgoing_l2_fingerprint_and_deformation_algebra_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage194_outgoing_l2_fingerprint_and_deformation_algebra_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage194_outgoing_l2_fingerprint_and_deformation_algebra_mathematica_audit.txt`

## What the paper claims

`\stagefield{Output}` (card line 15): "Derives the exact outgoing spherical DtN fingerprint and fixes \(\chi_Q=1\) on the canonical compact branch." The notes elaborate five distinct deliverables: (1) the exact outgoing spherical `l=2` DtN operator `Λ₂^out(z)=z d/dz ln h₂⁽¹⁾(z)` with the small-`z` fingerprint `-3 + z²/3 + z⁴/9 + i z⁵/9 - 2z⁶/27 - i z⁷/27 + …` and normalized branch `Ŷ₂^out = 1 + z²/9 + 4z⁴/81 + i z⁵/27 - 11z⁶/729 - i z⁷/243 + …`; (2) matching the retarded grouped-`P2` one-pole module (with `Ω_Q=3c_s/2a`, `σ_Q^can=9/(8Ω_Q⁵)=4a⁵/27c_s⁵`) to that fingerprint, fixing `\chi_Q=1`; (3) the isotropic DtN deformation family `Λ₂^def = S Λ₂^out(βz)+Σ₀+Σ₂z²+Σ₄z⁴+iΣ₅z⁵`, its canonical-even constraints (`Σ₂`, `Σ₄` forced) and the outgoing-normalization map `χ_Q = 3(Sβ⁵+9Σ₅)/(3S-Σ₀)`; (4) the localization statement that only `(β,Σ₀,Σ₅)` can move `χ_Q`; (5) the carry-forward invariant tuple `K̄₀=54Gc_s⁵/5a⁵c⁵`, `K̄₂=6Gc_s³/5a³c⁵`, `K̄₄=8Gc_s/15ac⁵`, `Γ̄₅=2G/5c⁵`. Card is `checkpoint: False`, not status-only — both engines required.

## What the script claims to verify

Both scripts verify all five deliverables. Section I derives the Hankel `l=2` fingerprint and its normalized form and checks the series against the boxed targets and the static slot `Λ(0)=-3`. Section II constructs/derives the retarded one-pole module, confirms `Ω_Q`, `σ_Q^can`, the low-frequency `Ŷ_Q^ret` form, that it equals the outgoing branch at `χ_Q=1`, and that odd-coefficient matching forces `χ_Q-1=0`. Section III builds the deformation `L₀..L₅`, the normalized compiler, solves the canonical-even constraints for `Σ₂,Σ₄`, and verifies the `χ_Q` and `χ_Q-1` deformation laws. Section IV reproduces the invariant tuple. All `expect_zero`/`expectZero` checks pass (per saved outputs).

## Paper ↔ script cross-check

| Paper deliverable | Script check | Status |
|---|---|---|
| (1) Λ₂^out / Ŷ₂^out fingerprint series | py:76–78, wl:91–93 (series vs boxed target + static slot) | match |
| (2) Ω_Q, σ_Q^can, χ_Q=1 from matching | py:112–125, wl:139–158 | match |
| (3) deformation L₀..L₅, Σ₂/Σ₄ constraints, χ_Q map | py:154–183, wl:192–228 | match |
| (4) only (β,Σ₀,Σ₅) move χ_Q | encoded in χ_Q-1 law py:183, wl:228 (depends only on β,Σ₀,Σ₅) | match |
| (5) invariant tuple K̄₀,K̄₂,K̄₄,Γ̄₅ | py:204–206, wl:248–251 | match |

`paper_alignment: aligned`. Every paper-side deliverable has a faithful, non-tautological script-side check on both engines.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored? |
|---|---|---|---|---|---|
| A1 | sympy | 76 | `expect_zero(Λ_out series − boxed target)` | claim 1 | yes |
| A2 | sympy | 77 | `expect_zero(Ŷ_out series − boxed target)` | claim 1 | yes |
| A3 | sympy | 78 | `expect_zero(Λ(0)+3)` | claim 1 (static slot) | yes |
| A4 | sympy | 112 | `expect_zero(σ_Q^can − 4a⁵/27c_s⁵)` | claim 2 | yes |
| A5 | sympy | 113 | `expect_zero(Yret series − expected)` | claim 2 | yes |
| A6 | sympy | 114–117 | `expect_zero(Yret\|_{χ=1} − Yout)` | claim 2 | yes |
| A7 | sympy | 125 | `expect_zero(odd-coeff match − (χ_Q−1))` | claim 2 | yes |
| A8 | sympy | 154 | `expect_zero(Ydef series − compiler)` | claim 3 | yes |
| A9 | sympy | 172–173 | `expect_zero(Σ₂/Σ₄ formula)` | claim 3 | yes |
| A10 | sympy | 182–183 | `expect_zero(χ_Q / χ_Q−1 law)` | claims 3,4 | yes |
| A11 | sympy | 204–206 | `expect_zero(K̄₂,K̄₄,Γ̄₅ − literals)` | claim 5 | yes |
| B1 | mathematica | 91–93 | `expectZero(series − target; slot+3)` | claim 1 | yes |
| B2 | mathematica | 139–141 | `expectZero(Solve-derived pole/σ − target)` | claim 2 | yes |
| B3 | mathematica | 142–146 | `expectZero(Yret − target; chi=1 match)` | claim 2 | yes |
| B4 | mathematica | 155–158 | `expectZero(odd-coeff match − (χ_Q−1))` | claim 2 | yes |
| B5 | mathematica | 192–196 | `expectZero(L₀..L₅, Ydef − compiler)` | claim 3 | yes |
| B6 | mathematica | 215–216 | `expectZero(Σ₂/Σ₄ from Solve − target)` | claim 3 | yes |
| B7 | mathematica | 227–228 | `expectZero(χ_Q / χ_Q−1 law)` | claims 3,4 | yes |
| B8 | mathematica | 248–251 | `expectZero(K̄₀.. from Solve − literals)` | claim 5 | yes |

No tautological rows. A4–A7/B2–B4 each compare a derived quantity against an *independently written literal* or a *Solve* result, so they can fail if the physics is wrong. A6/B3 cross-check two independently-constructed objects (retarded kernel vs Hankel branch) — substantive.

## Findings

### F1 — paper_misalignment (subtype: paper_missing_script_claim — card-text lag, paper-side)

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_194.tex:11`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage194_outgoing_l2_fingerprint_and_deformation_algebra_mathematica_audit.wl:1` (file exists and passes)

**What's wrong:** The card's `\stagefield{Verification}` (line 11) states: "SymPy audit: \StageFile{...sympy_audit.py}.  Mathematica audit: none yet." But a complete, passing Mathematica audit exists (`...mathematica_audit.wl`, saved output all `PASS`). The card text lags the actual artifact set.

**Why this matters:** A reader/auditor would conclude the second engine is absent and the stage single-engine, when in fact it is dual-engine and independently verified. This is the paper-cleanup class explicitly called out in the V.3 context.

**Required change:** Paper-side prose edit — route to user (Codex does not edit paper/). The card line should be updated to reference the existing `.wl`. NOT a Codex/script fix.

**Verification:** Card line 11 names the `.wl` audit instead of "none yet."

### F2 — stale_output

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage194_outgoing_l2_fingerprint_and_deformation_algebra_sympy_audit.txt` (mtime 2026-06-01 11:43:38)
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage194_outgoing_l2_fingerprint_and_deformation_algebra_sympy_audit.py` (mtime 2026-06-03 15:59:11)

**What's wrong:** The SymPy `.py` is newer than its saved `.txt` by ~2 days. The captured `.txt` banner reads `STAGE 177 — …` (lines 3, 143) and `STAGE 177 LEDGER`, whereas the current `.py` prints `STAGE 194` (lines 35, 208). So the committed output predates the 177→194 banner renumber and does not reflect the current script's banner. (The 177 == 194 − 17 offset confirms this is the known pre-renumber stale capture, not a content error: all numeric/symbolic *results* in the stale `.txt` still match the current script.) The Mathematica output is fresh (banner `STAGE 194`, line 3).

**Why this matters:** Informational. The result lines in the stale `.txt` are still correct, but the banner mismatch means the captured transcript is from a prior revision. A fresh re-run regenerates the correct `STAGE 194` banner.

**Required change:** Re-run `python3 scripts/moving_throat_pde_stage194_outgoing_l2_fingerprint_and_deformation_algebra_sympy_audit.py` and recommit the refreshed `.txt`. No source edit needed.

**Verification:** New `.txt` line 3 reads `STAGE 194 — …` and `.txt` mtime ≥ `.py` mtime.

## Independent-derivation check (Mathematica)

GENUINELY INDEPENDENT. Three corresponding sections show materially different routes, not transliteration:

1. **Hankel fingerprint (claim 1).** SymPy builds `h₂⁽¹⁾` from explicit trig closed forms — `j2 = (3/z³ − 1/z) sin z − 3 cos z/z²`, `y2 = −((3/z³ − 1/z) cos z + 3 sin z/z²)`, `h2 = j2 + I*y2` (`py:44–46`). Mathematica uses the **native special function** `hankel2 = FunctionExpand[SphericalHankelH1[2, x]]` (`wl:63`), which the output shows as `(I E^(I x)(-3 + 3I x + x²))/x³` — a wholly different representation reached without re-typing the trig. This is exactly the native-vs-closed-form independence signal flagged for the 105-family.

2. **Pole / σ extraction (claim 2).** SymPy **pins** the answers: `Omega_Q = 3 c_s/(2 a)` and `sigma_Q_can = 9/(8 Omega_Q⁵)` (`py:88–89`), then checks consistency against `4a⁵/27c_s⁵`. Mathematica instead **Solves the matching equations**: `poleCandidates = poleScale /. Solve[Coefficient[retardedSeriesRaw,freq,2] == Coefficient[outgoingOmegaSeries,freq,2], poleScale, Reals]` with `positiveRoot` branch selection (`wl:107–113`), and likewise `sigmaCandidates = … Solve[…freq,5…]` (`wl:115–121`). Pin-and-verify vs Solve-with-positive-root-selection is a genuinely different operation on the same physical premise.

3. **Invariant tuple (claim 5).** SymPy pins `Kbar0 = 54 G c_s⁵/(5a⁵c⁵)` and derives the rest by ratio (`py:191–194`). Mathematica **Solves** for `baseK0` from the `Γ̄₅` relation `9 baseK0/(32 poleFromEven⁵) == 2 G/(5 c⁵)` (`wl:233–238`), then derives the rest from the Solve-derived pole. Different anchor (Γ̄₅ inversion vs literal pin).

Section III (deformation algebra) is the most parallel — both assemble `L₀..L₅` and the normalized compiler whose *form* is fixed by the paper — but Mathematica extracts coefficients via `Series`+`Coefficient` and uses `Solve` for the even matching (`wl:198`), while SymPy uses `sp.series`+`sp.diff` and `sp.solve` (`py:156`). The shared structure here is dictated by the paper's definitions, not by echoing the other script's variable choreography; the operations differ. Verdict: **independent**, no `mathematica_transliteration` finding.

## Engine cross-check

Both engines emit identical results. Λ₂^out series, Ŷ₂^out series, `Ω_Q=3c_s/2a`, `σ_Q^can=4a⁵/27c_s⁵`, the `χ_Q=1` / `χ_Q-1` laws, the deformation `Σ₂/Σ₄` and `χ_Q` map, and the full invariant tuple all agree (SymPy `.txt` lines 34–140; Mathematica `.txt` lines 22–84). Mathematica additionally surfaces `chi_Q condition from odd coefficient = chiQ == 1` (`.txt:35`) consistent with SymPy's `chi_Q − 1 = 0`. No disagreement.

## Verdict justification

`findings` — two low-severity, non-math findings: a paper-side card-text lag ("Mathematica audit: none yet" while a passing `.wl` exists; F1, routes to user) and a stale SymPy `.txt` whose banner still reads STAGE 177 (F2, refresh-only). All physics holds: I attacked the Hankel series (verified `j₂/y₂` closed forms and the native `SphericalHankelH1` reduction agree), checked the static slot `Λ(0)=-3`, hand-verified the invariant tuple ratios (`K̄₂=6Gc_s³/5a³c⁵`, `K̄₄=8Gc_s/15ac⁵`, `Γ̄₅=2G/5c⁵` all reproduce exactly from `Ω_Q=3c_s/2a` and `K̄₀`), and confirmed the deformation `χ_Q` law depends only on `(β,Σ₀,Σ₅)`. No tautology, no hardcoded answer-against-itself, no symbol-domain error (`a,c_s,Ω_Q>0`, `3S−Σ₀≠0` justified). The Mathematica retrofit is a genuine independent re-derivation (native special function + Solve-based pole/σ/K̄₀ extraction), NOT a port. Paper alignment is exact on all five deliverables. Read paper card + notes + appendix before the scripts.

## Self-test notes

Variable-independence: the only `D[…]`/`diff` are `d/dz ln h₂⁽¹⁾` (h₂ depends on z — nonzero) and the `sp.diff(series, omega, 5)` coefficient extractions (the series genuinely carry `ω⁵` terms) — no identically-zero-derivative trap. Symmetry/parity: no unbounded-domain integrals in this stage (series-coefficient algebra only). Trivial-case: at `χ_Q=1` the retarded series collapses onto the outgoing branch (A6/B3 = 0, confirmed in both `.txt`s), and the `χ_Q−1` law vanishes when `β=1, Σ₀=0, Σ₅=0` — both consistent with the paper. Paths: F2 names `scripts/...sympy_audit.txt` (output dir) correctly; F1 is paper-side (no Codex path). No new `paper_misalignment` introduced by the F2 refresh.

## Value Reconciliation (pass-2 augmentation)

Authoritative record: script source (`.py`/`.wl`) + committed `.txt` outputs (SymPy `.txt` is stale-bannered but its *result* values are unchanged from current source; Mathematica `.txt` is fresh). Nothing was executed.

| value | source (py/wl + output line) | .tex/.md location | status |
|---|---|---|---|
| `Λ₂^out = -3 + z²/3 + z⁴/9 + i z⁵/9 - 2z⁶/27 - i z⁷/27` | py:53–60 / wl:70–77; sympy.txt:16/wl.txt:10 | notes:111–121 (boxed); appendix has Λ form (eq part05-Lambda-out) | MATCH |
| `Ŷ₂^out = 1 + z²/9 + 4z⁴/81 + i z⁵/27 - 11z⁶/729 - i z⁷/243` | py:61–68 / wl:78–85; both `.txt` | notes:136–146 (boxed); appendix:1280 (low orders) | MATCH |
| `Λ₂^out(0) = -3` (static slot) | py:78 / wl:93 | notes:124 | MATCH |
| `Ω_Q = 3 c_s/(2a)` | py:88 / wl:139; sympy.txt:35, wl.txt:22 | notes:72; card via Inputs | MATCH |
| `σ_Q^can = 9/(8Ω_Q⁵) = 4a⁵/27c_s⁵` | py:89,112 / wl:140–141; both `.txt` | notes:74 | MATCH |
| `χ_Q = 1` | py:115 / wl:144; sympy.txt:52, wl.txt:35 | card:15; notes:39,184; appendix:1285 | MATCH |
| `Σ₂ = -(3Sβ²-3S+Σ₀)/9` | py:165 / wl:209; sympy.txt:93 | notes:267 | MATCH |
| `Σ₄ = -(3Sβ⁴-3S+Σ₀)/27` | py:166 / wl:210; sympy.txt:98 | notes:269 | MATCH |
| `χ_Q = 3(Sβ⁵+9Σ₅)/(3S-Σ₀)` | py:177 / wl:223; sympy.txt:105, wl.txt:64 | notes:280; appendix:1299 | MATCH |
| `χ_Q-1 = (3S(β⁵-1)+Σ₀+27Σ₅)/(3S-Σ₀)` | py:178 / wl:224 | notes:287; appendix:1302 | MATCH |
| `K̄₀ = 54Gc_s⁵/(5a⁵c⁵)` | py:191 / wl:248; sympy.txt:116, wl.txt:73 | notes:309 | MATCH |
| `K̄₂ = 6Gc_s³/(5a³c⁵)` | py:192,204 / wl:249; both `.txt` | notes:311 | MATCH |
| `K̄₄ = 8Gc_s/(15ac⁵)` | py:193,205 / wl:250; both `.txt` | notes:313 | MATCH |
| `Γ̄₅ = 2G/(5c⁵)` | py:194,206 / wl:251; both `.txt` | notes:315 | MATCH |

INTERNAL (scaffolding, no prose expected): banner/ledger print strings, `L0/L2/L4/L5` intermediate slots, `expect_zero`/`expectZero` residuals (all 0), pass/fail flags, `positiveRoot` candidate lists, the ω⁶/ω⁷ higher-order coefficients used only to pin the series window.

reconciliation: complete; 14 deliverable values checked, 0 misaligned. The terse `.tex` card legitimately omits the intermediate symbolic results; all 14 live correctly in the notes (the natural carrier), several mirrored in the Part V appendix. No MISMATCH and no MISSING-DELIVERABLE — so no `value_mismatch` or `script_missing_paper_claim` finding is raised from reconciliation.

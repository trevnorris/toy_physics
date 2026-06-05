---
unit_id: 052
batch: III.2
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-06-05T00:00:00Z
verdict: findings
stop_cold: null
findings_count: 1
paper_alignment: aligned
scripts_checked:
  sympy: present
  mathematica: present
  engines_agree: true
  outputs_fresh: false
docs_read:
  paper_stage_tex: present
  notes_stage_files: [moving_throat_pde_stage052_nontwin_asymmetry_threshold.md]
  paper_appendix: present
---

# Audit unit 052 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_052.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage052_nontwin_asymmetry_threshold.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part03.tex` (row 82 + input/diagram lines)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage052_nontwin_asymmetry_threshold_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage052_nontwin_asymmetry_threshold_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage052_nontwin_asymmetry_threshold_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage052_nontwin_asymmetry_threshold_mathematica_audit.txt`

## What the paper claims

Stage 052 derives, in branch-product variables, the exact requirement for non-twin lowest-lane asymmetry once the symmetric twin lane fails. `\stagefield{Output}` reads verbatim: *"Regime split \eqref{eq:app-stage052-regimes}, excess formula \eqref{eq:app-stage052-Dzeta}, and resource factorization \eqref{eq:app-stage052-zeta0-general}."* The distinct deliverables are: (1) the closed form `zeta_req = (Pi_tr - C_mix)/[C_mix - eps(2 C_mix - Pi_tr)]` (boxed Eq. zeta-req); (2) its strict monotonicity `dzeta_req/dPi_tr = C_mix(1-eps)/[...]^2 > 0` (Eq. zeta-derivative); (3) the three-regime split keyed to `Pi_tr` vs `C_mix`/`2 C_mix` (boxed Eq. regimes); (4) the excess `Delta_zeta = zeta_req - 1 = (1-eps)(Pi_tr - 2 C_mix)/[C_mix - eps(2 C_mix - Pi_tr)]` (boxed Eq. Dzeta); (5) the resource factorization `zeta_0^phys = (K_W^eff/K_phi0^eff) Omega_0^2` (boxed Eq. zeta0-general). The notes additionally state the lowest-lane rescue thresholds (`Omega_0^2 >= zeta_req K_phi0/K_W`, `K_phi0 <= K_W Omega_0^2/zeta_req`), the twin diagnostic `zeta_0^twin = 1`, the equal-stiffness overlap threshold `Omega_0 >= sqrt(zeta_req)`, and the equal-overlap softening fraction `1 - 1/zeta_req = (1-eps)(Pi_tr - 2 C_mix)/(Pi_tr - C_mix)`. All deliverables are symbolic; there are no numeric constants.

## What the script claims to verify

Both scripts build `S_req = Pi_tr/C_mix` and `zeta_req = (S_req-1)/(1+eps(S_req-2))`, then assert: `zeta_req = 0` at `Pi_tr = C_mix`; `zeta_req = 1` at `Pi_tr = 2 C_mix` (the two regime-boundary anchors); `d zeta_req/dPi_tr` equals the paper's closed form (SymPy by direct `diff`, Mathematica additionally via an independent log-derivative path); `Delta_zeta` equals the paper's closed form; the resource factorization `zeta_0^phys = K_W Omega_0^2/K_phi0`; that `solve(zeta_phys = zeta_req)` for `Omega_0^2` and for `K_phi0` reproduces the hand-written rescue thresholds; the twin diagnostic `zeta_0^twin - 1 = 0`; and the equal-overlap softening fraction equals `(1-eps)(Pi_tr - 2 C_mix)/(Pi_tr - C_mix)` (Mathematica additionally via an independent `Together[1 - 1/zeta_req]` path). The regime-split inequalities themselves are printed as a ledger (SymPy lines 124-135) rather than asserted, but are anchored by the two boundary-value assertions.

## Paper ↔ script cross-check

| Paper deliverable | Script-side check | Status |
|---|---|---|
| (1) `zeta_req` closed form (Eq. zeta-req) | SymPy 45 / Math 42-44 build it; boundary anchors at 52-53 / 48-49 | match |
| (2) monotonicity `dzeta/dPi = …^2 > 0` (Eq. zeta-derivative) | SymPy 61 / Math 60-62 assert derivative = closed form; Math 56-59 independent path | match (form); `>0` sign is manifest from positive symbols, not separately asserted — partial on the strict-positivity word |
| (3) three-regime split (Eq. regimes) | SymPy 52-53 / Math 48-49 assert the two boundary values (`0` at `C_mix`, `1` at `2 C_mix`); inequalities printed 124-135 / ledger | match (boundaries anchored; inequalities are prose) |
| (4) excess `Delta_zeta` (Eq. Dzeta) | SymPy 69 / Math 68 | match |
| (5) resource factorization `zeta_0^phys` (Eq. zeta0-general) | SymPy 72 (defined+printed); 77/86-93 use it for thresholds; twin limit 105-106 | match |
| notes: rescue thresholds | SymPy 86-102 / Math 77-83 via `solve` | match |
| notes: twin diagnostic `zeta_0^twin=1` | SymPy 106 / Math 94 | match |
| notes: equal-overlap softening fraction | SymPy 122 / Math 91-98 | match |

`paper_alignment: aligned` — every paper-side deliverable has a faithful script-side check, in both engines.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 52 | `expect_zero(zeta_req|Pi_tr=C_mix)` | (3) lower boundary | yes |
| A2 | sympy | 53 | `expect_zero(zeta_req|Pi_tr=2C_mix - 1)` | (3) upper boundary | yes |
| A3 | sympy | 61 | `expect_zero(diff(zeta_req) - dz_expected)` | (2) | yes |
| A4 | sympy | 69 | `expect_zero(Delta_zeta - Delta_expected)` | (4) | yes |
| A5 | sympy | 89-93 | `solve(...)=Omega0_req_sq` | (5)/notes threshold | yes |
| A6 | sympy | 98-102 | `solve(...)=Kphi0_req` | (5)/notes threshold | yes |
| A7 | sympy | 106 | `expect_zero(zeta_twin - 1)` | notes twin diagnostic | partial (KW/KW-1; near-trivial but maps to claim) |
| A8 | sympy | 122 | `expect_zero(soft_frac - soft_expected)` | notes softening fraction | yes |
| B1 | math | 48-49 | `expectZero` boundary anchors | (3) | yes |
| B2 | math | 59 | `expectZero(dZdPi - dZdPiAlt)` indep path | (2) | yes |
| B3 | math | 62 | `expectZero(dZdPi - dZExpected)` | (2) | yes |
| B4 | math | 68 | `expectZero(deltaZeta - deltaExpected)` | (4) | yes |
| B5 | math | 78-83 | `solve` thresholds | (5)/notes | yes |
| B6 | math | 92 | `softFrac vs Together[1-1/zetaReq]` indep path | notes softening | yes |
| B7 | math | 94 | `expectZero(zetaTwin - 1)` | notes twin diagnostic | partial (trivial) |
| B8 | math | 98 | `expectZero(softFrac - softExpected)` | notes softening | yes |

## Findings

### F1 — stale_output

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage052_nontwin_asymmetry_threshold_sympy_audit.txt:3,77`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage052_nontwin_asymmetry_threshold_mathematica_audit.txt:3`
- (associated stale self-labels) `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage052_nontwin_asymmetry_threshold_sympy_audit.py:4-5`

**What's wrong:**
Both saved transcripts predate their scripts (output mtime 2026-05-22; script mtime 2026-06-03), and the staleness is content-bearing for the stage-number banners. The committed SymPy `.txt` prints `STAGE 35 — EXACT BRANCH-PRODUCT REGIME CLASSIFICATION` (line 3) and `STAGE 35 THEOREM LEDGER` (line 77), whereas the current `.py` banners (lines 42, 124) say `STAGE 52`. The committed Mathematica `.txt` prints `STAGE 035 — NON-TWIN ASYMMETRY THRESHOLD` (line 3) while the current `.wl` (line 32) says `STAGE 052`; the Mathematica transcript's closing line 36 ("Stage 052 Mathematica audit passed.") is already current, so that file is mixed-stale. Separately, the SymPy docstring carries a stale self-label: line 4 `moving_throat_pde_stage35_nontwin_asymmetry_threshold_sympy_audit.py` (filename label says 35, while prose line 6 says "Stage 52"). The *mathematical* content of both transcripts otherwise matches what the current scripts produce (every printed expression was traced and agrees).

**Why this matters:**
Only the displayed stage label is wrong; no numeric/symbolic result is affected. This is the known low-severity numbering/stale-banner class. It is informational and non-blocking; a fresh re-run refreshes the transcripts.

**Required change:**
Re-run both scripts to refresh the committed `.txt` (the orchestrator's independent exec re-run does this). For the stale self-label in the SymPy docstring, fix `stage35`→`stage052` on line 4 (and "Stage 52"→"Stage 052" on lines 6, plus "STAGE 52" banners on lines 42, 124 to the zero-padded canonical form if the orchestrator's label policy requires three-digit padding). Per the standing policy these self-labels are resolved by the orchestrator's label scope, not blocking.

**Verification:**
After re-run, `.txt` line 3 of each transcript should read `STAGE 052` (or the canonical banner), and no math line changes.

## Independent-derivation check (Mathematica)

The `.wl` is not a transliteration. It derives `zetaReq` by `Solve[]`-ing the lowest-lane support equation (lines 42-44) rather than copying SymPy's `(Sreq-1)/(1+eps(Sreq-2))` literal, and it adds two genuinely independent verification paths that have no SymPy counterpart: (i) logarithmic differentiation of the rational form to cross-check `dZdPi` (lines 53-59), and (ii) `Together[1 - 1/zetaReq]` to cross-check the softening fraction (lines 90, 92). The threshold equations are obtained by `Solve[]` on the physical equation `(kW omega0^2/kPhi0) - zetaReq == 0`. This is an independent re-derivation, not an echo of the SymPy algebra.

## Engine cross-check

The two engines agree on every shared result. SymPy `zeta_req` output `(C_mix - Pi_tr)/(-C_mix + eps(2 C_mix - Pi_tr))` equals (multiply num/denom by -1) the Mathematica `-((cMix - piTr)/(cMix - 2 cMix eps + eps piTr))` and the paper's `(Pi_tr - C_mix)/[C_mix - eps(2 C_mix - Pi_tr)]`. Both report `Delta_zeta`, both threshold solves, the twin diagnostic, and the softening fraction as residual-0 / PASS. No sign, factor, or domain disagreement.

## Verdict justification

`findings`, with the single finding being a low-severity, non-blocking `stale_output` (banner stage-number labels in both committed transcripts predate a 35→52 banner edit, plus a stale `stage35` self-label in the SymPy docstring). I attacked: (a) tautology — the boundary anchors, derivative check, excess check, threshold solves, and softening fraction are all `expr_A - expr_B == 0` where one side is computed (`diff`/`solve`/`Together`) and the other independently typed, so none is guaranteed by construction; the only near-trivial check is `zeta_twin - 1` (`KW/KW - 1`), which still maps to a real paper claim and is harmless; (b) symbol domains — SymPy declares `eps` only `positive`, omitting `eps<1`, but every SymPy assertion is a pure rational-function identity valid for all `eps`, so the missing upper bound invalidates nothing (Mathematica correctly carries `0<eps<1`); (c) sign claim — the paper's `>0` monotonicity is manifest from the positive closed form rather than separately asserted (partial, not a defect); (d) paper alignment — every `\stagefield{Output}` deliverable and every notes-level threshold has a faithful, anchored, dual-engine check. The math is sound and aligned; only the cosmetic transcript/label staleness is flagged.

## Value Reconciliation (pass-2 augmentation)

All stage-052 deliverables are symbolic closed forms; there are no computed numeric constants. Each emitted symbolic deliverable is reconciled against the `.tex` card and `.md` notes below.

| value | source (py/wl + output line) | .tex/.md location | status |
|---|---|---|---|
| `zeta_req = (Pi_tr - C_mix)/[C_mix - eps(2 C_mix - Pi_tr)]` | py 45 / wl 42-44; sympy out 10-13, math out 6 | tex 16-20 (Eq. zeta-req); md 34,58 | MATCH |
| boundary `zeta_req=0` at `Pi_tr=C_mix` | py 52; out 14 | md 71 | MATCH |
| boundary `zeta_req=1` at `Pi_tr=2 C_mix` | py 53; out 15 | md 73 | MATCH |
| `dzeta_req/dPi_tr = C_mix(1-eps)/[C_mix - eps(2 C_mix - Pi_tr)]^2` | py 60-61 / wl 60-62; out 17-21, math out 13 | tex 22-26 (Eq. zeta-derivative); md 64-65 | MATCH |
| regime split (`<=C_mix`; `C_mix<…<=2C_mix`; `>2C_mix`) | py 126-128 / printed; out 79-82 | tex 28-36 (Eq. regimes); md 77-79 | MATCH |
| `Delta_zeta = (1-eps)(Pi_tr - 2 C_mix)/[C_mix - eps(2 C_mix - Pi_tr)]` | py 68-69 / wl 66-68; out 27-30, math out 16 | tex 38-44 (Eq. Dzeta); md 93-95 | MATCH |
| `zeta_0^phys = (K_W/K_phi0) Omega_0^2` | py 72 / wl 70; out 36-39, math out 19 | tex 46-49 (Eq. zeta0-general); md 113 | MATCH |
| `Omega_0,req^2 = zeta_req K_phi0/K_W` | py 77 / wl 71; out 42-45, math out 20 | md 125-126,205 | MATCH |
| `K_phi0^req = K_W Omega_0^2/zeta_req` | py 78 / wl 72; out 47-51, math out 21 | md 132-133,207 | MATCH |
| `zeta_0^twin = 1` | py 105-106 / wl 85,94; out 58, math out 28 | md 153 | MATCH |
| `Omega_0 >= sqrt(zeta_req)` (equal stiffness) | py 108 / wl 86; out 60-63, math out 30 | md 163 | MATCH |
| equal-overlap softening fraction `(1-eps)(Pi_tr-2C_mix)/(Pi_tr-C_mix)` | py 121-122 / wl 91-98; out 70-73, math out 32 | md 177-179 | MATCH |

Internal scaffolding (no prose expected): `S_req = Pi_tr/C_mix` (intermediate ratio), `expect_zero`/`expectZero` residual-0 prints, PASS flags, the `sol_*`/`*Sol` solve handles, `dZdPiAlt`/`softFracDerived` cross-check intermediates.

reconciliation: complete; 12 values checked, 0 misaligned

## Self-test notes

Variable independence: every `diff`/`D[...,piTr]` is taken of `zeta_req`, which genuinely depends on `piTr`, so no derivative is identically zero. Symmetry/parity: no integrals over unbounded domains in this stage. Trivial-case pre-check: substituting `Pi_tr=C_mix` gives `zeta_req=0` and `Pi_tr=2C_mix` gives `zeta_req=1` (boundary asserts reduce to 0); `zeta_twin` reduces to `KW/KW=1` (so its assert is near-trivial but true and on-claim). Domains: SymPy's omission of `eps<1` is harmless because all its checks are pure identities; Mathematica carries `0<eps<1`. No script-side math fix is prescribed, so no paper round-trip risk — the only finding is cosmetic stale_output.

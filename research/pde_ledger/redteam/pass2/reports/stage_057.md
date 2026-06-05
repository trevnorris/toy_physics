---
unit_id: 057
batch: III.2
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-06-05T00:00:00Z
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
  notes_stage_files: [/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage057_physical_parameter_map.md]
  paper_appendix: present
---

# Audit unit 057 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_057.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage057_physical_parameter_map.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part03.tex` (row 92; `\input` at line 232)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage057_physical_parameter_map_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage057_physical_parameter_map_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage057_physical_parameter_map_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage057_physical_parameter_map_mathematica_audit.txt`

## What the paper claims

Stage 057 rewrites the abstract lowest-lane support family entirely in physical
throat-operator variables `(Pe, kappa, eta)`. The card's `\stagefield{Output}` is
"The physical support map \eqref{eq:app-stage057-zeta}", namely the boxed
`zeta_0(Pe,eta;kappa) = Omega_Pe^2 (kappa+pi^2/4)/(kappa+y(eta)^2)` with
`y(eta) tan y(eta) = eta` (eq:app-stage057-zeta), plus its zero-bias value
`zeta_0(0,eta;kappa) = (kappa+pi^2/4)/(kappa+y^2)` (eq:app-stage057-zero-bias).
The notes add a fuller enumeration: (1) `x = pi^2/(kappa+pi^2/4)`;
(2) softening `A_K = (kappa+pi^2/4)/(kappa+y^2)`; (3) the physical family
`zeta_0^(Pe+R)`; (4) monotonicity (increasing in `Pe`, decreasing in `eta` and
`kappa`); (5) ceiling `zeta_max(kappa) = (pi^2/4)(kappa+pi^2/4)/kappa`;
(6) stiffness ceiling `kappa_max = pi^4/[4(4 zeta_req - pi^2)]`; and (7) three
threshold formulas `Omega_req^2`, `y_req^2`, `kappa_req`. The appendix row 92
summarizes it as "Physical support ratio `zeta_0(Pe,eta;kappa)`" with status
ExactClosure.

## What the script claims to verify

Both engines: (a) reduce the abstract `A_K`-in-`x` form under `x=pi^2/(kappa+pi^2/4)`
to `(kappa+pi^2/4)/(kappa+y^2)` (the `.wl` ALSO re-derives `A_K` independently from
the physical stiffnesses `K_W^eff`, `K_phi0^eff`); (b) assert the full family
`zeta_phys = Omega_Pe^2 (kappa+pi^2/4)/(kappa+y^2)`; (c) assert the closed forms of
`partial_kappa zeta` and `partial_y zeta`, plus numerical sign sweeps establishing
`partial_kappa < 0` (5 y-points) and `partial_Pe > 0` (6 Pe-points); (d) assert
`zeta_max = (pi^2/4)(kappa+pi^2/4)/kappa` equals the iterated limit
`Pe->inf, y->0+`; (e) assert `kappa_max = pi^4/[4(4 zeta_req - pi^2)]` from solving
`zeta_req = zeta_max`; and (f) assert the threshold formulas `kappa_req`, `y_req^2`
against independent `solve`s. `Omega_req^2` is printed only (definitional).

## Paper ↔ script cross-check

| Paper/notes deliverable | Script-side check | Status |
|---|---|---|
| (1) `x = pi^2/(kappa+pi^2/4)` | `x_sub` printed (py:47 / wl:63), used in A_K | match (printed, used as input to asserted A_K) |
| (2) `A_K = (kappa+pi^2/4)/(kappa+y^2)` | `expect_zero` py:49 / wl:54 (phys) + wl:66 | match |
| Output box `zeta_0 = Omega_Pe^2 A_K` | `expect_zero` py:56 / wl:67 | match |
| zero-bias `zeta_0(0,..)` = A_K | implied (Omega_Pe(0)→1); not separately asserted | partial (not separately checked; Omega_Pe→1 is a Stage 056 carry-forward) |
| (4a) increasing in `Pe` | `partial_Pe` sign sweep py:87-91 / wl:99-107 | partial (6-point numeric sweep, no symbolic proof; honestly flagged) |
| (4b) decreasing in `eta` | `partial_y identity` py:72 + sweep-free; `dy/deta>0` not in-script | partial (verifies `partial_y<0`; eta-direction via external IFT sign) |
| (4c) decreasing in `kappa` | `partial_kappa identity` py:66 + sign sweep py:77-81 | match (identity gives `Omega^2(y^2-pi^2/4)/(...)^2 < 0` on 0<y<pi/2) |
| (5) `zeta_max(kappa)` | `expect_zero` vs iterated limit py:96 / wl:117 | match |
| (6) `kappa_max(zeta_req)` | `expect_zero` py:104 / wl:119 | match |
| (7) `Omega_req^2` | printed only py:118 / wl:127 | partial (definitional restatement; not asserted — acceptable) |
| (7) `y_req^2` | `y_req identity` vs independent solve py:129 / wl:140 | match |
| (7) `kappa_req` | `kappa_req identity` vs solve py:121 / wl:135 | match |

`paper_alignment: aligned` — every load-bearing deliverable is exercised by a
non-tautological assertion; the `partial` rows are either definitional prints or
genuinely-true carry-forwards honestly annotated in the script.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | mathematica | 54 | `expectZero[A_K(phys) - (kappa+pi^2/4)/(kappa+y^2)]` | deliverable 2 (independent route) | yes |
| A2 | sympy | 49 | `expect_zero(A_K - (kappa+pi^2/4)/(kappa+y^2))` | deliverable 2 | yes |
| A3 | both | py56 / wl67 | `expect_zero(zeta_phys - Omega^2 A_K)` | Output box | yes |
| A4 | both | py66 / wl76 | `partial_kappa identity == 0` | deliverable 4c | yes |
| A5 | both | py72 / wl80 | `partial_y identity == 0` | deliverable 4b | yes |
| A6 | both | py77-81 / wl87-94 | `partial_kappa < 0` sign sweep | deliverable 4c | yes (corroborates A4) |
| A7 | both | py87-91 / wl99-107 | `partial_Pe > 0` sign sweep | deliverable 4a | partial (numeric only) |
| A8 | both | py96 / wl117 | `zeta_max - iterated limit == 0` | deliverable 5 | yes |
| A9 | both | py104 / wl119 | `kappa_max identity == 0` | deliverable 6 | yes |
| A10 | mathematica | 120 | `zeta_max(kappa_max) - zeta_req == 0` | deliverable 6 (round-trip) | yes |
| A11 | both | py121 / wl135 | `kappa_req identity == 0` (vs solve) | deliverable 7 | yes |
| A12 | mathematica | 131 | `kappa_req defining equation == 0` | deliverable 7 | yes |
| A13 | both | py129 / wl140 | `y_req identity == 0` (vs solve) | deliverable 7 | yes |

## Findings

### F1 — stale_output

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage057_physical_parameter_map_sympy_audit.txt` (mtime 2026-05-26 03:06:16)
- `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage057_physical_parameter_map_mathematica_audit.txt` (mtime 2026-05-26 03:06:25)

**What's wrong:**
Both saved transcripts are older than their scripts (both scripts mtime
2026-06-03 15:59:11). The content also lags the current source: the SymPy output
banner reads `STAGE 40 — ...` (txt:3) while the current `.py` banner emits
`STAGE 57` (py:33); the Mathematica output banner reads `STAGE 040 — ...` (txt:3)
while the current `.wl` banner emits `STAGE 057` (wl:32). The mathematica output's
closing line already says "Stage 057 ... audit passed." (txt:44) but its banner is
stale, confirming the transcript predates the banner edit. The numeric/symbolic
result lines themselves still match what the current scripts would produce (I
hand-verified the A_K reduction, the `zeta_max` and `kappa_max` closed forms, and
the partial-derivative forms), so this is a freshness/label discrepancy, not a
result discrepancy.

**Why this matters:**
The committed transcript is the authoritative record for downstream readers; a
stale banner mislabels which stage the run belongs to.

**Required change:**
Re-run both scripts to refresh the committed `.txt` transcripts after F2 lands.

**Verification:**
After F2's self-label fixes, the re-run SymPy transcript banner should read
`STAGE 057 — ...` (or `STAGE 57`, matching whatever the corrected `.py` banner
emits) and the Mathematica transcript banner `STAGE 057 — ...`; both `.txt`
mtimes must be newer than their script mtimes.

### F2 — stale_output (numbering self-labels in scripts)

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage057_physical_parameter_map_sympy_audit.py:1` (docstring `"Stage 40 SymPy audit."`)
- `.../moving_throat_pde_stage057_physical_parameter_map_sympy_audit.py:33` (banner `"STAGE 57 — ..."` — un-padded)
- `.../moving_throat_pde_stage057_physical_parameter_map_sympy_audit.py:134` (`"Stage 40 audit passed."`)

**What's wrong:**
The SymPy script carries stale self-labels naming the wrong stage. The module
docstring (py:1) says "Stage 40 SymPy audit"; the closing print (py:134) says
"Stage 40 audit passed." — both should be 057. The banner (py:33) says
`"STAGE 57"`, which is the right number but un-zero-padded relative to the `.wl`
banner `"STAGE 057 — ..."` (wl:32). The Mathematica `.wl` source self-labels are
already correct (banner wl:32 `STAGE 057`, closing wl:143 `Stage 057`). This is the
known low-severity `stale_output` numbering-label class (incomplete EM-extension
renumber, +17 chain: 40→057).

**Why this matters:**
Non-blocking, but the docstring/closing-line mislabel the stage and the banner is
inconsistently padded; they propagate into the refreshed transcript (F1).

**Required change:**
- py:1: `"Moving-Throat PDE — Stage 40 SymPy audit."` → `Stage 057`.
- py:33: `banner("STAGE 57 — PHYSICAL ...")` → `STAGE 057`.
- py:134: `print("\nStage 40 audit passed.")` → `Stage 057 audit passed.`

**Verification:**
After edit, `grep -n "Stage 40\|STAGE 57\b\|Stage 057\|STAGE 057"` on the `.py`
shows only `057`-form labels; the refreshed SymPy `.txt` banner reads
`STAGE 057 — ...` and closing line `Stage 057 audit passed.`

## Independent-derivation check (Mathematica)

The `.wl` is NOT a transliteration. For the central `A_K` result the two engines
take genuinely different routes: SymPy reduces the abstract `1/(1 - x/4 + x y^2/pi^2)`
form under `x = pi^2/(kappa+pi^2/4)` (py:43-52), whereas Mathematica first
*independently* derives `A_K = K_W^eff/K_phi0^eff = (kappa+pi^2/4)/(kappa+y^2)`
from the physical operator stiffnesses `KW = KX + Pi^2 TX/(4 LL^2)`,
`Kphi0 = KX + TX y^2/LL^2` with `KX -> kappa TX/LL^2` (wl:47-56), and only THEN also
reproduces the x-substitution route (wl:58-66) as a redundant cross-check. The
limit, solve, and threshold blocks are structurally parallel between engines (same
target identities), but that is expected for a results-matching audit and the
load-bearing `A_K` derivation is independent, so no `mathematica_transliteration`
finding.

## Engine cross-check

The engines agree. SymPy `zeta_max = pi**2/4 + pi**4/(16*kappa)` (sympy txt:16)
equals Mathematica `(4 kappa Pi^2 + Pi^4)/(16 kappa)` (math txt:22) — same closed
form. `kappa_max`: SymPy `pi**4/(4*(4*zeta_req - pi**2))` (txt:18) =
Mathematica `Pi^4/(-4 Pi^2 + 16 zetaReq)` (txt:23). `partial_kappa`,
`partial_y`, `Omega_req^2`, `y_req^2`, `kappa_req` all match (the Mathematica
`kappa_req` carries a stripped `ConditionalExpression` wrapper, math txt:36, whose
condition is exactly the existence window for a positive root; `expectZero` strips
it per the documented idiom). Every `expectZero`/`expect_zero` reports residual `0`
and the `.wl` prints PASS on each (math txt:5-42). Both transcripts end with
"audit passed".

## Verdict justification

`verdict: findings` driven solely by two low-severity `stale_output` items
(F1 stale transcripts, F2 stale SymPy self-labels). The substantive physics holds
up under attack: I checked the A_K reduction by hand (denominator collapses to
`4(kappa+y^2)/(4kappa+pi^2)`, inverse gives `(kappa+pi^2/4)/(kappa+y^2)` — exactly
the printed form); the `partial_kappa` identity makes the sign claim a genuine
consequence (`y^2 - pi^2/4 < 0` on `0<y<pi/2`), not a hardcoded assertion; the
`kappa_max` and `kappa_req`/`y_req` identities are compared against independent
`solve`s, so they are non-tautological; and the `zeta_max` check is against an
iterated symbolic limit, not a pre-baked literal. The two `partial` rows I did NOT
upgrade to findings: `partial_Pe > 0` rests on a 6-point numeric sweep (the script
honestly annotates this as a Stage-056 carry-forward whose sign was not upstream-
proven), and `eta`-monotonicity is mediated by the standard `dy/deta>0` IFT sign
not re-derived in-script — both are acceptable scaffolding for a placement-map
stage and are openly documented in the comments, not silent gaps.
`material_change: false`. No `stop_cold`. I read the paper card, the notes, and the
appendix row before the scripts, and the script's verified claim matches the card's
`\stagefield{Output}`.

## Self-test notes

Variable independence: `zeta_phys` depends on `Pe, kappa, y` (via `Omega_Pe(Pe)`
and the `A_K` factor), so `diff(.,kappa)`, `diff(.,y)`, `diff(.,Pe)` are all
genuinely nonzero — no identically-zero-derivative trap. Parity/integral traps:
none — there are no integrals over unbounded domains here. Trivial-case check: at
`Pe=1, kappa=1, y=pi/4`, `partial_kappa zeta = Omega^2(pi^2/16 - pi^2/4)/(1+pi^2/16)^2
< 0` (numerator negative) ✓, consistent with the sweep PASS. Paper round-trip: the
F2 self-label fix is label-only and introduces no new constant or paper
disagreement.

## Value Reconciliation (pass-2 augmentation)

reconciliation: complete; 7 values checked, 0 misaligned

| value | source (py/wl + output line) | .tex/.md location | status |
|---|---|---|---|
| `x(kappa) = pi^2/(kappa+pi^2/4)` (printed `4 pi^2/(4 kappa+pi^2)`) | py:47 / wl:63; sympy txt:5, math txt:7 | notes:78-79 (`x = pi^2/(kappa+pi^2/4)`); notes:294 | MATCH |
| `A_K(eta;kappa) = (kappa+pi^2/4)/(kappa+y^2)` | py:48,49 / wl:64,66; sympy txt:6, math txt:8 | tex:26 (zero-bias eq); notes:102 | MATCH |
| `zeta_0^(Pe+R) = Omega_Pe^2 (kappa+pi^2/4)/(kappa+y^2)` | py:55,56 / wl:65,67; sympy txt:8, math txt:9 | tex:17-21 (boxed Output); notes:112-113,300-301 | MATCH |
| `partial_kappa zeta = Omega^2(y^2-pi^2/4)/(kappa+y^2)^2` | py:64,66 / wl:74,76; sympy txt:10, math txt:14 | notes:166-167 | MATCH |
| `zeta_max(kappa) = (pi^2/4)(kappa+pi^2/4)/kappa` | py:94 / wl:109; sympy txt:16, math txt:22 | tex(notes):202,306; notes:42 | MATCH |
| `kappa_max(zeta_req) = pi^4/[4(4 zeta_req - pi^2)]` | py:103 / wl:111; sympy txt:18, math txt:23 | notes:49,227 | MATCH |
| `kappa_req = (Omega^2 pi^2/4 - zeta_req y^2)/(zeta_req - Omega^2)` | py:120 / wl:129; sympy txt:22, math txt:36 | notes:278-279 | MATCH |

Also reconciled (deliverable-level, prose-carried): `Omega_req^2 = zeta_req(kappa+y^2)/(kappa+pi^2/4)`
(py:118 / sympy txt:20 ↔ notes:246) MATCH; `y_req^2 = (Omega^2/zeta_req)(kappa+pi^2/4) - kappa`
(py:119 / sympy txt:21 ↔ notes:264) MATCH; `partial_y zeta` closed form
(sympy txt:11 ↔ notes:158 sign claim) MATCH.

INTERNAL (scaffolding, no finding): `Omega_Pe` intermediate form,
`partial_kappa/partial_Pe` sign-sweep PASS flags, the `kappa_req` and `y_req`
`solve`-vs-closed-form residual-0 cross-checks, and the iterated-limit residual.

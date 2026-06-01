---
unit_id: 191
batch: V.3
verifier_model: claude-opus-4-8[1m]
verify_date: 2026-06-01T12:00:00-06:00
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 2
findings_total: 2
material_change: false
---

# Verification — unit 191

## Per-finding outcomes

### F1 — paper_misalignment (stale stage-number banner labels)

**Classification:** resolved

**What changed:**
The captured diff (`redteam/exec_logs/stage_191_diff.patch`) shows exactly two single-line edits to `scripts/moving_throat_pde_stage191_minimal_pde_data_packet_sympy_audit.py`:
- line 65: `banner("STAGE 174 — MINIMAL PDE DATA PACKET AND THE EXACT HOME-STRETCH THEOREM")` → `banner("STAGE 191 — ...")`
- line 284: `banner("STAGE 174 LEDGER")` → `banner("STAGE 191 LEDGER")`

No other hunk appears in the diff. The regenerated SymPy output header (`scripts/output/...sympy_audit.txt:3`) now reads `STAGE 191 — MINIMAL PDE DATA PACKET AND THE EXACT HOME-STRETCH THEOREM` and the section-V banner (`:333`) reads `STAGE 191 LEDGER`.

**Assessment:**
String-only as authorized under the settled canonical-stage-number convention (direction (a) in the directive). Only the two `print` banner string literals changed; no symbol, assertion, constant, series order, substitution, or `expect_zero` argument was touched. The change cannot affect any residual, and indeed every SymPy `expect_zero` still prints `= 0` / a zero matrix in the refreshed output. The notes-prose legacy `242`/`174` references were correctly left untouched (out of red-team scope; logged to PAPER_CLEANUP_TRACKER). Resolved.

### F2 — missing_mathematica (independent second-engine route)

**Classification:** resolved

**What changed:**
New file `mathematica/moving_throat_pde_stage191_minimal_pde_data_packet_mathematica_audit.wl` (+ regenerated `output/...mathematica_audit.txt`). It re-verifies every load-bearing SymPy assertion and asserts cross-engine agreement (each Mathematica-native result is subtracted from the SymPy closed form, expecting 0).

**Assessment — genuine independent route (NOT a transliteration):**
The `.wl` uses Mathematica-native primitives through a different decomposition than the `.py`, not a line-by-line port:
- **Taylor compilers (I):** SymPy does `sp.series(...).removeO()` minus a hand-built `1 + nu2 ω² + nu4 ω⁴` polynomial. The `.wl` instead *extracts individual coefficients* with `Coefficient[Normal[Series[...]], w, k]` AND independently cross-checks each via the Taylor-derivative route `D[singleResponse,{w,2}]/2 /. w->0` and `D[...,{w,4}]/24` (lines 190–202, 207–208). Two genuinely distinct computations agree before either is compared to the SymPy closed form.
- **Grouped coordinates (II):** SymPy hand-codes `grouped_trace_anomaly`/`grouped_inverse`. The `.wl` *solves the linear system* `Solve[Thread[triple == xMean*traceVec + xAnom*anomVec + xSplit*splitVec], {xMean,xAnom,xSplit}]` (lines 159–168) and reconstructs the triple, then compares to the SymPy formulas — a Solve-based route, not a re-typed inverse.
- **Projectors (II):** SymPy hardcodes `Pbar,Pa,Pb`. The `.wl` *constructs* them from the metric via `weightedProjector[v] = Outer[Times, v, v.groupMetric]/(v.groupMetric.v)` over `DiagonalMatrix[{1,2,2}]` (lines 178–180, 224–226) and then checks the constructed projectors against the SymPy targets plus full completeness/idempotency/orthogonality algebra.
- **Packet B / logs (IV):** SymPy hand-writes `q_from_m`/`m_from_q`. The `.wl` derives `qFromMNative = Log /@ rFromMTarget` and obtains `m_from_q` by *solving the log-linear system* `Solve[Thread[logToQ == qVars], {ellT,ellK,ellMu}]` then `Exp /@` the solution (lines 332–350) — a solve-then-exponentiate route distinct from the asserted closed form.
- **Home-stretch theorem (V):** the `.wl` adds an explicit propositional check absent from SymPy — `LogicalExpand[(And@@(Join[branch,orbit]==0)) ⟺ ((And@@(branch==0)) && (And@@(orbit==0)))]` over abstract slot symbols, evaluating to `True` (lines 370–380).

House idioms only (`expectZero`/`expectZeroLog`/`expectTrue`, `stripCE`, `$Assumptions` positivity, the `math -script` + `Exit[0]` convention) follow the established stage-187-style template. This is a genuinely independent derivation, not a transliteration.

**Assessment — substantive (no harmful X−X / vacuous checks):**
40 PASS lines, 0 FAIL. Load-bearing content is all present: Packet A → Delta_branch (section III, native coefficient+basis compilers vs SymPy packet, plus isotropic and one-pole-normalized vanishing), Packet B → Delta_orbit (section IV interconversion round-trips), two-packet home-stretch theorem (section V), Taylor compilers (I), projector algebra over diag(1,2,2) (II). The cross-engine checks subtract a Mathematica-native quantity from the SymPy closed form, so they are not trivial X−X. The two individually-trivial checks are the same ones the original auditor already flagged for SymPy: `Delta_orbit solve zero-set` (solving q==0 for q) and the orbit-lock substitutions — harmless, and collectively corroborated by the eight non-trivial interconversion/log round-trips in section IV. Non-blocking.

## Exec log assessment

**SymPy:** exit=0 (`stage_191_sympy.log:348` → `# exit_code: 0`). All `expect_zero` lines print `= 0` or a zero matrix; header now reads STAGE 191.

**Mathematica:** exit=0 (`stage_191_mathematica.log:126` → `# exit_code: 0`). 40 PASS, 0 FAIL. Representative lines:
- `response coefficients - SymPy compiler = {0, 0}` → PASS
- `Delta_branch native route - SymPy packet = {0, 0, 0, 0, 0, 0, 0, 0}` → PASS
- `Delta_branch on isotropic one-pole normalized Packet A = {0,...}` → PASS
- `formal closure zero-set iff both packet zero-sets vanish = True` → PASS

**Output freshness:** confirmed. SymPy output (11:35:23) is newer than its script (11:30:24); Mathematica output (11:35:23) is newer than the `.wl` (11:35:04). Both exec logs are post-fix (11:35).

## Material-change assessment

`material_change`: false. F1 is a pure print-string relabel with zero mathematical propagation. F2 adds a new independent verifier that confirms (does not alter) the existing SymPy results. No derived quantity changed, so no downstream unit's inputs are affected.

## Side observations (non-blocking)

- The `Delta_orbit solve zero-set - quotient lock` check (`.wl:366`) is `Solve[q==0 for q]` then subtract — near-vacuous on its own, but it is the literal orbit-lock zero-set statement and the substantive orbit math is carried by the eight surrounding round-trip checks. This mirrors the auditor's note about the SymPy A29–A32 triviality; not a regression.
- Notes-prose legacy stage numbers (`242` home-stretch theorem, `..._stage242_...` references) remain to be reconciled in PAPER_CLEANUP_TRACKER per the directive — out of scope here.

## Verdict justification

Both findings are resolved. F1 is confirmed string-only from the captured diff (two banner literals, STAGE 174 → STAGE 191) with the refreshed SymPy header reading STAGE 191 and all residuals still zero. F2's new `.wl` is a genuinely independent Mathematica route — Series+Coefficient with a derivative cross-check, Solve-based basis-coordinate recovery, metric-constructed projectors, solved log-linear interconversion, and a LogicalExpand home-stretch-theorem check — covering every load-bearing claim and asserting cross-engine agreement (40 PASS, 0 FAIL, exit 0). Both engines exit 0 and outputs are fresh. No regressions, no blocking vacuity. Verdict: verified.

---
unit_id: 222
batch: VII.1
verifier_model: claude-opus-4-8[1m]
verify_date: 2026-06-09T19:30:00Z
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 0
findings_total: 0
material_change: false
---

# Verification — unit 222

## Per-finding outcomes

The original auditor report (`reports/stage_222.md`) returned `verdict: clean` with
`findings_count: 0` and wrote no directive (the directive path does not exist — expected
for a zero-finding stage). There are therefore no `Fn` to classify. This verification
instead confirms the `clean` verdict holds: that both scripts still pass under their own
assertions, the load-bearing checks are non-tautological, and the `.wl` is a genuine
independent recomputation rather than a transliteration of the `.py`.

### (zero findings) — auditor verdict re-confirmation

**Classification:** resolved (verdict holds; nothing to fix)

**What changed:**
Nothing. Codex applied zero edits. `redteam/pass2/exec_logs/stage_222_diff.patch` is
0 bytes (confirmed via `wc -c`), and `git status --porcelain` on both
`scripts/...stage222..._sympy_audit.py` and
`mathematica/...stage222..._mathematica_audit.wl` is empty — both files are byte-identical
to HEAD. This is the legitimately-empty diff of the first zero-correction pass-2 batch, not
a capture failure.

**Assessment:**
Re-reading both scripts confirms the auditor's non-tautology and independence calls:

- **A2 / M2 (exact quartic pole polynomial).** The `.py` builds `F_y` as a hand-written
  product form (`((K-M y)(varpi²-y)-C²)((Ω_U²-y)(Ω_W²-y)-R²) - (varpi²-y)(...)`, py
  L92-95) and asserts `cancel(D - F(omega²)/((varpi²-omega²)Δ)) == 0` (L97-100). The `.wl`
  derives the quartic from D's *own* cleared denominator — `nativeNumerator =
  Numerator[Together[dBranch]]` then `quarticFromD = Collect[Expand[PowerExpand[
  nativeNumerator /. omega->Sqrt[y]]], y, FullSimplify]` (wl L138-148) — and only then
  cross-checks it three independent ways: against `fProductY` (L157), by reconstructing
  `Together[D]` from it (L153-156), and by `Exponent`/degree (L162-163). Two
  independently-constructed objects are equated, not a define-then-assert `x==x`. The
  native-numerator route and the coefficient-list extraction have no `.py` counterpart.

- **A4 / M3 (residue/linewidth cancellation + survival threshold).** The `.py` constructs
  `Aqq_star = 1/D0prime` and `gamma_star = Gamma5·omega_star⁵·N_star/D0prime` *separately*
  (py L114-115) then divides, so the `D0prime` cancellation is exercised by construction,
  not assumed (L116-117); the result is checked against the closed form
  `27 c_s⁵/(a⁵ omega_star⁵ N_star)`. The `.wl` mirrors this with `rqDerived =
  aqqStar/gammaStar` (wl L168-172). For the survival threshold the routes genuinely
  diverge: the `.py` writes the closed form directly (py L122-124) whereas the `.wl`
  *inverts* it — `thresholdFromPeak = rqRatio /. First[Solve[stagePeak == deltaVreq,
  rqRatio]]` (wl L179) — then checks equality (L183). Inversion vs. direct assertion.

- **Root finding (A5-A9 / M4-M6).** `.py` uses `sp.nroots(expand(poly), n=30)` (py L9);
  `.wl` uses `NSolve[poly==0, y, WorkingPrecision->50]` (wl L67). Different solvers; the
  `.wl` residuals (M5 ~1e-17, M5 pure-Q ~1e-13, M6 scan ~1e-14) are independently computed
  in-engine, not echoed from the `.py`.

- **A10 / M6 (monotone tension).** Genuine sign tests on numerically-distinct values — P0
  rises 0.00595 → 0.10847 while upper-wall R_Q,* falls 145.48 → 4.827 across the lambda_W
  scan. Not a trivial equality.

- **Independence signature.** Variable naming diverges throughout (`lam_B`/`K`/`M`/
  `omega_star`/`a`/`c_s` in py vs `lamB`/`kWall`/`mass`/`omegaStar`/`aScale`/`cs` in wl),
  the `.wl` is organized into M1-M6 banner blocks vs the `.py`'s flat flow, and the
  overlap integral is computed natively in each engine (`sp.integrate` vs `Integrate`).
  No shared-op port signature. Independence call: **independent** — confirmed.

## Exec log assessment

**SymPy:** exit=0. Notable lines:
- `kappa = 2*sqrt(2)/pi` (overlap constant verified, A1).
- `with D(omega) = F(omega^2)/((varpi^2-omega^2) Delta(omega))` — quartic identity (A2)
  printed after its assertion passed; `deg F = 4` (A3) implicit in the F(y) print.
- `R_Q,* = 27*c_s**5/(N_star*a**5*omega_star**5)` (A4 cancellation).
- `(0.2, 0.005947405317693074, 2.8272344215844973, 2.04402272302752, 145.4838586578641)`
  — the corrected lambda_W=0.2 upper-wall R_Q,* = 145.48… (NOT the historical 213.x typo).
- `All Stage 222 symbolic and numerical checks passed.` then `# exit_code: 0`.

**Mathematica:** exit=0. Notable lines:
- `PASS: M1 overlap constant`, `M1 overlap constant residual = 0`.
- `native quartic degree = 4`; `PASS: M2 numerator-over-denominator reconstructs
  Together[D]`; `PASS: M2 native numerator equals product F[y]` — the native-denominator
  route (residual 0) confirms the independent quartic derivation.
- `PASS: M3 residue/linewidth cancellation` (residual 0); `PASS: M3 survival threshold
  form` (Solve-inverted, residual 0).
- `lambda_W=0.2 upper-wall R_Q,* = 145.48385865786412765610…` (corrected value, NOT 213.x);
  `PASS: M6 P0 increases with lambda_W`; `PASS: M6 upper-wall R_Q,* decreases with
  lambda_W`.
- `All Stage 222 Mathematica checks passed.` then `# exit_code: 0`.

Every PASS line cited in the audit's assertion inventory (M1-M6, A1-A10) is present; no
FAIL lines; no warnings.

**Output freshness:** Both exec logs carry a fresh top-of-file timestamp
`# date: 2026-06-09T19:21:18-06:00` from the orchestrator's direct re-run (post-fix in the
batch sense, though there was no fix). The orchestrator reports the run was deterministic —
committed `.txt` was byte-identical to the fresh run — so the fresh exec logs are the
freshness evidence; the committed `.txt` mtime is not the criterion here. Both engines exit
0 on the current (HEAD) scripts.

## Material-change assessment

`material_change`: false. No script was edited; both `.py` and `.wl` are unchanged from
HEAD. No derived result moved, so no downstream unit (> 222) can be affected by this
verification. The R_Q,* = 145.48… value already in the scripts and notes is the
previously-corrected figure and is reconfirmed here, not changed.

## Side observations (non-blocking)

- The auditor's report already flagged that the string `213.483858657863` survives only in
  out-of-scope process trackers (STAGE_VERIFICATION_COVERAGE / CHECKPOINT_CONSTANT_PROVENANCE
  / PAPER_CLEANUP_TRACKER / MATHEMATICA_MIRROR_POLICY) as historical-record prose describing
  the typo correction itself. That is consistent with scripts-only scope; no action and not a
  finding. Noting only for completeness.
- `cleanExpr`/`stripConditional` in the `.wl` correctly guard against
  `ConditionalExpression[0, …]` leaking from `Solve`/`Together` (the known Mathematica
  idiom), so the M3 threshold inversion's zero-residual check is robust.

## Verdict justification

The auditor returned a clean, zero-finding verdict and wrote no directive; there were no
edits to verify. I confirmed the empty diff is legitimate (0-byte patch, both scripts
identical to HEAD), that both exec logs are freshly timestamped and exit 0 with every
audited PASS line present and no FAILs, and — the load-bearing check for a clean verdict — I
re-read both scripts and confirmed the load-bearing assertions are non-tautological (A2/M2
equate two independently-built objects; A4/M3 exercise the D0prime cancellation by explicit
construction and invert the threshold via Solve; A10/M6 are genuine sign tests on distinct
values) and that the `.wl` is a genuine independent recomputation (native Integrate, native
Numerator[Together[D]] quartic extraction, NSolve root-finding, Solve-based threshold
inversion, divergent variable naming) rather than a transliteration of the `.py`. The
`clean` verdict holds. **verified.**

---
unit_id: 228
batch: VII.1
verifier_model: claude-opus-4-8[1m]
verify_date: 2026-06-09T19:30:00-06:00
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 0
findings_total: 0
material_change: false
---

# Verification — unit 228

Audit unit 228 was reported `clean` with **zero findings** (report front-matter `verdict: clean`, `findings_count: 0`). There is no directive file for 228 (expected — no findings to direct). Codex applied no edits: the captured diff patch `stage_228_diff.patch` is genuinely empty (0 bytes), and `git status --porcelain` on both script files is empty (unchanged from HEAD). This verification therefore confirms that the `clean` verdict still holds under a fresh re-read of both scripts and both exec logs — it is not a fix-closure check.

## Per-finding outcomes

No findings were raised. Nothing to classify per-finding. The confirmation criteria for a zero-correction unit are:

**(a) Scripts re-read, unchanged from HEAD.** Both `.py` (469 lines) and `.wl` (379 lines) re-read in full. `git status --porcelain` empty for both; diff patch 0 bytes. No surviving stale values — grep for `247`, `251+215`, `215π`, `247π` across both scripts and both `.txt` outputs returned no match.

**(b) Both exec logs pass.** SymPy log ends `# exit_code: 0` with the final `Stage 228 audit completed successfully.` banner. Mathematica log ends `# exit_code: 0` with all M1–M9 `PASS:` lines present (no `FAIL:`). Both logs dated 2026-06-09T19:21:18-06:00 (fresh).

**(c) Load-bearing asserts are non-tautological and the `.wl` is a genuine independent recomputation.** Confirmed (detail below).

## Non-tautology spot-check

- **Row coefficients / determinant (A2/A4/M2/M3/M5).** The `\pi_1` and `\delta_1` rows and the reduced determinant are DERIVED in-script from the dressed-primitive linearization (`.py:136-186` symbolic `sp.diff(...).subs(eps,0)`; `.wl:97-151` `Coefficient[Normal[Series[...]],eps,1]`) and only THEN compared to the expected literals (`expected_delta` at `.py:194`, `expectedDelta` at `.wl:164`; `expected_det` at `.py:250-253`, `expectedDet` at `.wl:200-203`). If the derivation produced any other coefficient the asserts would fail — derivation-then-compare, not define-then-assert. Non-tautological.
- **Generators (A5/M6).** `v_num`/`v_den` come from `M_num.nullspace()[0]` / `First[numNull]` (`.py:273-274`, `.wl:214-215`) — computed from the matrix, then both the literal vector AND the derived identities `\Xi_1=-2\delta_1` / `\Xi_1=2\pi_1` are asserted (`.py:300-301`, `.wl:233-234`). The identity check is independent of the literal check.
- **Ceiling comparison (A8/M9).** The strict `>` comparisons (`.py:454-456`, `.wl:373-375`) compare independently-computed dynamic vs transported-static ceilings, not a value against itself.

## Independence check (Mathematica)

Confirmed the `.wl` is a genuine independent re-derivation, not a transliteration, on the decisive route — the dynamic slopes:

- **`.py:383-393`** computes log-slopes by symmetric finite difference (`step=1e-8`, evaluate at `±step`, `(log(p)-log(m))/(2*step)`).
- **`.wl:288-301`** computes them analytically via the implicit-function theorem: `fEps=D[fBranch,eps]`, `fY=D[fBranch,y]`, `dy=-fEps/fY`, `rqSlopes=(rqEps+rqY·dy)` evaluated symbolically at `eps=0`. Different mathematics, not the same algebra in two syntaxes.

The numerical fingerprint confirms it: the `R_{Q,-}` numerator-rigid slope is `0.7135849466877175` (SymPy, finite-diff; log line 37) vs `0.71358483642759108` (`.wl`, analytic; log line 76) — agreement to ~7 digits with divergence in the 8th. That is the finite-difference truncation-error signature a port cannot fake. Linearization (symbolic-diff vs Series-extraction) and pole census (`numpy.roots` on numeric coeffs vs `NSolve` at `WorkingPrecision->60`) are likewise distinct routes. Both scripts are self-contained (no `Get`/`Needs`/`Import`). **INDEPENDENCE CONFIRMED.**

## Deliverable confirmation

Corrected values emitted in BOTH engines, no stale residue:
- `\delta_1` coeff `196π²/(98π²−25)` on `\Omega_U,\Omega_W` — py-out l.6, wl-out l.5 (and M3 squared-difference test = 0, log line 17).
- Reduced determinant `196(200+147π²)(80000+343225π²+43218π⁴)/(475(8670000+14894275π²+2117682π⁴))` — py-out l.13, wl-out l.31, both byte-equivalent; M5 = 0 (log line 37).
- Grep for `247`, `251+215π²`, `215π` across both scripts and both `.txt` outputs: no match.

## Exec log assessment

**SymPy:** exit=0. Notable lines: `delta_1= -50*xLR/(-25 + 98*pi**2) + 196*pi**2*xOU/(-25 + 98*pi**2) + 196*pi**2*xOW/(-25 + 98*pi**2)` (l.11); `det[...] = 196*(200 + 147*pi**2)*(...)/(475*(...))` (l.18); `Stage 228 audit completed successfully.` (l.47). All `assert`s passed (Python aborts on the first failed assert; clean banner ⇒ all passed).

**Mathematica:** exit=0. Notable lines: `M5 reduced determinant = 0` → `PASS` (l.37-38); `M8 numerator dln R_Q,- ... diff = 3.57e-9` → `PASS` (l.87-88); all of M1–M9 carry `PASS:` lines, zero `FAIL:`; `Stage 228 Mathematica audit completed successfully.` (l.132).

**Output freshness:** orchestrator re-ran both engines directly (logs dated 2026-06-09T19:21:18). Per batch context the committed `.txt` outputs are byte-identical to the fresh deterministic run; the corrected deliverables grep-confirmed present in the committed `.txt`. No staleness concern; per instructions did not fail on committed `.txt` mtime.

## Material-change assessment

`material_change`: **false**. Codex made no edits; the scripts are byte-identical to HEAD. No derived result changed, so no downstream unit is affected by this verification.

## Side observations (non-blocking)

None. The scripts were already clean at audit time and remain so.

## Verdict justification

Audit unit 228 was a zero-finding `clean` unit; the empty diff patch and empty `git status` confirm Codex correctly applied no edits — this is the first pass-2 batch member needing zero script corrections, and the empty patch is legitimate, not a failure. Both exec logs exit 0 with all PASS lines and clean banners. The load-bearing assertions (row coefficients, reduced determinant, generators, ceiling comparison) are derivation-then-compare and non-tautological. The `.wl` is a genuine independent recomputation, decisively so on the dynamic slopes (analytic implicit-differentiation vs SymPy finite-difference), corroborated by the ~7-digit-then-diverge truncation signature on `R_{Q,-}`. The corrected `\delta_1` coefficient and reduced determinant are emitted by both engines with no surviving stale `247`/`251+215π²` residue. Verdict: **verified**, material_change false.

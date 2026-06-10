---
unit_id: 220
batch: VII.1
verifier_model: claude-opus-4-8[1m]
verify_date: 2026-06-09T00:00:00Z
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 1
findings_total: 1
material_change: false
---

# Verification — unit 220

## Per-finding outcomes

### F1 — paper_misalignment (paper_missing_script_claim; card-text-lag)

**Classification:** resolved (deferred-correctly; no script change required or possible)

**What changed:**
Nothing in any script. The diff patch
`redteam/pass2/exec_logs/stage_220_diff.patch` is legitimately empty (0 content
lines); both `scripts/...stage220...sympy_audit.py` and
`mathematica/...stage220...mathematica_audit.wl` are byte-identical to HEAD. The
directive's body explicitly instructs Codex to do **nothing** and contains no
`## Applied:` block — by design, because F1 is a doc-side STATUS annotation
(`stage_220.tex:11` `\stagefield{Verification}` reads "Mathematica audit: none
yet" while a committed, passing `.wl` exists).

**Assessment:**
Correct disposition. F1 is OUTSIDE scripts-only scope: it is a stale
verification-coverage label on the paper card (and the notes Section 10 wording),
not a script defect. The math is dual-engine green and there is no possible
script edit that resolves it. Per the standing user decision the doc relabel is
deferred to PAPER_CLEANUP P4-51 (the card understates coverage; it does not
misstate any value). Treating this as the sole, deferred finding — and NOT
rolling the verdict to needs_rework on its account — is the correct call. The
scripts-only verifier's job here is to confirm the scripts are clean, which they
are (below).

## Exec log assessment

**SymPy:** exit=0. Notable lines:
- "Verified determinant identity: det(K_dyn) = Delta_Pi * D_Pi"
- "Verified static reduction back to the Stage 219 one-port bundle"
- "Verified outgoing-port derivative identity: dV_mix/dPi = -1/2 * T_J^2"
- "All Stage 220 symbolic checks passed." then `# exit_code: 0`

**Mathematica:** exit=0. Notable lines:
- "PASS: M1 Det[Kdyn] - DeltaPi DPi" … through all of M1–M9
- "M6 Laurent support = {{-6, 0}, {-4, 1}, {-2, 2}}" + "PASS: M6 no additional
  x-y monomial families" (independent closure proof)
- "PASS: M7 self-test VmixExpected depends on Pi" (variable-independence trap guard)
- "PASS: M9 Pabs perfect-square form" + "All Stage 220 Mathematica checks passed."
  then `# exit_code: 0`

Both logs PASS every block and exit 0.

**Output freshness:** the orchestrator re-ran both engines directly (logs dated
2026-06-09T19:21:18-06:00, exit 0, deterministic; committed `.txt` reported
byte-identical to the fresh run). No staleness concern; I did not fail on
committed `.txt` mtime per the batch close-out note.

## Spot-check: load-bearing asserts non-tautological, `.wl` is a genuine recomputation

- **P_abs perfect-square (load-bearing; first-pass F2 fix).** In the `.py`
  (193–198) `real_part, imag_part = expand_complex(deltaV_linear_expected).as_real_imag()`
  DERIVES the imaginary part from the full complex `δV^(1) = -i/2·Γ·T_J0²`, then
  asserts `simplify(P_abs - omega*Gamma/2 * T_J0**2) == 0` where
  `P_abs = factor(-omega*imag_part)`. This is a derivation, not a posited
  identity: it holds only because `T_J0` is real on the conservative branch (all
  base symbols `real=True`), so `Im[...]` collapses to `-1/2·Γ·T_J0²`. The `.wl`
  (250–256) independently extracts `phaseImag = ComplexExpand[Im[dV1]]` and
  asserts `Pabs - omega·gammaOut/2·TJ0^2 === 0` — a route-distinct confirmation of
  the same perfect-square form. Non-tautological in both engines.
- **`.wl` is independent, not a transliteration.** It carries machinery absent
  from the `.py`: (1) M6 Laurent-support closure via `CoefficientRules` on
  `-2·VprimY·DeltaPi·DPi·x^6` (wl 223–229) — proves the spatial-family support set
  is *closed* `{{-6,0},{-4,1},{-2,2}}`, strictly stronger than the `.py`'s
  form-match; (2) M6b coefficient re-extraction via `Coefficient` (wl 216–221);
  (3) M7a self-test `(dVdPi /. sampleConservativeRules) != 0` (wl 235–238) guarding
  the identically-zero-derivative trap; (4) M9c dissipative-sample guard (wl
  257–264). The shared signature (native `Det`/`Inverse` compared to the same
  published closed forms) is correct second-engine behaviour, not a port.
- No tautological rows anywhere: every assert subtracts an independently
  constructed form from a natively computed quantity (native `det`/`inv`,
  `expand_complex`/`ComplexExpand`, `diff`/`D`), or is a genuine boolean guard.

## Material-change assessment

`material_change`: false. No script bytes changed; no derived result moved.
Downstream units depending on Stage 220 deliverables are unaffected.

## Side observations (non-blocking)

None. The card-text-lag is the single (deferred, out-of-scope) finding; nothing
else in the scripts warrants attention.

## Verdict justification

VERIFIED. This is the expected zero-script-correction close-out: the diff patch
is legitimately empty, Codex correctly applied nothing, and both engines pass all
blocks with exit 0 on a fresh deterministic re-run. The lone finding (F1) is a
stale verification-coverage label on the paper card — outside scripts-only scope
and deferred to PAPER_CLEANUP P4-51 by standing user decision — and correctly
does not force needs_rework. The load-bearing P_abs perfect-square assert is a
real derivation in both engines, and the `.wl` is a genuine independent
recomputation (Laurent-support closure, coefficient re-extraction, Pi-dependence
and dissipative self-tests) rather than a transliteration. Scripts clean.

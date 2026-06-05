---
unit_id: 047
batch: III.1
verifier_model: claude-opus-4-8[1m]
verify_date: 2026-06-05T00:00:00Z
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 1
findings_total: 1
material_change: false
---

# Verification — unit 047

## Per-finding outcomes

### F1 — stale_output (self-label, label-only)

**Classification:** resolved

**What changed:**
Codex edited only the two UNAMBIGUOUS stale self-labels in the SymPy source
`scripts/moving_throat_pde_stage047_coherent_kernel_map_sympy_audit.py`:
- line 3 (docstring header): `Stage 30 SymPy audit.` → `Stage 47 SymPy audit.`
- line 158 (closing pass-line): `All Stage-30 symbolic checks passed.` → `All Stage-47 symbolic checks passed.`

The captured diff (`exec_logs/stage_047_diff.patch`) is exactly these two hunks and nothing
else: two single-token `30`→`47` substitutions with all surrounding text preserved. The `.wl`
does not appear in the diff (untouched, as directed — its labels were already canonical).

**Assessment:**
The change is strictly LABEL-ONLY. I re-read the live `.py`:
- the only mutated tokens are the stage number in the docstring (line 3) and the closing
  print (line 158); the banner already read `STAGE 47` (line 32) and was left as-is.
- No equation, value, variable name, or assertion was touched. `expect_zero` body (lines
  25–29), all closed-form constructions, the nine identity asserts, the negative control
  (lines 121–156), and the spoiled-probe `.subs` dict are byte-for-byte as before the fix.
- The two INTENTIONAL DEFERRALS are correctly UNTOUCHED, confirmed by grep:
  - `.py` line 44 `... established upstream at Stage 28 ...` — cross-ref to upstream stage 045, owned by the dedicated numbering pass. Present, unchanged.
  - `.py` line 121 `... perturbed off the exact Stage-30 support-loading coefficient ...` — an ambiguous self-vs-cross comment that does not reach the transcript. Present, unchanged.
  - `.wl` line 40 `... established upstream at Stage 28 ...` — present, unchanged (the entire `.wl` was untouched).

Refreshed outputs confirm the labels propagated through a real re-run:
- SymPy txt: banner line 8 `STAGE 47 — COHERENT KERNEL MAP AUDIT`; close line 58
  `All Stage-47 symbolic checks passed.`
- Mathematica txt: banner line 8 `STAGE 047 — COHERENT KERNEL MAP`; close line 45
  `Stage 047 Mathematica audit passed.` (already canonical; `.wl` not in scope of this fix).

Every prior PASS / `= 0` residual survives unchanged in both transcripts (see exec assessment
below), and the spoiled-control probe still fires nonzero (`-44.8664…`, sympy txt line 56),
so the math content is provably identical to pre-fix. Finding fully resolved.

## Exec log assessment

**SymPy:** exit=0. Notable lines:
- `eps_phi - zeta_def*eps_W = 0` / `Z_phi - zeta_def*Z_W = 0` (txt 25–26)
- `M_tr - M_mix*S = 0` (txt 43), `product law = 0` (txt 50), `support-loaded R_target reconstruction = 0` (txt 51)
- `dR_target_loaded/dzeta = 0`, `dS/dzeta - (1-eps)/(1-zeta eps)^2 = 0`, `S(zeta=0)-1 = 0` (txt 52,54,55)
- `spoiled dR_target/dzeta at a probe point = -44.8664510213650197722046979701` (txt 56) — negative control fires
- `All Stage-47 symbolic checks passed.` (txt 58), `# exit_code: 0` (txt 59)

**Mathematica:** exit=0. Notable lines:
- `PASS: eps_phi - zeta_def eps_W`, `PASS: Z_phi - zeta_def Z_W` (txt 17,19)
- `PASS: M_tr - M_mix S`, `PASS: product law`, `PASS: support-loaded R_target reconstruction` (txt 29,33,35)
- `PASS: dR_target_loaded/dzeta`, `PASS: R_target_loaded(zeta) - R_target_loaded(0)`, `PASS: dS/dzeta - (1-eps)/(1-zeta eps)^2`, `PASS: S(zeta=0)-1` (txt 37,39,41,43)
- `Stage 047 Mathematica audit passed.` (txt 45), `# exit_code: 0` (txt 46)

All nine identities still hold under both engines; counts and PASS set match the original
auditor's inventory (A1–A9 / B1–B9 + A10 control). No new failures or warnings.

**Output freshness:** confirmed re-generated post-fix. mtimes:
- `.py` 2026-06-05 07:54:25; sympy `.txt` 2026-06-05 08:12:32 (newer ✓)
- `.wl` 2026-06-03 15:59:11 (untouched); mathematica `.txt` 2026-06-05 08:12:32 (newer ✓)
Both transcripts were captured after the edit; the stale `STAGE 30`/`STAGE 030` banners noted
in F1 are gone.

## Material-change assessment

`material_change`: false. The edit touched only two comment/print stage-number tokens in the
SymPy source. No equation, derived value, variable, or assertion changed; both engines emit the
identical `= 0` residuals and PASS set as before, and the `.wl` was not modified. Nothing
downstream of 047 depends on a docstring or pass-line label, so no downstream unit is affected.

## Side observations (non-blocking)

- The two deferred numbering refs (`.py:44` `Stage 28` cross-ref, `.py:121` `Stage-30`
  ambiguous comment, plus `.wl:40` `Stage 28`) remain and are correctly owned by the
  dedicated `NUMBERING_SCRIPT_OUTPUT_BAND_PLAN.md` pass — not defects here.
- The Mathematica transcript close line reads `Stage 047 Mathematica audit passed.` while the
  SymPy one reads `All Stage-47 ...`; both are canonical (047/47), just differing wording per
  engine. Not a defect; the `.wl` was out of scope for this label-only fix.

## Verdict justification

The sole finding (F1, low-severity stale_output / self-label) is fully resolved. The captured
diff is strictly label-only — two `30`→`47` token substitutions in the SymPy docstring and
closing print, with no equation, value, variable, or assertion altered, and the intentional
`Stage 28` cross-ref and `Stage-30` ambiguous comment correctly left for the dedicated
numbering pass. Both engines re-ran to exit 0 with every prior PASS/`= 0` residual intact and
the negative control still firing; both `.txt` transcripts are freshly regenerated (newer
mtimes) and now show canonical STAGE 047/47 banners. No regressions, no new findings.
Verdict: verified, material_change: false.

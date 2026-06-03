---
unit_id: 244
batch: VIII.1
verifier_model: claude-opus-4-8[1m]
verify_date: 2026-06-03T10:05:00-06:00
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 3
findings_total: 3
material_change: false
---

# Verification — unit 244

## Per-finding outcomes

### F0 — notes correction (USER-APPROVED 2026-06-03)

**Classification:** resolved

**What changed:**
`notes/stages/moving_throat_pde_stage244_selected_branch_leakage_and_scalar_photon_work_compiler_sympy_audit.md:366`: `196\sqrt2` → `128\sqrt2`. The git diff against HEAD shows exactly one changed line (the `\varrho`-form bulk-work coefficient inside the second boxed equation); no other notes line moved.

**Assessment:**
Correct and minimal. The script was NOT touched on this value: `scripts/...sympy_audit.py:106` still carries `W_bulk_pull_varrho_expected = 128 * sqrt(2) * ...` and the SymPy parity output line 60 prints `128*sqrt(2)*...`. The published card `paper/stages/stage_244.tex` is untouched (clean git status). The new `.wl` independently recomputes the `(1-eps)`-form coefficient `512√2/9` (M6 "Wbulk pullback epsilon" = 0), which equals `128√2` under `(1-eps)²=(9/4)ϱ²` — genuine cross-engine corroboration of the corrected value. Codex edited only notes:366 as authorized.

### F1 — insufficient_verification (variable-independence self-test trap)

**Classification:** resolved

**What changed:**
`scripts/...sympy_audit.py` Section 4 (diff lines 8-29): the three derivative computations and the `assert d_Rtr == 0 / d_Rtarget == 0 / d_epseta == 0` lines are GONE. Replaced with `orbit_syms = {R_tr, R_target, eps_eta}`, `support_syms = {Lam, varrho, eta_leak}`, and per-observable `assert orbit_syms.isdisjoint(free)` plus the positive control `assert support_syms.issubset(free)`.

**Assessment:**
Non-vacuous and load-bearing. The `isdisjoint` test fires the instant any orbit symbol leaks into a compiled observable; the `issubset` positive control rules out the degenerate "constant is free of everything" pass and would fail on a regression that collapsed an observable. The SymPy log confirms both are exercised against real expressions: `S_leak free symbols = [Lam, eta_leak, lam, mu_w, rho0, varrho]`, `orbit overlap = []`, `support coverage = [Lam, eta_leak, varrho]` (works add `q`). The old `assert d_* == 0` lines do not appear anywhere in the post-fix script. Matches the directive's prescribed structure verbatim.

### F2 — missing_verification_script (missing_mathematica)

**Classification:** resolved

**What changed:**
New `mathematica/...mathematica_audit.wl` created (240 lines), implementing the M1–M9 claim manifest with `expectZero`/`expectTrue` helpers (both `Exit[1]` on failure). `etaleak` declared `Element[..., Reals]` (not positive), enabling the parity sign-flip; positive physical scales (`lam,E0,muw,rho0,q,Lam>0`, `0<eps<1`) asserted separately.

**Assessment:**
Genuinely independent route, not a `.py` transliteration. The leakage and bulk integrals are derived by native `Integrate[..., {w,-Infinity,Infinity}]` (2 `Integrate` calls, M2/M3); M4 builds the work-leak relation from those integrated quantities (`bulkIntegral - q E0 leakageIntegral`), and M5 builds the session/quadratic relations from `bulkIntegral`/`leakageIntegral` — genuine engine-internal cross-checks, not restatements of closed forms. M7 uses `FreeQ` + `Not[FreeQ[...]]` positive controls, NOT `D[...]` — the F1 trap is not reproduced. The single `D[` in the file is `D[projector, w]` at line 77, the legitimate Gaussian-derivative integrand for the leakage moment, not a derivative-w.r.t.-orbit-symbol check. M8 parity substitutes `etaleak -> -etaleak`; M9 substitutes `etaleak -> 0`. The Mathematica log shows all 24 checks PASS and exit 0, with M6 reproducing the `512√2/9`, `2048/9`, and `32/(3π²)` constants that match the SymPy script and the card — the engine cross-check is now live.

## Exec log assessment

**SymPy:** exit=0. Notable lines:
- `S_leak orbit overlap = []` and `S_leak support coverage = [Lam, eta_leak, varrho]` (F1 non-vacuous, both asserts exercised against real expressions).
- `W_bulk pullback = 512*sqrt(2)*...*(eps - 1)**2/(9*pi**(9/2)*lam**3)` and `W_bulk(-eta) = 128*sqrt(2)*...varrho**2/(pi**(9/2)*lam**3)` (corroborates F0 128√2).
- `All symbolic checks passed.`

**Mathematica:** exit=0. Notable lines:
- `PASS: M2 leakage compiler`, `PASS: M3 bulk work`, `PASS: M4 work-leak relation` (native Integrate-derived, cross-relation).
- `PASS: Spull orbit-free` + `PASS: Spull support-dependent` (FreeQ split, non-vacuous, no D[] on orbit syms).
- `PASS: M8 Spull odd in etaleak` / `M8 Wbulkpull even` (parity). 24 PASS lines total; `All Stage 244 Mathematica checks passed.`

**Output freshness:** confirmed. `.py` mtime 1780501554 and `.wl` mtime 1780501771 both precede both `.txt` output mtimes (1780502062) — outputs were regenerated post-fix.

## Material-change assessment

`material_change`: false. No derived value or closed form changed: the SymPy script's math (all constants, integrals, closed forms) is byte-identical except for the verification surface in Section 4 (derivative → free-symbol set test). F0 fixed a notes-only numeric typo, not a script/card value. F2 is purely additive (a new independent engine). No downstream unit's inputs change. Units > 244 do not become stale on account of this fix.

## Side observations (non-blocking)

- The diff patch at `redteam/exec_logs/stage_244_diff.patch` captured only the `.py` hunk; the notes:366 edit and the new untracked `.wl` are not in that patch. I verified both directly against the working tree (git diff for notes, file read + git status for the untracked `.wl`), so coverage is complete — but the orchestrator may want the diff-capture step to also include notes edits and new untracked scripts for future units.

## Verdict justification

All three applied findings are genuinely resolved. F0 is a minimal, authorized one-line notes typo fix with the script value and published card both confirmed untouched and the corrected coefficient cross-corroborated by the new `.wl`. F1's vacuous orbit-derivative asserts are removed and replaced by a non-vacuous free-symbol disjointness test guarded by a passing positive control. F2 adds a genuinely independent Mathematica route (native Gaussian Integrate, FreeQ-based split, no `D[]` orbit trap, integral-derived cross-relations) that exits 0 with 24 passing checks and reproduces the SymPy constants. Both engines pass, outputs are fresh, and no derived result changed (material_change=false).

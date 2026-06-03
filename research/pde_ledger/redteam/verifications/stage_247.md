---
unit_id: 247
batch: VIII.1
verifier_model: claude-opus-4-8[1m]
verify_date: 2026-06-03T12:40:00-06:00
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 6
findings_total: 6
material_change: false
---

# Verification — unit 247

## Per-finding outcomes

### F1 — paper_misalignment (USER-RESOLVED, notes-only correct-to-script)

**Classification:** resolved

**What changed:**
`notes/stages/...stage247...sympy_audit.md:406` now reads `\Delta=142.17750000,` (was `210.17750000`). Confirmed via grep: the disputed `210.1775` no longer appears anywhere in the notes file. The adjacent notes value `D_0=3.76481862` (line 407) is consistent only with Δ=142.1775.

**Assessment:**
Correct and scoped exactly to the one authorized line.
- Published card `paper/stages/stage_247.tex` is UNAFFECTED: it is 93 lines, contains NO occurrence of `210.1775` OR `142.1775` (grep returned none), and is NOT in `git status` (untouched). The auditor's `stage_247.tex:407` citation was a misattribution; the value only ever lived in the notes.
- The script literal is unchanged: py:266 still asserts `abs(Delta_num - 142.1775) < 1e-8`. No collateral edit. Both engines independently compute Δ=142.1775 from `OmU2*OmW2 - Rmix^2` (sympy log line 62; .wl M6).

### F2 — tautological_check (Lvar inverse round-trip)

**Classification:** resolved

**What changed:**
The symbolic round-trip `assert sp.simplify(Lvar_from_W - Lvar) == 0` remains at py:122 (left as documentation, explicitly permitted by the directive), but it is no longer load-bearing. New numeric asserts added at py:234-236: `Wsess_from_Lvar` (W_sess re-evaluated at the inverted `Lvar_soft`) reproduces `Wsess_obs` to 1e-7, and `Lvar_soft` matches the paper figure `20.01677473` to 1e-6.

**Assessment:**
Non-tautological. `Lvar_soft` (py:227-230) is solved from `Wsess_obs` by an explicit closed-form sqrt; the `Lvar_soft ~ 20.01677473` pin is falsifiable — a wrong inversion constant breaks it (printed value 20.01677472685125). The sqrt-vs-square round-trip is no longer the verification of record. `.wl` M7 independently pins Lvar to the same figure (log line 72, PASS).

### F3 — script_missing_paper_claim (D_UV>=0, g>=2/pi, M_sigma>=0)

**Classification:** resolved

**What changed:**
- py:132-140: drain nonnegativity — `sp.ask(Q.nonnegative(D_UV)) in (True, None)` plus a numeric probe `D_UV_probe >= 0` (probe yields a positive literal 0.01867).
- py:221-225: branch positivity at `r_soft` — `g_soft >= 2/pi`, `g_soft < rF1`, `M_sigma_num >= 0`.

**Assessment:**
Genuine, non-vacuous. The probe assert is the load-bearing half (the `ask` returns a permissive set so it never blocks, but the probe pins a real positive literal). `g_soft` and `M_sigma_num` are concrete numbers over variables the expressions genuinely contain (g_r depends on r/r_sigma/a0/b0; M_sigma printed 0.18386120). `.wl` M4 strengthens this with a native `Reduce[... Duv<0] === False` (log line 38, no negative branch) plus a positive probe; M7 checks the branch bounds (PASS).

### F4 — insufficient_verification (V_eff <= V_short never checked)

**Classification:** resolved

**What changed:**
Identity assert kept at py:181 (valid structural check). New numeric inequality at py:244-247 using the paper literal `lambda_L_paper = sp.Float("0.26971918")`: `Veff_session` built forward, then `assert float(Vshort_num - Veff_session) >= 0`.

**Assessment:**
Genuine. The gap equals `lambda_L*S + Wsess + UVdrop + M_sigma`, a sum of independently-pinned positive numbers (~1.9946); it would fail if any packet term were wrongly signed. Not vacuous. `.wl` M8 forward decomposition corroborates.

### F5 — tautological_check (benchmark lambda_L self-solve)

**Classification:** resolved

**What changed:**
The tautological closure `assert abs(Vrebuild_soft - Veff_obs) < 1e-10` is DELETED (diff line 62; absent from current py — `Vrebuild_soft` is still computed at py:251 for the print banner only, no longer asserted). Replaced with: four pins (py:240-243) `Vshort_num~3.74163698`, `M_sigma_num~0.18386120`, `S_soft~0.31069599`, `lambda_L_soft~0.26971918`; and a FORWARD closure (py:248-250) `Veff_forward` built from the paper literal `lambda_L=0.26971918` asserted against `Veff_obs=1.74701126` to 1e-6.

**Assessment:**
The solve-then-resubstitute tautology is gone as a check. The forward decomposition genuinely reconstructs 1.74701126 from independently-pinned quantities + the paper's literal lambda — it would fail on a wrong packet term (hand-check: 0.26971918*0.31069599 + 1.51632107 + 0.21064278 + 0.18386120 subtracted from 3.74163698 = 1.7470113). The `lambda_L_soft~0.26971918` pin (the solved value) is a legitimate independent cross-check that the solved and paper lambda agree, not the closure of record. `.wl` M8 reproduces VeffForward = 1.74701126 forward (log line 84, PASS).

### F6 — missing_verification_script (missing_mathematica)

**Classification:** resolved

**What changed:**
New `mathematica/...stage247...mathematica_audit.wl` (untracked, mtime 10:16 < output 12:33). Covers M1-M8 with `expectZero`/`expectTrue`/`expectApprox` helpers (each `Exit[1]` on fail). Exit 0; all 16 expect-checks PASS in the log.

**Assessment:**
Genuinely independent route, not a `.py` transliteration:
- M1: native `Inverse[Kred]` + `FullSimplify`, compared to paper closed forms `chiClosed` derived independently from qMix/pMix/dZero (the directive prescribed this exact route).
- M3: native `Limit[gMoment[rLimit], rLimit->Infinity]` (= 2/Pi).
- M4: native `Reduce[... Duv<0, Reals] === False` (a stronger, structurally different check than the .py probe-only).
- Session uses exact rationals (e.g. Rmix->27/20, rF1->88899676773749/50000000000000), distinct from the .py float dict.
- M6/M8 use the formula-correct Δ target `1421775/10000 = 142.1775` and lambda literal `13485959/50000000 = 0.26971918` — does NOT bake in the disputed 210.1775, and uses the forward (not solved) lambda. Confirmed by grep: no 210.1775 in the .wl.
Expected PASS lines counted: M1, M2(×2), M3, M4(×2), M5, M6(×3), M7(×6), M8(×2) = 17 PASS lines; log shows all PASS, exit 0.

## Exec log assessment

**SymPy:** exit=0. Notable lines:
- `Delta(session) = 142.1775000000000`, `D0(session) = 3.764818624606566`, `V_short(session) = 3.741636979285524`
- `Lvar(session) = 20.01677472685125`, `S_leak(session) = 0.3106959887767410`, `M_sigma(session) = 0.1838612007190837`
- `lambda_L(session) = 0.2697191840048420`; `V_eff rebuilt = 1.747011260000000`; "All symbolic checks passed."

**Mathematica:** exit=0. Notable lines:
- `PASS: M1 inverse-entry closed forms` (residual matrix all-zero), `PASS: M3 g infinity`
- `M4 drain negative branch Reduce = False` → `PASS: M4 no negative drain branch`
- `M6 Delta session` PASS (142.1775), `M8 forward V_eff` PASS (VeffForward = 1.7470112612...), `PASS: M8 lambda_L positive`
- "All Mathematica checks passed."

**Output freshness:** confirmed. SymPy script mtime 10:13, output .txt 12:33 (newer). `.wl` mtime 10:16, output .txt 12:33 (newer). Both outputs regenerated post-fix.

## Material-change assessment

`material_change`: **false**. No derived script value changed: the SymPy script's numeric literals (142.1775, 3.76481862, 20.01677473, 0.31069599, 0.18386120, 0.26971918, 1.74701126) are unchanged — the edits only WIDEN verification (new asserts) and DELETE one vacuous assert. The only prose edit is a single notes typo correction (210→142) bringing the notes into agreement with the already-correct script. The new `.wl` is an independent corroborating engine, not a result change. No downstream unit's inputs are altered.

## Side observations (non-blocking)

- The `sp.ask(Q.nonnegative(D_UV)) in (True, None)` assert is intentionally permissive (passes on `None`); the load-bearing nonnegativity check is the numeric probe and the `.wl`'s `Reduce` False-branch. Acceptable as written; the directive anticipated `None` returns.
- `Vrebuild_soft` (py:251) is still computed but only printed, not asserted — harmless residual scaffolding, no longer load-bearing.

## Verdict justification

All six findings are resolved. F1 is a clean notes-only typo fix (210→142) with the published card confirmed untouched and value-free, and the script literal unchanged. F2 and F5 tautologies are no longer load-bearing — F5's solve-then-resubstitute assert is deleted and replaced with a falsifiable forward decomposition from independently-pinned quantities plus the paper's literal lambda. F3/F4 add genuine, non-vacuous numeric positivity and lowering-inequality probes. F6 adds a genuinely independent Mathematica engine (native Inverse/Limit/Reduce, exact-rational session, formula-correct Δ=142.1775, forward lambda) that passes all checks. Both engines exit 0 with all in-file checks PASS, outputs are fresh, and no derived result changed (material_change false).

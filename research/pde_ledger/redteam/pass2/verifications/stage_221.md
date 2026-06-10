---
unit_id: 221
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

# Verification — unit 221

This is the VII.1 batch checkpoint and the FIRST pass-2 zero-correction unit. The auditor returned `clean` with 0 findings, so there was no directive and nothing for Codex to apply; the git diff patch (`exec_logs/stage_221_diff.patch`) is legitimately 0 bytes and both scripts are unchanged from HEAD (confirmed via `git status --short` → no output). My job is to confirm the auditor's `clean` verdict holds: both engines exit 0 with substantive PASS lines, the load-bearing asserts are non-tautological, and the re-authored `.wl` is genuinely independent of the `.py`.

## Per-finding outcomes

No findings in the original report (`findings_count: 0`). There is nothing to classify as resolved/partial/regressed. The verification below confirms the auditor's affirmative verdict rather than re-checking applied fixes.

## Independence confirmation (checkpoint bar)

The auditor's INDEPENDENCE CALL = `independent`, with the decisive evidence in §II (opposite information flow on the Stage-220 derivative identity). I read the cited lines and confirm:

- **§II derivative identity `dD_Pi/dPi = -N` — opposite information flow.** `.py` L50-52 POSITS the perfect square as an INPUT: `Nfun = sp.together((A*G_W + R*G_U)**2/Delta_Pi**2)`, then asserts `sp.diff(D_Pi,Pi)+Nfun==0`. `.wl` L96-102 DERIVES `NfunDerived = Together[D[QPi/DeltaPi, portPi]]` natively (an OUTPUT), then checks BOTH (a) `NfunDerived - (Afun GW + R GU)^2/DeltaPi^2 == 0` (proves the derived N equals the perfect square) AND (b) `D[DPi, portPi] + NfunDerived == 0`. The perfect-square form is an output on the Mathematica side and an input on the SymPy side — strictly more independent, not a shared posited form. Confirmed.
- **`.wl`-only `Residue[]` extraction.** `.wl` L78-81: `Residue[chiPassive, {delta, I gammaStar}] - Astar` — native pole operation at the explicit pole `delta = I gamma_*`; the `.py` (L19-28) only does `together`-algebra equality and never extracts a residue. Confirmed present.
- **Line-shape split mechanism + order.** `.wl` L121-125 splits the GENERIC `delta`-dependent form FIRST via `ComplexExpand[Re/Im[...]]` then substitutes `delta -> r gamma`, with extra generic-form checks L130-131 the `.py` lacks. `.py` L72-73 substitutes `delta = r*gamma` FIRST then `sp.expand_complex(...).as_real_imag()`. Different mechanism and opposite order. Confirmed.

This is not the precedent shared-operation transliteration failure mode; no `mathematica_transliteration` defect. The re-author-from-transliteration is sufficient/independent (like VI.1-218).

## Non-tautology spot-check (load-bearing asserts)

Every load-bearing check is a derivation chain compared against an independently-written closed form, never `x==x`:
- `.py` L24/L28 `chi_lin - chi_bw`, passive form; L52 `diff(D_Pi,Pi)+Nfun`; L77-79 `re_r - re_expected`, `im_r - im_expected`, `re_r/im_r - r`; L91/L96 max + low-loss factorization; L109 `ratio_barrier - 1/r`; L143-147 survival round-trips.
- `.wl` mirrors via `expectZero[name, expr]` (L29-33), which `FullSimplify`s the residual and calls `fail[...]→Exit[1]` if nonzero — so a wrong derivation produces a nonzero residual and a nonzero exit. None of the LHS expressions is identically the RHS. Confirmed non-tautological.

## Exec log assessment

**SymPy:** exit=0 (`# exit_code: 0`, log L51). Notable lines:
- `Verified exact Stage 220 derivative identity: dD_Pi/dPi = -N(omega)` (L13-14)
- `|Re|/|Im| = r` (L23), `P_abs / (omega |U_disp|) = 1/r` (L32)
- probe slice deterministic: `Re chi = 1.05`, `|U_disp| = 13.125`, `required |A|/gamma = 0.204` (L43-50)

**Mathematica:** exit=0 (`# exit_code: 0`, log L94). Notable lines:
- `PASS: A_* is the residue at delta=I gamma_*` (L19) — the `.wl`-only Residue check
- `PASS: N is a perfect square` (L28), `PASS: dD_Pi/dPi + N(omega)` (L30) — the derived-N §II checks
- all 17 PASS lines present; final `All Stage 221 identities ... and survival-window relations verified.` (L93)
- banner correctly prints `STAGE 221` (L8), i.e. the stale `STAGE 204` banner the auditor flagged has been refreshed by the orchestrator re-run.

Both logs are freshly timestamped `2026-06-09T19:21:18-06:00` (sympy L4, mma L4) and ran via `timeout 600`.

**Output freshness:** the orchestrator re-ran both engines directly post-audit; the two exec logs are the freshness evidence (both exit 0, deterministic). The auditor's only blemish — the committed `.txt` STALE `STAGE 204` banner — is now resolved in the fresh Mathematica log (banner reads `STAGE 221`). Per the prompt, the committed `.txt` mtime is NOT a failure condition here. No content disagreement between engines.

## Material-change assessment

`material_change`: false. No script was edited (empty diff, both files unchanged from HEAD), so no derived result changed; the only refresh was the cosmetic banner in the committed Mathematica `.txt`, which carries no downstream-consumable value. Downstream units (>221) are unaffected by this verification.

## Side observations (non-blocking)

- `.wl` `$Assumptions` adds `eta <= 1` (L58) beyond the `.py` (which only declares `eta` positive). This is a tighter — not contradictory — slice consistent with the low-loss regime (`eta` small), and does not weaken any load-bearing identity; the shared identities are polynomial/rational and hold for all positive `eta`. Non-blocking.
- The auditor's assertion-inventory line citations are approximate (e.g. cites `.py` L143,144-147 for the survival window; actual asserts are L143 and L144-147). Cosmetic; the checks themselves are present and substantive.

## Verdict justification

`verified`. The report carries 0 findings, so there were no fixes to confirm — instead I confirmed the auditor's affirmative `clean` verdict holds at the checkpoint bar: both engines exit 0 with the full set of substantive PASS lines and the final all-verified line; every load-bearing assert compares a derivation chain to an independently-written closed form (non-tautological, `Exit[1]` on nonzero); and the re-authored `.wl` is genuinely independent of the `.py` on all three load-bearing quantities — most decisively the §II derivative identity, where `N` is a native-`D[]` OUTPUT on the Mathematica side and a posited INPUT on the SymPy side, plus the `.wl`-only `Residue[]` residue extraction and the different-mechanism/different-order line-shape split. The empty diff and unchanged scripts are the legitimate signature of a zero-correction batch, and the stale `STAGE 204` banner the auditor flagged is resolved in the fresh Mathematica exec log (`STAGE 221`). No regressions, no transliteration missed, no tautology.

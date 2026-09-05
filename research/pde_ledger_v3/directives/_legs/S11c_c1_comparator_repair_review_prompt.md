# Independent review — S11c-c1 comparator REPAIR (SCRIPT; ablate it)

You are an independent, adversarial reviewer of a Codex-written **repair** to a measurement instrument that was
already re-reviewed SOUND. The repair (diff `7141e6ad` → working tree) makes exactly two changes to
`research/pde_ledger_v3/scripts/S11c_c1_cross_engine_comparator.py`:
- **R1** — added `"cS0": "c_s0"` to the `C1_BARE_SYMBOL` spelling map (a bare→bare mechanical fold);
- **R2** — added a `hold_multi_range` rewrite in the **c1** `preprocess_wl` that turns a MULTI-range
  `Inactive[Integrate][…]` into `HeldInactiveIntegrate[…]` before the inherited one-range parser sees it.
Plus test-file changes (a replaced sound-speed test + new tests).

Your job: confirm the repair does exactly R1+R2, changes NO soundness behaviour, and introduces NO
false-agreement / hidden-disagreement / broken-seal path. ⛔ A prose judgement is worth nothing — ABLATE a
`/tmp` copy and report LITERAL stdout. Save every ablation script + output to named /tmp paths and report them.
Review against the FROZEN CONTRACT `research/pde_ledger_v3/directives/S11c_c1_SHARED_PHYSICS.md:580-587`, not the
build/repair directive.

Working dir: `/var/projects/toy_physics` (branch `ledger-v3-rebuild`). ⚠ The `.out` are large (~82/91 MB) —
sample with grep/awk/cut, never full-parse; `datalad get` both first; wrap comparator runs in `timeout 900` and
use `--family` to bound memory; ⛔ never two heavy runs at once; a full-run OOM is an operational deferral to
report, not a defect. First run `python3 -m pytest research/pde_ledger_v3/scripts/test_S11c_c1_cross_engine_comparator.py -q`
and report the literal result. ⛔ Copy the script to /tmp and ablate the COPY; never touch the working tree.

## Checks (CLEAN / DEFECT, each with the ablation command + literal stdout)
1. **Diff is exactly R1+R2.** `git diff 7141e6ad -- <comparator>` — confirm the ONLY comparator changes are the
   `cS0` map entry and the `hold_multi_range` preprocessor addition. Grep the diff for any change to: the 5
   seals (two-momentum `qOut`, μ_θ, ω-assumption, background density, regime/parity), the `raw_control_case`
   whitelist, `make_key`/axis typing, the three-valued residual objects, the accounting/inventory/scope
   emitters, exit semantics. Any change beyond R1+R2 is a DEFECT.
2. **R1 activation + injectivity + no-arg-strip.** Confirm (a) with the full 12-symbol required PY vocabulary
   present, a WL `cS0` folds to `c_s0` so a synthetic `c_s0`/`cS0` pair residuals to `Integer(0)`; (b) with the
   required vocabulary ABSENT, `cS0` does NOT fold (stays `Symbol('cS0')` — the all-or-none gate); (c) injectivity
   reports 0 collisions with the 12th entry; (d) adding `cS0` opened NO applied-head strip — the seal tests
   (`qOut[...]`, `rhoBrBgRho4Constant[...]`, `w1Profile/w1Jet[...]` stay applied) still pass. ⛔ Verify `cS0` and
   `c_s0` are the SAME physical quantity (the bulk sound speed) on both streams — if `cS0` ever named something
   else, folding it would be a false-agreement DEFECT.
3. **R2 held-parse of the REAL form + inherited path preserved.** Confirm (a) a 4-arg triple-range
   `Inactive[Integrate][f[xOne,xTwo,xThree],{xOne,-Infinity,Infinity},{xTwo,-Infinity,Infinity},{xThree,-Infinity,Infinity}]`
   and the nested `Inactive[Limit][<that>, outwardDistance -> Infinity]` parse to HELD `AppliedUndef` heads
   (not evaluated, not a native value, not `PARSE_FAILED`), keeping all args; (b) ⛔ the inherited ONE-range
   `Inactive[Integrate][f,{w,0,2}]` STILL canonicalizes to `BOUND_INTEGRAL` (binder+bounds) — R2 did not widen or
   break it; (c) `Inactive[Equal]→HeldEqual` still works; (d) on the real streams the ENERGY families are still
   UNPAIRED (join=0, SCENARIO axis mismatch) so R2 moved NO residual — it only converts `<PARSE_FAILED>` energy
   leaves to displayed held operands (the parse-failure COUNTS drop; no residual changes). Run the 3 ENERGY
   families (`--family ENERGY_RESIDUAL` etc.) and show the before/after parse_failed accounting on your /tmp copy.
4. **Seals still load-bearing (regression).** Pick one previously-surfaced seal residual (e.g. DTN_KERNEL's
   `qOut` two-momentum, or FACE_RESPONSE's μ_θ) and confirm it is STILL nonzero/surfaced after the repair — the
   repair must not have accidentally reconciled a seal. A name-map ablation on that seal must still move the
   residual.
5. **No new false-agreement.** Is there ANY family where the repair now prints a joined `Integer(0)` that the
   baseline did not, other than the intended `cS0`/`c_s0` collapse in a genuinely-same-quantity leaf? Hunt for a
   spurious new zero.

## Output
Per check CLEAN/DEFECT with the ablation command, literal stdout, and your /tmp script/output paths. Rank real
defects most-severe first; a clean pass with citations is equally useful. ⛔ Do NOT propose an expected measured
residual — a physics-bearing difference stays a SURFACED residual for post-run adjudication.

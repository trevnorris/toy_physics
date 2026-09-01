# Independent build review — S11c-b #89a WL §3a basis repair (Codex-written engine)

You are one of two independent legs reviewing a **Codex-written** modification to a Wolfram engine. Derive
the physics yourself from first principles, then judge the engine against your derivation. A prose derivation
is worth nothing here — **every physics claim you make must carry a script you ran and its literal stdout,
saved to a named absolute path you report.** Without the script+stdout, the claim is discarded.

## Artifact (working tree, uncommitted)
- Engine: `research/pde_ledger_v3/mathematica/S11c_b_brane_operator_mathematica_audit.wl` (just modified).
- Its regenerated output: `research/pde_ledger_v3/mathematica/out/S11c_b_brane_operator_mathematica_audit.out`
  (~156 MB — ⛔ NEVER `cat` it; `grep -a -o 'TAG…'` specific tags only).
- The build directive it must satisfy: `research/pde_ledger_v3/directives/S11c_b_89_wl_3a_repair_directive.md`.
  Read it — but judge the ENGINE against the PHYSICS and the SPEC, not against the directive's prose.
- Spec: `research/pde_ledger_v3/directives/S11c_b_SHARED_PHYSICS.md` (§1c symmetry group ~109–122, §3a ~245–265).

## What the build was supposed to do (its claim)
Replace the engine's hand-coded 8-per-source "spurion" new-invariant family with the **complete**
symmetry-allowed family, impose the exact thickness map on the new invariants, and let the existing
incremental-`MatrixRank` selection compute the count. Scope is the §3a **basis only**: the build should
**not** repair the operator freeze (the slab operator still profiles before the EL variation — a separate,
deferred defect) but should **emit a diagnostic** documenting it. The completed count was **withheld** from
the builder; the only public target was the restrict-to-original-8 regression `26`.

## Derive it yourself first (SymPy or any CAS — NOT Mathematica; save script+stdout)
1. **The complete family.** Build the complete set of scalar field-bilinears: **one** spurion factor (a
   vector `s`) × a quadratic in the DOF data `{u, ∇u (the full rank-2 shear ∂u_i/∂x_j; its trace is divU),
   θ, ∇θ, e_W, ∇e_W}`, contracting **all** index slots to a scalar using the metric δ only — **Kronecker
   contractions, O(3)** (rotations AND parity). ⛔ Exclude Levi-Civita (ε) contractions — they are
   parity-odd pseudoscalars and S11b parity / `w→−w` are unbroken. Impose the exact thickness map
   `e_W → (W_0/W_bg) e_W` before the rank test. Compute the raw independent **rank per source**, and the
   two-source total with the 10 carried uniform invariants. Save the family and the ranks.
2. **Confirm the frozen subset.** Transliterate the engine's original 8 forms (`newInvariantExpressions`,
   the committed pre-build version — recover it with `git show HEAD:…mathematica_audit.wl` and read
   `newInvariantExpressions`) and confirm they lie **within the span** of your complete family and span a
   strict subspace. Compute the restricted rank (with the 10 uniform) — this is the public `26`.
3. **The operator freeze (why the split).** In an abstract EL-signature space, confirm that freezing the
   background second jet in the EL (Hessian → 0) drops the rank relative to the live EL — i.e. the Hessian
   is genuine, non-absorbable operator content. This is what the engine must **document, not repair**.

## Then check the engine against your derivation (grep the .out tags; ablate a /tmp COPY of the .wl)
4. **Count matches your derivation AND spans the same space.** `grep -a -o 'ENERGY_BASIS_COUNT[^}]*}'`
   → does `COUNT_OPERAND` equal your complete two-source total, for BOTH anchorings? ⛔ Verify SPANS not just
   the integer: extract the engine's emitted new-invariant forms (`ENERGY_BASIS_NEW_INVARIANTS`) and confirm
   they span the **same** space as your complete family (rank of the union = rank of each). A matching
   integer over a different span is a defect (`26=26` was two undercomplete mechanisms).
5. **Count is COMPUTED, not hard-coded.** Grep the engine source for any hard-coded completed count or pasted
   invariant list; confirm the count is `Length` of the rank-selected indices from `MatrixRank`. Report any
   `Integer(<completed count>)` literal.
6. **⛔⛔ FORM ABLATION (MANDATORY — the only thing that has ever caught the worst defect).** On a `/tmp`
   COPY of the `.wl`, change the STRUCTURE of the enumeration and report the literal rank change:
   (a) restrict to the original 8 forms → count must move to `26` (and the emitted family back to 8/source);
   (b) ADD one Levi-Civita (ε) pseudoscalar candidate → the raw rank must RISE (proving the parity exclusion
   is load-bearing, not vacuous); (c) DROP one independent Kronecker contraction the shear admits → the rank
   must FALL (proving completeness is load-bearing). A count byte-identical under (a)/(c) means the repair
   did nothing. ⭐ Extract just the enumeration + `independentRepresentativeIndices` + profile functions into
   a MINIMAL harness — ⛔ do NOT re-run the full 156 MB engine.
7. **Thickness map imposed.** Confirm the new invariants use the mapped local field, not raw `ew`. Ablate:
   build the new family with raw `ew` vs mapped `localEw` and report the nonzero residual / grade change.
8. **Operator-freeze diagnostic is honest.** Confirm the engine emits the operator-row **frozen-EL vs
   live-EL rank + their difference** (a nonzero difference documenting the freeze), and that it does NOT
   assert them equal and does NOT claim the operator is repaired. Confirm the production operator
   (`evaluatedModel`) indeed still profiles before the variation (scope honored).
9. **Controls non-tautological + able to fail.** Ablate the load-bearing controls (per-anchoring independent
   rank; dimension derived-then-compared; coefficient control) — each must move when you break its input.
   Report any control whose residual is zero for any input (corollary 3), any emit conditioned on a payload
   value (corollary 4), any tag naming its value. Confirm the existing controls (`REP_INVARIANCE`,
   `INDEPENDENCE`, `FORM`, `UNIFORM`, `HOMOGENEITY`) still emit and pass.

## Physics filter
Report a finding only if it catches a way the physics/count/operator could be wrong, or a control that
cannot fail. Do not report "wrong on a different input."

## Operational constraints (MANDATORY — this artifact spawns Mathematica kernels)
- ⛔ Wrap EVERY Mathematica kernel run in `timeout 600`. A 600 s hit is a FAILED ablation — report it and
  move on. ⛔ NEVER raise the timeout. ⛔ Run at most ONE kernel at a time (the licence has TWO seats and the
  sibling leg may hold one).
- ⛔ Copy the `.wl` to `/tmp` and ablate the COPY. ⛔ Never modify the working tree. Use a MINIMAL harness
  (extract the enumeration/rank/profile functions); ⛔ do not regenerate the full `.out`.
- ⭐ Save every derivation and ablation script AND its literal stdout to named absolute paths under
  `~/.s11_build/S11c_b_89a_wl_buildleg_<yourname>/`, and report those paths with each finding.
- If a kernel orphans or memory runs low, `ps -eo pid,rss,pcpu,comm --sort=-rss | head` and `free -h`; do not
  relaunch blindly.

## Report
A numbered list: `severity — file:line — problem — the script+stdout path that shows it`. Severity ∈
{BLOCKER, SHOULD-FIX, NIT}. Then `VERDICT: CLEAR` or `VERDICT: N findings (B blockers)`. A leg that finds
nothing is weak evidence — name the two or three structural things you actively tried to break (with the
ablation that failed to break them) so the orchestrator can weigh the clearance.

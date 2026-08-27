# Independent build review — WL density-grounding fix (Codex-written)

The blind Wolfram engine of S11c-a was just modified so its **shifted-face density trace (T-e)** grounds
the density background in the supplied §2b representative instead of an undefined free-premise function.
Your job is to find any way the change is physically wrong, incomplete, or breaks something else.

## Artifact under review
The working-tree engine `research/pde_ledger_v3/mathematica/S11c_a_interface_geometry_mathematica_audit.wl`
and its regenerated transcript `…/mathematica/out/S11c_a_interface_geometry_mathematica_audit.out`.
The change is `git diff` against HEAD (`git show HEAD:<path>` is the pre-fix baseline).

## What the change is supposed to achieve
- The T-e density trace's background face value and **normal derivative** must be obtained by an actual
  differentiation of the supplied §2b density representative (in-plane, `w`-independent,
  representative-aware) — not from any free-premise or `w`-argumented abstract function.
- Because that background is representative-parameterised, the shifted-trace family must be keyed by the
  density representative (a DENSITY axis), in **every** place the engine builds a shifted trace: the
  primary `S11CA_FACE_SHIFT`, the form-ablation control's shifted-trace slice, and the uniform-limit
  regression's shifted-trace slice.
- The **perturbation** part of every traced field must be unchanged (evaluated at the background face
  `s·W_bg/2`, retaining first-shape-order dependence). Pressure/velocity/current backgrounds (already
  zero) must be unaffected.

## What to check — each with a script and its literal stdout
1. **Derive the correct §3c grounded density trace yourself** from the spec (`S11c_a_SHARED_PHYSICS.md`
   §1b/§2b/§2d/§3c). What must the background normal derivative be, and does it *emerge* from
   differentiating the representative rather than being asserted? Write a runnable snippet (SymPy is
   fine), save it and its stdout to absolute paths, and report both. A prose derivation is discarded.
2. **FORM ABLATION (mandatory).** Take the density background the engine now uses and replace it with an
   explicitly `w`-dependent test function; show that the emitted density-trace **shift term reappears**
   (i.e. the shift is genuinely wired through a real `D[·,w]`, not a hardcoded result). A **focused**
   ablation of the trace-construction logic is sufficient and preferred — you do **not** need to re-run
   the whole engine. Report the literal before/after of one density operand. If a coefficient rescale (not
   a form change) is all you do, that only tests arithmetic — you must change the *form* of the background.
3. **Independence / not-a-hardcode.** Confirm the grounded result is not a hand-written `0`: if you change
   *which* representative is used, does the emitted background face value change accordingly? Show it.
4. **Collateral.** Diff the regenerated transcript against the committed baseline (`git show HEAD:<path>`).
   Confirm the only families that changed are the shifted-trace-bearing ones (primary FACE_SHIFT + its
   form-control slice + its uniform-limit slice) — and that they gained the DENSITY axis. Report any other
   tag family whose payload moved; that would be a collateral defect.
5. **All three sites.** Confirm the density axis + grounded background are present in all three
   shifted-trace emissions, not just the primary one. Name the loop/keys for each.
6. **Blindness.** Confirm the engine still imports nothing and re-derives from the spec; the fix must not
   transcribe the sibling (SymPy) engine.

## Mathematica operational constraints (if you spawn any kernel)
- Wrap EVERY kernel run in `timeout 600`. A 600s hit is a FAILED ablation — report it and move on.
- NEVER raise the timeout, and NEVER run more than one kernel at a time (the licence has TWO seats).
- Copy the artifact to /tmp and ablate the COPY. NEVER modify the working tree.
- Prefer a focused snippet over a full-engine re-run.

## Physics filter
Report a finding only if it catches a way the shifted density trace is physically wrong, §3c-noncompliant,
incompletely applied (a site left ungrounded/unkeyed), or a way the change corrupts another family.
"It could be tidier" is not a finding.

## Return
A short verdict: is the grounded density trace §3c-correct and complete across all three sites, does the
form control bite, is anything else disturbed? Cite file:line and paste your scripts' literal stdout.

# Independent physics review — S11c-a WL background-current fix (a Codex-written engine change)

## Artifact under review
The change to the Wolfram engine
`research/pde_ledger_v3/mathematica/S11c_a_interface_geometry_mathematica_audit.wl` at the bulk-current field
definitions (~lines 435-448). Review the CHANGE (diff vs the committed baseline `a7459cb8`) and the engine's
behaviour under it. Do NOT review the directive; form your own view from the spec.

## Read the source of truth first, form your own view, then the change
- Spec `research/pde_ledger_v3/directives/S11c_a_SHARED_PHYSICS.md`: §1 (drain `v_bulk_normal_0`), §1b
  (`j = ρ_4D v_bulk`, continuity), §2d (`𝔅⁰`; `V_s⁰=J_s⁰=0`), §3c ("the rest-frame background current
  `ρ_4D⁰v_bulk⁰` vanishes"; "none may be introduced as a free premise").
- The change: the background bulk current was previously the free symbols `currentWBackground` /
  `currentXBackground{i}`; it is now built as `bulkCurrentZero = rhoBulkZero · bulkVelocityZero`, with
  `currentWZero = Last[bulkCurrentZero]`, `currentXZero[i_Integer] = bulkCurrentZero[[i]]`. The engine's
  `bulkVelocityZero` is the zero 4-vector (lines 115-119); `rhoBulkZero` is the nonzero background density.

## What to check — and the ablations to RUN (this is a SCRIPT artifact: derive + ablate, don't just read)
1. **Physics (derive it yourself).** From §1b/§2d/§3c, is the rest-frame background bulk current
   `j⁰ = ρ_4D⁰·v_bulk⁰` required to be zero in this scope, given `v_bulk⁰=0`? Is it correct that the density
   background stays nonzero (`rhoBulkZero`) while the current background vanishes? Quote the spec lines. If you
   make a derivation claim, back it with a CAS script + literal stdout at named /tmp paths.
2. **⛔ MANDATORY FORM ABLATION — the BACKGROUND-VELOCITY PROBE.** This is the load-bearing check and it must
   be RUN, not reasoned. On a COPY, replace the background bulk velocity `bulkVelocityZero` with a formal
   NONZERO probe (e.g. `{vp1,vp2,vp3,vpW}` or `{v[1],v[2],v[3],v[4]}` — symbolic, nonzero) and evaluate the
   background current `currentWZero`/`currentXZero[i]`. A FAITHFUL `j⁰ = ρ⁰·v⁰` construction makes the
   background current TRACK the probe: `currentWZero → rhoBulkZero·vpW`, `currentXZero[i] → rhoBulkZero·vp_i`
   (nonzero). If instead it stays zero under the probe, the construction is a disguised `:= 0` assertion (a
   post-hoc zeroing), which would FAIL this check. Report the literal probe output. You can do this on the
   ISOLATED current definitions (you do NOT need to run the whole 47MB engine for this ablation).
   ⭐ Also do the reverse: with the real `bulkVelocityZero` (zero), confirm the background current evaluates to
   0 and no `Part::pkspec1` message appears (concrete indices AND a symbolic-index call staying inert).
3. **No regression / correct scope.** Confirm the change touches ONLY the background current: the perturbation
   current (`currentWWave`/`currentXWave`) is unchanged, the nonzero background density is preserved, and the
   pure face-geometry objects are untouched. If you run the full engine output, confirm the non-current-
   consuming tags are unchanged vs baseline and that ZERO `currentWBackground`/`currentXBackground` symbols
   remain. (A clean full run `/tmp/s11ca_wl_bgcurrent_full.out` may be provided — you may use it, but verify
   claims about it yourself.)

## Mathematica ablation constraints — put these in effect (both legs identical)
⛔ Do NOT run the full audit engine (`wolframscript -file …audit.wl`): it takes ~14 min and has ballooned to
14 GB when mis-built — it EXCEEDS the 600 s ablation budget. Do the velocity-probe on the ISOLATED current
definitions instead (a few-second `wolframscript -code '…'` reproduction, pasting the engine's edited
`bulkCurrentZero`/`currentWZero`/`currentXZero` forms — read them from the `.wl` yourself). ⛔ Wrap any kernel
run in `timeout 600`; a hit is a FAILED ablation — report and move on. ⛔ ONE kernel at a time (2-seat licence).
⛔ Copy anything you execute to /tmp and ablate the COPY; never modify the working tree. ⭐ Save every ablation
script AND its literal stdout to named absolute paths and report them.
For OUTPUT-LEVEL claims (all 40 tags present, 0 `current*Background`, unrelated tags unchanged), a clean full
run is already available at `/tmp/s11ca_wl_bgcurrent_full.out` (and the committed baseline at
`research/pde_ledger_v3/mathematica/out/S11c_a_interface_geometry_mathematica_audit.out`) — grep/compare those
rather than regenerating; verify any claim you rely on yourself.

## Physics filter
Report a finding only if it catches a way the physics could be wrong (a background current that should/should
not be there; a broken perturbation current or density; a regression in another object). Not style, not "wrong
on a different input."

## Return
Your independent physics derivation (with script+stdout), the velocity-probe ablation result (does the
background current track the probe? does it vanish with the real zero velocity? any `Part::pkspec1`?), the
no-regression check, and a plain verdict: is the change faithful to §3c/§1b, or not?

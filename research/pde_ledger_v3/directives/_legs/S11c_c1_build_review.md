# Independent physics review — the S11c-c1 SymPy engine (Codex-written SCRIPT)

## Artifact
`/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_c1_bulk_closure_sympy_audit.py`
(Codex-written per the migrated build directive; it runs to exit 0 and writes the own-rows-delta export
`/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_c1_exports.py`, 44 rows / ~906 KB.)

## What to check
S11c-c1 solves the perturbed curved two-face outgoing bulk problem for the nonlocal DtN/impedance operator and
closes the permeable face response. This engine is one of two blind engines (the other is a from-scratch Wolfram
engine, not built yet); its job is to COMPUTE and PRINT objects, deciding nothing. Verify the physics is right
and the checks are real — AND that the export is a correct bind-closure own-rows delta (the migration's point).

## What you are handed (read these; ⛔ you are NOT given the build directive — derive from the spec)
- The engine script above.
- The PHYSICS AUTHORITY (derive independently from this, not from the script):
  `/var/projects/toy_physics/research/pde_ledger_v3/directives/S11c_c1_SHARED_PHYSICS.md` (§§0–8).
- The export design + fold module (for the delta-export checks):
  `/var/projects/toy_physics/research/pde_ledger_v3/directives/export_ledger_bind_closure_design.md`,
  `/var/projects/toy_physics/research/pde_ledger_v3/scripts/ledger_fold.py`.
- The frozen base ledger `/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_b_exports.py` (plain git,
  2441 rows). The sibling substrate specs it cites (S11c-a, S11b) are in the same directives/ dir.

## Required method — this is a SCRIPT
Derive independently. **Write your own derivation script BEFORE opening the engine, and save both the script and
its literal stdout to named absolute paths** — a prose "I re-derived it and got X" is discarded (CLAUDE.md).
Ablate every load-bearing check and report its literal output; code-reading alone has repeatedly missed real
defects. Probe for: a value verified with the definition that produced it (`c ≔ √x` then asserting `c²−x=0`); a
conclusion emitted as an unconditional literal; a check whose expected value lives inside the artifact; an
`assert` that precedes the value it guards (report it).

⛔⛔ **A FORM ABLATION IS MANDATORY — the only thing that has ever caught the worst defect.** Change the STRUCTURE
of a load-bearing object and re-run, reporting the literal diff. For THIS engine, the load-bearing structure is
the **two-momentum DtN kernel carrying both branch legs** `q_out(ω,k)` and `q_out(ω,k′)`: collapse the two legs
into one (`k′→k`, the WKB/local freeze the spec forbids under rule 17), or flip one branch sign, and confirm the
DtN kernel, the response resolvent `[I+(Λ_A/ρ_m²)Z]⁻¹`, and the §5e branch-liveness residuals MOVE. A byte-identical
output under this ablation = the object was never live (the rule-17 freeze, or a hand-typed payload). A COEFFICIENT
rescale tests arithmetic; only a FORM change tests physics.

## The physics to verify (derive each yourself; demand computation)
1. **The DtN is a two-momentum nonlocal OPERATOR** carrying BOTH legs `q_out(ω,k)`, `q_out(ω,k′)` explicitly
   (kernel `Z_s^{(1)}(ω;k,k′)` with `Ŵ_bg(k−k′)`). ⛔ NOT a single-`k` multiplier, ⛔ NOT a one-leg left-quantized
   `a(x,k)=W_bg(x)σ(q(k))` — both delete the propagating↔evanescent mode mixing (the rule-17 freeze). Confirm the
   engine's kernel is genuinely two-legged (ablate: does freezing one leg change it?). Check the rigid-shift
   (`k=k′`, `Ŵ(0)`) cancellation is emitted as a computed residual, not asserted.
2. **The permeable face response** `(δp_s,J_s,t_s)(V_s,μ_θ)` via the OPERATOR inverse `[I+(Λ_A/ρ_m²)Z]⁻¹`, not a
   scalar division. `Λ_X` appears ONLY in the traction `t_s`, not in `J_s`. `μ_θ` is carried OPAQUE (⛔ never
   expanded into slab DOFs — that is c2's work).
3. **The three-object dissipation audit, all distinct:** (a) bulk-radiation Hermitian part `H_a[Z]=(Z+Z^†)/2`
   under the true-area pairing; (b) the two-port permeable-port Hermitian form on `(V,μ_s)→(δp+Λ_X𝒜,J)`; (c) the
   INDEPENDENT energy balance whose face operand is the true-area **traction** pairing `½Re Σ∫a_s t_s·v*`
   (⛔ NOT `½Re(δp·V*)`, which equals the bulk flux and never sees `t_s`) vs the outgoing far-field Poynting flux.
   Confirm (c) is a genuine second route: one-sided corrupt the `t_s` sign and confirm ONLY the face operand
   moves (the residual goes nonzero). If it is `½Re(δp·V*)` in disguise, that is a defect.
4. **Rule-17 liveness:** `q_out`, the regime discriminant, `ω`, and the background `W_bg(x)` and its jets must be
   kept LIVE and differentiated — never frozen to a constant to proceed. The §5c uniform-limit and §5d zero-jet
   are regressions (§5d target is the UNMODIFIED S11b `z_impermeable`, ⛔ NOT B0b re-solved at a finite gap).
5. **No conclusion emitted as prose or as a hand-typed CAS object** with no data dependency (delete the derivation
   and see if the payload moves). Tag NAMES must not encode the answer (parity/sign/regime).

## The delta-export to verify (the migration's headline — run these)
- Load the export and confirm it is an **own-rows delta**: `python3 -c "import ledger_fold as lf; d,_=lf.load_model('S11c_c1_exports.py'); print(len(d), sorted(d))"` from the scripts/ dir. Expect ~44 rows — the 5
  model-level objects + only NEW `s11cc1_`-prefixed symbols — ⛔ NOT the 2441-row inherited model.
- Confirm the fold `load_model('S11c_b_exports.py','S11c_c1_exports.py')` adds the delta with ZERO overwrites, the
  bare `face_response` stays the S11b row (F9c predecessor preserved), and c1's own is `s11c_c1_face_response`.
- Confirm the engine's in-run §D3 guard is REAL, not theater: it calls `check_consumer` + `assert_lookups_equal_manifest` (the exact-44 smoke-test) + `assert_delta_is_minimal`. **Ablate it**: on a /tmp copy, add a bound
  fold lookup the manifest omits (or drop a declared key) and confirm the smoke-test RAISES — a guard that passes
  regardless is theater. Confirm IMPORT_KEYS is the exact 44 and the engine does not re-originate an inherited row.

## Physics filter
Report a finding only if it catches a way the physics could be wrong, a check that cannot fail, a rule-17 freeze,
or a broken delta/guard. Do not report "would be wrong on a different input" or style.

## Ablation sandbox
Copy the engine to /tmp and ablate the COPY; ⛔ never modify the working tree. This engine is pure SymPy (no
Mathematica kernel) — no seat/timeout constraints, but note it takes a few minutes to run; if you run it, capture
stdout to an absolute path and grep for the tags you need rather than holding it all.

# Premises behind `S11c_b_89_sympy_3a_repair_directive.md` (measurements)

Rule 2: every premise the directive asserts about the spec or the engine carries the command that produced
it, with literal output. This file verifies the directive's PREMISES only. ⛔ It does **not** state the
corrected basis count or the corrected operator (the withheld acceptance criterion, rule 5); those are the
build's computed outputs, diffed against the reference outside the build.

Run from repo root `/var/projects/toy_physics` at HEAD `7584c16a`.

---

## P1 — the spec explicitly requires retaining the generated higher background jet (§3a)

    $ sed -n '251,256p' research/pde_ledger_v3/directives/S11c_b_SHARED_PHYSICS.md

Literal: "A divergence or variation may generate a higher spatial derivative of a background profile, and
that term is retained at the background-bookkeeper order of its originating factor: a second spatial
derivative of `W_bg` is still first order in background amplitude, `O(η)` / `O(σ_W)`, and is not dropped.
Independence is judged as field bilinears with B1's constraint not applied; carry every independent
invariant with its own free symbolic coefficient."

⇒ Hessian retention is a spec requirement, not a design choice. The directive's §1 quote is faithful.

## P2 — §1d: the uniform quotient does not lift; the first-jet term is physics

    $ sed -n '163,168p' research/pde_ledger_v3/directives/S11c_b_SHARED_PHYSICS.md

Literal: "⛔ the uniform quotient does not lift trivially to variable coefficients: integrating by parts a
variable coefficient generates first-background-jet terms (`c∇·F ≡ −(∇c)·F` modulo a boundary term), so
representatives that were equivalent uniformly differ by a first-jet invariant that is physics in the
operator/kernel, not a representational identity."

## P3 — §3b: do not freeze a coefficient before differentiation

    $ sed -n '276,278p' research/pde_ledger_v3/directives/S11c_b_SHARED_PHYSICS.md

Literal: "Retain the full spatial dependence of every background coefficient (`μ_R,bg`, `W_bg`,
`ρ_4D,bg⁰`, `ρ_br,bg⁰`, and the `Σ_E⁰` map) and its first jet; do not freeze a coefficient at its constant
binding before differentiation."

⇒ Spec-confirm verdict: §1d/§3a/§3b already mandate the Hessian-live construction. #89 is a
spec-COMPLIANCE repair; no spec change is warranted.

## P4 — Path A basis quotient is DOF-only (spurion frozen)

    $ grep -nE 'def basis_euler_signatures|def quotient_independent_indices|^basis_fields' \
        research/pde_ledger_v3/scripts/S11c_b_brane_operator_sympy_audit.py

Literal:
    936:def basis_euler_signatures(
    972:def quotient_independent_indices(
    1025:basis_fields = tuple(

`basis_fields` (L1025–1028) lists only `bu/bG`, `btheta/bq`, `be/br` — the DOF abstract fields. The spurion
`bg` (L1014) is not among them, so `basis_dx` (L947–956) differentiates no spurion factor. Call site:

    $ sed -n '1270,1273p' research/pde_ledger_v3/scripts/S11c_b_brane_operator_sympy_audit.py

Literal:
    candidates = enumerate_new_candidates(g_vector)
    ...
    signatures = basis_euler_signatures(expressions, basis_fields)
    selected, omitted = quotient_independent_indices(expressions, signatures)

## P5 — Path B strong rows use the global frozen `dx`; `DERIVATIVE_MAP` stops at the first background jet

    $ sed -n '611,621p' research/pde_ledger_v3/scripts/S11c_b_brane_operator_sympy_audit.py

Literal:
    for i in DIRECTIONS:
        DERIVATIVE_MAP[i][W_bg] = grad_W[i]
        DERIVATIVE_MAP[i][mu_R_bg] = grad_mu[i]
    ...
    def dx(expression, direction):
        result = sp.Integer(0)
        for atom, derivative in DERIVATIVE_MAP[direction].items():
            if atom in expression.free_symbols:
                result += sp.diff(expression, atom) * derivative
        return sp.expand(result)

There is no `DERIVATIVE_MAP[i][grad_W[j]] = …` entry ⇒ `dx(grad_W[i]) = 0`. `operator_from_density`
(L1459) calls this global `dx` at L1473/L1482/L1489 for the U/θ/e_W rows.

## P6 — the correct Hessian-retaining pattern already exists (operator_dx / background_dx)

    $ sed -n '1850,1869p' research/pde_ledger_v3/scripts/S11c_b_brane_operator_sympy_audit.py

Literal (abridged to the added entries):
    def operator_dx(expression, direction):
        derivative_map = dict(DERIVATIVE_MAP[direction])
        for jet_direction in DIRECTIONS:
            ...
            derivative_map[grad_W[jet_direction]]  = sigma_W * w_profile_second / L_W
            derivative_map[grad_mu[jet_direction]] = mu_R * sigma_W * m_profile_second / (W0 * L_W)

`w_profile_second`/`m_profile_second` are the `w1_profile_d{i}d{j}`/`m1_profile_d{i}d{j}` second-jet atoms.
`background_dx` (L2122) repeats the pattern. So the retained-Hessian divergence is already implemented for
the coupling/admissibility routes; Path A and Path B are the inconsistent frozen ones.

## P7 — the frozen public regression target (SUPPLIED, committed, cross-engine agreed)

    $ grep -a S11CB_ENERGY_BASIS_COUNT \
        research/pde_ledger_v3/scripts/out/S11c_b_brane_operator_sympy_audit.out | head

Literal: both LAB_HELD and MATERIAL_ADVECTED report `('VALUE', Integer(26))`. This is the frozen value the
§4 frozen-switch regression must reproduce, and the only target value the directive exposes.

## P8 — partial scaffolding exists; the repair DOES need new scaffolding (corrected after the decision legs)

    $ grep -nE "w1_profile_d\{i\}d\{j\}|def ablated_jets|CONTROL_TASKS|EXPORT_PATH|def material_pullback" \
        research/pde_ledger_v3/scripts/S11c_b_brane_operator_sympy_audit.py

Literal (abridged):
    30:EXPORT_PATH = SCRIPT_DIR / "S11c_b_exports.py"
    39:CONTROL_TASKS = ("REP_INVARIANCE", "INDEPENDENCE", "FORM", "UNIFORM", "HOMOGENEITY")
    356:            f"w1_profile_d{i}d{j}",   # ledgered via bind_additional_inherited (≈L353–360)
    1218:def ablated_jets(source, direction):
    1348:def material_pullback(

Corrected reading (an earlier draft of this file over-claimed "no new scaffolding"): the `w1_profile_d{i}d{j}`
2nd-jet atoms are ledgered (≈L353–360), but the `m1_profile_d{i}d{j}` modulus 2nd jets are created **locally**
inside `operator_dx` (not ledgered), there is **no abstract-space `bg` 2nd-jet table**, and there are **no
3rd-jet atoms** anywhere. `ablated_jets` (L1218) zeroes only one direction of one source (does not collapse a
family). `EXPORT_PATH` (L30, written at L2930) is the tracked `S11c_b_exports.py` the run regenerates. So the
repair requires: a symmetric abstract `bg` Hessian table (divergence map only), shared ledgered `W`/`μ` jet
tables extended to the needed depth, a full-source jet-zero ablation, a dedicated Hessian-freeze switch, and
authorization to regenerate `S11c_b_exports.py`. These are folded into directive §0/§3/§4/§5.

## P9 — decision-leg fold (2 legs, one pass, rule 7; artifacts recorded)

The directive was reviewed by two independent legs before any builder:
`~/.s11_build/S11c_b_89_sympy_directive_codex.txt` (Codex, xhigh, 6 defects) and
`~/.s11_build/S11c_b_89_sympy_directive_grok.txt` (Grok, high, 2 defects; both a subset of Codex's, converging
on the same two). All eight were verified against the engine bytes (rule 13) and folded once into the revised
directive:
1. (Codex 1) `S11c_b_exports.py` is a tracked file the run rewrites ⇒ authorize its regeneration (§0).
2. (Codex 2) `material_pullback` (L1363–1364) is a THIRD frozen path (MATERIAL anchoring) ⇒ consumer C (§2/§3).
3. (Codex 3) the coupling cascade needs the 3rd jet `operator_dx` drops ⇒ consumer D + the jet tower (§2/§3).
4. (Codex 4 = Grok 2) `ablated_jets` is per-direction, does not collapse the family ⇒ full-source ablation (§5).
5. (Codex 5 = Grok's flagged trap) spurion must enter the divergence map, NOT the variational-field loop; need
   a symmetric abstract `bg` Hessian table + shared ledgered `W`/`μ` jet tables (§3).
6. (Codex 6) the live derivative makes `σ_W≥2` terms the projection removes ⇒ grade-project both sides; compare
   at the projected level; raw byte diff only as a diagnostic (§3/§4).
7. (Grok 1) the frozen-limit regression must use a dedicated Hessian-freeze switch (first jet LIVE), NOT
   `uniformize`/UNIFORM (which kills the first jet and cannot reach 26) (§4).
Per rule 7 this is one fold; the directive now goes to build (Codex), reviewed by its own two build legs. No
re-legging of the decision list (never iterate a decision list to green).

## Provenance / withholding note

The corrected basis count and corrected operator are the #86/#88 reference result
(`_measurements/S11c_b_86_reference_result.md`, `_measurements/S11c_b_88_blast_radius_result.md`). Those
records exist in the tree; per CLAUDE.md rule 12 blindness is enforced by ABSENCE of a matching acceptance
criterion in the directive and by the input-driven engine + form ablation, not by hiding files. The
directive states no live count; the frozen-switch `26` cannot leak the live value; a hard-coded count fails
the form ablation. That is the safety this repair relies on.

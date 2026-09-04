# Directive design review — S11c-c1 Wolfram engine build directive (decision-leg gate, rule 7)

## What you are reviewing

`research/pde_ledger_v3/directives/S11c_c1_wl_build_directive.md` — an **orchestrator-written** build directive
for the *blind Wolfram engine* of sub-step S11c-c1 (curved-interface bulk closure). No builder has run yet: this
is the rule-7 decision-leg gate that must pass **before** the Codex builder writes the `.wl`. Your job is to find
every way this directive would make the build wrong, make the cross-engine control vacuous, or leak an answer —
**now, while it is still cheap to fix.**

This directive governs the SECOND of two engines. The first (SymPy) engine is already built + reviewed + committed;
its directive is `directives/S11c_c1_sympy_build_directive.md`. The two engines exist so they can **disagree**: the
downstream T7 comparator measures the difference, and the disagreement is the measurement (`CLAUDE.md` rule 16).

## Sources of truth — read these and form your OWN view FIRST, then read the directive

⭐ Read the specs and form your own view of what the WL engine must compute and how it must stay blind **before**
you open the directive. Quote both sides (spec text vs directive text, with line numbers) for every finding.

- `directives/S11c_c1_SHARED_PHYSICS.md` — the **sole physics authority** both engines read (§§0–8: scope; the
  inherited setup; the curved-bulk problem; §3a two-momentum DtN operator; §3b permeable response + Fredholm loci
  + three-object dissipation audit; §4 the emit list; §5 controls; §6 method/locus-protocol/three-clauses; §7
  names + tag grammar; §8 supplied-vs-computed).
- `directives/S11c_a_SHARED_PHYSICS.md` — the S11c-a substrate the WL engine RE-DERIVES (its §§1–3 setup + §4
  T-a..T-i shape-derivative objects). `directives/S11b_SHARED_PHYSICS.md` — the bulk acoustics / radiation
  condition / branch object `q_out` (§1b, §2) and the flat B0b/B0c reference (§9 B0) the WL engine re-derives.
- `directives/S11c_c1_sympy_build_directive.md` — the SIBLING engine's directive. The two directives must wire the
  **same** spec obligations and pin the **same** naming/binding conventions (so the T7 join lands on the same
  objects), while differing correctly: SymPy imports via the fold + writes an own-rows delta; WL imports NOTHING
  and writes only a stdout tag stream.
- `research/pde_ledger_v3/directives/S9_export_chain_rebuild_directive.md:14-18` — the ONE blindness control
  ("each step's Wolfram engine imports nothing and re-derives independently"). `directives/S11c_decisions.md` `N8`
  (the frozen T7 comparator contract, inherited verbatim) and `N9` (the denylist stays CUT — blindness is
  absence, never a do-not-read list). `CLAUDE.md` rules 2, 5, 6, 7, 11, 12, 16, 17.
- Mechanical WL precedent (idioms only, NOT physics): `mathematica/S11c_b_brane_operator_mathematica_audit.wl`
  head; the master WL build directive `directives/S11c_b_wl_build_directive.md`.

## The questions this review must answer (report a finding for any "no")

Answer each with spec/directive quotes. A clean "yes" with the citation is as useful as a finding.

1. **Spec fidelity.** Does the directive faithfully wire every §0–§8 obligation without distorting one? In
   particular: §3a (the DtN is a two-momentum OPERATOR carrying BOTH legs `q_out(ω,k)`, `q_out(ω,k′)` — ⛔ not a
   single-`k` multiplier, ⛔ not a one-leg left-quantized `a(x,k)`); §3b (the response is an operator inverse
   `[I+(Λ_A/ρ_m²)Z]⁻¹`, ⛔ not scalar division; the three dissipation objects are DISTINCT, the energy route's
   face operand is the true-area traction pairing and its bulk operand is the far-field Poynting flux at
   `|w|→∞`, ⛔ not `½Re(δp·V*)`); the `Λ_X`-only-in-traction placement (§1d).
2. **Leak / rule 5.** Does the directive state, imply, or pre-register ANY computed value, sign, parity,
   regime-pair behaviour, locus, dissipation sign, or residual value? Scrutinize the sentence claiming `μ_θ`
   "cancels out of every c1 cross-engine residual" (directive §"blindness") — is that a legitimate naming
   consequence, or does it pre-register an agreement the comparator is supposed to MEASURE (rule 16 /
   "don't design away the disagreement")? Flag any acceptance criterion that references an expected value.
3. **Blindness soundness (rule 12 / N9).** Is "imports nothing" enforced by ABSENCE (the engine has no import
   machinery), not by a prohibition/denylist? Does the "re-derive the S11c-a substrate T-a..T-i by
   shape-differentiating the S11c-a setup" instruction actually keep the WL engine blind — i.e. can T-a..T-i be
   re-derived from the S11c-a spec's SUPPLIED setup (§§1–3) WITHOUT reading an S11c-a *result*? Or does the c1
   spec / this directive hand the WL engine a T-substrate *result* (e.g. the §2a outward-normal form) in a way
   that makes the "re-derivation" a transcription, collapsing the cross-engine control on the substrate to
   vacuous? Is the copy-to-empty-dir Acceptance check a sound structural proof of blindness?
4. **`μ_θ`-opaque handling.** The directive treats `μ_θ` as a single opaque symbol carried identically by both
   engines (so it is a supplied, unfalsifiable operand, not a measurement target). Is that faithful to spec
   §0/§3b (`μ_θ` kept symbolic, its expansion reserved for c2) and to the SymPy directive (which binds
   `LEDGER['mu_theta_operator']` opaque)? Is the assumption that the SymPy `mu_theta_operator` row is a SINGLE
   atom (so the two engines' `μ_θ` coincide under transliteration) correct — or could SymPy's imported `μ_θ` be a
   composite expression, so that c1 objects containing `μ_θ` would mismatch trivially cross-engine? If the latter,
   the directive has a real defect.
5. **Emit-list completeness.** Is the directive's emit list EXACTLY the spec §4 object list — every §3a/§3b
   primary AND every §5/§6 control operand+residual — with nothing dropped, nothing invented, and the shared
   grammar `WL_S11CC1_<QUANTITY>` (one tag per named object; cases as keyed Association)? Cross-check tag-for-tag
   against spec §4 and §3a/§3b's `⇒` blocks.
6. **T7-joinability.** Do the pinned naming/binding conventions (the two-momentum leg convention `k`=output /
   `k′`=input; the standard names for shared symbols so transliteration is name-for-name; `v_bulk_normal_0` never
   aliased to `v_0`) match the SymPy directive's, so the T7 axis-typed join and the §5e leg-corruptions land on
   the SAME object in both engines? Does the WL-symbols-cannot-contain-`_` handling correctly separate internal
   spellings from emitted standard names?
7. **Designed-to-agree (rule 16 / rule 6).** Does anything in the directive push the WL engine to reproduce the
   SymPy representation/form/values rather than construct independently? The §5c/§5d regression operand B for WL
   is an INDEPENDENT re-derivation of the flat S11b object (WL imports nothing) — is that stated correctly, and is
   §5d operand B the BARE half-space `Z` with no `W_0→W̄₀(1+η)` and no two-face re-solve (spec §5d)?
8. **WL idioms + run discipline.** Are the mandated idioms correct and sufficient (per-tag `Flush`;
   `stripConditional`; `Inactive[Equal]` for relationals/loci; the spool/streaming idiom for the heavy
   two-momentum kernel; pole test `1/q_out==0`; regime classification NOT forced through `_EQUATIONS`/`_SOLUTION`;
   the §6 locus protocol only for finite-dim algebraic loci)? Is the run discipline right (one kernel / 2-seat;
   `danger-full-access`; 600s-no-output + RSS kill recorded-not-narrowed; detached launch for the ~87s reap;
   defer-heavy-out-of-band, never drop)? Is the typed-boolean → T7-rejects-native-boolean linkage correct?
9. **Structural rule + three clauses.** Are §6's three script clauses and the structural rule ("physical symbols
   combined by hand only in the supplied §§1–2 setup; every §3–§5 expression reached by computation; every control
   re-enters at an ansatz/map/level-set/operator/law, never at a result") stated exactly and bindingly?
10. **Anything that would make the builder resolve an ambiguity wrongly** — an under-specified object, a missing
    spec pointer, a conflict between the directive and the spec, an instruction the builder could satisfy while
    producing the wrong physics.

## Method and output

This is a DOCUMENT review, not a script ablation — but a prose claim is worth nothing without the spec quote that
grounds or refutes it (`CLAUDE.md` rule 2). For every finding: quote the directive text and the spec text it
violates or distorts, with file:line. State severity: **MUST-FIX** (the build would be wrong / the control
vacuous / an answer leaked) vs **NIT** (clarity, no physics consequence). If you believe the directive is safe to
build as-is, say so explicitly and name the two or three obligations you checked most closely. Do not invent work
for the builder that the spec does not require. Report only findings that catch a way the physics or the
cross-engine control could be wrong.

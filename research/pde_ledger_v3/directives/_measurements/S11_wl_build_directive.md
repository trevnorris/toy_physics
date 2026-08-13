# Measurements — S11_wl_build_directive.md

Every factual claim in the directive, with the command that produced it and its literal output.
Run from `/var/projects/toy_physics`. All at HEAD `cebd0e5e` (2026-08-13).

The directive was folded once after two review legs (Codex + Grok, reports at
`~/.s11_build/wl/directive_leg_codex.log` and `~/.s11_build/wl/directive_leg_grok.log`); the pre-fold
defects those legs found and this fold's changes are recorded in the leg logs, not restated here.

## The current `.wl` exists and predates this build

```
$ git log --oneline -2 -- research/pde_ledger_v3/mathematica/S11_stray_longitudinal_mathematica_audit.wl
fc920079 S11 engines rebuilt, spec repaired -- and the comparison method is changing
b037ce93 S11: both engines, the registry rows, and the two repairs their reviews forced

$ wc -l research/pde_ledger_v3/mathematica/S11_stray_longitudinal_mathematica_audit.wl
910 research/pde_ledger_v3/mathematica/S11_stray_longitudinal_mathematica_audit.wl
```

## The blindness-control citation (`S9_export_chain_rebuild_directive.md:16-18`)

```
$ sed -n '16,18p' research/pde_ledger_v3/directives/S9_export_chain_rebuild_directive.md
construction rather than reconciled against a second copy. Each step's Wolfram engine imports nothing and
re-derives independently — **that is the only blindness control in this design, and nothing else in it
should be built pretending to be one.**
```

## Spec citations

Engine-local scoping `:1097-1105` (the `_LOCAL_` infix and why parity needs the declaration):

```
$ sed -n '1097,1105p' research/pde_ledger_v3/directives/S11_SHARED_PHYSICS.md
### ⭐ Engine-local tags — declare them, so parity is meaningful

Some tags **cannot** exist in both engines: §Q6r's previous-step export-chain comparison has no counterpart
in the engine that does not import that chain's `LEDGER`, and each CAS emits its own solver-condition tags.
⛔ Those are **not** disagreements, and a parity checker that reports them as gaps trains its reader to
ignore it.

⭐ Give every such tag the infix `_LOCAL_` immediately after the engine prefix — `WL_S11_LOCAL_…`,
`PY_S11_LOCAL_…` — and ⭐ emit one tag per engine listing every `_LOCAL_` name it produced.
```

Streaming/flush rule `:1044-1047`:

```
$ sed -n '1044,1047p' research/pde_ledger_v3/directives/S11_SHARED_PHYSICS.md
⭐⭐ **Emit every tag as it is computed, and FLUSH.** ⛔ Do not buffer output to the end of a package, of a
dimension, or of the run. ⭐ A long run that is visibly completing cells is **acceptable**; a run producing
**no output for a long stretch** is the failure, because from outside nothing distinguishes it from a solve
that will never return.
```

Runtime rule header `:1042` (why the kill criterion is silence, not elapsed time):

```
$ sed -n '1042p' research/pde_ledger_v3/directives/S11_SHARED_PHYSICS.md
⚠⚠ **Runtime — the rule is about OBSERVABLE PROGRESS, ⛔ not total elapsed time.**
```

Every-package obligation `:301-302` (why acceptance item 2 demonstrates one cell per §7 package):

```
$ sed -n '301,302p' research/pde_ledger_v3/directives/S11_SHARED_PHYSICS.md
⛔ **This section says what to produce. It does not say what any of it equals.** Everything below is
required for **every package** defined in §7, at **every** `D` in that package's sweep.
```

`InputForm` payload convention `:1088`:

```
$ sed -n '1088p' research/pde_ledger_v3/directives/S11_SHARED_PHYSICS.md
(`InputForm` in Wolfram; `sympy.srepr`-safe string or `str()` of the expression in SymPy).
```

## The S9 `.wl` emit-site citation (`:25-30`)

```
$ sed -n '25,30p' research/pde_ledger_v3/mathematica/S9_light_requires_shear_mathematica_audit.wl
emittedTags = {};
emit[tag_, value_] := Module[{emittedTag = Lookup[standardEmissionNames, tag, tag]},
  If[MemberQ[emittedTags, emittedTag], Exit[91]];
  AppendTo[emittedTags, emittedTag];
  Print[emittedTag <> ": " <> ToString[value, InputForm]]
];
```

The claim that this site shows no explicit flush operation is read directly off the six lines above:
the only output call is `Print`.

## Sibling-directive precedent for the pointing and corollary-4 phrasings

```
$ rg -ni "point to rather than duplicate" research/pde_ledger_v3/directives/S11_sympy_build_directive.md
12:and wins every conflict. Implement its §§1–10 directly. In particular, point to rather than duplicate §4;

$ rg -n "corollary 4 applies as a property" research/pde_ledger_v3/directives/S11_sympy_build_directive.md
17:exit policy. Section 5 corollary 4 applies as a property, without named exceptions.
```

## Out-directory convention

The command establishes only where engine stdout lands by convention — existence, nothing more:

```
$ ls research/pde_ledger_v3/mathematica/out/
S10_brane_mode_spectrum_mathematica_audit.out
S11_stray_longitudinal_mathematica_audit.out
S9_light_requires_shear_mathematica_audit.out
```

## Not measurements — decisions and supplied constraints, marked as such

- The input enumeration (spec = sole physics input; directive, `CLAUDE.md` and the two precedent
  citations as the only other inputs): orchestrator decision, structured as a positive enumeration
  after both legs flagged the earlier "exactly two files" sentence as self-contradictory denylist prose.
- Two licence seats, one kernel at a time during a build: user-supplied constraint (my paraphrase;
  standing since 2026-08-03, restated in `.claude/skills/review-legs/SKILL.md`).
- Kill criteria (600 s of NO new output; 6 GB RSS), replacing an earlier elapsed-time `timeout 600`
  that both the spec (`:1042`, quoted above) and a leg showed to be a contradiction: orchestrator
  decision aligned to the spec's observable-progress rule.
- The argv contract (`[PACKAGE D]` positional, no other shape): orchestrator decision, so the
  orchestrator can reproduce every demonstration command from this directive alone.
- The build never running the declared sweep to completion; demonstration stdout only under the
  build's scratch paths; one demonstration cell per §7 package: orchestrator decisions.
- The product being the stdout stream and nothing else: orchestrator decision, mirrored from the
  closed sibling directive's product clause (its line 5-7).

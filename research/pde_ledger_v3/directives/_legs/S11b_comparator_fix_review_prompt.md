# Directive-design review — S11b comparator fix round 1

You are one of two independent legs reviewing an ORCHESTRATOR-written FIX directive BEFORE the builder runs.
The comparator's one cardinal sin is reporting a genuine difference as anything other than DISAGREE. Three
such holes were found by two prior legs and are to be fixed. Your job: does this fix directive close them
WITHOUT opening a new one, and is its added acceptance decisive? Find defects; a leg that finds nothing is
weak — show the checks you ran.

## Artifact under review
`/var/projects/toy_physics/research/pde_ledger_v3/directives/S11b_comparator_fix_round1_directive.md`

## Sources of truth — read first, form your own view, then read the directive
- The comparator being fixed (committed baseline): `research/pde_ledger_v3/scripts/S11b_cross_engine_comparator.py`
  — the relevant code: `transliteration_collisions` (~L505), `_transliterate_basic` (~L525),
  `authoritative_undecided_token` (~L691), `compare_records` UNDECIDED branch (~L800),
  `_convert_parsed_containers` / `_text_value` (~L405-465).
- The frozen contract: `research/pde_ledger_v3/directives/S11b_comparator_build_directive.md` (the 7-item
  delta + the "disagreements are the expected output" rule).

## Checks
1. **Do the three fixes actually close the holes?** For each: fix 1 (function-head collision guard), fix 2
   (residual siblings before UNDECIDED), fix 3 (promote only textual-atom keys, never bare Symbol). Is the
   prescription correct and complete, or could the hole persist (e.g. a function↔symbol collision, a nested
   status token, a mixed Str/Symbol key tuple)?
2. **Do any of the fixes OPEN a new hiding path?** Especially: could "DISAGREE dominates" ever be bypassed?
   Could the per-object residual guard (item 5) ever land on AGREE or silently drop instead of UNCOMPARED?
   Could restricting tuple promotion break a legitimate Str-keyed Association comparison?
3. **Is the added acceptance decisive?** For each new fixture, is there a bad fix that still passes it? Does
   the mixed-status fixture actually force DISAGREE (not UNDECIDED)? Does the collision fixture exercise a
   FUNCTION head, not just a symbol?
4. **Rule-5 / value-free**: does the directive leak a real value, a real count, or use the real `.out` as an
   acceptance gate?
5. **Regression**: does it adequately require the 7 confirmed-correct core rules to stay unchanged?

## Method
Document review; you may run small scripts against the baseline comparator to check a claim (save script +
literal stdout to named absolute paths, report them). ⛔ Do not edit the working tree.

## Report — under ~20 lines
Numbered findings (directive line, the failure it permits, the concrete fix), then a one-line verdict: safe
to build from, or fold first?

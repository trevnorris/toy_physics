# Independent review — round-4 census-instrument repair directive (decision list)

## Artifact

`/var/projects/toy_physics/research/pde_ledger_v3/directives/_legs/S11_census_instr_repair4_directive.md`
at commit `c37817a7`, with its measurements twin in `../_measurements/` (same name).

## What to check

This is an orchestrator-written decision list for a SINGLE defect class: a builder will
implement it verbatim against the instruments at `fd9a5835`. The defect was verified by both
round-3 review legs and orchestrator reproduction (commands in the twin); your job is the
DIRECTIVE:

1. **Is each Required semantics CORRECT and sufficient?** Check Required 1 (constant
   zero-status by certificate, UNDECIDED on failure) and Required 2 (sampled coverage refutes
   only on a certified-nonzero union product point) against the actual code at `fd9a5835`
   (`research/pde_ledger_v3/reduction/s11_census_math.py:544-578, 689-717, 829`) and real
   record payloads. Can a wrong verdict survive both rules — e.g. a numeric certificate at
   fixed precision on an algebraic number that is nonzero but smaller than the working
   precision (what separation bound is achievable, and does the directive's "rigorous
   separation bound" have a decidable reading)? Does UNDECIDED-on-constants propagate soundly
   into every consumer named (membership, spurious, probe, definedness) without flipping a
   certified verdict class?
2. **Is anything under-specified** — satisfiable two materially different ways with different
   census verdicts?
3. **Does any requirement weaken a certified direction?** The closing constraint lists them;
   check the twin's bounding note (sampled refutations upheld in round 3 must not change) is
   actually implied by the Required rules rather than merely asserted.
4. **Level-above misses**: an interaction with the round-3 classes the directive does not
   order (e.g. constant-UNDECIDED inside the witness conjunct classification of round-3
   defect 4 — does a nested-radical constant premise atom now become UNDECIDED where it was
   decided, and is that sound?); a calibration plant that cannot fail as specified; a plant
   whose specified behavior encodes an expected census outcome a builder could iterate toward.
5. **Twin claims**: verify each measured claim in the twin by computation (the union-product
   zero, the four pairing statuses, the D3/D4 cand[0]/cand[1] split).

## What you are handed

The directive + twin; the round-2/3 directives (`33babf8d`, `ef9085c6`) and the brief
(`fa8c58b3`); the instruments at `fd9a5835`; round-3 census outputs under
`/home/trevnorris/.s11_build/census_build3/`; both committed records (WL `a4cf6539`, PY
`19591194`); `DEFECT_REGISTER.md` (untouched, 7 entries).

## Required method

Where a directive claim is checkable by computation, check it by computation and show the
literal command and output; save scripts and stdout to named absolute paths under `/tmp/` and
report the paths. A prose objection to a computable claim will be discarded. If you claim a
Required semantics is wrong, construct a concrete counterexample payload (byte-shaped like the
records) — a hypothetical is not a finding.

## Physics filter

Report a finding only if it changes what the builder would build, what the census may verdict,
or what the round may claim. Style, tone, and formatting are out.

## Sandbox

Read-only on the working tree; write only under /tmp; never commit. No Wolfram kernel. Cap any
single script run at `timeout 600`.

## Report format

Numbered findings, most severe first: claim, the directive line it targets, evidence (command +
literal output or concrete counterexample payload). If nothing survives the filter, list what
you checked with literal outputs.

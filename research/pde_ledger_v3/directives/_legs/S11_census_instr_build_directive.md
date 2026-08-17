# Build directive — S11 locus-census instruments

You are building the two census instruments governed by
`/var/projects/toy_physics/research/pde_ledger_v3/directives/_legs/S11_locus_census_instrument_brief.md`
at commit `fa8c58b3`, in the working tree at `/var/projects/toy_physics`, with its generated
measurements twin in `../_measurements/`. The brief's eight numbered obligations and its census
population section ARE the acceptance surface; nothing below overrides them.

## The one structural rule

Scripts print computed objects; they never state conclusions. Per-record lines carry operands and
computed memberships/residuals, then the verdict token; nothing is asserted before it is printed.
The instruments read committed records; they author no physics expression of any kind.

## Order of work

1. Build the containment/soundness census and the undecided-record probe census as repo-homed
   tools under `/var/projects/toy_physics/research/pde_ledger_v3/reduction/`, taking a record
   file path as argv. Implement the brief's population enumeration and reconciliation first and
   report its literal reconciliation output for both committed records before implementing any
   verdict logic.
2. Build the planted-record calibration per brief obligation 2 (byte-shaped real lines, both
   dialects, one defect each; fed through the production entrypoint). Run it FIRST and report
   its literal output — a calibration that cannot fail is a broken instrument and the build
   stops there.
3. Run both instruments over both committed records (WL canonical at `a4cf6539`:
   `research/pde_ledger_v3/mathematica/out/S11_stray_longitudinal_mathematica_audit.out`; PY at
   `19591194`: `research/pde_ledger_v3/scripts/out/S11_stray_longitudinal_sympy_audit.out`),
   preserving complete literal stdout to files under your scratch directory.
4. Build and run the named acceptance reducer (brief obligation 6) over the census stdout;
   report its literal output.
5. Write the build report: for every measured claim, the command and its literal output. Census
   verdicts on real records are FINDINGS to report with their per-record lines — never repair
   anything the census finds, never re-run an engine, never modify a committed record.

## Operational constraints

- The containment census must not spawn a Wolfram kernel. If the probe census needs one: one
  kernel at a time, never in parallel with anything, and it reads records — it never runs
  engine cells. Prefer a kernel-free design.
- Absolute paths everywhere. Scratch and logs under `/home/trevnorris/.s11_build/census_build/`
  (create it), never under the repository and never under /tmp.
- Do not modify either engine, the spec, any committed record, or DEFECT_REGISTER.md. Do not
  commit; the orchestrator commits after independent review.
- Long-running census stages must be visibly emitting per-record lines; long + silent is the
  failure mode this program removes.

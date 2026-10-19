"""Package only the final round-11 captures after mechanical validation."""
from pathlib import Path
import gzip,hashlib,json,shlex
ROOT=Path(__file__).resolve().parents[4]
DEST=Path(__file__).resolve().parent
STATE=ROOT/'_scratch/s9b_wl_repair11'
WORK=STATE/'ablation-final'
PREFIX='S9b_repair_wolfram_round11_'
v=json.loads((DEST/(PREFIX+'validation.json')).read_text())
def write(suffix,text): (DEST/(PREFIX+suffix)).write_text(text)
def link(suffix,label): return f'[{label}]({PREFIX}{suffix})'
for source,suffix in [(ROOT/'_scratch/s11c/s9b-r11-full/stdout','harness.stdout.gz'),(WORK/'BASELINE.stdout','baseline.stdout.gz')]:
    with gzip.open(DEST/(PREFIX+suffix),'wb') as f: f.write(source.read_bytes())
captures={p.name:p.read_text() for p in sorted(WORK.iterdir()) if p.suffix in ('.wl','.stdout')}
with gzip.open(DEST/(PREFIX+'captures.json.gz'),'wt') as f: json.dump(captures,f)
hashes=json.loads((DEST/(PREFIX+'source_hashes.json')).read_text())
hashes[str(Path(__file__).resolve().relative_to(ROOT))]=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
write('source_hashes.json',json.dumps(hashes,indent=2)+'\n')
full=v['runs']['s9b-r11-full'];k10=v['witnessSearches']['K10']
commands=[]
for name in ['s9b-r11-checks-final','s9b-r11-full','s9b-r11-engine-repo','s9b-r11-engine-alone']:
    r=json.loads((STATE/(name+'.run.json')).read_text())
    commands.append(f"### {name}\n\nWorking directory: `{r['cwd']}`.\n\n```sh\n{shlex.join(r['command'])}\n```\n\nLiteral guard stdout:\n\n```json\n"+(STATE/(name+('.attempt'+str(r['admissionRetries']) if r['admissionRetries'] else '')+'.guard.stdout')).read_text()+"```\n")
knifes=[]
engine=(STATE/'engine-pinned-final.wl').read_text()
for name,mutation in v['mutations'].items():
    lineno=engine[:engine.index(mutation['construction'])].count('\n')+1
    knifes.append(f"| {name} | {lineno} | `{mutation['construction'].replace(chr(10),' ')}` |")
witness_rows=[]
for name,row in v['witnessSearches'].items():
    witness_rows.append(f"| {name} | {row['result'].strip(chr(34))} | `{row['evaluations']}` |")
checks=(DEST/(PREFIX+'structural_checks.txt')).read_text().splitlines()
warning_count=full['stderr'].count('Limit::alimv: Warning:')
report=f'''# S9b Wolfram round 11 — final K10-only repair

The final-code baseline/K1–K11 harness, checks and both engine portability runs completed with exit code 0 through `--pool s9b --memory-gib 8 --tasks-max 64`, one Wolfram kernel at a time, without execution deadlines. No guard kill, memory exhaustion, task-limit or item-12 event occurred. Guard admission refusals: **{len(v['guardRefusals'])}**, all caused only by other jobs' memory commitments; admission was retried at the same 8 GiB limit. **Commit IDs: none.**

K10 now replaces the computed along-ray branch-existence `ForAll` by its logical negation only when `rayGate` is assembled at engine line **266**. Ray construction, traversal, and all base/forward `BRANCH_EXISTENCE`, `PATH_TRAVERSAL` and `BRANCH_TYPE` emissions are unchanged. The inventory excludes those independent diagnostics and includes the observable/reduction domains and their dependent parts. K10's witness is the conditional domain of `A_NONRECIPROCAL_PATH_DEPENDENCE["OnDomain"]`; its evaluations are **`{k10['evaluations']}`** (corrupted, baseline), using the computed common full-ray solution. Both whole emitted objects are evaluated. The engine, other knives, structural comparison and K11 method are unchanged.

The engine emits **{v['tags']} tags = {v['shared']} shared + {v['local']} local**; no required shared name or item-10 output is missing. All **{len(v['witnessSearches'])}** searches record witnesses. **{v['baselineNotEstablished']}** baseline tags contain conditional `NOT_ESTABLISHED`; complete occurrence lists are linked below. Structural marks claim literal identity only.

K11 retains exact mechanical prefix extraction with identical configuration, untruncated local dispersion/Fermat and the supplied twelve-grade one-form ring: **{v['K11Reached']}** reached tags print both ordered pairs; **{v['K11Unreached']}** tags print not evaluated/not compared. The method residual supplies no witness.

Item 13 remains general radial functions/functionals, exterior Abel reduction with a GM anchor and signed velocity branches, and coordinate Cartesian mass divergence. No profile family, fold binding or premise is added. Supplied-input bindings/spec lines, restrictions and every knife site are listed in the measurements.

Validation includes **{len(checks)}** checks records, byte-identical nonempty repository/standalone/baseline engine streams, only the copied file in the standalone directory, observed capture growth while the kernel was alive, and empty production-out status. Full-run peak cgroup memory: **{full['peakMemoryBytes']} bytes**; maximum sampled tasks: **{full['maximumSampledTasks']}**; sampled swap: **{full['maximumSampledSwapBytes']} bytes**. **{warning_count}** printed limit-assumption warnings and their suppression notice are preserved.

Read: {link('measurements.md','commands and item-17 details')}, {link('validation.json','validation')}, {link('witness_searches.txt','attempts')}, {link('k11_not_evaluated_tags.txt','K11 unreached tags')}, {link('source_hashes.json','source hashes')}. No build question remains open; interpretation and review belong to the orchestrator. Nothing was committed or pushed; no other AI, reviewer or agent was launched.
'''
(DEST/'S9b_repair_wolfram_build_report.md').write_text(report)
measurements=f'''# S9b Wolfram round 11 — final-code measurements

This packet uses only the final round-11 sources and their current executions. The engine is unchanged. The only harness repair is K10's observable gate, site inventory and domain witness; the scratch output destination is isolated for this round. The checks exercise synthetic gate factors and conditional-object readings. Input snapshots are scope evidence, not a source of measured results. This is an execution record, not a review verdict or universal equality claim.

## Artifacts and reproduction

- {link('harness.stdout.gz','Complete harness stdout')}; {link('baseline.stdout.gz','baseline stdout')}; {link('captures.json.gz','all engine captures and generated mutation/compact sources')}.
- {link('witness_searches.txt','Every attempt and density-sign trial')}; {link('site_inventory.txt.gz','site inventories, exclusions and search plans')}; {link('k11_reached_pairs.txt','reached K11 ordered pairs')}.
- {link('structural_checks.txt','Literal checks stdout')}; {link('harness.stderr.txt','literal harness stderr')}; {link('validation.json','mechanical validation')}; {link('validation_stdout.txt','literal validator stdout')}.
- {link('validate.py','Mechanical validator')}; {link('package.py','packager')}; {link('source_hashes.json','source hashes')}; {link('input_sources.json.gz','handed-in source snapshots')}; {link('artifacts.json','artifact hashes and sizes')}.

The validator checks every printed operand against its actual capture, the vocabulary and local-name registry, unique single mutations, unchanged engine/other-knife source, K11 mechanical prefix extraction, pinned sources, all witness check flags and sign trials, fixture records, resources and portability. It invokes no CAS. Literal structural comparison in the harness is `SameQ` on `BinarySerialize` integer bytes of frozen parsed `HoldComplete` InputForm; no evaluation, simplification or canonicalization precedes a mark. Approximate-number encodings are printed. Marks are only **identical** or **not compared** and decide no mathematical equality or change.

```sh
python research/pde_ledger_v3/directives/_measurements/{PREFIX}validate.py s9b-r11-full
python research/pde_ledger_v3/directives/_measurements/{PREFIX}package.py
```

## Guarded commands and literal outputs

All commands use 8 GiB, zero swap, TasksMax 64 and pool s9b. Before each, the local driver polls `pgrep -x WolframKernel` until it is empty; only one kernel is launched. No wall-clock, CPU-time or inactivity deadline is imposed. The final harness/checks executions were awaited by polling until completion; no detached driver, watcher, queued resumption or recurring schedule was launched for them. Logs and final source pins remain under `_scratch/`.

{chr(10).join(commands)}

Guard admission refusals: **{len(v['guardRefusals'])}**, all solely because available memory did not cover outstanding commitments plus this job. Admission was retried after waits of about three minutes at the unchanged 8 GiB limit; no CAS workload started in those refused attempts. Exact refusal tracebacks are in validation JSON. Guard kills, Wolfram memory exhaustion, task-limit events and item-12 events: **none**. Full-run peak cgroup memory: **{full['peakMemoryBytes']} bytes**; maximum sampled tasks: **{full['maximumSampledTasks']}**; maximum sampled swap: **{full['maximumSampledSwapBytes']} bytes**. The cgroup and native address-space caps are **8589934592 bytes**. Per-run outcomes and final memory-event counters are recorded in validation JSON.

Harness stderr is preserved verbatim. It prints **{warning_count}** `Limit::alimv` warning messages and a suppression notice: those limit calls ignore global assumptions involving the limit variable. The checks and both standalone engine runs have empty stderr. No warning is reported as a guard event or a method failure.

## Tags and search outcomes

**{v['tags']} tags = {v['shared']} shared + {v['local']} local**: {link('tags.txt','complete list')}. No required shared name, local registry entry or item-10 output is missing. Every non-K11 copy emits all tags; the harness prints every baseline/corrupted pair and baseline-against-baseline control. K11 prints all full-baseline tags with reached/unreached status.

**{v['baselineNotEstablished']}** baseline tags carry conditional `NOT_ESTABLISHED`: {link('not_established_tags.txt','baseline occurrences')}, {link('harness_not_established_tags.txt','harness pair occurrences')}, {link('capture_not_established_tags.txt','all engine capture occurrences')}. These are literal conditional branches, not numerical decisions for an independently chosen ray.

| Knife | Search outcome | Evaluations (corrupted, baseline) |
|---|---|---|
{chr(10).join(witness_rows)}

Each search prints the union of all direct/transitive site parts across copies and every exclusion, independently of structural marks. K11 includes compact-unreached parts. Mathematical values/domains are distinguished from containers; branch sets and ordered first-match selections retain their engine readings. Attached conditions are carried through containing levels. Density labels partition the shared response object without imposing an extra condition. Absent entries are not emitted and receive no invented value. Additional exclusions are extracted symbol/function/row-name lists and literal replacement-source text. All readings and computations are in the inventories and attempts.

All searches stop at their recorded witness. Other enumerated parts receive no no-change decision. No sample, unresolved evaluation, agreeing point or method residual supplies a universal conclusion. K7 uses the engine's literal status text at `B_IMPLIED_JN_DEFLECTION["VelocityStatus"]`, with a printed CAS literal comparison. K6 evaluates the emitted local radar integral kernel and finite-endpoint source with its constructor conditions. Other numeric witnesses print all roles, exact profiles, points, fixing equations/solution intersections, attached conditions and finite evaluations.

K9a and K9b each try the shared object at concrete inputs in both density signs and print both trials before returning the first witness. The negative input has `f(r)=-9/r`, `rho0=1`, `r=3`, `n=3`; the positive input has `f(r)=1/(100r)`. Both have the supplied asymptotic density. The exact density evaluations, condition/value computations and profile certificates are in the attempt file and validation JSON. Test profiles enter no engine object and are not profile-family restrictions.

## K10 construction and actual witness

The single replacement is at `rayGate = propagationGate && traversalGate;` (engine **266**). `With` inserts the computed traversed-radius domain into a held pattern. In the already computed `propagationGate`, that pattern replaces the complete along-ray branch-existence `Inactive[ForAll]` by `Not` of itself. The quantifier body is retained, not negated pointwise. The enclosing ray-construction conjuncts and the independently computed traversal gate are retained. `propagationGate` itself is not assigned a new value.

The generated source has exactly this one replacement. All local tags, ray-domain construction and all base/forward `BRANCH_EXISTENCE`, `PATH_TRAVERSAL`, `BRANCH_TYPE` payloads are byte-identical to baseline. The synthetic gate-factor fixture applies the actual manifest replacement and prints unchanged diagnostic factors and the selected negated quantifier. These unchanged diagnostics are excluded from K10's site inventory. Observable/reduction ray-gate domains and values depending on their attached domains remain included.

The witness site is **`WL_S9B_A_NONRECIPROCAL_PATH_DEPENDENCE`**, payload **`{{"OnDomain","conditional domain"}}`**, reached through engine **266 → 271 → 287**. Its domain evaluations are **`{k10['evaluations']}`**. The full traversal evaluations are **`{{{', '.join(k10['traversalEvaluations'])}}}`**. Whole emitted object evaluations, corrupted then baseline:

```wl
{chr(10).join(k10['emittedObjectEvaluations'])}
```

Whole-object evaluation first resolves each scoped conditional domain, then evaluates its value only on a true domain. `OnDomain`'s false condition supplies a decided absence (`Undefined`); its complementary `OutsideDomain` supplies `NOT_ESTABLISHED`. The association introduces no additional condition. This is a witness in an emitted observable domain, not in the unchanged branch diagnostic or a local pre-gate predicate.

Witness-only profiles are flat: delta, velocity, xi, bulk fraction and exchange zero; brane density one, stiffness four and light speed two. The printed independent point has `GM=0`, `b=2`, `bFar=1`, `ZE=ZR=2`, with every remaining input and role shown. Each copy's actual ray construction is independently solved. Turning candidates are filtered by their full domains; ordinary endpoint integrals are evaluated on the derived admissible interior and the closed upper endpoint is checked separately. Active quantifiers retain the full traversed ranges. Complete real solution sets and their intersection are computed; turning radii and radar impact are fixed by those equations, not independently assigned. The intersection is:

```wl
{k10['solutionIntersection']}
```

Both copies use that common solution. Every other attached condition and the traversal gate is computed true, and both emitted objects are evaluated. The attempt includes scoped bound variables, equations, solution selection and all CAS operations.

## Item 13, bindings and restrictions

General live radial functions and integral functionals represent the profiles. Every-b reduction uses exterior-half-line Abel inversion, a constant-GM anchor and signed squared-velocity branches carrying integration-constant domains. Coordinate Cartesian component differentiation computes mass divergence before radial evaluation. No fold/LEDGER binding, flux-to-GM relation or premise is introduced.

There is no engine profile-family restriction. Item-10 flow-only sets delta/xi to zero; speed-only sets velocity/xi to zero; tilt-only sets delta/velocity to zero; bulk-only sets velocity/xi to zero. The forward case uses its labelled `jn=0` premise and flux-defined velocity, first with delta/xi zero, then live. These additional restricted outputs do not replace general-profile objects. Part C applies its supplied responses and first-order f counting, retaining the exact shared response/domain at each response site.

Supplied-input/spec-line bindings in `S9b_SHARED_PHYSICS.md`: dispersion/metric **45**; optical ratio **52**; polarizations **59**; advection **64**; mass balance **69**, coordinate measure **72–74**; embedding **87**; historical uniform anchor **95**; LAB_HELD/background anchoring **103–104**; reference speed **116**; delta/nonnegative root **124**; retained grades **128**, optical counting **135,138–145**; independent GM **202–204**; references **208–210**; pressure, sound square and bulk fraction/asymptote **219–224**; Part-C responses **244–246**. Directive item 9 supplies root, zero-power and density-sign conventions; item 10 supplies the labelled forward premise. Part D is not built.

| Knife | Engine line | Original construction |
|---|---:|---|
{chr(10).join(knifes)}

Every construction is present. Exact replacements and unique occurrence counts are preserved in the full SITE records and validation JSON.

## K11 and executable checks

K11 uses exact mechanical extraction through `COMPACT_K11_BOUNDARY`, with identical configuration for corrupted and unmutated copies. It retains untruncated local dispersion/Fermat and the twelve-grade one-form ring, with no added truncation, samples, seeds or inexact arithmetic. Coverage is not an attached condition or a universal equality claim.

**{v['K11Reached']} reached tags**: {link('k11_evaluated_tags.txt','complete list')}. Each prints the full baseline, compact-corrupted/compact-unmutated difference pair and compact-unmutated/full-baseline method residual pair, with structural marks and coverage. The method residual is not a witness. The full corrupted K11 engine is not run.

**{v['K11Unreached']} not-evaluated tags**: {link('k11_not_evaluated_tags.txt','complete list')}. Both pairs are explicitly not evaluated/not compared; no no-change decision is made. The inventory includes their downstream site parts.

The checks script prints **{len(checks)}** records before validation: syntax and structural fixtures, exact approximate-number encodings, domains, binders, ordered cases, common/disjoint solution sets, semantic branch identities, literal text/list readings, density-sign/asymptote bookkeeping, whole-quantifier versus pointwise negation, and the actual K10 replacement on synthetic gate factors. There are no added expected physics values.

The repository and copied-alone engine runs exit 0, produce nonempty byte-identical streams, and also match the full harness baseline byte for byte. Only the copied `.wl` remains in the standalone directory. Capture growth was observed while the kernel was alive; observations are in validation JSON. Production-out `git status --porcelain` is empty. No other checkout or orchestrator memory was read. No commit, push, other AI, reviewer or agent was launched. Build/report work stops here; interpretation/review remains with the orchestrator.
'''
write('measurements.md',measurements)
artifacts={p.name:{'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()} for p in sorted(DEST.glob(PREFIX+'*')) if p.name!=PREFIX+'artifacts.json'}
write('artifacts.json',json.dumps(artifacts,indent=2)+'\n')
print(json.dumps({'packagedArtifacts':len(artifacts),'report':'S9b_repair_wolfram_build_report.md'},indent=2))

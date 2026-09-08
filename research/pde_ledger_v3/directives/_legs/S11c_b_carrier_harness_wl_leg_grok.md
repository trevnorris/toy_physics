# Independent review leg — S11c-b WL carrier ablation HARNESS (astra-written), WL EXECUTION-ABLATION ONLY

You are the **second** independent review leg (Grok) for a code engine's (gpt-6-astra) WL ablation harness. A
fresh-Claude leg already reviewed the full pair and the SymPy harness is separately two-leg-closed; **your job is
the one piece not yet independently done: the FORM ablation of the WOLFRAM harness by execution.** Working dir
`/var/projects/toy_physics`. Repo read access; write/run check scripts under `/tmp` ONLY; ⛔ do NOT modify the
working tree (copy to /tmp, ablate the copy). Prose is discounted unless grounded in a cited line or a check
script whose absolute path + literal stdout you report.

## ⛔⛔ BOUNDS — read first (a prior leg WEDGED by ignoring these)
- ⛔ **Write your report and exit.** ⛔ Do NOT spawn agents, and ⛔ do NOT build any `run_all` / `watch` /
  supervisor / loop orchestration — iterating is not your job. Run each command yourself, in the foreground, one
  at a time, read its output, move on.
- ⛔ **Do NOT run the committed driver `scripts/S11c_b_carrier_ablation_harness_wl.py`** — it serializes all 10
  modes (~50 min) and holds a Mathematica seat. You invoke `wolframscript` on individual MODES directly, per the
  recipe below.
- ⛔ **Wrap EVERY kernel run in `timeout 600`** and run **only one WL kernel at a time** (the licence has two
  seats; an orphaned kernel leaks memory). A run that hits 600s is a **failed** step: report it as a finding
  ("WL <mode> exceeded budget") and MOVE ON — ⛔ never raise the timeout, ⛔ never block waiting, ⛔ never retry
  more than once. Each canonical build is ~4–5 min, so 600s is real headroom.
- If WL execution is unavailable/times out repeatedly, fall back to a **static** read of the harness + astra's
  committed transcript and say so as a finding — ⛔ do NOT wedge.

## Artifacts
- WL harness under review: `research/pde_ledger_v3/mathematica/S11c_b_carrier_ablation_harness.wl`
- The engine it WRAPS (must NOT be re-derived): `research/pde_ledger_v3/mathematica/S11c_b_brane_operator_mathematica_audit.wl`
- astra's committed run transcript (baseline evidence): `research/pde_ledger_v3/_measurements/S11c_b_carrier_ablation_harness_wl.md`
- The knife-list it implements: `research/pde_ledger_v3/directives/S11c_b_carrier_ablation_harness_directive.md`

## What the harness must satisfy (the physics)
The carrier is `∂(operator row)/∂(native pressure atom) |_{atoms→0}` for every scalar row × pressure slot, computed
on the LIVE engine's own emitted operator (`evaluatedModel["EULERIAN","MATERIAL_ADVECTED","RHO4_CONSTANT"]["OPERATOR"]`),
never a hand-typed carrier. Each knife patches ONE construction site by exact-string `StringReplace` of the engine
source and recomputes the carrier from the patched definitions:
- **K_A** = drop `lambdaAResponse affinity` from `faceSources` flux (closure addend);
- **K_T** = set the whole `faceSources` `virtualWork` product to 0 (traction channel);
- **K_W** = redefine `pressureField[-1]` (lower face) from `pressureUpper` (native-pressure face collapse).

## ⛔⛔ YOUR CORE CHECK — ABLATE THE HARNESS, do NOT trust its self-report
A harness that PRINTS a "bite" whether or not the knife did anything is worthless. Prove, by execution, that (a)
each knife's printed carrier diff is caused by the patch (**no-op ⇒ diff identically 0**), and (b) a
**coefficient** change is not misread as a **structural** (FORM) one.

### Setup (once)
```
mkdir -p /tmp/gwl && cd /tmp/gwl
cp /var/projects/toy_physics/research/pde_ledger_v3/mathematica/S11c_b_carrier_ablation_harness.wl harness_live.wl
cp harness_live.wl harness_noop.wl
ENG=/var/projects/toy_physics/research/pde_ledger_v3/mathematica/S11c_b_brane_operator_mathematica_audit.wl
mkdir -p wd
```
In **`/tmp/gwl/harness_noop.wl`** make all three knives NO-OPs by replacing the THIRD element of each `patches`
entry with its already-defined baseline string (so `StringReplace` changes nothing):
- `patches["K_A"]` third element `"  flux = lambdaVResponse normalVelocity;"` → **`fluxBefore`**
- `patches["K_T"]` third element `"  virtualWork = 0;"` → **`workBefore`**
- `patches["K_W"]` third element `"pressureField[-1] := pressureUpper[xOne, xTwo, xThree, time];"` → **`pressureBefore`**
(These symbols are defined at the top of the file and equal the baseline source; `validateSites` still passes
because the `spec[[2]]` before-strings are untouched.)

### Runs (each its OWN foreground command, `timeout 600`, one at a time)
Invocation shape (worker MODEs export `<MODE>.wxf` into the workdir AND print `COMPUTED_OBJECT` with the carrier):
```
S11CB_CARRIER_MODE=<MODE> S11CB_CARRIER_WORKDIR=/tmp/gwl/wd S11CB_CARRIER_ENGINE=$ENG \
  timeout 600 wolframscript -file /tmp/gwl/<HARNESS> 2>&1 | tee /tmp/gwl/<MODE>.out
```
1. `MODE=BASELINE`, `<HARNESS>=harness_live.wl` → `wd/BASELINE.wxf` (also carries `ZERO_FIRST_CARRIER` for the
   extractor-order self-test).
2. `MODE=K_A`, `MODE=K_T`, `MODE=K_W`, `<HARNESS>=harness_noop.wl` → `wd/{K_A,K_T,K_W}.wxf` — these are the NO-OP
   runs; each patched engine equals baseline.
3. `MODE=RESCALE_K_A`, `<HARNESS>=harness_live.wl` → `wd/RESCALE_K_A.wxf` (the built-in ×2 coefficient variant;
   add `RESCALE_K_T`/`RESCALE_K_W` if you have budget).
4. `MODE=DEAD_PATH`, `<HARNESS>=harness_live.wl` → `wd/DEAD_PATH.wxf` (deletes a **pressure-free** term
   `kineticEwLive` from `THICKNESS_ROW`).

### Adjudicate (write your OWN tiny comparison wolframscript — cite its path + literal stdout)
Import the `.wxf` records (each has key `"CARRIER"`; BASELINE also has `"ZERO_FIRST_CARRIER"`) and print, per
mode, the slot-wise `Expand[modeCarrier[row][slot] - baselineCarrier[row][slot]]`. Report the LITERAL output and
whether:
- **NO-OP K_A/K_T/K_W diffs are ALL identically zero.** Any nonzero entry ⇒ the harness **fabricates** that
  knife's bite ⇒ NOT-SOUND (name the row/slot).
- **RESCALE_K_A diff is a coefficient multiple of the SAME live term** that K_A removes (the affinity carrier ×1
  extra), ⛔ NOT a vanished/new structural coupling. (Cross-check: the *live* K_A bite in astra's committed
  transcript is that same term with coefficient −1; the rescale must be +1× it.) A rescale showing as a
  structural change ⇒ NOT-SOUND.
- **DEAD_PATH diff is identically zero** (a pressure-free change must not move the carrier — proves the extractor
  responds only to pressure-bearing structure). Nonzero ⇒ NOT-SOUND.
- **EXTRACTOR_ORDER:** report the literal `CARRIER − ZERO_FIRST_CARRIER` diff and state whether the harness's
  chosen order (∂ then →0) is the physically correct one for a carrier defined as `∂/∂p |_{p→0}`.

## Also confirm by READING (no execution)
- The harness computes on the engine's emitted `evaluatedModel[...]["OPERATOR"]`, ⛔ not a re-implemented
  `faceSources`/`pressureField` or a hand-typed carrier (cite the lines).
- Each knife is exactly the directive's ONE site × ONE FORM, applied by exact-string patch of the production
  source (cite `patches`/`validateSites`), and `validateSites` genuinely pins each site to a single occurrence.
- ⛔ No `PASS`/`FAIL`/verdict/"bite" payload and ⛔ no asserted-zero/nonzero diff anywhere (PRINT-not-PASS).

## Physics filter
Report a finding only if it makes the harness certify the WRONG carrier or fabricate/suppress a bite. ⛔ Do not
report "it would be wrong on a different input".

## Output
Per-check verdict with cited lines / your check-script absolute path + its literal stdout, then a final line:
**SOUND** (the WL harness correctly certifies the carrier and cannot self-report a fabricated bite) or
**NOT-SOUND** (exact file + line + the fix for every blocking issue).

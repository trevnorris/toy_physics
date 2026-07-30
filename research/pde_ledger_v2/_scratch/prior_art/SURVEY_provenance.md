# PRIOR-ART SURVEY — artifact→producer binding

**Date:** 2026-07-28 · **Mode:** read-only survey · **Installed:** nothing · **Files changed in the repo:** none but this one.

---

## 0. THE HOLE, RESTATED PRECISELY

`ledger_dimensions.emit_dimension_sidecar` (`scripts/ledger_dimensions.py:184-222`) writes
`<stage>.dimensions.txt` whose header carries `source_sha256` (SHA-256 of the stage `.py` bytes,
`:214`) and `ledger_dimensions_sha256` (SHA-256 of the shared module bytes, `:215`). The comparator
recomputes both from disk and rejects a mismatch
(`scripts/compare_dimension_artifacts.py:204-230`). The module pin adds an *external* expected value
(`scripts/ledger_dimensions.accepted.sha256:1`, checked at `scripts/check_ledger_dimensions_pin.py:75-98`).

Everything above is a statement about **files on disk**. None of it is a statement about **execution**.
`emit_dimension_sidecar` never reads the destination and the comparator never invokes the stage, so a
hand-authored sidecar carrying the two correct digests passes every gate — as §9 of
`manifests/DIMENSION_REWRITE.md` records by demonstration ("the forged one says `PASS mismatches=0`,
exit 0, while the `.py` on disk still declares `a_dim: Dim(2,0,0)`").

Two distinct properties are in play, and most tooling conflates them:

| property | current status |
|---|---|
| **P1 staleness** — "the artifact names the source bytes now on disk" | ✅ closed by the header digests |
| **P2 origin** — "the artifact is what running that source emits" | ❌ **open on the `.py` side** |
| **P3 correctness** — "what the source emits is right" | out of scope here (comparator §11 + model read (c2)) |

The Mathematica side already closes P2 by the only mechanism that closes it: §4-**(g2)**, "run
`math -script <the .wl>` … and confirm it reproduces the committed `.out` byte-for-byte." That is
**re-execution by a party other than the producer, followed by byte comparison.**

---

## 1. ⭐ THE STRUCTURAL FINDING (this decides the whole survey)

> **No attestation scheme on the list binds an artifact to its producer. Every tool that closes P2
> closes it by re-executing the producer and comparing bytes. They differ only in (i) how hermetic the
> re-execution is, and (ii) whether the executor is outside the producer's reach.**

Reason, stated as a mechanism rather than an opinion: a digest, a signature, or a link-metadata record
is computed **over bytes that already exist**. Nothing in the act of hashing or signing can distinguish
bytes a process wrote from bytes a hand wrote. The only two families that can are:

- **(A) Re-derive and compare** — recompute the artifact and diff. Bazel, Nix, Make/CI, and the
  project's own (g2) are all this, with different amounts of hermeticity bolted on.
- **(B) Observe the write** — syscall/interpreter tracing that records *which process opened the file
  for write during the run* (ReproZip, Sciunit, noWorkflow, strace-class). This is the only family that
  binds P2 without re-execution — but the trace it emits is itself a file, so it only helps if the
  observer is outside the producer's reach, which lands you back at key custody or an outside executor.

Consequence for this repo: **the missing control is not a tool. It is (g2) for Python.** And unlike
Mathematica — barred from verification agents because of the licence seat — Python needs no seat, so
the Python-side (g2) can be **mechanised** instead of left to "the orchestrator regenerates the sidecar
itself before the (i) commit."

---

## 2. PER-TOOL FINDINGS

### 2.1 in-toto — **PARTIAL** (closes artifact *flow between recorded steps*, not artifact→producer)

- **Solves:** a signed record per supply-chain step. `in-toto-run` records `materials` (inputs) and
  `products` (outputs) as `{PATH: HASH}` dictionaries, plus the command line, and signs the resulting
  *link*. `in-toto-verify` checks links against a signed *layout* using artifact rules
  (MATCH / CREATE / MODIFY / DISALLOW).
- **Mechanism, and why it is not this problem — cited, not inferred.** The spec's own text on
  `expected_command`: *"during verification, clients should only show a warning if the expected command
  field does not match its counterpart in the link metadata"*, with the rationale that the field *"can
  easily be forged (e.g. by changing the PATH environment variable in a host) and thus it should not be
  trusted for security checks."* The link is produced by **hashing declared paths before and after the
  wrapped command** — there is no tracing of *which* process wrote the product. So
  `in-toto-run --products sidecar.txt -- python3 stage.py`, executed in a tree where `sidecar.txt` was
  hand-written and `stage.py` never touches it, yields a link that verifies. What in-toto genuinely
  gives is the MATCH-rule chain: *product of step N equals material of step N+1, by hash, and a
  DISALLOW rule fails on unaccounted artifacts.* That is tamper-evidence **between** steps.
- **By construction?** No. It is another self-attested digest, with the trust root moved to a
  functionary key. If the key is held by whoever runs the scripts, this is the current situation with
  extra ceremony.
- **Local:** absent (`in-toto-run`/`in-toto-verify` not on PATH; `in_toto` and `securesystemslib` not
  importable).
- **Cost/blast radius:** a layout, a key hierarchy, and every ad-hoc `python3 scripts/x.py` becomes
  `in-toto-run --step-name … -- python3 scripts/x.py`. High ceremony for zero P2 gain.
- **Verdict: PARTIAL.** Closes P1 more thoroughly than today (materials cover the declared dep set, not
  just one module). Leaves P2 entirely open.

### 2.2 SLSA / `slsa-github-generator` — **PARTIAL, and it is a different threat model**

- **Solves:** provenance about the **build platform**: builder ID, source repo + commit, artifact
  digest, signed by the platform. Build L2 = *"builds run on a hosted platform that generates and signs
  the provenance."* Build L3 adds isolation that *"prevent[s] secret material used to sign the
  provenance from being accessible to the user-defined build steps."*
- **That L3 clause is a genuine by-construction property — of the signature, not of the artifact.** It
  makes provenance unforgeable *by the build steps*. It says nothing about whether the artifact is the
  output of a particular script: a workflow whose only step is `cp committed_sidecar.txt out/` produces
  perfectly valid L3 provenance. SLSA is explicit that producer-submitted bad code is out of scope
  (*"SLSA v1.0 does not address source threats"*); its threat model targets a compromised builder or
  dependency, not a lying build script.
- **By construction?** For "who built it and on what platform," yes. For P2, **no** — a SLSA setup only
  closes P2 if the workflow *regenerates the artifact and fails on diff*, at which point the closing
  mechanism is §1(A) and the signing is decoration on top.
- **Local:** absent (`slsa-verifier` not on PATH). Repo has **no `.github/`** directory. Repo *is* on
  GitHub (`git@github.com:trevnorris/toy_physics.git`) and is **public**, so Actions minutes are free.
- **Cost/blast radius:** requires adopting GitHub Actions and expressing "the build" as a workflow.
  ⛔ **Mathematica cannot run on a hosted runner** (licence seat); a self-hosted runner restores
  Mathematica but discards the ephemeral-isolation property that L3 rests on.
- **Verdict: PARTIAL / mostly DIFFERENT PROBLEM.** The *useful* part for this repo is the plain
  Actions runner underneath it — see §4.

### 2.3 sigstore / cosign — **DIFFERENT PROBLEM**

- **Solves:** keyless signing (short-lived Fulcio certs bound to an OIDC identity) + a public
  transparency log (Rekor); `cosign attest` wraps an in-toto predicate.
- **Mechanism:** signs whatever bytes it is handed. It answers *who vouched* and *is this tamper-evident
  after the fact*. It cannot distinguish a produced artifact from a written one.
- **Local:** `cosign` absent; `sigstore` Python package absent.
- **Repo-specific note:** `commit.gpgsign=true` with `user.signingkey=251CA676820DC7F3` is **already
  configured**, so the "who vouched, attributably" layer exists today. The incremental value of
  sigstore/gitsign here is ≈ zero. Also: Rekor is a **public** log; this repo is public, so that is not
  a disclosure blocker, but it is worth knowing.
- **Verdict: DIFFERENT PROBLEM.** (One narrow place it *would* help is §9's complaint that
  `--accept` *"moves the baseline without a reason field, signature, or second witness"* — but signing
  is already on; what is missing there is a reason field and a second witness, which is process, not
  tooling.)

### 2.4 DVC — **DIFFERENT PROBLEM (dependency pinning / data versioning)**

- **Solves:** content-addressed storage for data/model files outside git, plus `dvc.yaml` stages with
  declared `deps`/`outs` and a `dvc.lock` recording hashes; `dvc repro` re-runs stages whose dep hashes
  changed.
- **Mechanism that decides it — cited.** `dvc commit` *"stores the current contents of files and
  directories tracked by DVC in the cache, and updates dvc.lock or .dvc files if/as needed"*, explicitly
  **without re-executing stage commands**. Its documented use cases include *"after executing stage
  commands by hand (outside `dvc repro`), use `dvc commit` to register the outputs"* and *"when you
  modify code files in ways that don't affect output … use `dvc commit` instead of reproducing the
  pipeline."* That is the forgery affordance, blessed and documented. `dvc repro` also *skips* stages it
  believes unchanged and, when it does run, **overwrites** the output rather than diffing it.
- **By construction?** No — it is P1 with better ergonomics. It *would* be a real P1 upgrade: today the
  pin covers exactly one module (`ledger_dimensions.py`) and nothing else — not sympy 1.14.0, not
  Python 3.10.12, not any other import — whereas DVC hashes the declared dep set.
- **Local:** absent (no `dvc` binary, no `dvc` module).
- **Cost/blast radius:** `dvc init`, a `dvc.yaml` stage per script (~30 for the rewrite scope, 91 files
  tracked in `scripts/`), a cache directory, and invocation changes from `python3 scripts/x.py` to
  `dvc repro x`. Moderate. Mathematica works fine as a `cmd:` string (it runs on the host, no sandbox),
  so no licence conflict — but it buys nothing for P2.
- **Verdict: DIFFERENT PROBLEM.**

### 2.5 Snakemake — **DIFFERENT PROBLEM (execution orchestration)**

- **Solves:** rule-based DAG construction and incremental re-execution.
- **Mechanism that decides it — cited.** `--touch` / `-t`: *"Touch output files (mark them up to date
  without really changing them) instead of running their commands. This is used to pretend that the
  rules were executed, in order to fool future invocations of snakemake. … Note however that you lose
  the provenance information when the files have been created in reality. Hence, this should be used
  only as a last resort."* Snakemake ships a documented "pretend it ran" switch. Its up-to-date decision
  uses `--rerun-triggers` (`code`, `input`, `mtime`, `params`, `software-env`) — i.e. it decides
  *whether to run*, and never asserts that an existing output equals what a run would produce.
- **By construction?** No. It is a useful **re-run harness** (§1(A) vehicle), not a guarantee.
- **Local:** absent (no binary, no module).
- **Cost:** a `Snakefile` with a rule per stage; invocation becomes `snakemake <target>`. Moderate.
- **Verdict: DIFFERENT PROBLEM.**

### 2.6 Nextflow — **DIFFERENT PROBLEM / weak PARTIAL**

- **Solves:** dataflow workflow execution with per-task isolated work directories, container execution,
  and a `-resume` cache keyed on script text + inputs + container.
- **Mechanism:** a task's outputs *are* produced inside its own work dir, which is closer to
  by-construction than Snakemake — but `publishDir` then copies/links results out, and nothing binds the
  **repo-committed** file to the published one. Editing the published copy afterwards is unpoliced.
- **By construction?** No, for the committed artifact.
- **Local:** absent. Requires a JVM.
- **Cost/blast radius:** highest ergonomic delta of the workflow engines — DSL2 rewrite of every
  invocation, channel semantics, and containers. ⛔ Running Mathematica in a Nextflow container is a
  licence-activation problem, not a packaging problem.
- **Verdict: DIFFERENT PROBLEM.**

### 2.7 Bazel — **CLOSES IT (Python side), at very high cost; PARTIAL for Mathematica**

- **Solves:** hermetic, sandboxed, content-addressed builds. Actions declare inputs/outputs; the action
  cache keys on input digests + command line; outputs land in `bazel-out/`, **not** the source tree.
- **Mechanism, precisely.** `linux-sandbox` makes *"the entire filesystem read-only except for the
  sandbox directory"*, and the sandbox strategy *"moves the known output artifacts out of the sandbox
  into the execroot."* So a file in `bazel-out/` was written by the declared action, and a hand-written
  file in the source tree is an **input**, never mistakable for an output. ⚠ Note a documented
  limitation that bears on the claim: the sandbox *"does not hide undeclared inputs"* — processes *"can
  freely access all files on the file system"*; the guarantee is about **modification**, not visibility.
  And *"you cannot easily run `linux-sandbox` inside a Docker container, unless you use
  `docker run --privileged`."*
- **How it actually closes P2 for a checked-in generated file:** the `diff_test` /
  `write_source_files` pattern. Per its docs, *"a `diff_test` target (`{name}_test`) is generated that
  ensures the source tree file or directory to be written to is up to date"* and *"will fail if the file
  is out of date and print out instructions on how to update it."* ⭐ Note what that means: the closing
  mechanism is **regenerate + diff** (§1(A)); the sandbox's contribution is only that the *freshly built*
  side is guaranteed genuine.
- **Local:** absent (`bazel`, `bazelisk` not on PATH).
- **Cost/blast radius: very high, and it changes how everything is invoked.** `MODULE.bazel`,
  `rules_python`, a Python toolchain, a `py_binary` + `genrule`/`write_source_files` per script for 91
  files, and every ad-hoc run becomes `bazel run //…`. For Mathematica: the licence/`mathpass` is an
  undeclared input; declaring it puts licence material into the action key and the sandbox, and a
  floating licence server needs network the sandbox is designed to remove.
- **Verdict: CLOSES IT (Python) / PARTIAL (Mathematica).** Disproportionate for a collection of
  standalone scripts run ad hoc.

### 2.8 Nix — **CLOSES IT (Python side); PARTIAL-BLOCKED for Mathematica; not installable here**

- **Solves:** purely functional derivations. With `sandbox = true` (default on Linux), builds are
  *"isolated from the normal file system hierarchy and will only see their dependencies in the Nix
  store, the temporary build directory, private versions of `/proc`, `/dev`, `/dev/shm` and
  `/dev/pts`"*, run *"in private PID, mount, network, IPC and UTS namespaces"*, and non-fixed-output
  derivations get **no network**. Output paths are store paths the builder alone writes.
- **By construction?** For the store path, yes — the same shape as Bazel and with a stronger
  environment guarantee (this would also neutralise the `sitecustomize.py`/`PYTHONPATH` `hashlib`
  interposition that §9 measured, because the build environment is constructed from the derivation, not
  inherited from a developer shell). Policing the **committed** copy still requires
  `nix build .#sidecar && cmp result <committed>` — §1(A) again.
- **Local:** absent, and installing Nix is a system-level change (daemon, `/nix` store) — out of scope
  by instruction.
- **Cost/blast radius:** lower than Bazel if used narrowly ("run this one script in a sandbox"): a flake
  with pinned nixpkgs providing python3 + sympy, and one derivation per stage. The **scripts need not
  change**. ⛔ Mathematica: `sandbox-paths` could bind-mount the licence in, but that either exposes
  licence material or breaks purity, and a network licence server is exactly what the sandbox removes.
- **Verdict: CLOSES IT (Python) / PARTIAL-BLOCKED (Mathematica).**

### 2.9 Observed-execution provenance — ReproZip / Sciunit / noWorkflow / strace-class — **PARTIAL**

- **Solves:** the one thing nothing else on the list does — recording *which process actually accessed
  or wrote a file during a run.* ReproZip traces OS syscalls to capture the files, binaries and
  dependencies a command touched; noWorkflow captures Python-level provenance (function activations,
  file accesses, dynamic slicing) without a workflow definition; Sciunit versions and replays execution
  histories.
- **By construction?** This is family §1(B): it can bind P2 **without** re-execution. But the trace is
  itself a file the producer's operator can fabricate unless it is signed or captured by an outside
  party — so it does not escape the key-custody question, it only relocates it. It is also strictly
  heavier to *verify* than simply re-running a 0.5–82 s script.
- **Local:** none installed.
- **Verdict: PARTIAL** — real prior art for the exact question, and the right thing to cite if
  re-execution ever becomes too expensive. It is not, here (§4).

### 2.10 reprotest / diffoscope — **DIFFERENT PROBLEM (complementary)**

- **Solves:** reproducibility testing — run a build twice under deliberately varied environments and
  diff the results (`diffoscope` renders the diff of arbitrary binary/structured files).
- **Relation to this hole:** these do not attest anything, but they are the correct instrument to prove
  that a regenerate-and-compare gate **cannot false-positive** (i.e. that the producer is deterministic
  across time, path, locale, `PYTHONHASHSEED`, umask). §4's measurement is a hand-rolled one-variable
  version of what reprotest does systematically.
- **Local:** both absent.
- **Verdict: DIFFERENT PROBLEM, complementary.**

### 2.11 Signed git commits / gitsign — **DIFFERENT PROBLEM**

Attests *who committed*. `commit.gpgsign=true` is already set here, so this layer exists. Says nothing
about P2.

### 2.12 Make / pytest / a plain CI job — **CLOSES IT** (see §4)

The "checked-in generated file is up to date" pattern: regenerate, then `git diff --exit-code` (or an
explicit `cmp`). This is the same mechanism the Bazel `write_source_files` rule implements, minus
hermeticity and caching. `make`, `git`, `gh` and a reachable `docker` daemon are all present locally.

---

## 3. LOCAL AVAILABILITY — MEASURED

Probed with `command -v` and `importlib.util.find_spec`; nothing installed.

| tool | status |
|---|---|
| `in-toto-run`, `in-toto-verify`, `in_toto`, `securesystemslib` | **absent** |
| `slsa-verifier`, `attest` | **absent** |
| `cosign`, `sigstore` (py) | **absent** |
| `dvc` (bin + module) | **absent** |
| `snakemake` (bin + module) | **absent** |
| `nextflow` | **absent** |
| `bazel`, `bazelisk` | **absent** |
| `nix`, `nix-build`, `guix` | **absent** |
| `diffoscope`, `reprotest` | **absent** |
| `podman`, `bubblewrap`, `firejail` | **absent** |
| `make` | **PRESENT** `/usr/bin/make` |
| `docker` | **PRESENT** `/usr/bin/docker` — **daemon reachable** (`docker info` exit 0) |
| `git`, `gh` | **PRESENT**; `gh` authenticated (repo query succeeded) |
| `gpg`, `ssh-keygen` | **PRESENT** |
| `python3` | **PRESENT** 3.10.12, sympy 1.14.0 |
| `math`, `wolframscript` | **PRESENT** `/usr/local/bin/math`, `/usr/bin/wolframscript` |

Repo context: `origin = git@github.com:trevnorris/toy_physics.git`, **visibility PUBLIC** (so GitHub
Actions minutes are free), **no `.github/` directory exists**. `commit.gpgsign=true`,
`user.signingkey=251CA676820DC7F3`. 6 `*.dimensions.txt` and 44 `mathematica/out/*.out` are git-tracked.

⚠ **Incidental hazard found while checking tracking:** `.gitignore:20` contains a bare `*.out`. The 44
existing `.out` files are tracked and therefore unaffected, but **a newly produced `.out` for a
newly-converted stage will be silently ignored** unless force-added. That is a provenance hazard in the
same family as this survey (the reference half of the universal gate can go missing without a signal).

---

## 4. ⭐ THE HONEST CHEAP OPTION

**Yes, there is a materially better cheap control than "a human re-runs it," and it is the same
mechanism the expensive tools use.**

### 4.1 The control

`scripts/check_sidecar_regeneration.py <stage>`:

1. copy `scripts/` to a temp directory outside the repo (or a `git worktree` of the commit);
2. **delete the target `.dimensions.txt` in the copy** — the committed bytes must not be an input;
3. clear `__pycache__` (§9's stale-`.pyc` landmine applies), run `python3 <stage>.py` in the copy;
4. `cmp` the emitted sidecar against the committed one; exit 1 on any difference or on a missing
   emission;
5. emit a structured marker (`SIDECAR_REGEN_OK|stage=…|sha256=…` / `SIDECAR_REGEN_MISMATCH|…`) in the
   same style as `MODULE_PIN_OK`.

Wire it into `compare_dimension_artifacts.py` alongside the existing pin call, and into
`run_all_audits.sh`. ⚠ **`run_all_audits.sh` must be fixed at the same time** — confirmed by reading
`:115-140`: it tallies `fail` and then reaches EOF with no `exit` derived from `$fail`, so it reports a
non-zero `Fail:` count and still exits 0. Adding a gate to a runner that swallows failures adds nothing.

### 4.2 Measured feasibility — this is not a hypothetical

I copied `scripts/` to the scratchpad, removed `__pycache__`, deleted each target sidecar in the copy,
and ran all six converted stages there (different directory, different time, committed sidecar absent):

| stage | exit | wall | sidecar vs committed |
|---|---|---|---|
| 004 | 0 | 0.5 s | **BYTE-IDENTICAL** |
| 011 | 0 | 2.0 s | **BYTE-IDENTICAL** |
| 012 | 0 | 15.8 s | **BYTE-IDENTICAL** |
| 013 | 0 | 82.3 s | **BYTE-IDENTICAL** |
| 016 | 0 | 11.3 s | **BYTE-IDENTICAL** |
| 018 | 0 | 10.6 s | **BYTE-IDENTICAL** |

**6 / 6 byte-identical; ≈ 123 s total for the whole converted set.** The emission is deterministic by
construction as well as by measurement: the header is stem + axes + two digests, records iterate an
insertion-ordered `dict`, and exponents render through `sp.sstr` of exact `Rational`s — no set
iteration, no float formatting, no timestamps, no absolute paths.

⚠ **Budget caveat:** stage013 alone is 82 s and ~30 more stages are in scope. Make this a **per-stage**
gate in the §4 loop (where it belongs — one stage at a time), not an always-on whole-corpus gate.

### 4.3 ⭐ Judged by the required standard: can the producing script forge it?

- **The comparison itself: no.** The checker computes it; the committed bytes are never fed to the run.
- **Self-echo (script reads its own committed sidecar and reprints it): defeated**, because step 2
  deletes it in the copy. ⭐ This is precisely why the clean-copy form matters and why the naive
  in-place `regenerate && git diff --exit-code` form is weaker — in-place, the committed file is present
  during the run.
- **Producer emits wrong values honestly: not this gate's job.** That is P3, owned by the cross-engine
  comparator (§11: 18/21 mismatches under an L↔M relabel) and the model read (c2). ⛔ Do not let anyone
  upgrade a green regeneration into a correctness claim — it is exactly the same category error §11
  already warns about for a bare `python3 …` run.
- **Environment interposition (a `sitecustomize.py` on `PYTHONPATH`, as §9 measured against the pin):
  still possible**, inherent to any in-process check. Cheap mitigation available *today*: run step 3
  inside a `python:3.10-slim` container with pinned sympy — the docker daemon is reachable here — which
  builds the interpreter environment from an image rather than inheriting the developer shell.
- **Deleting the checker's invocation: still possible.** Unchanged from §9's "the pin moves the single
  point of failure; it does not abolish it." ⭐ **Only an executor outside the repo fixes this**, and the
  cheapest one is already available: a GitHub Actions workflow (repo is public → free minutes) that
  checks out the commit under review, regenerates, and `git diff --exit-code`s the sidecars. The
  workflow file is itself at the commit, so deleting it is visible in the diff being reviewed.
  ⛔ Python-only — Mathematica cannot run on a hosted runner — but the Python side is exactly the side
  that lacks (g2).

### 4.4 Recommended ladder (cheapest first, each strictly stronger)

1. **`check_sidecar_regeneration.py` + fix `run_all_audits.sh`'s exit code.** ~60 lines. Converts an
   unverifiable human assertion into a mechanical, re-runnable, able-to-fail check. **Biggest single
   step; do this one.**
2. **Run step 1 inside the pinned container.** Removes the interpreter-environment interposition
   surface. ~5 extra lines, docker already reachable.
3. **Run step 1 in GitHub Actions on push.** Moves the executor outside the repo — the only thing that
   addresses "someone deletes the invocation."
4. *(Optional, unrelated to P2)* if P1 breadth ever matters — pinning sympy/Python and transitive
   imports rather than one module — **DVC or a lockfile** is the honest answer, not in-toto.

⛔ **Do not adopt Bazel, Nix, Nextflow, in-toto or SLSA for this.** Each either (a) closes P2 only by
performing step 1 with extra machinery, or (b) does not close P2 at all — and every one of them
requires restructuring how scripts are invoked, which is disproportionate for a directory of standalone
scripts run ad hoc, and each collides with the Mathematica licence seat.

---

## 5. VERDICT TABLE

| tool | one-line problem it solves | closes artifact→producer **by construction**? | local | verdict |
|---|---|---|---|---|
| **in-toto** | signed per-step record of declared materials/products | **No** — hashes paths before/after; `expected_command` mismatch is a *warning*, spec says the field "should not be trusted for security checks" | absent | **PARTIAL** (closes inter-step flow; P2 open) |
| **SLSA / slsa-github-generator** | platform-signed provenance: builder, source, digest | **No for the artifact** — L3 makes the *signature* unforgeable by build steps; a workflow that `cat`s a file still gets valid provenance; producer bad code explicitly out of scope | absent; no `.github/` | **PARTIAL / DIFFERENT** |
| **sigstore / cosign** | keyless signing + transparency log | **No** — signs whatever bytes it is given | absent | **DIFFERENT PROBLEM** |
| **DVC** | data versioning + dep-hashed pipeline stages | **No** — `dvc commit` records on-disk output hashes *without re-executing*, documented for "you ran it by hand" | absent | **DIFFERENT PROBLEM** (P1 upgrade only) |
| **Snakemake** | rule-based incremental re-execution | **No** — `--touch` documented as "pretend that the rules were executed, in order to fool future invocations" | absent | **DIFFERENT PROBLEM** (re-run harness) |
| **Nextflow** | dataflow engine, isolated work dirs, containers | **No** for the committed copy; publishes then leaves it unpoliced | absent | **DIFFERENT PROBLEM** |
| **Bazel** | hermetic sandboxed build + action cache | **Yes (Python)** via sandbox + `diff_test`/`write_source_files` — but the closing act is still regenerate-and-diff | absent | **CLOSES IT (py) / PARTIAL (wl: licence vs hermeticity)** — cost disproportionate |
| **Nix** | pure sandboxed derivations, no network, store outputs | **Yes (Python)** via `nix build` + `cmp`; also kills env interposition | absent (install = system change) | **CLOSES IT (py) / PARTIAL-BLOCKED (wl)** |
| **ReproZip / Sciunit / noWorkflow** | syscall/interpreter-traced execution provenance | **The only §1(B) family** — binds P2 without re-running, but the trace needs an outside signer | absent | **PARTIAL** |
| **reprotest / diffoscope** | determinism testing across varied environments | **No** — but validates that a regen gate can't false-positive | absent | **DIFFERENT PROBLEM** (complementary) |
| **signed commits / gitsign** | who vouched | **No** | present & already enabled | **DIFFERENT PROBLEM** |
| **⭐ regenerate-in-a-clean-copy + `cmp` (make / pytest / CI)** | "the checked-in generated file is what the producer emits" | **Yes** — and it is the mechanism inside Bazel/Nix/(g2) | `make`/`git`/`gh`/`docker` present | **CLOSES IT** (Python side) |

---

## 6. ⛔ WHAT I COULD NOT VERIFY

- **Every verdict on an absent tool rests on primary documentation, not on execution here.** I installed
  nothing. Specifically unverified by execution: in-toto's warning-only `expected_command` and its
  before/after hashing; `dvc commit`'s no-re-run semantics; Snakemake `--touch`; Bazel sandbox
  behaviour and `write_source_files`'s diff test; Nix sandbox behaviour; SLSA L3's key isolation.
- **I did not run Mathematica at all**, and did not test whether a licence seat can be provisioned
  inside a container, a Nix sandbox, a Bazel action, or a self-hosted Actions runner. The Mathematica
  columns above are reasoned from documented sandbox/hermeticity semantics plus the project's own
  standing constraint, not measured.
- **Determinism is measured on 6 stages only** — the 6 that currently emit sidecars — on this machine,
  in one locale/timezone/umask, under one `PYTHONHASHSEED` regime, with Python 3.10.12 / sympy 1.14.0. I
  did not vary those axes (that is what `reprotest` would do). Unconverted stages have no sidecar and
  could not be tested; long-running ones could change the cost picture.
- **I did not test the proposed control** — no script was written. The measurement in §4.2 exercises its
  core step (clean copy, absent committed artifact, byte compare) but is not the control itself.
- **Repo/CI facts I read but did not exercise:** Actions availability, runner minutes, and whether any
  branch protection exists were not checked beyond repo visibility.

---

## SOURCES

- [SLSA v1.0 threats](https://slsa.dev/spec/v1.0/threats) · [SLSA v1.0 levels](https://slsa.dev/spec/v1.0/levels)
- [in-toto specification](https://github.com/in-toto/docs/blob/master/in-toto-spec.md) · [in-toto spec v0.9](https://github.com/in-toto/docs/blob/v0.9/in-toto-spec.md)
- [DVC `commit`](https://doc.dvc.org/command-reference/commit)
- [Snakemake CLI](https://snakemake.readthedocs.io/en/stable/executing/cli.html)
- [Bazel sandboxing](https://bazel.build/docs/sandboxing)
- [Nix configuration `sandbox`](https://nix.dev/manual/nix/2.24/command-ref/conf-file.html)
- [Aspect `write_source_files`](https://docs.aspect.build/rulesets/aspect_bazel_lib/docs/write_source_files/) · [Bazel can write to the source folder](https://blog.aspect.build/bazel-can-write-to-the-source-folder)
- [ReproZip packing docs](https://docs.reprozip.org/en/1.1/packing.html) · [noWorkflow (VLDB)](https://dl.acm.org/doi/abs/10.14778/3137765.3137789) · [user-space provenance collectors benchmark](https://mir.cs.illinois.edu/marinov/publications/GraysonETAL24ProvenanceCollectorsSurvey.pdf)

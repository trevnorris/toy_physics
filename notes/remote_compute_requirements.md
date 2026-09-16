# Remote compute — requirements (for the research pass)

**What this is:** requirements for renting on-demand compute so the orchestrator (Claude) and Codex can
**boot a box → run a memory-bound CAS job → tear it down**, with our installs preserved between runs so we
never re-provision from scratch. **Status:** requirements only — hand this to a web-enabled research pass to
pick a provider + concrete setup. **Non-physics infra doc:** changes no computation, premise, check, or claim,
so no review legs (CLAUDE.md "Other" row). Author: orchestrator, 2026-09-15.
**Rev 2026-09-15:** §8 Wolfram licensing **resolved** → on-demand Wolfram Engine paid from existing Service
Credits. §3, §6, §7, §9, and §10 updated to match.
**Rev 2026-09-15 (b):** folded in the external review of the research plan:
- §1: memory sizing at a stated concurrency.
- §3: lifecycle rewrite plus layered billing guard.
- §4: detached-launch fix.
- §7: exact Wolfram pin.
- §8: spend ceiling and env-var injection.
- §9: failure-mode demos.
- §10.1: research errata.

---

## 1. Why (the concrete driving jobs — so research sizes it right)

The dev box **cannot** run our heaviest CAS jobs; they OOM. The binding constraint is **RAM**, not GPU (the
dev machine's GPU is disabled and these are FP64 symbolic/numeric, not GPU work). The jobs waiting on this:

1. **Near-term / the acceptance job:** the deferred **cross-engine Wolfram (WL) comparison** for S11c-b / c1 /
   c2 against the repaired exports. Scripts already exist (`research/pde_ledger_v3/mathematica/S11c_b_brane_operator_mathematica_audit.wl`
   and siblings); est. **≥64 GB**; never run because the dev box can't hold it. This is the instrument that
   catches the S11c-b sign/orientation defect class — high value, blocked purely on RAM.
2. **Deferred:** c1 "giants" + full self-energy residual (S11c-e), also **≥64 GB**, possibly more.
3. **General:** any heavy SymPy / Mathematica batch run that OOMs locally (e.g. the source/profile quadrature
   Codex is currently parallelizing).

**Job shape:** memory-bound, long-running (minutes → hours), **batch/non-interactive**, emits `.out`
transcripts that belong in GIN. No GPU. No human in the loop during the run.

**⚠ Memory sizing.** The "≥64 GB" estimates don't say at what concurrency they apply.
- Parallel workers each hold their own large symbolic expressions, so memory scales with the worker count. A 64 GB
  box is **not** automatically enough for a "~64 GB" job run with 4 workers.
- **Declare the worker count per job.**
- **Measure peak memory at that concurrency** by sampling system-wide used memory during the run.
- Don't rely on per-process max RSS (`/usr/bin/time -v`): it reports only the largest single process and
  undercounts parallel jobs.

## 2. Hard requirements

| Axis | Requirement |
|---|---|
| **RAM** | **≥64 GB floor** (stated threshold for the WL comparison). Must be able to select **128 GB** for the deferred giants. Size per job. |
| **CPU** | Multi-core — we parallelize quadratures (Codex just split one into 4 workers). **≥8 vCPU, prefer 16.** |
| **Disk** | ≥100 GB. Repo clone + `datalad get` of needed inputs + heavy WL working space (tens of GB). |
| **GPU** | **None** — do not pay for it. |
| **Arch** | **x86-64** (Wolfram + our pinned Python wheels; avoid ARM unless verified). |
| **OS** | Linux, Ubuntu/Debian 22.04-class base (matches dev box). |
| **Billing** | **Per-hour or per-second; cheap when torn down.** Idle cost = only the preserved image/volume storage, not a running instance. |
| **Automation** | **Fully scriptable, non-interactive** create/start/stop/destroy (CLI or API the agent calls from bash). ⛔ No web-console clicking — the agents can't click. |

## 3. Lifecycle the agents must be able to drive (fire up → run → bring down)

All steps must run from a **non-interactive script** an agent invokes. Every run gets a **run ID**
(`<UTC timestamp>-<job>`), which is used in labels, paths, and the run manifest.

1. **Pre-flight.**
   - Check for existing ephemeral resources (§3.1).
   - Create the per-run Wolfram entitlement (§8.3).
   - Fix the job's instance size **and worker count** (§1).
2. **Create** an instance from the **preserved image** (§5), sized for the job (64 or 128 GB).
   - Labels: `role=ephemeral`, `run=<id>`, `expires=<epoch>`.
   - Cloud-init must install our login key into the **`cas` user's** `authorized_keys`. The image disables root
     login, and provider key injection typically targets root (verify on the chosen provider).
3. **Wait for readiness.**
   - Get the IP, then retry an **authenticated** `ssh -o BatchMode=yes cas@<ip> true` until it succeeds.
   - Enforce a hard deadline (e.g. 5 min); on timeout, tear down.
   - An open port 22 is **not** readiness.
   - Use a **per-run `known_hosts`** file with `StrictHostKeyChecking=accept-new`. Every box has new host keys and
     IPs get reused, so the global `known_hosts` would fail non-interactively.
4. **Provision secrets** (§6) over SSH into a `cas`-owned `0600` file on tmpfs (`/run/cas/secrets.env`), sourced
   by the job wrapper.
   - ⚠ An `export` in one SSH session does **not** persist into later sessions or into the detached job.
5. **Sync.**
   - `git clone` at a **pinned commit**.
   - **Configure and verify the GIN sibling + annex access on the box.** The dev box's git/annex config does not
     travel with a fresh clone.
   - `datalad get` the inputs.
   - Write a **run manifest** recording:
     - commit SHA and input annex keys;
     - Wolfram `$VersionNumber` and SymPy version;
     - worker count and instance type.
6. **Run detached** through a wrapper that:
   - writes a **status file containing the exit code** (not a bare `DONE`);
   - fully detaches stdio (§4).

   Prefer a transient systemd unit (`systemd-run`), which also surfaces OOM kills. Sample system-wide memory use
   during the run (§1). Poll the status file.
7. **Validate.** Success requires **both** exit code 0 **and** the audit's own checks passing.
   - Audit scripts must exit non-zero on any failed check (`Exit[1]`).
   - A process that finished is not an audit that passed.
8. **Retrieve — on success *and* failure, before any teardown.** `rsync` the log, status file, manifest, memory
   samples, and outputs back to the dev box. This is the recovery copy: never delete the box before this step has
   been attempted.
9. **Publish** (success path only), then **verify** that content is retrievable from GIN (`git annex whereis`
   lists gin, or `datalad get` from a fresh clone):
   1. `datalad save` the outputs. `push` only publishes already-saved state.
   2. `datalad push --to gin` (content + pointers).
   3. `git push origin`.

   If publishing fails, the step-8 copy on the dev box is the recovery path.
   - *Option:* publish **from the dev box** after step 8 instead. The box then needs only read access to GIN and
     no git push credentials.
10. **Teardown.**
    - **Delete** the instance; do not power it off.
    - Delete any billable IP that isn't auto-deleted.
    - Delete the per-run entitlement and record `CreditsSpent` in the manifest.
    - **Confirm** that no `role=ephemeral` servers, primary IPs, or volumes remain. An empty server list alone does
      not prove billing stopped.

### 3.1 ⚠ Orphaned-billing guard (mandatory, layered)

- **Labels.** Ephemeral resources carry `role=ephemeral`. The preserved image carries `role=image` **plus provider
  delete protection**. `down`, `nuke`, and the reaper select **only** `role=ephemeral`, never images.
- **Pre-flight.** Refuse to create a box if an ephemeral resource for the same job already exists. List and
  destroy must be idempotent.
- **Controller trap** (`trap … EXIT`) — best effort only. It cannot run if the dev box loses power or the
  controller is `kill -9`'d.
- **Independent reaper — required.**
  - A scheduled job on **neither** the dev box nor the compute box that deletes any `role=ephemeral` resource past
    its `expires` label.
  - It is the only layer that covers controller loss **and** a hung VM.
  - Candidate: a scheduled CI workflow holding a provider token. Scheduled CI runs can be delayed, and some
    platforms disable schedules on inactive repos, so choose one whose schedule we can monitor.
  - Alternative: a provider with enforced max-lifetime deletion (e.g. GCP `--max-run-duration` +
    `--instance-termination-action=DELETE`).
- *Optional* **on-box self-delete** at `expires` via the provider API.
  - It must **delete**, not `shutdown`: a powered-off Hetzner server still bills.
  - It puts a provider token on the box, so run this compute in a **dedicated provider project** to limit what
    that token can touch.
- **Wolfram credits.** Kernels keep charging until they are terminated. The per-run `EntitlementExpiration` bounds
  that exposure (§8.2–8.3).
- A provider spending alert, where available, as the final backstop.

## 4. Access model — recommended: box is a **dumb compute target**

Keep the agents **local**; the box just runs jobs.

- **Claude + Codex reach the box over SSH** as `cas` (key-based, non-interactive) from the dev machine. The
  local detached-launch pattern needs two fixes on the box:
  - **Detach all stdio.** `setsid` alone leaves the inherited stdio attached, so ssh can hang. Redirect the
    **whole** wrapper: `ssh box 'setsid bash -c "…" </dev/null >/dev/null 2>&1 &'`, or launch via `systemd-run`.
  - **Record the exit code, not a bare marker.** The wrapper ends by writing its exit code to a status file.
    Poll that file and rsync the log (§3 steps 6–8).
- Codex needs `--sandbox danger-full-access` **locally** to run ssh/rsync; **on** the box, sandboxing is moot
  (ephemeral, isolated) but run jobs as a **non-root** user for hygiene.
- ⛔ We do **not** need to install Codex/Claude on the box (simpler, no agent auth to manage remotely).
  *Alternative if we later want a resident agent: install the CLIs on the box — more setup, note it but don't
  build it now.*

## 5. Persistence — "preserve our core files + installs"

Split by what changes at what rate; **GIN + git already ARE our data/code persistence**, so the box is nearly
stateless — only the **installs** need preserving:

- **Installs → a custom image / snapshot, baked once** (OS + full toolchain, §7). Boot from it every time;
  rebake only when the toolchain changes. This is what removes the "reinstall everything" cost.
- **Code → `git clone`** at boot (already versioned; fast).
- **Data → GIN via datalad** — `datalad get` inputs, `datalad push --to gin` outputs. The `.out` store is
  ~370 MB total; the box needs no persistent data disk.
- **Secrets → injected at boot** (§6), never in the image.
- *Optional:* a small persistent volume only if boot-from-clean proves too slow (pip/datalad caches).
  Prefer to avoid — image + GIN should cover it.

Research should confirm the chosen provider supports **custom images/snapshots** cleanly (vs. Docker container
vs. attached volume for installs) and how cheap idle image storage is.

## 6. Secrets & security

- **SSH key-only**, no passwords; restrict inbound to SSH, ideally from our IP. Outbound internet must stay
  open: the Wolfram kernel contacts Wolfram's license server during runs (§8.4), and the box also talks to GIN
  and git.
- Secrets are **injected at boot** (into `/run/cas/secrets.env`, §3 step 4) and **never baked** into the image: **GIN token** (datalad keyring credential name
  `gin`), **Wolfram on-demand entitlement ID** (§8, preferably per-run), git creds, SSH keys. Box is destroyed at teardown, so they don't
  linger.
- Non-root job user.

## 7. The toolchain to bake (the "core set of installs") — pinned to the dev box

| Component | Pin (dev box, 2026-09-15) | Notes |
|---|---|---|
| Python | **3.10.12** | + venv |
| SymPy | **1.14.0** | ⭐ must match — the exports chain is built against it |
| mpmath | **1.3.0** | |
| numpy | **1.26.4** | |
| scipy | **1.15.3** | |
| **Wolfram Engine** | dev box: full Mathematica at `/usr/local/Wolfram/Wolfram/`, driven via **`wolframscript` / kernel**. Box: headless **Wolfram Engine** + `wolframscript`, **exact same version + build as the dev box** (`$VersionNumber` / `$ReleaseNumber`: TBD — record them) | ⭐ **No activation in the image.** Licensed at runtime by on-demand entitlement (§8) |
| datalad | **0.19.5** | GIN `.out` I/O |
| git-annex | **10.20221003** | |
| git, rsync, ssh, openssh-client | present | |
| build tools | gcc/make | native wheels |

Don't prune Wolfram paclets from the image until the real audit scripts have passed on it.

We have **no `requirements.txt`** in the repo — the bake list above is the source of truth; research/setup
should turn it into a pinned lockfile.

## 8. ⭐ Wolfram licensing — RESOLVED (2026-09-15): on-demand Wolfram Engine, paid from Service Credits

**Decision:** on the box, run the `.wl` audits under the **headless Wolfram Engine** via `wolframscript`, licensed
by **on-demand licensing** (pay-as-you-go against the Service Credits already on our Wolfram Personal Plan: **200
credits available** as of 2026-09-15).
- Nothing is activated inside the image.
- No seat transfer, no license server.

Sources: Wolfram docs for `CreateLicenseEntitlement`, `LicenseEntitlementObject`, and the WolframScript reference;
Wolfram's `WL-FunctionCompile-CI-Template` repo; the `wolframresearch/wolframengine` Docker image docs.

### 8.1 Routes considered

| Route | Verdict | Reason |
|---|---|---|
| **On-demand Engine (Service Credits)** | ✅ **chosen** | Built for ephemeral cloud batch. The entitlement is injected at runtime, and the credits are already paid for. |
| Free Wolfram Engine (node-locked) | ⛔ rejected | Activation is limited to ~2 machines per account, so destroy-every-run would burn activations. Its terms (pre-production dev / personal non-commercial) are *arguable* for this work — not clearly prohibited — but the activation limit decides it regardless. |
| Move the 2-seat Mathematica license onto the box | ⛔ rejected | Single-machine activation needs a system transfer per move, which doesn't work with destroy-every-run. |
| MathLM network license | fallback only | Clean, and it would enforce the 2-kernel cap, but it needs a network license plus a reachable license server. Revisit only if credit spend becomes material. |

### 8.2 Cost and spend ceiling — ✅ RATE VERIFIED LIVE (2026-09-15)

**Rate — verified against our own account** (created + inspected + deleted a real entitlement, 0 credits spent):
- **10 credits per kernel-hour** under the **Standard** policy (`WSDS`), billed in **900 s intervals**.
- **Parallel subkernels cost the SAME** (`KernelCosts -> <|Standard -> 10., Parallel -> 10.|> Credits/Hours`).
- ⛔ The **4-credit CI-template sample does NOT apply to us** — the research report's 10 was correct.

**What 200 credits buys at 10/kernel-hour**
- **~20 single-kernel hours** in total.
- Acceptance run (~3 h, **1 kernel**): **~30 credits**.
- 1 main kernel + 4 subkernels: **50 credits/h**.
- ⭐ **SymPy/Python parallelism costs 0 Wolfram credits** — only `LaunchKernels`/`Parallel*` *in a `.wl`* costs.
  Keep WL audits single-kernel (`ParallelKernelLimit -> 0`); parallelize SymPy freely.

**⭐ Credits are effectively unconstrained (2026-09-15).** Top-ups are ~**5000 credits for $25** ($0.005/credit):
acceptance run ≈ **$0.15**, worst-case parallel run ≈ **$0.75**. ⇒ **the real cost/risk envelope is the Hetzner box,
not Wolfram** — the scarcity concerns below (worst-case ceiling vs balance, zero-balance urgency) relax once topped
up; keep the per-run entitlement expiry only as a hung-kernel guard.

**How billing works**
- Creating an entitlement is free.
- A charge is applied when a kernel starts, and charging continues **until the kernel terminates**.
- A hung kernel keeps billing.

**⚠ `CheckCreditsBalance` is not a ceiling**
- It only refuses to create an entitlement if the balance can't cover the configured kernel limits for **1 hour**.
- It doesn't reserve credits for the whole job and doesn't cap spend. Keep it on anyway.

**The real ceiling is the entitlement itself.** Worst case per run is roughly:

> (`StandardKernelLimit` + `ParallelKernelLimit`) × rate × `EntitlementExpiration` (h),
> **plus** one billing interval per kernel start (a restart loop costs extra).

Size the limits so this fits the budget. Example: 1 + 4 kernels × 4 × 10 h = **200 credits**, which is the entire
current balance.

**❓ Open**
- Do the Personal Plan's credits renew (monthly or annually), or are they a one-time balance? Do unused credits roll
  over?
- What happens to a running kernel when the balance hits zero? Ask Wolfram; don't test by draining the balance.

### 8.3 Entitlement settings (critical) — one entitlement per run

```wolfram
e = CreateLicenseEntitlement[<|
  "StandardKernelLimit"   -> 1,
  "ParallelKernelLimit"   -> 0,                      (* raise ONLY if the script uses Parallel*/LaunchKernels *)
  "LicenseExpiration"     -> Quantity[8, "Hours"],   (* per-kernel max runtime *)
  "EntitlementExpiration" -> Quantity[10, "Hours"]   (* per-run; worst case here: 1 × 4 × 10 = 40 credits *)
|>];
e["EntitlementID"]
```

**Expiration settings**
- ⚠ **`LicenseExpiration` caps each kernel's runtime.** Wolfram's template uses 1 h; copying that would kill a
  multi-hour audit mid-run.
  - Set it above the longest expected job.
  - 8 h covers the near-term targets; raise it for the S11c-e giants.
  - Recompute the §8.2 worst case whenever you raise it.
- Kernels **cannot outlive `EntitlementExpiration`**, whatever `LicenseExpiration` says.

**Kernel limits**
- `ParallelKernelLimit` defaults to **0**.
- Set it to the job's declared worker count (§1), and include the subkernels in the §8.2 worst case.

**Lifecycle, driven by the wrapper (§3)**
1. At pre-flight: **create** the entitlement on the dev box. This briefly uses a dev-box seat kernel via
   `wolframscript`.
2. At teardown: **record** `e["CreditsSpent"]` in the run manifest, then **delete** the entitlement with
   `DeleteObject[e]`.

A long-lived entitlement kept in the local secret store is simpler. The cost is that a leaked ID or an orphaned
kernel isn't time-bounded.

### 8.4 Injection & runtime

**Treat the entitlement ID as a license key.** Handle it like the GIN token (§6):
- write it into the box's secrets file (§3 step 4);
- never put it in the image or the repo.

**Pass it via the environment, not argv.**
- The job wrapper sources the secrets file, which sets **`WOLFRAMSCRIPT_ENTITLEMENTID`** (documented WolframScript
  env var), then runs `wolframscript -file <script>.wl`.
- Prefer this over `-entitlement <ID>`: command-line arguments are visible to other users via `ps`.
- `LicensingSettings` configures `RemoteBatchSubmit` jobs; it does **not** license our SSH-launched process.

**Build time is unlicensed.**
- The image bake has **no** entitlement, so bake-time checks are limited to "binary present, version matches".
- The licensed smoke test runs only on a booted box (§8.5).
- Never leave a build-time entitlement in the snapshot.

**The box needs outbound internet.**
- On-demand kernels regularly contact Wolfram's license server.
- Inbound access stays SSH-only.
- ⚠ New failure mode: a license-server or network blip mid-run. The status file (§3 step 6) must capture it; watch
  the first long runs.

**Seat cap**
- On-demand kernels count against the entitlement's own limits, not the 2 Mathematica seats. This is our reading of
  the docs; confirm it.
- The "never exceed 2 concurrent Mathematica kernels" rule still governs **seat-licensed** kernels on the dev box,
  including the brief entitlement-creation kernel.

**Version pin**
- Install the Engine at the **exact version and build** of the dev box's Mathematica (§7).
- This keeps a WL version change from confounding results against dev-box WL runs.

### 8.5 Verify before the first real run (cheap, falsifiable)

1. **Rate check (free).** ✅ **DONE 2026-09-15** — created/inspected/deleted a real entitlement, 0 credits spent,
   `LicenseEntitlements[]` clean after. Result: **10 credits/kernel-hour (Standard `WSDS`), 900 s interval, parallel
   = same rate.** §8.2 updated. `$CloudConnected = True`.
2. **Licensed smoke test on the first booted box.** ⚠ still pending (needs a booted box).
   - Run `WOLFRAMSCRIPT_ENTITLEMENTID=<ID> wolframscript -code '1+1'`.
   - Check `e["CreditsSpent"]` on the dev box. Expect ~one billing interval (~1 credit at the §8.2 rate).
   - If any job will use parallel workers, also confirm that `LaunchKernels[n]` gets `n` subkernels licensed under
     the entitlement.
3. **Credit terms.** ⚠ still open (owner/Wolfram question) — renewal/rollover and zero-balance behavior
   (§8.2 open items). ⚠ **now more material at 10/hr:** a worst-case parallel entitlement can exceed the 200 balance.
4. **Front-end check.** ✅ **DONE 2026-09-15 — both false positives.** `S11b_interface_coupling_law` line 485–486
   matched the *variable* `omegaDynamic[...]`; `S11c_a_interface_geometry` line 1737 matched `rawDynamic[...]` and
   the string `"DYNAMIC"` — not the `Dynamic[` primitive. **All `.wl` audits are headless-clean.**

## 9. Acceptance test (how we know the whole setup works)

**Prerequisites — small, cheap, and done first.** Each one tests one of the most consequential failure modes:

1. **§8.5 checks.** Credit rate, licensed smoke test on a booted box, credit terms, and the front-end grep.
2. **Failure demo.**
   - Run a deliberately failing job (a script that `Exit[1]`s).
   - The status file and the validation step must report it as **failed**.
   - Its logs must come back to the dev box.
3. **Controller-loss demo.**
   - Start a short run with a near `expires`, then `kill -9` the controller.
   - The **independent reaper** (§3.1) must delete the box.
   - No `role=ephemeral` resources may remain.
4. **Publish-failure demo.**
   - Break the GIN push (e.g. a bad credential).
   - Results must be preserved on the dev box and recoverable.
   - Teardown must still happen.

**End-to-end — one fully scripted run, zero interactive steps.** Each arrow is a step; every step must succeed
before the next starts:

> 1. Pre-flight + per-run entitlement
> 2. Boot from the preserved image — **128 GB if peak memory at the chosen worker count is still unmeasured**,
>    otherwise size from the measurement
> 3. Authenticated SSH readiness
> 4. Secrets file
> 5. `git clone` at the pinned commit
> 6. `datalad get` the S11c-b/c1/c2 exports
> 7. Run `S11c_b_brane_operator_mathematica_audit.wl` (+ c1/c2) at a **stated worker count**
> 8. Status file shows exit 0 **and** the audit checks pass
> 9. Retrieve the bounded `.out`, the manifest, and the peak-memory record
> 10. `datalad save`
> 11. `datalad push --to gin` + `git push origin`
> 12. Verify retrieval
> 13. **Destroy the box** and delete the entitlement
> 14. Confirm:
>     - no `role=ephemeral` resources remain;
>     - **zero reinstalls** were needed;
>     - `CreditsSpent` is consistent with §8.2.

That job is also item #1 on our real work queue, so the acceptance test does double duty.
The recorded peak memory sizes all later runs of these scripts.

## 10. Open questions for the research pass

1. **Provider** best fitting: 64–128 GB RAM CPU instances, per-hour/second billing, custom image/snapshot
   support, a **mature scriptable CLI/API**, cheap idle image storage, region near GIN (EU, `gin.g-node.org`).
   Candidates to compare on $/hr @64 GB and @128 GB, snapshot $/mo, CLI maturity: **Hetzner** (cheap, EU),
   **DigitalOcean**, **AWS EC2 r-series**, **GCP**, **Vultr**, **Linode/Akamai**, **Fly.io**. (RunPod is in
   our notes but GPU-focused — check if it has big-RAM CPU instances worth it.)
2. ~~**Wolfram licensing** (§8) — the highest-leverage unknown.~~ **Resolved 2026-09-15** → on-demand Engine
   (§8). Remaining: confirm the credit rate, the renewal terms, and the front-end grep (§8.5).
3. **Install-persistence mechanism** the provider supports best: custom image/snapshot vs. Docker vs. volume.
4. **Data locality** — is `datalad get`/`push` to GIN (EU) fast enough from the provider's region? Pick region
   accordingly.
5. **Automation surface** — provider CLI (`hcloud`/`doctl`/`aws`/`gcloud`) vs. Terraform vs. a thin wrapper.
   Whichever the agent can drive most reliably non-interactively.
6. **Cost envelope** — rough $/run for a ~64 GB, few-hour WL job, plus monthly idle image storage, so we know
   the standing cost of keeping this capability warm.

### 10.1 Research outcome & review errata (2026-09-15)

**Research pick.**
- Primary: Hetzner Cloud CCX (CCX43 = 64 GB, CCX53 = 128 GB), German regions.
- Fallback: AWS r7i or GCP n2-highmem in Frankfurt.
- The corrected plan is **`remote_compute_plan.md`**. It replaces the original research report and is consistent
  with this doc; if the two ever conflict, this doc wins.
- The corrections below are already applied in that plan.

Corrections from the external review:

- **Pricing method.**
  - Hetzner bills per **started** hour and publishes its own USD tariff. Use the account-currency tariff, not a
    EUR→USD conversion.
  - The review reports **CCX43 $0.5216/h** and **CCX53 $1.0088/h**; a third-party tracker lists CCX43 from
    $0.5138/h. Confirm against the live price list.
  - A 3 h job plus provisioning and publishing usually means **4 billable hours**: ≈ **$2.09 (CCX43) /
    $4.04 (CCX53)**. That excludes Wolfram credits, IPv4, and VAT.
  - The report's per-run totals were too optimistic.
- **Linode/Akamai** was wrongly excluded. The review reports official High Memory plans up to 300 GB. Re-check them
  against §2 (vCPU, disk, x86-64, EU region).
- **Credits.**
  - The report's 10 credits/kernel-hour figure is unverified, as is the 4-credit sample (§8.2).
  - Its "buy credits" step ignores the existing 200-credit balance.
- **Watchdog.** The report's `shutdown` option leaves a Hetzner box billing. See §3.1.
- **`nuke`.** The report's version deletes the preserved snapshot, which contradicts §5. See §3.1 (labels +
  delete protection).
- **Image bake.** The report's Packer smoke test runs `wolframscript` unlicensed. See §8.4–8.5 (licensed smoke test
  only after boot).
- **Launch/publish.** The report's steps omit the exit-code status file, `datalad save`, the GIN sibling setup, and
  failure-path retrieval. See §3.
- **Licensing language.** The "organizational use ⇒ prohibited" claim is softened. See §8.1.
- **Front end.** "Front-end features are irrelevant to every batch script" is too broad. See §8.5 step 4.
- **Paclets.** Don't prune them before the real scripts pass (§7).

## Related
Implementation plan: `remote_compute_plan.md` (same directory).
GIN/datalad `.out` policy: root `CLAUDE.md` (annex + GIN section) + the `project-datalad-gin-out-storage`
memory. GPU-disabled dev box: `project_gpu_disabled_machine.md`. Mathematica 2-seat: the mathematica-seat
memory. Detached-launch pattern: the background-process-launch / `codex-exec-hangs-on-stdin` memories.

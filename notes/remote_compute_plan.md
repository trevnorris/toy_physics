# Remote compute — research plan (revised 2026-09-15)

**What this is:** the provider and setup plan for renting memory-bound, CPU-only boxes that agents boot, use, and
destroy for CAS jobs.
- It **replaces the original research report** and folds in everything from the external review.
- Companion doc: `remote_compute_requirements.md` (rev b). The two are consistent; if they ever conflict, **the
  requirements doc wins.**

**For the implementing agent:**
- Build in the order given in §7.
- Anything marked ❓ is unverified or undecided. Resolve it or ask; **don't** hard-code a guess.
- Anything marked ⚠ is a known way to lose money or results silently.

---

## 1. Summary

| Decision | Choice | Status |
|---|---|---|
| Provider | **Hetzner Cloud, dedicated-vCPU CCX**, German location (`nbg1`/`fsn1`, near GIN) | Chosen |
| Sizes | **CCX43** (16 vCPU / 64 GB / 360 GB), **CCX53** (32 vCPU / 128 GB / 600 GB) | Chosen |
| Fallback | AWS EC2 r7i (eu-central-1) or GCP n2-highmem (europe-west3); Linode/Akamai High Memory re-opened (§3) | Not needed yet |
| Wolfram | **Headless Wolfram Engine + on-demand licensing**, paid from the existing **200 Service Credits** (Personal Plan); **one entitlement per run** | Chosen; rate unverified ❓ |
| Install persistence | **Packer-baked Hetzner snapshot**, labeled `role=image`, delete-protected | Chosen |
| Automation | **Thin bash wrapper (`cas`) over the `hcloud` CLI**; no Terraform | Chosen |
| Billing safety | Layered: controller trap + **independent reaper (required)** + optional on-box self-delete + per-run entitlement expiry | Reaper host undecided ❓ |
| Publishing | **Recommended: publish from the dev box** after retrieval (requirements §3 allows either) | Recommended |

**Rough cost of a ~3 h job:**
- Compute: ≈ **$2.09** on CCX43 or **$4.04** on CCX53, since Hetzner bills 4 started hours.
- Wolfram: ~12 credits from the existing balance (at the unverified sample rate).
- Extra: IPv4 and tax.
- Idle cost: snapshot storage only.

---

## 2. Job profile (drives sizing)

**Job shape**
- Memory-bound, CPU-only (FP64 symbolic/numeric work), batch, minutes to hours, no human in the loop.
- The first job is the cross-engine WL audit for S11c-b/c1/c2, estimated at **≥64 GB**.
- Later jobs (S11c-e giants) may need 128 GB or more.

**⚠ Memory scales with the worker count.** Parallel workers each hold their own large expressions, so a "~64 GB"
estimate doesn't make a 64 GB box enough for a 4-worker run.

**Rules**
1. Every job declares its **worker count**.
2. Measure peak memory **at that worker count** by sampling system-wide used memory (`MemTotal − MemAvailable`)
   every ~10 s.
   - Per-process max RSS (`/usr/bin/time -v`) reports only the largest single process, so it undercounts parallel
     jobs.
   - 10 s sampling can miss short spikes. After a failure, check `dmesg` / `journalctl -k` for OOM kills.
3. If peak memory at the chosen worker count is unknown, **run first on CCX53 (128 GB)**, record the peak, and size
   later runs from it.

---

## 3. Provider comparison

Prices are for Linux on-demand instances, excluding tax. "Research" means the figure comes from the original
research pass and was **not re-verified** in review.

| Provider / shape | vCPU / RAM | $/hr @64 GB | $/hr @128 GB | Billing | Bills when stopped? | EU near GIN | CLI | Notes |
|---|---|---|---|---|---|---|---|---|
| **Hetzner CCX43 / CCX53** | 16/64; 32/128 | **$0.5216** (review) / €0.4423 (tracker) | **$1.0088** (review) / €0.8550 (tracker) | Per **started** hour, monthly cap | **Yes — must delete** | nbg1, fsn1 | `hcloud` | Dedicated AMD EPYC, x86-64. A tracker lists CCX43 from $0.5138/h. Use the **account-currency tariff**, not an FX conversion. ❓ Confirm on the live price list. |
| AWS r7i.2xlarge / 4xlarge | 8/64; 16/128 | ~$0.59 Frankfurt (research) | ~$1.18 (research) | Per second (60 s min) | No (EBS only) | eu-central-1 | `aws` | Large instances may need a vCPU quota raise. Spot is cheaper but interruptible. |
| GCP n2-highmem-8 / -16 | 8/64; 16/128 | ~$0.52 on-demand (research) | ~$1.05 (research) | Per second (60 s min) | No (PD only) | europe-west3 | `gcloud` | **Provider-enforced max lifetime** (`--max-run-duration` + `--instance-termination-action=DELETE`). |
| DigitalOcean Memory-Optimized | 8/64; 16/128 | ~$0.50 (research) | ~$1.00 (research) | Per second (research) | **Yes** | FRA1 | `doctl` | Pricier per GB than Hetzner. |
| Vultr | 16/64; 32/128 | ~$0.48 (research) | ~$0.96 (research) | Hourly | **Yes** | Frankfurt | `vultr-cli` | Check big-RAM stock. |
| **Linode/Akamai High Memory** | up to **300 GB** (review, official plans) | ❓ | ❓ | Hourly, monthly cap | Yes | Frankfurt | `linode-cli` | **The original report wrongly excluded this.** Check vCPU, disk, and region against requirements §2. |
| Fly.io / RunPod | ❓ | ❓ | ❓ | Per second | varies | fra / EU | `flyctl` / API | Not verified for CPU-only ≥64 GB. RunPod is GPU-oriented. Excluded for lack of evidence, not disproven. |

**Why Hetzner**
- It is the cheapest dedicated-vCPU option in Germany at these sizes.
- It is physically close to GIN.
- `hcloud` is scriptable end to end (JSON output, labels, firewalls, image protection, official Packer builder).

**Caveats**
- Hetzner charges for powered-off servers. That doesn't matter as long as every path **deletes** the server.
- Prices went up several times in 2026, and rescales are billed at current prices.
- ❓ New accounts may have low server limits. Request a raise before the first CCX53 run.

---

## 4. Wolfram licensing

### 4.1 Decision: on-demand Wolfram Engine

The `.wl` audits run through `wolframscript` with no notebook front end. On the box they run under the **headless
Wolfram Engine**, licensed per run by an **on-demand entitlement**. Usage is billed against the Service Credits
already on the Personal Plan (**200 available** as of 2026-09-15).

| Route | Verdict | Reason |
|---|---|---|
| **On-demand Engine** | ✅ | Designed for ephemeral cloud batch. Nothing is activated in the image; the ID is injected at runtime; the credits are already paid for. |
| Free Wolfram Engine | ⛔ | Node-locked activation is limited to about two machines per account, so destroy-every-run would burn activations. Its terms (pre-production development, personal non-commercial projects) are **arguable** for this personal research, not clearly prohibited. The activation limit decides it regardless. |
| Move the 2-seat Mathematica license | ⛔ | A single-machine license needs a system transfer per move. |
| MathLM network license | Fallback | It needs a network license and a reachable server. Revisit only if credit spend becomes material. |

### 4.2 Cost and spend ceiling

**Rate — ✅ VERIFIED LIVE 2026-09-15:** **10 credits/kernel-hour** (Standard `WSDS`; standard *and* parallel same),
**900 s** interval. (The CI-template 4-credit sample does not apply to us; the original report's 10 was right.)

**⭐ Credits are effectively unconstrained.** Top-ups are ~**5000 credits for $25** ($0.005/credit):
- A 3 h single-kernel run (~30 credits) ≈ **$0.15**; 1 main + 4 subkernels for 3 h (150 credits) ≈ **$0.75**.
- ⇒ **The real cost/risk envelope is the Hetzner box** (§5), not Wolfram. Credit-scarcity constraints below are
  relaxed: size entitlements generously; the worst-case-exceeds-balance concern is moot once topped up.
- Preferring single-kernel WL is now about *not wasting 5×*, not a budget limit.

**How billing works** (unchanged, still worth the hung-kernel guard)
- Creating an entitlement is free.
- Each kernel start triggers a charge, and charging continues **until the kernel terminates**, so a hung kernel keeps
  billing — but at $0.05/kernel-hour the dollar exposure is now trivial; keep `EntitlementExpiration` bounded anyway.

**⚠ `CheckCreditsBalance` is not a ceiling.**
- It only refuses to create an entitlement if the balance can't cover the kernel limits for one hour.
- It neither reserves credits nor caps spend. (Low stakes now, but keep it on.)

**The ceiling is the entitlement itself.** Worst case per run is roughly:

> (`StandardKernelLimit` + `ParallelKernelLimit`) × 10 credits/h × `EntitlementExpiration` hours,
> **plus** one billing interval per kernel start (restart loops cost extra).

Example: (1 + 4) × 10 × 10 h = **500 credits ≈ $2.50** — trivial against a topped-up 5000 balance.

**❓ Open**
- Do Personal Plan credits renew or roll over?
- What happens to a running kernel when the balance reaches zero? Ask Wolfram; don't test by draining the balance.

### 4.3 Per-run entitlement

The dev box creates the entitlement at pre-flight; this briefly uses one dev seat.

```bash
ENT_ID=$(wolframscript -code '
  CreateLicenseEntitlement[<|
    "StandardKernelLimit"   -> 1,
    "ParallelKernelLimit"   -> 0,
    "LicenseExpiration"     -> Quantity[8, "Hours"],
    "EntitlementExpiration" -> Quantity[10, "Hours"]
  |>]["EntitlementID"]' | tr -d '[:space:]')
```

**Expiration settings**
- ⚠ **`LicenseExpiration` caps each kernel's runtime.** Wolfram's template uses 1 h, which would kill a multi-hour
  audit.
  - Set it above the job's deadline.
  - Recompute the §4.2 worst case whenever you change it.
- Kernels can't outlive `EntitlementExpiration`.

**Kernel limits**
- `ParallelKernelLimit` defaults to 0.
- Set it to the job's declared worker count, and **only** if the script calls `Parallel*` / `LaunchKernels`.

**Recording and deletion**
- Record `ENT_ID` in the local run manifest; it is a secret, so keep the manifest local and `0600`.
- At teardown, first record `LicenseEntitlementObject["<ID>"]["CreditsSpent"]`, then run
  `DeleteObject[LicenseEntitlementObject["<ID>"]]`.

### 4.4 Runtime rules

**Injection**
- Put the ID in the box's secrets file as **`WOLFRAMSCRIPT_ENTITLEMENTID`**, a documented WolframScript environment
  variable.
- Don't use `-entitlement <ID>`: command-line arguments are visible to other users via `ps`.
- `LicensingSettings` applies to `RemoteBatchSubmit` jobs, **not** to our SSH-launched process.

**Build time is unlicensed.**
- The image bake runs **no licensed kernel**.
- The licensed smoke test happens only on a booted box.
- Never snapshot an image that contains an entitlement or activation.

**Network**
- The box needs **outbound** internet, because kernels contact Wolfram's license server regularly.
- Inbound access stays SSH-only.
- ⚠ A license-server or network blip mid-run is a new failure mode. The status file must capture it.

**Seat cap**
- On-demand kernels count against the entitlement's limits, not the 2 Mathematica seats. This is our reading of the
  docs ❓ confirm.
- The 2-kernel rule still governs seat-licensed kernels on the dev box, including the entitlement-creation kernel.

**Version pin**
- Install the Engine at the **exact version and build** of the dev box's Mathematica.
- ❓ Record the dev box's `$VersionNumber` and `$ReleaseNumber`.

**Front end**
- The Engine has no notebook front end.
- ❓ Confirm whether the `Dynamic[` / `Manipulate[` hits in `S11b_interface_coupling_law_*.wl` and
  `S11c_a_interface_geometry_*.wl` are real code.
- Don't assume that "batch" means "front-end-free" for scripts that haven't been checked.

---

## 5. Cost envelope

| Item | CCX43 (64 GB) | CCX53 (128 GB) |
|---|---|---|
| ~3 h job + provisioning + retrieval → **4 billable hours** | ≈ **$2.09** | ≈ **$4.04** |
| Wolfram, 1 kernel, ~3 h (sample rate ❓) | ~12 credits (existing balance) | same |
| Wolfram worst case with the §4.3 entitlement | ≤ 40 credits + restart intervals | same |
| Primary IPv4, tax | small ❓ | small ❓ |

**Standing cost**
- Snapshot storage only: compressed snapshot size × Hetzner's per-GB-month snapshot rate.
- ❓ Measure the size after the first bake and check the live rate. The research estimate was well under €1/month.
- The reaper adds roughly nothing if it runs on an existing CI plan (§6.6).

**Earlier totals were too optimistic.** The original report's per-run totals ($1.50–$3.50) ignored hour rounding and
used an FX-converted rate.

---

## 6. Design

### 6.1 Provider project and naming

- Use a **dedicated Hetzner project** (e.g. `cas-compute`) with its own API token.
  - Hetzner tokens are project-scoped, so any token that leaves the dev box can touch only this compute.
- Run ID format: `<UTC yyyymmddThhmmZ>-<job>`, lowercase.
  - Server name: `cas-<run-id>`.
- Labels:

| Resource | Labels |
|---|---|
| Ephemeral server | `role=ephemeral`, `run=<run-id>`, `job=<job>`, `expires=<unix epoch>` |
| Preserved image | `role=image`, `toolchain=<version tag>`, plus **delete protection** (`hcloud image enable-protection <id> delete`) |
| Firewall | `role=infra`. Reused across runs; inbound TCP 22 from our IP only. |

**⚠ Selection rule.** Every cleanup path (`down`, `nuke`, reaper) selects **only** `role=ephemeral`. Nothing
automated ever deletes `role=image` resources.

### 6.2 Image bake (Packer, `hcloud` builder)

**Base and toolchain**
- Base image: Ubuntu 22.04.
- Toolchain: pinned per requirements §7. Turn it into a lockfile:
  1. Write a `requirements.in` with exact pins (`sympy==1.14.0`, `mpmath==1.3.0`, `numpy==1.26.4`,
     `scipy==1.15.3`, `datalad==0.19.5`).
  2. Run `pip-compile --generate-hashes` to produce `requirements.lock`.
  3. Commit both files.
- Python: create a venv at `/opt/cas/venv` from the system `python3.10`.
  - ❓ Assert that the version is exactly 3.10.12 during the bake; jammy updates currently ship 3.10.12.
- **git-annex 10.20221003:**
  - ❓ Ubuntu 22.04's own archive ships an **older** git-annex.
  - Match however the dev box got 10.20221003 (check `apt-cache policy git-annex`, NeuroDebian, standalone build,
    etc.).
  - Assert the exact version during the bake.
- **Wolfram Engine:**
  - Download the Linux installer once to the dev box; Packer uploads it.
  - Install non-interactively at the pinned version/build. ❓ Verify the installer's unattended flags from its own
    help output.
  - **Don't prune paclets** until the real audit scripts have passed on this image.

**Users and SSH**
- Create user `cas` (non-root).
- Harden sshd: `PasswordAuthentication no`, `PermitRootLogin no`.
- Add `/etc/tmpfiles.d/cas.conf` containing `d /run/cas 0700 cas cas -`. This gives the job a tmpfs secrets
  directory on every boot.

**Scripts**
- Install the on-box job wrapper at `/usr/local/bin/cas-job` (§6.4).

**Bake-time checks (unlicensed only)**
- Python and package versions match the lockfile, and SymPy is exactly `1.14.0`.
- git-annex and datalad versions match.
- The Wolfram binaries are present and the installed version matches the pin.
- ⚠ **These checks must not launch a kernel.**

**Snapshot**
- Use Packer `snapshot_labels` to set `role=image,toolchain=<tag>`.
- Enable delete protection after the build.
- Re-bake only when the toolchain changes.
- ❓ Record the compressed snapshot size.

### 6.3 Boot, SSH, secrets

**Create** (in `cas up`):

```bash
hcloud server create --name "cas-$RUN" --type ccx43 --location nbg1 \
  --image "$IMAGE_ID" --ssh-key orchestrator --firewall cas-ssh \
  --user-data-from-file "$RUNDIR/cloud-init.yaml" \
  --label role=ephemeral --label run="$RUN" --label job="$JOB" --label expires="$EXPIRES" \
  -o json
```

- Pass `--ssh-key` even though root login is disabled. Otherwise Hetzner generates and emails a root password.
- **cloud-init must install our key for `cas`.** Provider key injection targets root, which the image locks out.
  - Use `write_files` (with `defer: true`, owner `cas:cas`, mode `0600`) to write `/home/cas/.ssh/authorized_keys`.
  - User-data is stored by the provider and readable from the box's metadata service, so it may carry **only public
    keys, never secrets**.
  - The readiness check below proves the key landed.

**Readiness**
- Get the IPv4 from the JSON output.
- Loop `ssh -o BatchMode=yes -o ConnectTimeout=5 -o UserKnownHostsFile="$RUNDIR/known_hosts"
  -o StrictHostKeyChecking=accept-new cas@$IP true` until it succeeds.
- Hard deadline: 5 min. On timeout, run teardown.
- An open port 22 alone doesn't count as ready.
- Use a per-run `known_hosts` because every box has new host keys and Hetzner reuses IPs.

**Secrets**
- Stream them over SSH, e.g. `ssh … 'umask 077; cat > /run/cas/secrets.env' < <(generate_secrets)`.
- Contents:
  - `WOLFRAMSCRIPT_ENTITLEMENTID`
  - GIN read credential (`DATALAD_CREDENTIAL_GIN_TOKEN` or a key, per §6.5)
  - any git read credential
- `cas-job` sources the file with `set -a`.
- Exporting variables in an SSH session doesn't carry over to later sessions or to the detached job.

### 6.4 Sync, run, validate, retrieve

**Sync** (as `cas`)
- `git clone` the repo, then `git checkout <pinned SHA>`.
- **Configure and verify GIN access:**
  - A fresh clone doesn't carry the dev box's remotes/annex config.
  - Replicate the dev box's `datalad siblings` setup for `gin` (❓ exact command depends on how the dataset records
    it).
  - Confirm with a `datalad get` of one input.
- `datalad get` all the job's inputs.
- Write `manifest.json` containing:
  - run ID and commit SHA;
  - input annex keys (`git annex find --format='${key}\n' <paths>`);
  - Wolfram `$VersionNumber`/`$ReleaseNumber`, run via the licensed kernel at job start inside the job;
  - SymPy version, worker count, instance type.

**On-box job wrapper** (`/usr/local/bin/cas-job <rundir> <cmd…>`):

```bash
#!/usr/bin/env bash
set -uo pipefail
rundir=$1; shift
set -a; . /run/cas/secrets.env; set +a
cd "$rundir" || exit 90
echo $$ > pid
( while :; do
    awk -v t="$(date +%s)" '/^MemTotal/{m=$2} /^MemAvailable/{a=$2} END{print t, m-a}' /proc/meminfo
    sleep 10
  done ) > mem.samples &
sampler=$!
"$@" > run.log 2>&1
rc=$?
kill "$sampler" 2>/dev/null
echo "$rc" > status.tmp && mv status.tmp status
```

**Launch** (detach **all** stdio, or ssh can hang):

```bash
ssh … cas@$IP "setsid /usr/local/bin/cas-job $RUNDIR wolframscript -file $SCRIPT </dev/null >/dev/null 2>&1 &"
```

A `sudo systemd-run --uid=cas --unit=cas-$RUN …` transient unit is an acceptable alternative. It also reports
`Result=oom-kill`.

**Poll** (controller), until the job deadline:

| Condition | Outcome |
|---|---|
| `status` exists | Finished — read the exit code. |
| `status` missing, `pid` not alive | **Failed** (wrapper killed, e.g. OOM). |
| Past the deadline | **Failed** — kill the job. |

**Validate.** Success requires **all** of:
1. `status` is `0`;
2. the audit's own pass criterion is met;
3. no OOM entries appear in `journalctl -k` for the run window.

**⚠ Audit scripts must signal failure.** They need to call `Exit[1]` on any failed check and print a final
machine-readable result line (e.g. `AUDIT-RESULT: PASS 12/12`).
- ❓ If the existing `.wl` scripts don't already do this, adding it **changes the audit instrument**. Route that
  through the repo's normal review legs; it is not an infra change.
- The failure demo (§8) is what proves that `wolframscript` propagates the exit code.

**Retrieve — always, before any teardown:**
- `rsync` the whole run directory back to the dev box: `run.log`, `status`, `pid`, `mem.samples`, `manifest.json`,
  and outputs.
- Also fetch the kernel log (`journalctl -k` for the run window).
- Record the peak from `mem.samples` in the manifest.

### 6.5 Publish (recommended: from the dev box)

**Why the dev box**
- The box then needs only **read** access to GIN and git: no GIN write key, no push credential.
- Results are already local before publishing starts, which makes publish failures recoverable by construction.

**Steps**
1. Place the retrieved outputs at their dataset paths in the dev box clone. The clone must be at the pinned commit or
   a descendant.
2. `datalad save -m "<run-id>: <job>" <paths>`. `push` publishes only saved state.
3. `datalad push --to gin` (content + pointers).
4. `git push origin`.
5. **Verify:** `git annex whereis <outputs>` lists gin, or `datalad get` succeeds from a throwaway fresh clone.

On any failure, **keep the local copy and the run directory**, report it, and continue to teardown. The box holds
nothing that isn't already on the dev box.

**If publishing from the box instead** (requirements §3 allows it):
- The box needs GIN write access. GIN annex content needs SSH, which means injecting a GIN-registered private key.
- Retrieval (§6.4) must still happen first.

### 6.6 Teardown and billing guard

**Teardown** (`cas down <run>`, also called from the trap):
1. Retrieve (§6.4), best effort.
2. `hcloud server delete cas-$RUN`. **Delete — never `poweroff`/`shutdown`.**
3. Primary IPs created with the server are normally auto-deleted. Verify anyway:
   - `hcloud primary-ip list -o json`;
   - flag any **unassigned** IP, since unassigned IPs keep billing.
4. Record `CreditsSpent`, then delete the entitlement (§4.3).
5. **Confirm:**
   - `hcloud server list -l role=ephemeral` is empty for this run;
   - no unassigned primary IPs remain;
   - no volumes are labeled `role=ephemeral`.

   An empty server list alone does not prove billing stopped.

**Layers of the guard**

| Layer | Covers | Status |
|---|---|---|
| Pre-flight: refuse if a `job=<job>` ephemeral server exists; list/delete are idempotent | duplicates | required |
| Controller `trap … EXIT` → `cas down` | normal errors | required, **best effort only** (no protection against power loss or `kill -9`) |
| **Independent reaper** on a schedule, running on **neither** machine | controller loss, hung VM, forgotten runs | **required** |
| On-box self-delete at `expires` via the API (must delete, not shut down) | controller loss | optional; puts a project token on the box, which is acceptable only with the dedicated project (§6.1) |
| Per-run `EntitlementExpiration` | hung kernels burning credits | required |
| Provider spending alert, if available | anything else | nice to have |

**Reaper**
- Every N minutes (e.g. 15):
  1. Run `hcloud server list -l role=ephemeral -o json`.
  2. Delete every server whose `expires` label is earlier than now; compare in `jq`, since label selectors can't
     compare numbers.
  3. Report unassigned primary IPs.
- The reaper must **never** match `role=image`.
- ❓ **Host undecided.** The default candidate is a scheduled CI workflow holding the project token.
  - Scheduled CI runs can be delayed.
  - Some platforms disable schedules on inactive repos (GitHub does this for public repos after 60 days without
    activity).
- **Health gate:** `cas up` must refuse to start if the reaper's last successful run is older than 2× its interval
  (e.g. check via `gh run list`). This makes a dead reaper visible instead of silent.
- Alternative: a provider with enforced max lifetime (GCP).

**`cas nuke`**
- Deletes all `role=ephemeral` servers.
- Reports unassigned primary IPs.
- Deletes the entitlements recorded in **local run manifests**, never "all entitlements on the account".
- Never touches images.

### 6.7 `cas` command surface

| Command | Does |
|---|---|
| `cas up <job> --size 64\|128 --workers N --deadline H` | reaper health gate → pre-flight → entitlement → create → readiness → secrets → sync |
| `cas run <run>` | launch detached (§6.4) |
| `cas status <run>` | poll: running / done(rc) / failed(reason) |
| `cas pull <run>` | retrieve run directory + kernel log |
| `cas publish <run>` | §6.5 from the dev box |
| `cas down <run>` | §6.6 teardown + confirmation |
| `cas list` | ephemeral servers, unassigned IPs, local manifests with live entitlements |
| `cas nuke` | §6.6 |

**CLI conventions**
- Always JSON (`-o json` + `jq`).
- Non-zero exit on any failure.
- Every command is idempotent and safe to re-run.
- `cas up` installs `trap 'cas down "$RUN"' EXIT` until `down` succeeds.

---

## 7. Build order (for the implementing agent)

### Phase 0 — verify (no compute spend)

**✅ RESOLVED by dev-box lookups (orchestrator, 2026-09-15)** — three Phase-0 items closed with local
fact-lookups (no compute spend, no account action):
- **Dev-box Wolfram pin:** `$Version = "15.0.1 for Linux x86 (64-bit) (July 2, 2026)"`, `$VersionNumber = 15.`,
  `$ReleaseNumber = 1`, `$SystemID = Linux-x86-64`. ⇒ bake the Engine at **exactly Wolfram 15.0.1, Linux-x86-64**.
- **git-annex source:** installed `10.20221003-2~nd22.04+1` from **NeuroDebian** (`http://neuro.debian.net/debian
  jammy/main`), not the Ubuntu archive (which ships 8.x). ⇒ the bake must add the NeuroDebian apt repo, then pin.
- **Front-end grep → BOTH FALSE POSITIVES.** `S11b_interface_coupling_law` line 485–486 matched the *variable*
  `omegaDynamic[w,x,t]`; `S11c_a_interface_geometry` line 1737 matched the *variable* `rawDynamic[...]` and the
  string `"DYNAMIC"` — neither is the `Dynamic[` front-end primitive. ⇒ **all** our `.wl` audits are kernel-clean
  and run headless under the Engine. Front-end concern **closed**.

**Credit rate — ✅ VERIFIED LIVE (orchestrator, 2026-09-15; owner-approved account check):** created a real
entitlement on the dev box, read its cost properties, deleted it — **0 credits spent, no dangling entitlement**
(`LicenseEntitlements[]` → `{}` after). Account policy = **Standard (`WSDS`)**:
- **Rate = 10 credits / kernel-hour** — for **both** `Standard` *and* `Parallel` kernels
  (`KernelCosts -> <|Standard -> 10. Credits/Hours, Parallel -> 10. Credits/Hours|>`). ⛔ **The 4-credit CI-template
  sample does NOT apply to our account** — the research report's 10 was right. **§4.2 / §5 / §8.2 must use 10.**
- **BillingInterval = 900 s** (confirmed). `$CloudConnected = True`.
- ⇒ Recompute: 200 credits ≈ **~20 single-kernel hours**; acceptance run (~3 h, **single kernel**) ≈ **~30
  credits**; **1 main + 4 subkernels = 50 credits/h** (a 3 h such run ≈ 150 credits = ¾ of balance).
- ⭐ **Mitigation:** SymPy/Python parallelism costs **zero** Wolfram credits (just CPU); only `LaunchKernels`/
  `Parallel*` *inside a `.wl`* costs. Our WL audits are single-kernel → keep `ParallelKernelLimit -> 0` and
  parallelize SymPy freely. Near-term WL comparison stays ~30 credits, not 150.

**Things the owner must supply or confirm**
- Credit renewal/rollover terms.
- Zero-balance behavior (ask Wolfram).
- ~~Dev box `$VersionNumber`/`$ReleaseNumber`.~~ ✅ resolved above (15.0.1 / 15. / 1).
- ~~How git-annex 10.20221003 was installed.~~ ✅ resolved above (NeuroDebian).

**Script checks**
- ~~Front-end grep on the two flagged scripts.~~ ✅ resolved above (both false positives; all `.wl` headless-clean).
- Whether the audit scripts already `Exit[1]` and print a result line. ⚠ **If not, adding it CHANGES the audit
  instrument → it is a physics-bearing script change requiring the repo's normal two-leg review (G1), NOT an infra
  edit.** See the success/verdict note below.

**⛔ Physics-discipline note the plan must honor (orchestrator):** the WL audit is a **blind cross-engine
comparison** — its whole purpose is that a **disagreement with SymPy is the FINDING**, not a job failure.
So "job success" (infra layer) = *ran cleanly, bounded output, no OOM, retrieved* — it must **NOT** be gated on
"the engines agree." A WL run that disagrees is a **successful run that produced a result to preserve and
adjudicate**, ⛔ never a failure to retry or a reason to withhold the output. The `Exit[1]`/result-line convention
(§6.4) is for *internal self-consistency* checks the audit already owns; keep it distinct from the cross-engine
verdict, which the orchestrator adjudicates (G4) — ⛔ never auto-decided by the runner.

**Decisions**
- Reaper host.
- Publish location (dev box recommended).

### Phase 1 — provider

- Create the dedicated Hetzner project and token.
- Request a server-limit raise to cover CCX53.
- Create the SSH key and the `cas-ssh` firewall.
- Confirm the live CCX43/CCX53 tariffs in the account's currency.

### Phase 2 — image

- Packer template, lockfile, bake, and unlicensed checks (§6.2).
- Label the image and enable delete protection.
- Record the snapshot size.

### Phase 3 — tooling

- Build `cas` (§6.7) and `cas-job` (§6.4).
- Build the reaper (§6.6) and its health gate.

### Phase 4 — demos (§8, prerequisites)

### Phase 5 — acceptance run (§8, end to end)

---

## 8. Acceptance

**Prerequisites (small and cheap; do these first)**

1. **Licensed smoke test** on a booted box.
   - Run `wolframscript -code '1+1'` with the env var set.
   - `CreditsSpent` should be about one interval.
   - If any job will use workers, confirm that `LaunchKernels[n]` gets `n` subkernels licensed under the entitlement.
2. **Failure demo.**
   - Run a script that calls `Exit[1]`.
   - The status file shows non-zero, `cas status` reports **failed**, and the logs come back.
3. **Controller-loss demo.**
   - Start a short run with a near `expires`, then `kill -9` the controller.
   - The **reaper** deletes the box, and no ephemeral resources remain.
4. **Publish-failure demo.**
   - Break the GIN push (e.g. a bad credential).
   - Results are intact on the dev box, teardown still happens, and a later `cas publish` succeeds.

**End to end: one scripted run, zero interactive steps**

> 1. Reaper health gate → pre-flight → per-run entitlement
> 2. Boot from the image: **CCX53** if peak memory at the stated worker count is unmeasured
> 3. Authenticated readiness → secrets file
> 4. Clone at the pinned SHA → GIN access verified → `datalad get` the S11c-b/c1/c2 exports → manifest
> 5. Run `S11c_b_brane_operator_mathematica_audit.wl` (+ c1/c2) at the stated worker count
> 6. Status 0 **and** the audit result line passes **and** no OOM
> 7. Retrieve the bounded `.out`, manifest, and memory samples
> 8. `datalad save` → `datalad push --to gin` → `git push origin` → retrieval verified
> 9. Delete the box → record `CreditsSpent` → delete the entitlement
> 10. Confirm:
>     - no ephemeral servers, unassigned IPs, or volumes remain;
>     - **zero reinstalls** were needed;
>     - credits spent are consistent with §4.2;
>     - the peak-memory record is stored for future sizing.

---

## 9. Risks and unknowns

**Money**
- ⚠ **Orphaned boxes.** Hetzner bills stopped servers. Mitigated only if the reaper actually runs, which is why the
  health gate exists.
- ⚠ **Credit exhaustion.** The sample rate is unverified; `CheckCreditsBalance` isn't a cap; a crash loop costs one
  interval per start; zero-balance behavior is unknown. Mitigate by sizing the per-run entitlement ceiling (§4.2).

**Correctness**
- ⚠ **False pass.** Mitigated by exit code + result line + OOM check. This depends on the audit scripts signaling
  failure (possible review-gated script change).
- ⚠ **Lost results.** Mitigated by retrieving before teardown and publishing from the dev box.
- **Memory.** Parallel workers multiply memory. Measure at concurrency, and do the first run on 128 GB.
- **Version drift.** Pin the Wolfram build and git-annex exactly; record both in every manifest.

**Environment**
- **License-server dependency.** A network blip mid-run can kill the kernel; the status file and logs must show it.
- **Hetzner pricing volatility.** There were three 2026 increases, and rescales reprice. Re-check tariffs before
  budgeting.
- **cloud-init behavior on a custom snapshot** (key for an existing user) is unproven. The authenticated readiness
  check catches a failure immediately.
- **Front-end dependence** in the two flagged scripts is unconfirmed.

**Fallbacks**
- Linode, Fly.io, and RunPod are not evaluated in depth.
- AWS/GCP/DO/Vultr figures come from the research pass and were not re-verified.

---

## 10. Changes from the original research report

1. **Watchdog.** Removed the `shutdown` self-destruct: it leaves Hetzner billing. Added the required independent
   reaper with a health gate. Stated plainly that the shell trap is best effort.
2. **`nuke`.** It no longer deletes snapshots. Images are `role=image` and delete-protected; cleanup targets only
   `role=ephemeral`.
3. **Run launch.** `DONE` is replaced by an exit-code status file. All stdio is detached. Added pid liveness and
   deadline handling, plus audit-level validation and an OOM check.
4. **Publishing.** Added `datalad save`, GIN sibling setup on fresh clones, retrieval verification, and
   retrieve-before-teardown. Publishing now defaults to the dev box.
5. **Licensing lifecycle in the executable steps.**
   - Per-run entitlement with explicit expirations and kernel limits.
   - Env-var injection.
   - `CreditsSpent` recorded.
   - Entitlement deleted at teardown.
   - Licensed smoke test only after boot.
   - Removed the "buy credits" step (the existing 200-credit balance is used).
6. **Spend ceiling.** `CheckCreditsBalance` is correctly described as a one-hour affordability check. The worst-case
   formula has been added.
7. **SSH bootstrap.** cloud-init installs the `cas` key. Readiness is an authenticated SSH connection with a deadline.
   Uses a per-run `known_hosts`. Secrets go in a tmpfs file, not session exports.
8. **Reproducibility.**
   - Exact Wolfram build pin; git-annex source flagged.
   - Commit SHA and input keys recorded in a manifest.
   - Per-job worker counts; peak memory measured at concurrency.
   - Paclet pruning deferred.
   - Headless compatibility is checked, not assumed.
9. **Pricing.**
   - Hetzner USD tariffs used instead of an FX conversion.
   - Hour rounding applied, so ~3 h jobs cost 4 billable hours ($2.09 / $4.04).
   - The narrow per-run totals were withdrawn.
10. **Linode.** The wrong "tops out below 64 GB" exclusion was replaced; the review reports High Memory plans up to
    300 GB.
11. **Credit rate.** The 10 credits/hour figure is marked unverified alongside the 4-credit sample; the real rate
    comes from our own entitlement.
12. **Licensing language.** The categorical "organizational use ⇒ prohibited" claim is softened to "arguable"; the
    decision rests on the activation limit.

---

## 11. Sources

**Wolfram**
- Free Wolfram Engine terms: https://www.wolfram.com/legal/terms/wolfram-engine.html
- Wolfram Engine FAQ: https://www.wolfram.com/engine/faq/
- `CreateLicenseEntitlement`: https://reference.wolfram.com/language/ref/CreateLicenseEntitlement.html
- `LicenseEntitlementObject`: https://reference.wolfram.com/language/ref/LicenseEntitlementObject.html
- `LicensingSettings`: https://reference.wolfram.com/language/ref/LicensingSettings.html
- WolframScript (incl. `WOLFRAMSCRIPT_ENTITLEMENTID`): https://reference.wolfram.com/language/ref/program/wolframscript.html
- Wolfram Engine Docker image (on-demand usage, charging until termination): https://hub.docker.com/r/wolframresearch/wolframengine
- Sample entitlement output (4 credits/h, 900 s): https://github.com/WolframResearch/WL-FunctionCompile-CI-Template and https://github.com/arnoudbuzing/function-compile-test
- Service Credits: https://www.wolfram.com/service-credits/

**Hetzner and pricing**
- Hetzner pricing tracker (EUR, hourly rounding): https://costgoat.com/pricing/hetzner
- Spare Cores CCX43 (USD from $0.5138/h): https://sparecores.com/server/hcloud/ccx43
- 2026 Hetzner increases: https://northflank.com/blog/hetzner-cloud-server-price-increases
- Packer `hcloud` builder: https://developer.hashicorp.com/packer/integrations/hetznercloud/hcloud/latest/components/builder/hcloud

**DataLad / GIN**
- DataLad GIN siblings: https://docs.datalad.org/en/maint/generated/man/datalad-create-sibling-gin.html

**Not re-checked here**
- The review's Hetzner USD tariffs and Linode 300 GB figure come from the reviewer's check of official pages.

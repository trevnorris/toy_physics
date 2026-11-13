"""Retrieve the accepted S15 results into the step record's durable evidence index."""
from collections import Counter, defaultdict
from hashlib import sha256
import gzip
import json
from pathlib import Path
import re
import subprocess
import sys
import tarfile

BASE = Path(__file__).resolve().parents[1]
ROOT = BASE.parents[1]
MEAS = BASE / "directives/_measurements"
ARCHIVE = BASE / "scripts/out/S15_reading_comparator_details.out"
FIELDS = ("Q", "E_geo", "E_resp", "separation_derivative", "departure_partition/residual",
          "item3/family_integral", "item3/object", "item3/Q_closed_operand",
          "item3/exponent_branch_inventory", "velocity_part_1/residual", "velocity_part_2/residual")
OUTPUT = BASE / "steps/_measurements/S15_record_measurements.md"


def data(name, revision="41adb784"):
    path = (MEAS / name).relative_to(ROOT)
    return json.loads(subprocess.check_output(["git", "show", f"{revision}:{path}"], cwd=ROOT))


def archived(archive, name):
    with gzip.open(archive.extractfile(name)) as stream:
        return json.load(stream)


def item4_conditions(manifest):
    """Stream the accepted Wolfram output, keeping final conditions from completed passes."""
    stream_hash = sha256()
    groups = {}
    current = defaultdict(list)
    pending = None
    base_pass = True
    header = re.compile(rb"^reading/(.+?)/item4/(condition(?:_\d+)?|Q_under_condition(?:_\d+)?|status|unresolved_parameter_branches): ")
    with (ROOT / manifest["transcript"]).open("rb") as stream:
        for line in stream:
            stream_hash.update(line)
            if line.startswith(b"END_OF_BASE_PASS: "):
                base_pass = False
                current = defaultdict(list)
                pending = None
            elif line.startswith(b"END_OF_ITEM4_PASS: "):
                name = json.loads(line.partition(b": ")[2])
                if pending != name or not current["status"]:
                    raise ValueError(f"Incomplete final pass for {name}")
                groups[name] = dict(current)
                current = defaultdict(list)
                pending = None
            elif line.startswith(b"reading/"):
                match = header.match(line[:600])
                if match:
                    name, field = [s.decode() for s in match.groups()]
                    if name != pending:
                        current = defaultdict(list)
                        pending = name
                    value = line[match.end():].decode().rstrip("\n")
                    current[field].append(value)
                    if base_pass and field == "status" and value == manifest["item4_statuses"].get(name) and "Not applicable" not in value:
                        groups[name] = dict(current)
    if stream_hash.hexdigest() != manifest["transcript_sha256"]:
        raise ValueError("Full Wolfram transcript fingerprint mismatch")
    for name, status in manifest["item4_statuses"].items():
        if "Not applicable" in status:
            continue
        group = groups[name]
        if group["status"][-1] != status:
            raise ValueError(f"Item-4 status mismatch for {name}")
        for key, conditions in group.items():
            if re.fullmatch(r"condition(?:_\d+)?", key):
                if len(conditions) != len(group[key.replace("condition", "Q_under_condition", 1)]):
                    raise ValueError(f"Unpaired item-4 condition for {name}")
    return groups


def main():
    run = data("S15_reading_comparator_run_manifest.json")
    summary = BASE / "scripts/out/S15_reading_comparator.out"
    coverage = MEAS / "S15_reading_comparator_coverage.json"
    for key, path in (("transcript", summary), ("detail_archive", ARCHIVE), ("final_coverage", coverage)):
        if sha256(path.read_bytes()).hexdigest() != run[key]["sha256"]:
            raise ValueError(f"Wrong accepted input {path}")
    lines = summary.read_text().splitlines()
    outcomes = [json.loads(line) for line in lines]
    wl = data("S15_separation_departure_wolfram_run_manifest.json")
    conditions = item4_conditions(wl)
    sections = ["# S15 record measurements", "Reproduce from the repository root with "
                "`python3 research/pde_ledger_v3/scripts/S15_record_measurements.py`. "
                "The sections below are the generator’s output. The published summary, detail archive, "
                "branch inventory and full native item-4 transcript are checked against the accepted manifests.",
                "## M1. Accepted comparator output", "```text"]
    kinds = {"reading_to_original", "final_field", "final_family_result", "signed_crossing",
             "all_dimension_records", "observed_argument_checks", "item3_exponent_coverage",
             "additional_branch_samples", "precision_rechecks", "unpaired_internal_kernel_occurrences"}
    sections.extend(line for line, record in zip(lines, outcomes) if record["record"] in kinds)
    sections += ["```", "## M2. Native examples and sufficient existence conditions", "```text"]
    with tarfile.open(ARCHIVE) as archive:
        members = set(archive.getnames())
        config = json.load(archive.extractfile("configuration.json"))
        spans = {entry["id"]: entry for entry in archived(archive, "coverage.json.gz")}
        roster = config["required_objects"]
        inventory = []
        counts = defaultdict(Counter)
        for i, name in enumerate(roster):
            member = f"objects/{i:03d}.json.gz"
            obj = archived(archive, member)
            if (obj["object"], obj["index"], obj["build_id"]) != (name, i, config["build_id"]):
                raise ValueError(member)
            allrows = obj["rows"]
            rows = {r["field"]:r for r in allrows if r["role"] == "cross"}
            final = {key:rows[key]["status"] for key in FIELDS if key in rows}
            refinement = f"refinement/{i:03d}.json.gz"
            inventory.append(f"### {i:03d}. {name}")
            detail = []
            if refinement in members:
                ref = archived(archive, refinement)
                if ref["branch"]:
                    if not ref["branch"]["same_elementary_boundaries"]:
                        raise ValueError(f"Unsupported branch disposition in {member}")
                    final["item3/exponent_branch_inventory"] = "same_elementary_loci"
                logname = f"log_branches/{i:03d}.json.gz"
                logs = {}
                if logname in members:
                    logs = {(r["field"],tuple(r["case"])): r for r in archived(archive,logname)["rechecks"]}
                for field in ref["fields"]:
                    key = field["field"]
                    status_counts = Counter()
                    for point in field["additional_samples"]:
                        effective = logs.pop((key,tuple(point["case"])), point)
                        status_counts[effective["status"]] += 1
                    status_counts.update(r["status"] for r in field["rechecks"])
                    if set(status_counts) - {"within_estimate","within_precision_estimate"}:
                        raise ValueError(f"Unexpected unresolved refinement at {member} {key}")
                    if final[key] not in ("exact_identity", "sample_coverage_incomplete") and status_counts:
                        final[key] = "sample_agreement"
                    detail.append(f"{key} {dict(status_counts)}")
                if logs:
                    raise ValueError(logname)
                inventory.append(f"Refinement `{refinement}`; additional log source "
                                 f"`{logname}`." if logname in members
                                 else f"Refinement `{refinement}`.")
                inventory.extend(detail)
            wl_name = spans[rows["Q"]["right"][0]]["tag"][len("reading/"):-len("/Q")]
            inventory.append(f"Archive `{member}`. Single-engine item-4 status for `{wl_name}` "
                             f"{wl['item4_statuses'][wl_name]}")
            for row in allrows:
                if row["role"] == "cross" and row["field"] not in final:
                    continue
                key = row["field"]
                status = row["status"]
                end = final[key] if row["role"] == "cross" else status
                if row["role"] == "cross":
                    counts[key][end] += 1
                inventory.append(f"- `{row['role']} / {key}` {status} → {end}. "
                                 f"{', '.join(row['left'])} ↔ {', '.join(row['right'])}")
            if name.startswith("both_surfaces/sum_crossing"):
                sections.append(name+" E_resp "+str(rows["E_resp"]["sample_counts"]))
                if i == 0:
                    sections.extend(rows["departure_existence_condition"]["residual"])
                    sections.extend(rows["item3/status"]["residual"])
                    sections.append("E_resp first stored sample "+json.dumps(rows["E_resp"]["samples"][0]))
            native_fields = ()
            if name == "inner/transport_centre/component_1/labels_drain2":
                native_fields = ("item3/Q_closed_operand", "item3/object")
            elif name.startswith("both_surfaces/sum_crossing/") or (
                name.startswith("region/") and name.split("/")[-1] in
                ("momentum", "integrated_force", "momentum_rate")
            ):
                native_fields = ("velocity_part_1/residual", "velocity_part_2/residual")
            if native_fields:
                sections.append(name)
                for key in native_fields:
                    r = rows[key]
                    for locator in r["left"]+r["right"]:
                        entry = spans[locator]
                        engine = locator.split(":")[0]
                        with (ROOT/config["inputs"][engine]["path"]).open("rb") as stream:
                            stream.seek(entry["span"]["offset"])
                            value = stream.read(entry["span"]["length"]).removesuffix(b"\n")
                        if sha256(value).hexdigest() != entry["payload"]:
                            raise ValueError(locator)
                        sections.append(f"{key} {locator} {value.decode()}")
        assert len(roster) == 336
        for row in outcomes:
            if row["record"] == "final_field":
                assert counts[row["field"]] == row["counts"], row
    sections += ["```", "## M3. Per-object comparison index",
                 "Each locator gives the native stream identifier and byte offset. Fields show initial → final "
                 "disposition. The associated refinement lists its observed checks. Archive members retain the other "
                 "operand and dimension diagnostics. M4 retains Wolfram’s separate one-engine output. "
                 "The tracked [exception inventory](../../directives/_measurements/S15_reading_comparator_coverage.json) "
                 "gives the tested and untested simultaneous exponent equalities object by object."]
    sections.extend(inventory)
    sections += ["## M4. Wolfram’s printed conditions and restricted quantities",
                 "The keys and values below are retrieved verbatim from the final completed pass of each eligible object. "
                 "Separately numbered condition and restricted-quantity lists pair in occurrence order. "
                 "Their standing is the native single-engine status in M3. A zero restricted quantity carries zero amplitude."]
    for name, group in sorted(conditions.items()):
        sections += [f"### {name}", "```text"]
        for tag, values in group.items():
            sections.extend(f"{tag}: {value}" for value in values)
        sections.append("```")
    register_path = "research/pde_ledger_v3/SUBSTRATE_REQUIREMENTS.md"
    old_register = subprocess.check_output(["git", "show", "83237f4a:" + register_path], cwd=ROOT, text=True)
    def census(source):
        rows = re.split(r"^### (R-[^ ]+) ", source.split("## Population passes")[0], flags=re.M)
        return {key: re.search(r"\*\*status\*\* (OPEN|RETIRED|DELIVERED)", entry)[1]
                for key, entry in zip(rows[1::2], rows[2::2])}
    before = census(old_register)
    after = census((ROOT/register_path).read_text())
    sections += ["## M5. Requirements-register census", "```text",
                 f"Before {dict(Counter(before.values()))}", f"After {dict(Counter(after.values()))}",
                 f"Added {sorted(after.keys()-before.keys())}", f"Removed {sorted(before.keys()-after.keys())}",
                 f"Status changes {sorted(k for k in before.keys() & after.keys() if before[k] != after[k])}", "```"]
    return "\n\n".join(sections).replace("```text\n\n","```text\n")+"\n"


if __name__ == "__main__":
    result = main()
    if sys.argv[1:] == ["--check"]:
        if OUTPUT.read_text() != result:
            raise SystemExit("Measurements differ from retrieval")
    elif sys.argv[1:]:
        raise SystemExit("Use no arguments to generate or --check to verify")
    else:
        temp = OUTPUT.with_suffix(".md.tmp")
        temp.write_text(result)
        temp.replace(OUTPUT)
    print(f"{OUTPUT.name} {len(result.encode())} bytes")

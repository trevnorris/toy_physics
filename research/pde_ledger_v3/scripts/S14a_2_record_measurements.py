#!/usr/bin/env python3
"""Retrieve the filed S14a part 2 record evidence, without running a witness.

Only existence, verbatim text/git retrieval and literal-match counts are used.
No CAS, numerical work, engine import, witness execution or physics predicate.
Any working output belongs under _scratch/s14a2/rec/work/.
"""
from pathlib import Path
import argparse
import io
import shlex
import subprocess


ROOT = Path(__file__).resolve().parents[3]
V3 = "research/pde_ledger_v3/"
D = V3 + "directives/"
DM = D + "_measurements/"
WITNESSES = DM + "S14a_2_review_witnesses/"
REGISTER = V3 + "SUBSTRATE_REQUIREMENTS.md"
PART1 = V3 + "steps/S14a_1_projection.md"
SCRIPT = V3 + "scripts/S14a_2_record_measurements.py"
BASE = "978c359c"
OUTPUT = ROOT / V3 / "steps/_measurements/S14a_2_record_measurements.md"
WORK = ROOT / "_scratch/s14a2/rec/work"
TITLES = {
    "M01": "Record contract, its gate and commit provenance",
    "M02": "Part 2 decision list, revisions 1 and 2",
    "M03": "Both review dispositions, lookup existence and filed assembly audit",
    "M04": "Pass 1 Codex: moving sources and constant totals",
    "M05": "Pass 1 Codex: constant sources do not select a limit",
    "M06": "Pass 1 Codex: regular versus full face-source totals",
    "M07": "Pass 1 Grok: throughput families and order-budget assembly",
    "M08": "Pass 2 Codex: long-time average versus a slow family",
    "M09": "Pass 2 Codex: spreading and the order of limits",
    "M10": "Pass 2 Grok: projected balance and two volume-flux gates",
    "M11": "Accepted part 1, carried citation sources and S12 owners",
    "M12": "Register schema, both-route method, census and selection entry",
    "M13": "Plan, v2 and G0 comparison inventories, retrieved only",
    "M14": "Second-route candidates and existing objection handoffs",
    "M15": "Final pass/status and exact ancillary diff against the starting commit",
    "M16": "Record round-0 dispositions, V2 user decision and committed-generator provenance",
    "M17": "Record round-1 dispositions and preserved review provenance",
}
WITNESS_NAMES = {
    "M04": "p1-codex_check_steady_source_motion",
    "M05": "p1-codex_check_limit_not_selected",
    "M06": "p1-codex_check_face_total_scope",
    "M07": "p1-grok_steady_throughput_and_budget",
    "M08": "p2-codex_check_dc_vs_slow",
    "M09": "p2-codex_check_spreading_order",
    "M10": "p2-grok_dc_limit_objects",
}


def command(out, *argv):
    """Emit the retrieval command and its literal stdout; no shell or timeout."""
    print("$ " + shlex.join(argv), file=out)
    result = subprocess.run(argv, cwd=ROOT, check=True, stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE, text=True)
    if result.stderr:
        raise RuntimeError(result.stderr)
    out.write(result.stdout)
    if result.stdout and not result.stdout.endswith("\n"):
        out.write("\n")
    out.write("\n")
    return result.stdout


def cat(out, path, revision=None):
    if revision:
        return command(out, "git", "show", revision + ":" + path)
    return command(out, "cat", path)


def excerpt(out, path, ranges):
    command(out, "sed", "-n", ranges, path)


def stored_block(out, path, heading, revision=None):
    """Verbatim retrieval of a named Markdown block, without interpreting it."""
    value = (subprocess.check_output(["git", "show", revision + ":" + path],
                                    cwd=ROOT, text=True) if revision else
             (ROOT / path).read_text())
    lines = value.splitlines(keepends=True)
    first = next(i for i, line in enumerate(lines) if line.startswith(heading))
    last = next((i for i in range(first + 1, len(lines))
                 if lines[i].startswith("### ")), len(lines))
    print("STORED BLOCK " + path + (" AT " + revision if revision else "")
          + " / " + heading, file=out)
    out.write("".join(lines[first:last]))
    out.write("\n")


def retrieve(section):
    out = io.StringIO()
    if section == "M01":
        cat(out, D + "S14a_2_record_decision_list.md")
        for name in ["S14a_2_record_setup_lookups.md",
                     "S14a_2_record_decision_list_review_p1_disposition.md",
                     "S14a_2_record_decision_list_review_p1_lookups.md"]:
            cat(out, DM + name)
        command(out, "git", "show", "-s", "--format=%H %s", BASE,
                "d72bb551", "945a1cd2", "e5b20813", "e9eed614",
                "94fddf2f", "de81b9fc", "4a4dfaf5", "ac909aaf")
        cat(out, DM + "S14a_2_review_witnesses.sha256")
        command(out, "git", "ls-tree", "-r", "--name-only", BASE, "--",
                DM + "S14a_2_review_witnesses", DM + "S14a_2_review_witnesses.sha256")
    elif section == "M02":
        for revision in ["d72bb551", "945a1cd2"]:
            cat(out, D + "S14a_2_decision_list.md", revision)
        cat(out, D + "S14a_2_decision_list.md")
        cat(out, DM + "S14a_2_setup_lookups.md")
    elif section == "M03":
        for pass_number in [1, 2]:
            cat(out, DM + f"S14a_2_decision_list_review_p{pass_number}_disposition.md")
            command(out, "git", "ls-tree", "--name-only", BASE, "--",
                    DM + f"S14a_2_decision_list_review_p{pass_number}_lookups.md")
        # Some historical cut-truncated lookup lines contain incomplete UTF-8.
        # Leave them untouched. Retrieve the complete primary witness text in
        # M04-M10 and the specific additional assembly stdout needed here.
        excerpt(out, DM + "S14a_2_decision_list_review_p2_lookups.md", "224,239p")
    elif section in WITNESS_NAMES:
        name = WITNESS_NAMES[section]
        for extension in [".py", ".stdout"]:
            cat(out, WITNESSES + name + extension)
    elif section == "M11":
        cat(out, D + "S14a_1_decision_list.md")
        for round_number in [0, 1]:
            cat(out, DM + f"S14a_1_record_r{round_number}_review_disposition.md")
        cat(out, PART1, BASE)
        cat(out, PART1)
        excerpt(out, D + "S12_decision_list.md", "90,99p")
        excerpt(out, D + "S11_requirements_pass2_directive.md", "30,65p")
        cat(out, V3 + "steps/S12_order_source.md")
    elif section == "M12":
        # The two counts are exact literal occurrences, not a status classifier.
        for revision in [BASE, None]:
            value = (subprocess.check_output(["git", "show", BASE + ":" + REGISTER],
                                            cwd=ROOT, text=True) if revision else
                     (ROOT / REGISTER).read_text())
            print("LITERAL COUNTS " + REGISTER + " AT " + (revision or "working tree"), file=out)
            for literal in ["\n### R-", "**status** OPEN"]:
                print(repr(literal) + " = " + str(value.count(literal)), file=out)
            print(file=out)
            if revision:
                # Full baseline text keeps schema, all entries and prior passes retrievable.
                cat(out, REGISTER, revision)
            else:
                excerpt(out, REGISTER, "1,115p")
            stored_block(out, REGISTER, "### R-S12-02", revision)
            stored_block(out, REGISTER, "### Method for population passes", revision)
        stored_block(out, REGISTER, "### R-S14a-01")
        stored_block(out, REGISTER, "### S14a part 2 pass")
    elif section == "M13":
        excerpt(out, V3 + "V3_STEP_PLAN.md", "566,675p")
        cat(out, V3 + "V3_STEP_PLAN.md", "bcbc2716")
        for path, ranges in [
            ("software/stage1_solver/reports/pathA_21_emergent_G_mass_bridge.md", "16,23p;62,66p"),
            ("research/pde_ledger_v2/notes/stages/ledger_stage004_gnls_action_dimensional_foundation.md", "148,150p"),
            ("software/em_charge_attribute/g0_closure_card_v0.md", "291,322p;387,398p"),
            ("research/pde_ledger_v2/notes/stages/ledger_stage002_matter_stress_force_assembly.md", "51,55p;109,113p"),
            ("research/pde_ledger_v2/notes/stages/ledger_stage006_two_phase_chiB_ontology.md", "78,87p"),
        ]:
            excerpt(out, path, ranges)
    elif section == "M14":
        excerpt(out, D + "O2_steady_brane_balance_scoping.md", "182,208p")
        for heading in ["### O2 pass", "### S9b pass", "### Light-closing pass",
                        "### S14a part 1 pass"]:
            stored_block(out, REGISTER, heading, BASE)
        excerpt(out, V3 + "steps/S9b_flow_oracle_comparison.md", "151,178p")
    elif section == "M15":
        stored_block(out, REGISTER, "### S14a part 2 pass")
        cat(out, "STATUS.md")
        command(out, "git", "diff", BASE, "--name-only")
        print("The pre-existing .claude/skills/prep-for-compact/SKILL.md change is unrelated; "
              "its contents are not retrieved.", file=out)
        command(out, "git", "diff", "--no-ext-diff", "--no-textconv", BASE, "--",
                REGISTER, PART1, "STATUS.md")
        command(out, "git", "status", "--short", "--", REGISTER, PART1, "STATUS.md",
                V3 + "steps/S14a_2_slow_limit.md",
                V3 + "steps/_measurements/S14a_2_record_measurements.md", SCRIPT)
    elif section == "M16":
        cat(out, DM + "S14a_2_record_r0_review_disposition.md", "12aa65d5")
        command(out, "git", "show", "-s", "--format=%H %s", "2ae682e4", "12aa65d5")
        command(out, "git", "ls-tree", "4a4dfaf5", "--",
                V3 + "scripts/S14a_1_record_measurements.py")
        command(out, "git", "ls-tree", "2ae682e4", "--", SCRIPT,
                V3 + "steps/_measurements/S14a_2_record_measurements.md",
                DM + "S14a_2_record_r0_review_disposition.md",
                DM + "S14a_2_record_r0_review_lookups.md")
    elif section == "M17":
        cat(out, DM + "S14a_2_record_r1_review_disposition.md", "85904c4e")
        command(out, "git", "show", "-s", "--format=%H %s", "85904c4e")
        command(out, "git", "ls-tree", "85904c4e", "--",
                DM + "S14a_2_record_r1_review_disposition.md",
                DM + "S14a_2_record_r1_review_lookups.md")
    else:
        raise ValueError(section)
    return out.getvalue()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--section", choices=TITLES)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    args = parser.parse_args()
    if args.section:
        print(retrieve(args.section), end="")
        return
    destination = args.output.resolve()
    if destination != OUTPUT and not destination.is_relative_to(WORK):
        parser.error("working output must be under _scratch/s14a2/rec/work/")
    result = ["# S14a part 2 record measurements\n\n",
              f"Generated by `{SCRIPT}` from starting commit `{BASE}`.\n\n",
              "Commands run from the repository root. Each Mnn block is the literal stdout of\n",
              f"`python3 {SCRIPT} --section Mnn`; inner commands and their literal stdout are included.\n",
              "Named-block retrievals and census counts are text retrieval and exact literal matches only.\n",
              "The filed witness scripts are read as text, never imported or executed. Their stdout is\n",
              "saved review evidence, not a new run or an accepted part 2 build. The checksum manifest\n",
              "and historical rerun reports are retrieved as filed; no new verification run is claimed.\n\n"]
    for section, title in TITLES.items():
        result += [f"## {section} — {title}\n\n",
                   f"````````text\n$ python3 {SCRIPT} --section {section}\n",
                   retrieve(section), "````````\n\n"]
    WORK.mkdir(parents=True, exist_ok=True)
    # Keep the generator's working file out of every other job's scratch directory.
    temporary = WORK / "S14a_2_record_measurements.generated.md"
    temporary.write_text("".join(result))
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_bytes(temporary.read_bytes())


if __name__ == "__main__":
    main()

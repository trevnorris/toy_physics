"""Mechanical artifact checks; no evaluation of Wolfram physics objects."""
from collections import Counter
import hashlib
import json
from pathlib import Path
import re
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[4]
DEST = Path(__file__).resolve().parent
STATE = ROOT / '_scratch/s9b_wl_repair8'
RUN = sys.argv[1]
FULL = ROOT / '_scratch/s11c' / RUN
WORK = STATE / 'ablation-final'
PREFIX = 'S9b_repair_wolfram_round8_'
KNIVES = ['K1', 'K2', 'K3', 'K4a', 'K4b', 'K4c', 'K5', 'K6', 'K7', 'K8', 'K9a', 'K9b', 'K10', 'K11']


def fields(text):
    """Split top-level WL association/list operands without parsing/evaluating."""
    start, stack, quoted, escaped, parts = 0, [], False, False, []
    for i, char in enumerate(text):
        if quoted:
            if escaped:
                escaped = False
            elif char == '\\':
                escaped = True
            elif char == '"':
                quoted = False
            continue
        if char == '"':
            quoted = True
        elif char in '[{(':
            stack.append(char)
        elif char in ']})':
            assert stack and stack.pop() == {']': '[', '}': '{', ')': '('}[char]
        elif char == '<' and text[i:i+2] == '<|':
            stack.append('<')
        elif char == '>' and i and text[i-1] == '|':
            assert stack.pop() == '<'
        elif char == ',' and not stack:
            parts.append(text[start:i].strip())
            start = i+1
    assert not stack and not quoted
    parts.append(text[start:].strip())
    return parts


def association(text):
    assert text.startswith('<|') and text.endswith('|>')
    return dict(field.split(' -> ', 1) for field in fields(text[2:-2]))


def capture(path):
    lines = path.read_text().splitlines()
    result = dict(line.split(': ', 1) for line in lines)
    assert len(result) == len(lines) and result
    return result


baseline = capture(WORK / 'BASELINE.stdout')
shared = {x.removeprefix('WL_S9B_') for x in baseline if x.startswith('WL_S9B_')}
grades = [f'{a}{b}{c}' for a in range(2) for b in range(3) for c in range(2)]
first = ['100', '010', '020', '001']
observations = ['DEFLECTION', 'RADAR']
responses = ['CONSTANT', 'FIXED_RATIO', 'POWER']
gates = {'BRANCH_EXISTENCE', 'PATH_TRAVERSAL', 'BRANCH_TYPE'}
all_a = {f'A_{name}_G{g}' for name in ['DEFLECTION', 'ROUND_TRIP', 'ONE_WAY_ER', 'ONE_WAY_RE', 'NONRECIPROCAL'] for g in grades}
radar = {f'A_RADAR_LOG_G{g}' for g in first}
b = {f'B_{name}_{obs}' for name in ['GAMMA', 'RESIDUAL', 'CONDITION', 'IMPLIED_JN'] for obs in observations} | {'B_GAMMA_DIFFERENCE'}
c = {f'C_{resp}_{name}_{obs}' for resp in responses for name in ['CONDITION', 'IMPLIED_JN'] for obs in observations} | {f'C_{resp}_N_DEPENDENCE' for resp in responses}
expected = gates | all_a | radar | b | c | {'A_NONRECIPROCAL_PATH_DEPENDENCE', 'F_MASS_SOLUTION'}
for name in ['FLOW_ONLY', 'SPEED_ONLY', 'TILT_ONLY']:
    expected |= {f'R_{name}_'+x for x in radar | b | {f'A_DEFLECTION_G{g}' for g in grades}}
expected |= {'R_BULK_ONLY_'+x for x in c if any(x.startswith('C_'+resp+'_') for resp in ['FIXED_RATIO', 'POWER']) and '_N_DEPENDENCE' not in x}
for stage in ['FLOW', 'LIVE']:
    expected |= {'F_'+stage+'_'+x for x in gates | all_a | radar | {x for x in b if '_IMPLIED_JN_' not in x}}
expected |= {'F_LIVE_'+x for x in c if '_CONDITION_' in x}
assert shared == expected, {'missing': sorted(expected-shared), 'extra': sorted(shared-expected)}
local = {x for x in baseline if x.startswith('WL_LOCAL_')}
registry = set(re.findall(r'"(WL_LOCAL_S9B_[^"]+)"', baseline['WL_LOCAL_S9B_ENGINE_LOCAL_NAMES']))
assert local == registry

captures = {name: capture(WORK / (name+'.stdout')) for name in KNIVES if name != 'K11'}
compact_b = capture(WORK / 'K11_COMPACT_BASELINE.stdout')
compact_a = capture(WORK / 'K11_COMPACT_CORRUPTED.stdout')
assert compact_b.keys() == compact_a.keys()
counts = {name: Counter() for name in KNIVES+['HARNESS_DEAD_PATH']}
witnesses, seen, harness_ne, sites = {}, set(), [], {}
witness_lines = []
k11_lines = []
with (FULL / 'stdout').open() as stream:
    for line in stream:
        tag, payload = line.rstrip('\n').split(': ', 1)
        if tag.endswith('_SITE'):
            sites[tag.removeprefix('WL_LOCAL_S9B_ABLATION_').removesuffix('_SITE')] = association(payload)
        if tag.endswith('_WITNESS_SEARCH'):
            witness_lines.append(line)
            name = tag.removeprefix('WL_LOCAL_S9B_ABLATION_').removesuffix('_WITNESS_SEARCH')
            values = association(payload)
            witnesses[name] = {'result': values['"Result"'], 'site': values['"Site"']}
            if values['"Result"'] == '"witness"':
                checks = association(values['"Checks"'])
                assert all(x == 'True' for x in checks.values()), (name, checks)
                witnesses[name]['evaluations'] = values['"Evaluations"']
                witnesses[name]['conditions'] = values['"ConditionEvaluations"']
                witnesses[name]['solutionIntersection'] = association(values['"FixingAndSolutions"'])['"SolutionIntersection"']
            continue
        match = re.fullmatch(r'WL_LOCAL_S9B_ABLATION_(K\d+[abc]?|HARNESS_DEAD_PATH)_(WL_.+)', tag)
        if not match:
            continue
        name, engine_tag = match.groups()
        assert (name, engine_tag) not in seen
        seen.add((name, engine_tag))
        values = association(payload)
        assert values['"Baseline"'] == baseline[engine_tag]
        if 'NOT_ESTABLISHED' in payload:
            harness_ne.append(tag)
        if name == 'K11':
            if engine_tag not in compact_b:
                assert values['"DifferencePair"'] == values['"MethodResidualPair"'] == 'Missing["NOT_EVALUATED"]'
                assert values['"DifferencePairMark"'] == values['"MethodResidualPairMark"'] == '"not compared"'
                counts[name]['not evaluated'] += 1
                continue
            k11_lines.append(line)
            assert fields(values['"DifferencePair"'][1:-1]) == [compact_a[engine_tag], compact_b[engine_tag]]
            assert fields(values['"MethodResidualPair"'][1:-1]) == [compact_b[engine_tag], baseline[engine_tag]]
            for key, left, right in [('DifferencePairComparison', compact_a[engine_tag], compact_b[engine_tag]), ('MethodResidualPairComparison', compact_b[engine_tag], baseline[engine_tag])]:
                result = association(values['"'+key+'"'])['"Mark"']
                assert result == ('"identical"' if left == right else '"not compared"'), (name, engine_tag, key)
                counts[name][key+': '+result.strip('"')] += 1
        else:
            corrupted = baseline if name == 'HARNESS_DEAD_PATH' else captures[name]
            assert values['"Corrupted"'] == corrupted[engine_tag]
            result = association(values['"StructuralComparison"'])['"Mark"']
            assert result == ('"identical"' if corrupted[engine_tag] == baseline[engine_tag] else '"not compared"'), (name, engine_tag)
            counts[name][result.strip('"')] += 1
assert set(witnesses) == set(KNIVES)
assert seen == {(name, tag) for name in counts for tag in baseline}

engine_path = ROOT / 'research/pde_ledger_v3/mathematica/S9b_light_bending_mathematica_audit.wl'
engine_bytes = engine_path.read_bytes()
assert (ROOT / 'research/pde_ledger_v3/mathematica/S9b_light_bending_mathematica_ablation.wl').read_bytes() == (STATE / 'harness-pinned-final.wl').read_bytes()
assert engine_bytes == subprocess.check_output(['git', 'show', 'HEAD:'+str(engine_path.relative_to(ROOT))], cwd=ROOT)
source = engine_bytes.decode()
mutations = {}
for name in KNIVES:
    mutated = (WORK / (name+'.wl')).read_text()
    # Export[Text] strips a terminal newline; no executable extra difference.
    from difflib import SequenceMatcher
    changes = [(source[a:b], mutated[c:d]) for op, a, b, c, d in SequenceMatcher(None, source.rstrip('\n'), mutated.rstrip('\n'), autojunk=False).get_opcodes() if op != 'equal']
    mutations[name] = changes
    assert changes, name
    original_site = json.loads(sites[name]['"Construction"'])
    replacement_site = json.loads(sites[name]['"Mutation"'])
    assert source.count(original_site) == 1
    assert mutated.rstrip('\n') == source.replace(original_site, replacement_site).rstrip('\n')
boundary = '(* COMPACT_K11_BOUNDARY *)'
for suffix, src in [('BASELINE', source), ('CORRUPTED', (WORK / 'K11.wl').read_text())]:
    assert (WORK / ('K11_COMPACT_'+suffix+'.wl')).read_text() == src.split(boundary)[0]

repo = ROOT / '_scratch/s11c/s9b-r8-engine-repository/stdout'
standalone = ROOT / '_scratch/s11c/s9b-r8-engine-standalone/stdout'
assert repo.read_bytes() == standalone.read_bytes() == (WORK / 'BASELINE.stdout').read_bytes()
assert sorted(x.name for x in (STATE / 'standalone').iterdir()) == [engine_path.name]
out_status = subprocess.check_output(['git', 'status', '--porcelain', 'research/pde_ledger_v3/mathematica/out/'], cwd=ROOT).decode()
assert out_status == ''
run_names = [RUN, 's9b-r8-engine-repository', 's9b-r8-engine-standalone', 's9b-r8-structural-checks5']
runs = {}
for name in run_names:
    folder = ROOT / '_scratch/s11c' / name
    outcome = json.loads((folder / 'outcome.json').read_text())
    assert outcome['exitCode'] == 0 and outcome['childOutcome']['guardReason'] is None
    assert (folder / 'stderr').stat().st_size == 0
    runs[name] = {'outcome': outcome, 'invocation': json.loads((folder / 'invocation.json').read_text()), 'resourceSamples': str(folder / 'resource-samples.jsonl')}
growth = json.loads((STATE / 's9b-r8-engine-standalone.run.json').read_text())['growthWhileKernelAlive']
assert growth
checks = capture(ROOT / '_scratch/s11c/s9b-r8-structural-checks5/stdout')
assert checks['HARNESS_PARSE_HEAD'] == 'HoldComplete'
assert len(checks) == 13
assert association(association(checks['SAME_SYNTAX'])['"Comparison"'])['"Mark"'] == '"identical"'
for key in checks.keys()-{'HARNESS_PARSE_HEAD', 'SAME_SYNTAX', 'COMMON_FIXED_SOLUTION_SETS', 'DISJOINT_FIXED_SOLUTION_SETS', 'ORDERED_CASE_SELECTION'}:
    assert association(association(checks[key])['"Comparison"'])['"Mark"'] == '"not compared"'

def write(suffix, text):
    (DEST / (PREFIX+suffix)).write_text(text)

for suffix, tags in [('tags.txt', baseline), ('not_established_tags.txt', [x for x in baseline if 'NOT_ESTABLISHED' in baseline[x]]), ('harness_not_established_tags.txt', harness_ne), ('k11_evaluated_tags.txt', compact_b), ('k11_not_evaluated_tags.txt', baseline.keys()-compact_b.keys())]:
    write(suffix, '\n'.join(sorted(tags))+'\n')
hashes = {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in [engine_path, ROOT / 'research/pde_ledger_v3/mathematica/S9b_light_bending_mathematica_ablation.wl', ROOT / 'research/pde_ledger_v3/mathematica/S9b_light_bending_mathematica_ablation_checks.wl', Path(__file__).resolve(), ROOT / 'scripts/s11c_guarded_run.py']}
write('witness_searches.txt', ''.join(witness_lines))
write('k11_reached_pairs.txt', ''.join(k11_lines))
write('structural_checks.txt', (ROOT / '_scratch/s11c/s9b-r8-structural-checks5/stdout').read_text())
write('source_hashes.json', json.dumps(hashes, indent=2)+'\n')
summary = {'tags': len(baseline), 'shared': len(shared), 'local': len(local), 'missingShared': [], 'baselineNotEstablished': sum('NOT_ESTABLISHED' in x for x in baseline.values()), 'structuralCounts': counts, 'witnessSearches': witnesses, 'K11Reached': len(compact_b), 'K11Unreached': len(baseline)-len(compact_b), 'engineUnchanged': True, 'engineBytes': len(engine_bytes), 'standaloneStreamsByteIdentical': True, 'standaloneDirectoryOnlyEngine': True, 'productionOutStatus': out_status, 'growthWhileKernelAlive': growth, 'mutations': mutations, 'runs': runs}
write('validation.json', json.dumps(summary, indent=2)+'\n')
print(json.dumps({key: summary[key] for key in ['tags', 'shared', 'local', 'baselineNotEstablished', 'K11Reached', 'K11Unreached', 'engineUnchanged', 'standaloneStreamsByteIdentical', 'standaloneDirectoryOnlyEngine', 'productionOutStatus']}, indent=2))
print(json.dumps({'structuralCounts': counts, 'witnessSearches': witnesses}, indent=2))

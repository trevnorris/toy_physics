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
STATE = ROOT / '_scratch/s9b_wl_repair10'
RUN = sys.argv[1] if len(sys.argv)>1 else 's9b-r10-full'
FULL = ROOT / '_scratch/s11c' / RUN
WORK = STATE / 'ablation-final'
PREFIX = 'S9b_repair_wolfram_round10_'
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
inventories = {}
search_plans = []
sign_trials = {}
site_lines = []
k11_lines = []
with (FULL / 'stdout').open() as stream:
    for line in stream:
        tag, payload = line.rstrip('\n').split(': ', 1)
        if tag.endswith('_SITE'):
            sites[tag.removeprefix('WL_LOCAL_S9B_ABLATION_').removesuffix('_SITE')] = association(payload)
        if tag.endswith('_SITE_INVENTORY'):
            name = tag.removeprefix('WL_LOCAL_S9B_ABLATION_').removesuffix('_SITE_INVENTORY')
            values = association(payload)
            parts = [association(x) for x in fields(values['"SiteParts"'][1:-1]) if x]
            exclusions = [association(x) for x in fields(values['"Exclusions"'][1:-1]) if x]
            ids = [x['"ID"'] for x in parts+exclusions]
            assert len(ids) == len(set(ids))
            assert parts and exclusions
            assert all(x['"Reason"'] for x in exclusions)
            inventories[name] = {'parts': len(parts), 'exclusions': len(exclusions),
                'tags': sorted(set(x['"Tag"'].strip('"') for x in parts)),
                'compactUnreachedParts': sum(x['"CompactReached"']=='False' for x in parts),
                'oneCopyParts': sum('not emitted' in x['"Copies"'] for x in parts),
                'literalProfilesParts':sum(x['"Role"']=='"literal profile list"' for x in parts)}
            assert not any('\"Profiles\"' in x['"Location"'] and 'additional exclusion 1' in x['"Reason"'] for x in exclusions)
            site_lines.append(line)
            continue
        if tag.endswith('_SEARCH_PLAN'):
            search_plans.append(line)
            continue
        if tag.endswith('_DENSITY_SIGN_TRIALS'):
            name = tag.removeprefix('WL_LOCAL_S9B_ABLATION_').removesuffix('_DENSITY_SIGN_TRIALS')
            trials = [association(x) for x in fields(payload[1:-1])]
            assert len(trials) == 2
            assert [x['"DensitySignTrial"'] for x in trials] == ['"rho <= 0"','"rho > 0"']
            assert all(x['"SignDecision"'] == 'True' for x in trials)
            assert all(x['"AsymptoticDensityEvaluation"'] == '1' for x in trials)
            sign_trials[name] = [{'sign': x['"DensitySignTrial"'], 'densities': x['"DensityEvaluationsEachEmittingCopy"'],
                'evaluations': x['"Evaluations"'], 'result': x['"Result"']} for x in trials]
            witness_lines.append(line)
            continue
        if tag.endswith('_ATTEMPT'):
            witness_lines.append(line)
            continue
        if tag.endswith('_WITNESS_SEARCH'):
            witness_lines.append(line)
            name = tag.removeprefix('WL_LOCAL_S9B_ABLATION_').removesuffix('_WITNESS_SEARCH')
            values = association(payload)
            witnesses[name] = {'result': values['"Result"'], 'site': values['"Site"'], 'attempt': values.get('"Attempt"'), 'evaluations': values.get('"Evaluations"')}
            if values['"Result"'] == '"witness"':
                checks = association(values['"Checks"'])
                assert all(x == 'True' for x in checks.values()), (name, checks)
                witnesses[name]['evaluations'] = values['"Evaluations"']
                if name == 'K7':
                    assert values['"LiteralEqualityDecision"'] == 'False'
                    witnesses[name]['literalStatusTexts'] = values['"Texts"']
                else:
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
assert set(witnesses) == set(inventories) == set(KNIVES)
assert set(sign_trials) == {"K9a", "K9b"}
assert len(search_plans) == len(KNIVES)
assert seen == {(name, tag) for name in counts for tag in baseline}
assert not any(x.startswith('WL_LOCAL_') for x in inventories['K10']['tags'])
assert witnesses['K10']['site'] and 'WL_S9B_BRANCH_EXISTENCE' in witnesses['K10']['site']

engine_path = ROOT / 'research/pde_ledger_v3/mathematica/S9b_light_bending_mathematica_audit.wl'
engine_bytes = engine_path.read_bytes()
assert (ROOT / 'research/pde_ledger_v3/mathematica/S9b_light_bending_mathematica_ablation.wl').read_bytes() == (STATE / 'harness-pinned-final.wl').read_bytes()
assert engine_bytes == subprocess.check_output(['git', 'show', 'HEAD:'+str(engine_path.relative_to(ROOT))], cwd=ROOT)
source = engine_bytes.decode()

current_harness = (ROOT / 'research/pde_ledger_v3/mathematica/S9b_light_bending_mathematica_ablation.wl').read_text()
previous_harness = (STATE / 'before_harness.wl').read_text()
# Only the K10 manifest and scratch destination change outside the witness block.
def outside_search(text):
    left,rest=text.split('(* BEGIN_WITNESS_METHOD *)',1)
    _,right=rest.split('(* END_WITNESS_METHOD *)',1)
    left=left[:left.index('knives = ')] + left[left.index('(* BEGIN_LITERAL_METHOD *)'):]
    return (left+right).replace('"s9b_wl_repair9_amend12"','"s9b_wl_repair10"')
assert outside_search(current_harness)==outside_search(previous_harness)
assert (ROOT/'research/pde_ledger_v3/mathematica/S9b_light_bending_mathematica_ablation_checks.wl').read_bytes()==(STATE/'checks-pinned-final.wl').read_bytes()
assert (STATE/'checks-pinned-final.wl').read_text().startswith((STATE/'before_checks.wl').read_text())
# Restore only the named witness edits and compare every remaining byte.
def witness_block(text):
    return text.split('(* BEGIN_WITNESS_METHOD *)',1)[1].split('(* END_WITNESS_METHOD *)',1)[0]
old_w=witness_block(previous_harness);new_w=witness_block(current_harness)
a=new_w.index('(* Round 10: decide the engine');b=new_w.index('(* END_ROUND9_SITE_SEARCH *)',a)
oa=old_w.index('searchReachedSite["K10"');ob=old_w.index('(* END_ROUND9_SITE_SEARCH *)',oa)
new_w=new_w[:a]+old_w[oa:ob]+new_w[b:]
a=new_w.index('(* K10 is introduced only at propagationGate');b=new_w.index('publicPart[record_]',a)
new_w=new_w[:a]+new_w[b:]
new_w=new_w.replace('"K10"->{},','"K10"->{"Domain"},')
new_w=new_w.replace('"K10"->{"ray",{"WL_S9B_BRANCH_EXISTENCE"},"propagationGate = rayConstruction"}',
 '"K10"->{"ray",{"WL_S9B_A_NONRECIPROCAL_PATH_DEPENDENCE","OnDomain"},"physicalPropagation = branchGateLocal"}')
new_w=new_w.replace('extractedNameKeys={"FreeInputs","RemainingInputs"};','extractedNameKeys={"Profiles","FreeInputs","RemainingInputs"};')
new_w=new_w.replace('add[z,path,Which[Length[path]>0&&Last[path]==="Profiles","literal profile list",StringQ[z],','add[z,path,Which[StringQ[z],')
new_w=new_w.replace('Which[Length[path]>0&&Last[path]==="Profiles","literal descriptive list of profile symbols written directly at engine 401; not extracted from expressions, not a condition, and no profile-function substitution is applied to its elements",ListQ[z],','Which[ListQ[z],')
new_w=new_w.replace('!If[knife==="K10",k10PartUses[rec],siteUses[knife,rec["Tag"],rec["Location"]]]','!siteUses[knife,rec["Tag"],rec["Location"]]')
new_w=new_w.replace('If[AssociationQ[#1],If[#1["Role"]==="literal profile list",<|"Evaluation"->#1["Object"],"Computation"->"evaluate the source-written list of symbols literally; no profile replacement; not an extracted-name or text comparison"|>,siteEvaluate[#1["Object"],#2]],"not emitted"]','If[AssociationQ[#1],siteEvaluate[#1["Object"],#2],"not emitted"]')
assert new_w==old_w,'Change outside the named K10/Profiles witness edits'

# Compare all source manifest entries except the expressly repaired K10.
def manifest(text):
    raw=text.split('knives = ',1)[1].split('(* BEGIN_LITERAL_METHOD *)',1)[0].strip().removesuffix(';')
    return association(raw)
old_manifest=manifest(previous_harness);new_manifest=manifest(current_harness)
assert old_manifest.keys()==new_manifest.keys()
assert all(old_manifest[k]==new_manifest[k] for k in old_manifest if k!='"K10"')

mutations = {}
for name in KNIVES:
    mutated = (WORK / (name+'.wl')).read_text()
    original_site = json.loads(sites[name]['"Construction"'])
    replacement_site = json.loads(sites[name]['"Mutation"'])
    mutations[name] = {'construction': original_site, 'mutation': replacement_site}
    assert source.count(original_site) == 1
    assert mutated.rstrip('\n') == source.replace(original_site, replacement_site).rstrip('\n')
# The repaired knife wraps exactly the original completed assignment RHS.
k10_original=mutations['K10']['construction'];k10_mutation=mutations['K10']['mutation']
assert k10_original.startswith('propagationGate = ')
assert k10_mutation=='propagationGate = Not['+k10_original.removeprefix('propagationGate = ').removesuffix(';')+'];'
assert 'branchGateLocal = propagatingLocal;' in (WORK/'K10.wl').read_text()
assert inventories['K7']['literalProfilesParts']>0 and inventories['K7']['oneCopyParts']>0
boundary = '(* COMPACT_K11_BOUNDARY *)'
for suffix, src in [('BASELINE', source), ('CORRUPTED', (WORK / 'K11.wl').read_text())]:
    assert (WORK / ('K11_COMPACT_'+suffix+'.wl')).read_text() == src.split(boundary)[0]

out_status = subprocess.check_output(['git', 'status', '--porcelain', 'research/pde_ledger_v3/mathematica/out/'], cwd=ROOT).decode()
assert out_status == ''
run_names = [RUN, 's9b-r10-checks-final', 's9b-r10-engine-repo', 's9b-r10-engine-alone']
runs = {}
for name in run_names:
    folder = ROOT / '_scratch/s11c' / name
    outcome = json.loads((folder / 'outcome.json').read_text())
    assert outcome['exitCode'] == 0 and outcome['childOutcome']['guardReason'] is None
    stderr = (folder/'stderr').read_text()
    if name == RUN:
        warning = 'Limit::alimv: Warning: Assumptions that involve the limit variable are\n    ignored.'
        suppression='General::stop: Further output of Limit::alimv\n     will be suppressed during this calculation.'
        assert not stderr.replace(warning,'').replace(suppression,'').strip(), stderr
    else:
        assert stderr == ''
    limits=json.loads((folder/'effective-limits.json').read_text())
    assert limits['memory.max']=='8589934592' and limits['nativeAddressSpaceBytes']==[8589934592,8589934592]
    assert limits['pids.max']=='64' and limits['memory.swap.max']=='0'
    samples=[json.loads(x) for x in (folder/'resource-samples.jsonl').read_text().splitlines()]
    assert all('oom_kill 0' in x['memory.events'] for x in samples)
    runs[name] = {'outcome':outcome,'invocation':json.loads((folder/'invocation.json').read_text()),
        'peakMemoryBytes':max(int(x['memory.peak']) for x in samples),
        'maximumSampledTasks':max(int(x['pids.current']) for x in samples),
        'maximumSampledSwapBytes':max(int(x['memory.swap.current']) for x in samples),
        'memoryEvents':samples[-1]['memory.events'], 'memoryMaxBytes':int(limits['memory.max']),
        'nativeAddressSpaceBytes':limits['nativeAddressSpaceBytes'],'stderr':stderr}
growth = json.loads((STATE / (RUN+'.run.json')).read_text())['growthWhileKernelAlive']
assert growth
checks = capture(ROOT / '_scratch/s11c/s9b-r10-checks-final/stdout')
assert checks['HARNESS_PARSE_HEAD'] == 'HoldComplete'
assert 'propagationGate = Not[rayConstruction' in checks['WHOLE_BRANCH_CONDITION_NEGATION']
assert 'literal profile list' in checks['LITERAL_PROFILES_READING']
assert '"ExtractedNameKeys" -> {"FreeInputs", "RemainingInputs"}' in checks['LITERAL_PROFILES_READING']
assert len(checks) == 23
assert 'unique linear-GM solution' in checks['SEMANTIC_BRANCH_CATALOG']
assert 'negative square-root solution' in checks['SEMANTIC_BRANCH_CATALOG']
assert 'positive square-root solution' in checks['SEMANTIC_BRANCH_CATALOG']
assert 'rho <= 0' in checks['DENSITY_SIGN_CATALOG'] and 'rho > 0' in checks['DENSITY_SIGN_CATALOG']
assert association(checks['NEGATIVE_DENSITY_ASYMPTOTE'])['"AsymptoticDensity"'] == '1'
assert association(checks['NEGATIVE_DENSITY_ASYMPTOTE'])['"DensityAtPoint"'] == '-2' 
assert association(association(checks['SAME_SYNTAX'])['"Comparison"'])['"Mark"'] == '"identical"'
for key in checks.keys()-{'HARNESS_PARSE_HEAD', 'SAME_SYNTAX', 'COMMON_FIXED_SOLUTION_SETS', 'DISJOINT_FIXED_SOLUTION_SETS', 'ORDERED_CASE_SELECTION', 'DOMAIN_BOOLEAN_DECISION', 'DOMAIN_UNRESOLVED_DECISION', 'BOUND_VARIABLE_SCOPE', 'BOUND_ASSOCIATION_SCOPE', 'SEMANTIC_BRANCH_CATALOG', 'DENSITY_SIGN_CATALOG', 'TEXT_COMPARISON_DECISIONS', 'NEGATIVE_DENSITY_ASYMPTOTE','WHOLE_BRANCH_CONDITION_NEGATION','LITERAL_PROFILES_READING'}:
    assert association(association(checks[key])['"Comparison"'])['"Mark"'] == '"not compared"'

def write(suffix, text):
    (DEST / (PREFIX+suffix)).write_text(text)

for suffix, tags in [('tags.txt', baseline), ('not_established_tags.txt', [x for x in baseline if 'NOT_ESTABLISHED' in baseline[x]]), ('harness_not_established_tags.txt', harness_ne), ('k11_evaluated_tags.txt', compact_b), ('k11_not_evaluated_tags.txt', baseline.keys()-compact_b.keys())]:
    write(suffix, '\n'.join(sorted(tags))+'\n')
hashes = {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in [engine_path, ROOT / 'research/pde_ledger_v3/mathematica/S9b_light_bending_mathematica_ablation.wl', ROOT / 'research/pde_ledger_v3/mathematica/S9b_light_bending_mathematica_ablation_checks.wl', Path(__file__).resolve(), ROOT / 'scripts/s11c_guarded_run.py']}
write('witness_searches.txt', ''.join(witness_lines))
write('harness.stderr.txt',(FULL/'stderr').read_text())
capture_ne=[]
for name,data in {'BASELINE':baseline,**captures,'K11_COMPACT_BASELINE':compact_b,'K11_COMPACT_CORRUPTED':compact_a}.items():
    capture_ne.extend(name+':'+tag for tag,payload in data.items() if 'NOT_ESTABLISHED' in payload)
write('capture_not_established_tags.txt','\n'.join(capture_ne)+'\n')
import gzip
with gzip.open(DEST / (PREFIX+'site_inventory.txt.gz'), 'wt') as archive:
    archive.write(''.join(site_lines+search_plans))
repo_output = (ROOT / '_scratch/s11c/s9b-r10-engine-repo/stdout').read_bytes()
alone_output = (ROOT / '_scratch/s11c/s9b-r10-engine-alone/stdout').read_bytes()
assert repo_output == alone_output == (WORK/'BASELINE.stdout').read_bytes()
alone_folder=STATE/'engine-alone'
assert sorted(x.name for x in alone_folder.iterdir()) == ['S9b_light_bending_mathematica_audit.wl']
assert (alone_folder/'S9b_light_bending_mathematica_audit.wl').read_bytes() == engine_bytes
write('k11_reached_pairs.txt', ''.join(k11_lines))
write('structural_checks.txt', (ROOT / '_scratch/s11c/s9b-r10-checks-final/stdout').read_text())
write('source_hashes.json', json.dumps(hashes, indent=2)+'\n')
summary = {'status':'round-10 full harness and checks completed from final code at 8 GiB','tags': len(baseline), 'shared': len(shared), 'local': len(local), 'missingShared': [], 'baselineNotEstablished': sum('NOT_ESTABLISHED' in x for x in baseline.values()), 'structuralCounts': counts, 'witnessSearches': witnesses, 'K11Reached': len(compact_b), 'K11Unreached': len(baseline)-len(compact_b), 'engineUnchanged': True, 'engineBytes': len(engine_bytes), 'productionOutStatus': out_status, 'growthWhileKernelAlive': growth, 'mutations': mutations, 'runs': runs, 'siteInventories': inventories, 'densitySignTrials': sign_trials, 'standaloneByteIdentity':True,'guardRefusals':[],'guardKills':[],'item12Events':[],'commits':[],'sourceUnchangedFromFinalPins':True,'scope':'K10 construction/search/dependencies; literal Profiles reading; round-8 historical note'}
write('validation.json', json.dumps(summary, indent=2)+'\n')
print(json.dumps({key: summary[key] for key in ['tags', 'shared', 'local', 'baselineNotEstablished', 'K11Reached', 'K11Unreached', 'engineUnchanged', 'productionOutStatus']}, indent=2))
print(json.dumps({'structuralCounts': counts, 'witnessSearches': witnesses, 'siteInventories': inventories, 'densitySignTrials': sign_trials}, indent=2))

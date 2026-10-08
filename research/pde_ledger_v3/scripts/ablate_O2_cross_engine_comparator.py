#!/usr/bin/env python3
"""Run isolated comparator-logic ablations, under the caller's pooled guard.

No production transcripts or engine imports. Every variant runs the complete
serialized synthetic suite. Assertion failures are required; operational errors,
kills, and surviving variants stop the run. No retries or duration limits.
"""
import argparse
import difflib
import json
import os
from pathlib import Path
import re
import subprocess
import sys


# label, source replacement, required failing test. Each edit removes comparator
# logic, not a fixture, assertion, reader invocation or test-runner behaviour.
ABLATIONS = (
    ('same_role_join',
     "'EXCHANGE_MOMENTUM/MaterialIdentification/@1',378,212",
     "'MASS_INPUT/RHS',378,212", 'test_mass_loss_join_keeps_rhs_separate'),
    ('plain_sign', 'sign *= -1', 'sign *= 1', 'test_plain_sum_orientation'),
    ('held_aggregate_sign', '            role = None\n',
     '            linear_arguments = ()\n            role = None\n', 'test_held_aggregate_orientation'),
    ('sympy_derivative_body_sign',
     'signed_actions(child,sign if i == body else 1,dependency or i != body)',
     'signed_actions(child,1,dependency or i != body)', 'test_all_held_linear_orientations'),
    ('applied_functions_as_symbols',
     "return sp.Function(key)(*(algebra(x,bindings,engine) for x in a[1:]))",
     "return sp.Symbol(key)", 'test_exact_argument_content'),
    ('addition_instead_of_subtraction', 'delta = sp.cancel(left-right)',
     'delta = sp.cancel(left+right)', 'test_common_translation_preserves_residual'),
    ('constant_residual', 'delta = sp.cancel(left-right)',
     'delta = sp.Integer(0)', 'test_one_sided_operand'),
    ('operand_binding_removed',
     "('S12_reaction_system','OPENReactionSystem','§3.3,5 S12 additional momentum reaction system',160,209)",
     "('S12_reaction_system','DifferentReactionSystem','§3.3,5 S12 additional momentum reaction system',160,209)",
     'test_open_operand_binding_consistency'),
    ('closed_part_omitted', 'closed_terms.append(term)', 'closed_terms.extend(())',
     'test_closed_balance_entry'),
    ('evaluation_point_collapsed', '    return algebra(n,bindings,engine)\n\n\ndef layout',
     "    return sp.Symbol('fixed_point',real=True)\n\n\ndef layout", 'test_derivative_evaluation_points'),
    ('open_field_deltas_omitted',
     'field:bag_delta(left.get(field,{}),right.get(field,{}))',
     'field:[]', 'test_action_field_differences'),
    ('action_repoint_ignored',
     "return {'py':{r.py:r.py for r in actions},'wl':{r.wl:r.py for r in actions}}",
     "return {'py':{r.py:r.py for r in ACTION_TABLE},'wl':{r.wl:r.py for r in ACTION_TABLE}}",
     'test_each_action_binding_repoint'),
    ('transpose_dropped', "if row.layout == 'transpose' and a.head == 'Matrix':",
     'if False:', 'test_transpose_layout_matches_explicit_transpose'),
    ('canonical_derivative_key_removed',
     "return Node('ProfileDerivative',(seq(*(convert(v) for v in orders)),",
     "return Node('WolframDerivativeSyntax',(seq(*(convert(v) for v in orders)),",
     'test_profile_derivative_is_one_object_key'),
    ('canonical_sqrt_removed',
     "return convert(Node('Pow',(a[1],atom('Number','1/2'))))",
     "return Node('WolframSqrtSyntax',(convert(a[1]),))",
     'test_identical_open_objects_have_empty_deltas'),
    ('object_sets_as_occurrence_counts',
     'inventory.update({key:1 for key in values if key not in inventory})',
     'inventory.update(values)', 'test_object_inventory_ignores_repeated_spelling'),
    ('relation_left_only', 'for x,y in zip(a.args,b.args):',
     'for x,y in zip(a.args[:1],b.args[:1]):', 'test_each_relation_operand_both_engines'),
    ('relation_right_only', 'for x,y in zip(a.args,b.args):',
     'for x,y in zip(a.args[1:],b.args[1:]):', 'test_each_relation_operand_both_engines'),
    ('closed_part_first_term', "Node('Add',tuple(closed_terms)) if closed_terms",
     "Node('Add',tuple(closed_terms[:1])) if closed_terms", 'test_each_closed_term_in_multiterm_balance'),
    ('function_metadata_lost', "'Function': 'FunctionName'}[h], value=a[0].value, options=options)",
     "'Function': 'FunctionName'}[h], value=a[0].value, options=())",
     'test_lossless_function_metadata_survives'),
)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    if not os.environ.get('S11C_POOLED_GUARD_MANIFEST'):
        raise RuntimeError('run through s11c_guarded_run.py in pooled mode')
    args.output.mkdir(parents=True,exist_ok=False)
    source=Path(__file__).with_name('O2_cross_engine_comparator.py').read_text()
    tests=Path(__file__).with_name('test_O2_cross_engine_comparator.py').read_text()
    summaries=[]
    for label,before,after,required in ABLATIONS:
        if source.count(before)!=1:
            raise RuntimeError('ablation anchor is not unique: '+label)
        folder=args.output/label
        folder.mkdir()
        modified=source.replace(before,after)
        (folder/'O2_cross_engine_comparator.py').write_text(modified)
        (folder/'test_O2_cross_engine_comparator.py').write_text(tests)
        (folder/'change.diff').write_text(''.join(difflib.unified_diff(
            source.splitlines(True),modified.splitlines(True),fromfile='baseline',tofile=label)))
        command=['/usr/bin/time','-f','elapsed_seconds=%e\npeak_rss_kib=%M',
                 '-o',str((folder/'resources.txt').resolve()),sys.executable,
                 str((folder/'test_O2_cross_engine_comparator.py').resolve())]
        with (folder/'stdout').open('w') as out,(folder/'stderr').open('w') as err:
            result=subprocess.run(command,stdout=out,stderr=err,
                                  env={**os.environ,'PYTHONDONTWRITEBYTECODE':'1'})
        transcript=(folder/'stderr').read_text()
        failures=re.findall(r'^FAIL: (.+)$',transcript,re.M)
        errors=re.findall(r'^ERROR: (.+)$',transcript,re.M)
        summary={'ablation':label,'command':command,'exit_code':result.returncode,
                 'failing_tests':failures,'errors':errors,'required_test':required,
                 'resources':(folder/'resources.txt').read_text()}
        summaries.append(summary)
        (args.output/'summary.json').write_text(json.dumps(summaries,indent=2)+'\n')
        print(json.dumps(summary),flush=True)
        if result.returncode!=1 or errors or not any(required in test for test in failures):
            raise RuntimeError('ablation did not produce its required assertion failure: '+label)


if __name__=='__main__':
    main()

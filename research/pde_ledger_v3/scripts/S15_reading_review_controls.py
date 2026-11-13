"""Apply one-sided R0 mutations to accepted native records via the candidate readers."""
import argparse
import gzip
import hashlib
import json
from collections import defaultdict
from pathlib import Path
import sys

import S15_reading_comparator as c
import S15_reading_context_checks as k
import S15_reading_math as m


def main(root):
    cfg = json.loads((root/'configuration.json').read_text())
    index = json.load(gzip.open(root/'coverage.json.gz'))
    refs = defaultdict(lambda: defaultdict(list))
    lookup = {}
    for row in index:
        stream = row['id'].split(':')[0]
        if stream in ('py_reading','wl_reading') and 'payload' in row:
            refs[stream[:2]][row['tag']].append(row['id'])
            lookup[row['id']] = stream,row['span'],row['payload']

    def read(rid):
        stream, span, digest = lookup[rid]
        path = c.inputs.ROOT / cfg['inputs'][stream]['path']
        with path.open('rb') as f:
            payload = c.inputs.payload(f, span)
        if hashlib.sha256(payload).hexdigest() != digest:
            raise ValueError('native digest mismatch ' + rid)
        return payload.decode()

    def edit(text, old, new, count=1):
        if text.count(old) < count:
            raise ValueError('missing native mutation target '+old)
        return text.replace(old, new, count)

    family = list(c.family_context(refs['py'],refs['wl'],read))
    if any(row['status']!='same_values' for row in family):
        raise ValueError(family)
    for engine, input_tag, field, old, new in [
        ('py','family/centre_substitution','family/centre_substitution','X_1_1(t):','X_2_1(t):'),
        ('py','family/centre_substitution','family/centre_substitution',' + X_1_1(t),',' + X_2_1(t),'),
        ('wl','profile_family/domain','profile_family/chart_condition','Reals','Integers'),
        ('wl','profile_family/domain','profile_family/chart_condition',' | aa[2]','')]:
        target, = refs[engine][input_tag]
        def mutated(rid):
            return edit(read(rid),old,new) if rid == target else read(rid)
        row, = (r for r in c.family_context(refs['py'],refs['wl'],mutated) if r['field']==field)
        print(field, old, '->',new,row['status'],flush=True)
        assert row['status']=='value_difference'

    for tag,old,new,checker,prefix,expected in [
        ('distance_endpoints','-integration_radius + kernel_distance','integration_radius - kernel_distance',
         k.kernel_checks,'positive endpoint 0','remaining_symbolic_difference'),
        ('azimuth_operand','(az, 0, 2*pi)','(az, 0, pi)',
         k.angular_moments,'angular moment (3, 0, 0)','azimuth_domain_difference')]:
        for changed in (False,True):
            cat = k.Catalogue(root)
            if changed:
                key = cat.tags['py','item3_kernel/'+tag][0]['payload']
                cat.blobs[key]=edit(cat.blobs[key],old,new)
            checker(cat)
            row, = (r for r in cat.rows if r['field']==prefix)
            print(prefix, 'modified' if changed else 'original', row['status'],
                  row.get('domain_residual',row.get('py_bounds')), flush=True)
            assert row['status']==(expected if changed else 'exact_identity')

    tag='reading/both_surfaces/sum_crossing/component_1/E_resp'
    left=m.parse(read(refs['py'][tag][0]),'py')
    wl=read(refs['wl'][tag][0])
    changes=[(None,None,0),('change[oo[1, 2]]','change[oo[1, 1]]',1),
             ('change[oo[1, 2]]','change[oo[1, 3]]',1),
             ('change[oo[1, 1]]','change[oo[1, 2]]',2)]
    failures=[]
    for old,new,n in changes:
        right=m.parse(edit(wl,old,new,n) if n else wl,'wl')
        outcome=m.sample_values(left,right,())
        first=outcome['samples'][0]
        print('net response',old,new,n,outcome['sample_counts'],first['values'][-1],
              first['residual'],first['tolerance'],flush=True)
        passed = {tuple(row['case'][:2]) for row in outcome['samples'] if row['status']=='separated'}
        print('detected form/stratum pairs', sorted(passed), flush=True)
        if n and passed!={(f,t) for f in m.FORMS for t in m.STRATA}:
            failures.append((old,new,'missing detections',sorted(passed)))
        if not n and outcome['sample_counts']!={'within_estimate':30}:
            failures.append(('unchanged',outcome['sample_counts']))
    if failures:
        raise ValueError(failures)


if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('root',type=Path)
    sys.setrecursionlimit(10000)
    main(ap.parse_args().root)

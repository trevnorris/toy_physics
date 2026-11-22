"""Measure shear retention, compression and rearrangement in finite 4D receiver holds."""
from collections import defaultdict
import argparse
import hashlib
import json
from math import sqrt
from pathlib import Path
from statistics import covariance, mean, stdev, variance

from scipy.stats import t


def interval(values):
    mid = mean(values)
    err = t.ppf(.975, len(values)-1)*stdev(values)/sqrt(len(values))
    return {'means_by_preparation': values, 'mean':mid, 't95':[mid-err,mid+err]}


def fraction(early, late):
    n = len(early)
    q = t.ppf(.975,n-1)**2/n
    x,y = mean(early),mean(late)
    a,b,c = x*x-q*variance(early), x*y-q*covariance(early,late), y*y-q*variance(late)
    disc=b*b-a*c
    result = dict(early=interval(early), late=interval(late), ratio=y/x if x else None)
    if a>0 and disc>=0:
        result['bounded_ratio_t95']=[(b-sqrt(disc))/a,(b+sqrt(disc))/a]
    else:
        # An unbounded or empty set cannot certify a small retained fraction.
        result['bounded_ratio_t95']=None
    return result


def weighted(blocks,key):
    return sum(b['steps']*b['mean'][key] for b in blocks)/sum(b['steps'] for b in blocks)


def analyze(root):
    design=json.loads((root/'design.json').read_text());done=json.loads((root/'DONE.json').read_text())
    groups=defaultdict(list); per_run={}
    expected={case['name'] for cases in design['jobs'].values() for case in cases}
    if set(done['results'])!=expected: raise ValueError('Incomplete material comparison')
    for worker, cases in design['jobs'].items():
        rows={}
        for case in cases:
            name=case['name']; raw=(root/'runs'/(name+'.json')).read_bytes()
            if hashlib.sha256(raw).hexdigest()!=done['results'][name]:raise ValueError(name)
            data=json.loads(raw);b=data['blocks'];assert len(b)==16
            assert (data['strain'],data['pressure_external'],data['seed'],data['history'],data['axis_index'])==(case['strain'],case['p_ext'],case['seed'],case['history'],case['axis'])
            rows[case['label']]=data
        zero=rows['zero']; assert all(row['prepared']==zero['prepared'] for row in rows.values())
        report = {'pressure':weighted(zero['blocks'][-4:],'pressure'),'density':weighted(zero['blocks'][-4:],'density'),
                  'temperature':weighted(zero['blocks'][-4:],'temperature'), 'min_cluster':min(b['structure']['max_cluster'] for b in zero['blocks']),
                  'absent_fraction': zero['blocks'][-1]['structure']['absent']/zero['blocks'][-1]['structure']['original'],
                  'nonaffine_rms':zero['blocks'][-1]['structure']['nonaffine_rms']}
        for arm,label in [('g','full'),('h','half')]:
            p,m=rows[arm+'p'],rows[arm+'m']; den=p['strain']-m['strain']
            report[label+'_early']=(p['initial_reading']['shear']-m['initial_reading']['shear'])/den
            report[label+'_late']=(weighted(p['blocks'][-4:],'shear')-weighted(m['blocks'][-4:],'shear'))/den
            report[label+'_windows']=[(weighted(p['blocks'][i:i+4],'shear')-weighted(m['blocks'][i:i+4],'shear'))/den for i in range(0,16,4)]
        report['full_minus_half_late']=report['full_late']-report['half_late']
        p,m=rows['pp'],rows['pm'];den=p['pressure_external']-m['pressure_external']
        report['compliance']=(weighted(m['blocks'][-4:],'log_volume')-weighted(p['blocks'][-4:],'log_volume'))/den
        report['pressure_difference']=weighted(p['blocks'][-4:],'pressure')-weighted(m['blocks'][-4:],'pressure')
        per_run[worker]=report;groups[zero['history']].append(report)
    def prep(key):return [mean(d[key] for d in g) for _,g in sorted(groups.items())]
    responses={label:fraction(prep(label+'_early'),prep(label+'_late')) for label in ('full','half')}
    summary={k:interval(prep(k)) for k in ('density','pressure','temperature','absent_fraction','nonaffine_rms','full_minus_half_late','compliance','pressure_difference')}
    min_cluster=min(r['min_cluster'] for r in per_run.values()); summary['minimum_largest_cluster']=min_cluster
    bounded=[r['bounded_ratio_t95'] for r in responses.values()]
    small_bound=all(b is not None and b[0]>=-.25 and b[1]<=.25 for b in bounded)
    return dict(cases=len(expected),histories=per_run,ratio=responses,summary=summary,
                neutral_screen=dict(small_late_shear_at_both_amplitudes=small_bound,all_sites_connected=min_cluster==128,
                                    resolved_compression=bool(summary['compliance']['t95'][0]>0)))


if __name__=='__main__':
    cli=argparse.ArgumentParser(description=__doc__);cli.add_argument('root',type=Path);cli.add_argument('out',type=Path)
    args=cli.parse_args();r=analyze(args.root)
    args.out.write_text(json.dumps(r,indent=2,allow_nan=False)+'\n')
    print(json.dumps({k:v for k,v in r.items() if k!='histories'},allow_nan=False))

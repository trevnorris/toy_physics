"""Evaluate the additional joint-branch cases that contain derivatives of Li2."""
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from contextlib import ExitStack
import gzip
import json
import multiprocessing
import os
from pathlib import Path
import sys
import numpy as np
import sympy as s
import S15_reading_refinement as ref
from S15_reading_refinement import c,m


class LogBranchSample(m.Sample):
    def compute(self,v,environment):
        if v.func==s.polylog and v.args[0] in (0,1):
            z=self.value(v.args[1],environment)
            if not (-1<z<1):raise m.Unsupported('outside the printed real log chart')
            return -np.log1p(-z) if v.args[0]==1 else z/(1-z)
        return super().compute(v,environment)


def work(path):
    obj=json.load(gzip.open(path))
    with ExitStack() as files:
        handles={k:files.enter_context(p.open('rb')) for k,p in c.PATHS.items()}
        cmp=c.ObjectComparison(handles)
        out=[]
        for field in obj['fields']:
            for row in field['additional_samples']:
                if row['status']!='unresolved' or row.get('reason')!='polylog order outside Li2':continue
                a,b=map(cmp.operand,(field['left'][0],field['right'][0]))
                form,stratum,point,_=row['case']
                exponents=tuple(map(s.Rational,row['exponents']))
                vals=[]
                for order in m.RULES:
                    ctx=LogBranchSample(form,stratum,point,exponents,order)
                    vals.append([float(ctx.value(v)) for v in (a,b)])
                if not np.isfinite(vals).all():raise ValueError('nonfinite log branch')
                x,y=vals[-1]
                scale=max(1.,abs(x),abs(y))
                delta=[max(abs(vals[j][i]-vals[j-1][i]) for i in (0,1)) for j in (1,2)]
                tol=8*max(delta)+1.e-9*scale
                result=dict(field=field['field'],case=row['case'],exponents=row['exponents'],values=vals,residual=x-y,
                            tolerance=tol,status='within_estimate' if abs(x-y)<=tol and delta[-1]<=1.e-6*scale else 'unresolved')
                if result['status']=='unresolved' and not (a.has(s.Integral) or b.has(s.Integral)):
                    result['precision_check']=ref.refine_closed(a,b,result)
                    result['status']=result['precision_check']['status']
                out.append(result)
    written=dict(object=obj['object'],build_id=obj['build_id'],rechecks=out)
    c.save(folder/path.name,written)
    return Counter(r['status'] for r in out)


if __name__=='__main__':
    sys.setrecursionlimit(10000)
    base=Path(sys.argv[1])
    config=json.loads((base/'configuration.json').read_text())
    folder=base/'log_branches'
    folder.mkdir(exist_ok=False)
    c.PATHS={k:c.inputs.ROOT/v['path'] for k,v in config['inputs'].items()}
    c.LOCATOR={r['id']:(r['id'].split(':')[0],r['span'],r['payload'])
               for r in json.load(gzip.open(base/'coverage.json.gz')) if 'payload' in r}
    jobs=[path for path in sorted((base/'refinement').glob('*.gz')) if json.load(gzip.open(path))['branch']]
    with ProcessPoolExecutor(len(os.sched_getaffinity(0)),mp_context=multiprocessing.get_context('fork')) as pool:
        total=sum(pool.map(work,jobs),Counter())
    print('previously unexercised log cases',dict(total))
    c.save(folder/'summary.json',dict(objects=len(jobs),counts=total))

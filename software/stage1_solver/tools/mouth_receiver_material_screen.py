"""Screen finite, pressure-held 4D preparations for density, order and site exchange."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from mouth_receiver_structure import grid, intensity, independent_peaks


def signature(state, known=None):
    box=np.array(state['box']); y=np.array(state['pos'])@np.linalg.inv(box)
    h=grid(box)
    powers=np.concatenate([intensity(h[i:i+256],y) for i in range(0,len(h),256)])
    idx=independent_peaks(h,powers)
    result={'peak':float(max(powers)),'four':powers[idx].tolist(),'indices':h[idx].tolist()}
    if known is not None:
        result['previous_four']=intensity(np.asarray(known),y).tolist()
    return result


def read_one(case, state):
    assert state['source_temperature']==1.024 and state['temperature']==case['temperature']
    assert state['strain']==state['pressure_external']==0 and (state['seed'],state['history'])==(case['seed'],case['history'])
    blocks=state['blocks']; last=blocks[-4:]
    def avg(b,k):return sum(x['steps']*x['mean'][k] for x in b)/sum(x['steps'] for x in b)
    density=avg(last,'density');previous=avg(blocks[-8:-4],'density'); ratio=density/previous
    cluster=min(b['structure']['max_cluster'] for b in last);end=last[-1]['structure'];lost=end['absent']/end['original'] if end['original'] else None
    dense=density>=1.2 and cluster>=116 and abs(ratio-1)<=.10
    pre=post=None
    # Evaporated cases already fail the density test. Their empty volume need not be surveyed for peaks.
    if dense:
        pre=signature(state['prepared']);post=signature(state['final'],pre['indices'])
    return dict(history=case['history'],temperature=state['temperature'],rho=density,rho_last_over_previous=ratio,
                pressure=avg(last,'pressure'),minimum_late_cluster=cluster,original_missing=lost,rms=end['nonaffine_rms'],
                prepared_fourier=pre,ending_fourier=post,retains_dense_sample=dense,
                qualifies_for_shear_probe=(dense and pre['peak']<=.2 and post['peak']<=.2 and lost is not None and lost>=.5 and end['nonaffine_rms']>=1.))


def main(root, output):
    d=json.loads((root/'design.json').read_text());manifest=json.loads((root/'DONE.json').read_text());result={}
    for group in d['jobs'].values():
        for c in group:
            raw=(root/'runs'/(c['name']+'.json')).read_bytes();assert hashlib.sha256(raw).hexdigest()==manifest['results'][c['name']]
            result[c['name']]=read_one(c,json.loads(raw))
    qualified=[]
    for temp in sorted({v['temperature'] for v in result.values()}):
        sample=[v for v in result.values() if v['temperature']==temp];assert len({s['history'] for s in sample})==3
        if all(s['qualifies_for_shear_probe'] for s in sample):qualified.append(temp)
        print(json.dumps(dict(temperature=temp,summary=[{k:z[k] for k in ('history','rho','pressure','minimum_late_cluster','original_missing','rms','qualifies_for_shear_probe')}|{'end_peak':z['ending_fourier']['peak'] if z['ending_fourier'] else None} for z in sample]),allow_nan=False))
    output.write_text(json.dumps(dict(cases=result,screened_temperatures=qualified,selected=min(qualified) if qualified else None),indent=2,allow_nan=False)+'\n')

if __name__=='__main__':
    cli=argparse.ArgumentParser(description=__doc__);cli.add_argument('root',type=Path);cli.add_argument('output',type=Path)
    args=cli.parse_args();main(args.root,args.output)

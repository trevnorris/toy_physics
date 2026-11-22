"""Distinguish how an independently chosen preparation affects structure and shear."""
import argparse
from hashlib import sha256
import json
from pathlib import Path
from statistics import mean

from mouth_receiver_pressure_readout import fraction, interval, weighted
from mouth_receiver_material_screen import signature


def main(base,out):
    read=lambda p:json.loads(p.read_bytes());design=read(base/'design.json');done=read(base/'DONE.json');summaries={}
    for group, cases in design['jobs'].items():
        rows={}
        for c in cases:
            f=base/'runs'/(c['name']+'.json');payload=f.read_bytes();assert sha256(payload).hexdigest()==done['results'][c['name']]
            v=json.loads(payload);assert (v['temperature'],v['pressure_external'],v['strain'],v['seed'],v['axis_index'],v['history'])==(1.088,c['p_ext'],c['strain'],c['seed'],1,c['history']);rows[c['label']]=v
        zero=rows['zero'];assert all(r['prepared']==zero['prepared'] for r in rows.values());pre=signature(zero['prepared'])
        arm={}
        for k,v in rows.items():
            blocks=v['blocks']; end=blocks[-1]['structure']; sig=signature(v['final'],pre['indices'])
            arm[k]=dict(initial=v['initial_reading']['shear'],late=weighted(blocks[-4:],'shear'),late_density=weighted(blocks[-4:],'density'),
              late_pressure=weighted(blocks[-4:],'pressure'),late_logv=weighted(blocks[-4:],'log_volume'),late_temperature=weighted(blocks[-4:],'temperature'),
              min_cluster=min(b['structure']['max_cluster'] for b in blocks),lost=end['absent']/end['original'],rms=end['nonaffine_rms'],
              peak=sig['peak'],four=sig['four'],preselected_peaks_at_end=sig['previous_four'])
        g=zero['history'];metric={}
        for k,a,lab in (('g',.02,'full'),('h',.01,'half')):
            metric[lab]={q:(arm[k+'p'][q]-arm[k+'m'][q])/(2*a) for q in ('initial','late')}
        metric['compliance']=(arm['pm']['late_logv']-arm['pp']['late_logv'])/.5; metric['pressure_response']=arm['pp']['late_pressure']-arm['pm']['late_pressure']
        summaries[g]=dict(source_temperature=zero['source_temperature'],prepared_density=zero['prepared_reading']['density'],prepared_peak=pre['peak'],prepared_four=pre['four'],
            arms=arm,measured=metric)
    fit={}
    for strength in ('full','half'):
        e=[summaries[h]['measured'][strength]['initial'] for h in range(3)];s=[summaries[h]['measured'][strength]['late'] for h in range(3)];fit[strength]=fraction(e,s)
    fit['compliance']=interval([summaries[h]['measured']['compliance'] for h in range(3)])
    result=dict(cases=done['groups'],history_ids={0:'cooled h0',1:'cooled h1',2:'cooled h2',3:'cooled independent cold seed 41',40:'continued ordered h0',41:'continued low-order h1',42:'continued ordered h2'},histories=summaries,training=fit)
    out.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    for k,h in summaries.items():
        z=h['arms']['zero'];m=h['measured'];print(json.dumps(dict(history=k,pre_peak=h['prepared_peak'],final_peak=z['peak'],driven_peak_range=[min(r['peak'] for r in h['arms'].values()),max(r['peak'] for r in h['arms'].values())],late_density=z['late_density'],min_cluster=z['min_cluster'],lost=z['lost'],rms=z['rms'],**m),allow_nan=False))
    print(json.dumps({'three_cooled_preparations':fit},allow_nan=False))
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('base',type=Path);p.add_argument('out',type=Path);a=p.parse_args();main(a.base,a.out)

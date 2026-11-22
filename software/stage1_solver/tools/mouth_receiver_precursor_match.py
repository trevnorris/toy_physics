"""Classify a fresh sample before any forces are read, at the same pressure and temperature."""
import argparse
from hashlib import sha256
import json
from pathlib import Path

from mouth_receiver_material_screen import signature
from mouth_receiver_pressure_readout import weighted


def main(root):
    plan=json.loads((root/'plan.json').read_bytes()); output={}
    for name,excursion in plan['excursions'].items():
        path=root/(name+'_back.json'); raw=path.read_bytes();d=json.loads(raw); g=root/(name+'_back_guard'); outcome=json.loads((g/'outcome.json').read_bytes());inv=json.loads((g/'invocation.json').read_bytes());assert outcome['childOutcome']['exitCode']==outcome['exitCode']==0 and outcome['limitsVerified'] and Path(inv['command'][2]).samefile(path) and d['temperature']==1.28 and d['source_temperature']==excursion and d['strain']==d['pressure_external']==0.
        b=d['blocks'];density=weighted(b[-4:],'density');drift=density/weighted(b[-8:-4],'density')-1.;cluster=min(z['structure']['max_cluster'] for z in b[-4:]);end=b[-1]['structure'];lost=end['absent']/end['original'] if end['original'] else None
        pre=signature(d['prepared']) if d['prepared_reading']['density']>.3 else None
        final=signature(d['final']) if density>.3 else None
        passes=(.95<=density<=1.2 and abs(drift)<=.10 and cluster>=116 and pre is not None and pre['peak']<=.2 and final is not None and final['peak']<=.2 and lost is not None and lost>=.5 and end['nonaffine_rms']>=1.)
        output[name]=dict(excursion=excursion,density=density,relative_density_change=drift,min_cluster=cluster,original_absent=lost,rms=end['nonaffine_rms'],pre_peak=None if pre is None else pre['peak'],end_peak=None if final is None else final['peak'],matches=passes,sha256=sha256(raw).hexdigest());print(json.dumps(dict(condition=name,**output[name]),allow_nan=False))
    passing=[k for k in output if output[k]['matches']];selected=min(passing,key=lambda k:output[k]['excursion']) if passing else None
    (root/'classification.json').write_text(json.dumps(dict(cases=output,selected=selected),indent=2,allow_nan=False)+'\n')
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('root',type=Path);main(p.parse_args().root)

#!/usr/bin/env python3
"""Focused independent current contraction and actual-census checks."""
import json
from pathlib import Path
import numpy as np
import S11c_d_continuum_currents as c


def main():
    rng=np.random.default_rng(11920)
    a={g:rng.normal(size=(3,2))+1j*rng.normal(size=(3,2)) for g in c.G}
    b={g:rng.normal(size=(3,3))+1j*rng.normal(size=(3,3)) for g in c.G}
    b={g:v+v.conj().T for g,v in b.items()};b[(0,0)]+=10*np.eye(3)
    value=c.quadratic(a,b);residual=c.direct_quadratic(a,b,value)
    c.f.require(c.boundary.norm(residual)<1e-11,'independent dense current contraction')
    changed={**b,(1,0):b[(1,0)]+b[(0,0)]};mutation=c.subtract(c.quadratic(a,changed),value)
    c.f.require(c.boundary.norm(mutation)>0,'one-sided current coefficient sensitivity')
    denominator={0:np.diag([2.,3.]),1:np.diag([.3,.8]),2:np.diag([.1,.4])}
    expected={i:rng.normal(size=2) for i in range(3)}
    numerator={i:np.diag(sum((np.diag(denominator[j])*expected[i-j] for j in denominator if j<=i),np.zeros(2))) for i in range(3)}
    quotient,qr=c.quotient(numerator,denominator)
    qdiff=max(float(np.max(abs(quotient[i]-expected[i]))) for i in expected)
    c.f.require(qdiff<1e-12,'independent quotient polynomial coefficient reconstruction')
    source,cp,rp=c.f.accepted_packet(c.f.M/'S11c_d_continuum_response_checkpoint.json','continuum-response.pickle')
    ends,bcp,bp=c.f.accepted_packet(c.f.M/'S11c_d_continuum_boundary_checkpoint.json','continuum-boundary.pickle')
    ref=next(Path(n) for n in source['inputPackets'] if Path(n).name=='modal.pickle');reference=c.f.unpickle(ref)[0]
    selection=c.selectors(source,ends,reference);metrics,metric_residuals=c.open_metrics(source,ends,selection)
    pencils={};input_hashes={str(p):c.f.digest(p) for p in (rp,bp,ref)}
    for end in ('LEFT','RIGHT'):
        pencils[end],_,p=c.f.accepted_packet(c.f.M/'S11c_d_continuum_boundary_checkpoint.json',end.lower()+'-pencil.pickle');input_hashes[str(p)]=c.f.digest(p)
    domains=c.bulk_domains(pencils)
    report={'scriptSha256':c.f.digest(Path(c.__file__)),'focusedSha256':c.f.digest(Path(__file__)),
        'independentCurrentResidual':c.boundary.norm(residual),'currentMutation':c.boundary.norm(mutation),
        'quotientDifference':qdiff,'quotientResidual':c.boundary.norm(qr),
        'actualOpenRanks':selection['ranks'],'actualMetricResidual':c.boundary.norm(metric_residuals),
        'actualCandidateRecords':sum(len(v['census']) for v in ends['ends'].values()),
        'actualSelectedClusters':sum(len(v) for v in selection['census'].values()),
        'actualSelectedDirections':sum(len(row['columns']) for rows in selection['census'].values() for row in rows),
        'bulkDepthSquared':{e:[str(v['depthSquared']) for v in rows] for e,rows in domains.items()},
        'bulkRealPropagationSets':{e:[str(v['realPropagationSet']) for v in rows] for e,rows in domains.items()},
        'inputPackets':input_hashes,'scope':'Arithmetic controls plus actual source/classifier and acoustic-domain joins; complete physical contractions remain in production.'}
    c.f.save(c.f.M/'S11c_d_continuum_currents_focused.json',report);print(json.dumps(report,indent=2))


if __name__=='__main__':main()

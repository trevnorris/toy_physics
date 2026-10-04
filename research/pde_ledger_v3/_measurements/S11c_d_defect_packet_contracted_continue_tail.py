"""Only the original unfinished geometry tail; no completed prefix calls."""
from fractions import Fraction
from S11c_d_defect_packet_contracted_continue_resume import require

def continue_geometry(raw,manifest,J,sp,C,G,N,entries,tails):
    get=raw.__getitem__
    old=lambda name:get("contraction/complete/"+name+".json")
    plans={};receipts={};kap=G.Quad(0,1,Fraction(119,20))
    for K,T in ((27,122),(29,124)):
        for carrier in (kap,kap*0):
            key=str(K)+'/'+str(carrier.b)
            label_alias='new-label-transport-'+('matching' if carrier==kap else 'zero')+'-K'+str(K)+'-square-full-arrangement'
            labels=old(label_alias)
            require(labels['originalPlan']==get('geometry/'+label_alias.removeprefix('new-label-transport-')+'.json'),'actual same window/carrier original labels')
            require({v['primitive'] for v in labels['transport']}=={'J','Dr','Dh','Dq'},'all four primitive label transports')
            for primitive in ('J','Dr','Dh','Dq'):
                J.emit('new-selected-label-transport-'+key.replace('/','-')+'-'+primitive,{'originalPlan':labels['originalPlan'],'transport':[v for v in labels['transport'] if v['primitive']==primitive],'sourceAlias':label_alias,'nativeMethod':'contracted variables; original external labels are retained as provenance, not an old quadrature request'})
            J.start('new-geometry-'+key.replace('/','-'),{'K':K,'T':T,'kappa':kap,'carrier':carrier,'width':8,'length':10,'inheritedLabels':labels})
            plan=G.mesh(kap,K,T,carrier,8,10);rec=J.emit('new-geometry-'+key.replace('/','-')+'-full',G.packed(plan));J.finish({'audit':True,'slabs':len(plan['slabs']),'fullWings':True})
            for side in (-1,1):J.emit('new-wing-control-'+key.replace('/','-')+'-'+str(side),G.packed(G.wing_control(plan,side)))
            plans[key]=plan;receipts[key]=rec
    J.emit('new-nested-allocation-lemma',{'epsilon':sp.Rational(1,80*10**11),'pointDensityTarget':'epsilon*(1+1/abs(q(m)))/(16*(3*M+4))','momentumNondimensionalizedInNativeUnit':True,'bound':'integral[-M,M]1/abs(q(m)) dm=pi+2 acosh(M/kappa)<4+M for kappa>2 and M>=149; hence integral w<3M+4','discreteGuardSeparate':'actual absolute K/G or GL weighted nested indicators <=epsilon/4','empiricalIndicatorsNotRigorous':True})
    return {'entries':entries,'geometryReceipts':receipts,'scope':manifest['scope'],'ruleReceipts':{r:manifest['savedInputs']['saved/rules/'+n+'.json'] for r,n in [('A24','A-GL24'),('A48','A-GL48'),('B50','B-G7-K15')]},'sourceManifestSha256':manifest['selfSourceIdentity'],'originalProofDependence':True,'inferredGammaUnits':True},plans,tails

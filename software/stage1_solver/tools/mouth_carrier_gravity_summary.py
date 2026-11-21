"""Collect the conditional gravity checks and plot only verified radius responses."""
import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

from mouth_flared_receiver import emit


BASE=Path('_scratch/s11c')


def rows(name,kind):
    folder=BASE/name
    for f in ('outcome','child-outcome'):
        if json.loads((folder/(f+'.json')).read_text())['exitCode']!=0:
            raise ValueError(f'Incomplete result: {folder}')
    if not json.loads((folder/'limit-validation.json').read_text())['verified']:
        raise ValueError(f'Unverified run: {folder}')
    return [r for line in (folder/'stdout').read_text().splitlines()
            if (r:=json.loads(line)).get('test')==kind]


def main(plot):
    choices={
        'depth scales with radius':['half_01','threequarters_02','one_02','threehalves_02','double_01'],
        'fixed depth':['fixedhalf_02','fixedthreequarters_02','one_02','fixedthreehalves_02','fixeddouble_02']}
    curves={}
    for name,names in choices.items():
        series=[]
        for folder in names:
            r=rows('carrier_radii_'+folder,'family_response')[-1]
            series.append((r['radius'],r['acceleration_per_incident_speed']))
        curves[name]=np.array(series)
    wide=np.array([(r['radius'],r['ratio']) for n in ('whalf','wone','wtwo','wfour')
                   for r in rows('carrier_action_'+n+'_01','matched_radius')[-1:]])
    weak=np.array([(r['radius'],r['weak_ratio_at_reference_bias']) for n in ('qfine','h','one','two','fourfine')
                   for r in rows('carrier_weak_wide_'+n+'_01','weak_geometry')])
    def relative(data):
        return data[:,1]/data[data[:,0]==1,1][0]
    emit(test='radii',finite_box14={k:v.tolist() for k,v in curves.items()},finite_box28=wide.tolist(),
         scaled_depth_max_over_min=np.ptp(curves['depth scales with radius'][:,1])/curves['depth scales with radius'][:,1].min()+1,
         fixed_depth_box14_max_over_min=curves['fixed depth'][:,1].max()/curves['fixed depth'][:,1].min(),
         fixed_depth_box28_R4_over_R1=relative(wide)[-1],weak_box28_R4_over_R1=relative(weak)[-1])
    free=rows('carrier_free_c34_01','free_response')
    prop=rows('carrier_inertia_c34_01','slow_acceleration')[0]['actuator_mass_volume']
    held=rows('carrier_inertia_c34_01','steady_balance')[0]['interior'][0]
    emit(test='independent_free',force_over_mass=held/prop,inverse_mass=1/prop,
         harmonic=[{'drive':r['drive'],'acceleration':r['acceleration']} for r in free if r['omega']==.01])
    pair=[r for r in rows('carrier_pair_energy_01','pair') if r['radii']==[1.,1.] and r['count']==28 and r['profile']=='uniform']
    for left,right in zip(pair[:-1],pair[1:]):
        emit(test='two_separations',distances=[left['d'],right['d']],
             force_at_smaller=left['force_toward'],force_at_larger=right['force_toward'],
             total_ratio=right['force_toward']/left['force_toward'],
             feeding_ratio=right['force_feed']/left['force_feed'])
    fig,(top,bottom)=plt.subplots(2,1,figsize=(7.6,6.5),sharex=True,layout='constrained',
                                  gridspec_kw={'height_ratios':[1.35,1]})
    colors=['#86532c','#155c82']
    for (label,data),color in zip(curves.items(),colors):
        top.plot(data[:,0],data[:,1],'-o',color=color,label=label)
    top.plot(wide[:,0],wide[:,1],'x--',color='#118269',label='fixed depth, wider reservoirs')
    top.set_yscale('log')
    top.set_ylabel('Acceleration per unit\nincident speed (model units)')
    top.legend(frameon=False,fontsize=9)
    top.set_title('What changes when the mouth radius changes?',loc='left',weight='bold')
    for data,label,style,color in ((curves['fixed depth'],'Finite bias, outer radius 14','-o',colors[1]),
                                   (wide,'Finite bias, outer radius 28','--x','#118269'),
                                   (weak,'Weak-bias intake / resting inertia, outer radius 28',':^','#744986')):
        bottom.plot(data[:,0],100*(relative(data)-1),style,label=label,color=color)
    bottom.axhline(0,color='.3',linewidth=.7)
    bottom.set_ylabel('Change relative to\nthe radius-one mouth (%)')
    bottom.set_xlabel('Mouth radius (units of the slab thickness)')
    bottom.legend(frameon=False,fontsize=8,loc='lower left')
    bottom.set_title('Fixed depth along the fourth dimension',loc='left',fontsize=11)
    for ax in (top,bottom):
        ax.set_xscale('log',base=2)
        ax.set_xticks([.25,.5,.75,1,1.5,2,4],['¼','½','¾','1','1½','2','4'])
        ax.grid(axis='y',alpha=.2)
        ax.spines[['top','right']].set_visible(False)
    fig.savefig(plot,dpi=170)
    plt.close(fig)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--plot',type=Path,required=True)
    main(p.parse_args().plot)

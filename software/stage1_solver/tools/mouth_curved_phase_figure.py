"""Draw the numerically solved stationary current and its crossing of the cap."""
import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

from mouth_curved_phase import Domain, quadratic, tensor


def draw(folder, output):
    outcome=json.loads((folder/'outcome.json').read_text())
    assert outcome['exitCode']==outcome['childOutcome']['exitCode']==0 and outcome['limitsVerified']
    state=json.loads((folder/'stdout').read_text().splitlines()[0])
    d=Domain(state['n'],state['far'])
    with np.load(folder/'field.npz') as f:
        np.testing.assert_array_equal(d.points,f['points'])
        psi,r,cap,q=(f[k] for k in ('psi','r','h_0.4','q_0.4'))
    rv=np.linspace(0,1.45,290); wv=np.linspace(-.48,.5,180)
    R,W=np.meshgrid(rv,wv); jr=np.full(R.shape,np.nan); jw=jr.copy()
    left=(R<=d.rc)&(np.hypot(d.rc-R,W-d.rounding)>=d.rounding)
    jr[left],jw[left]=d.current(psi,R[left],W[left])
    right=(R>d.rc)&(W<0)
    p,v=R[right],-W[right]; pg,vg,indices=d.blocks[1]
    i=np.searchsorted(pg,p)-1; j=np.searchsorted(vg,v)-1
    hp=np.diff(pg)[i]/2; hv=np.diff(vg)[j]/2
    ap,dp=quadratic((p-pg[i])/hp-1); av,dv=quadratic((v-vg[j])/hv-1)
    values=psi[indices[i*(len(vg)-1)+j]]
    def eval_shape(a,b):
        return (values*tensor(a,b)).sum(axis=1)
    f=eval_shape(ap,av)
    jr[right]=(f.conj()*eval_shape(dp,av)).imag/hp
    jw[right]=-(f.conj()*eval_shape(ap,dv)).imag/hv
    fig,(ax,bx)=plt.subplots(1,2,figsize=(9.0,3.65),gridspec_kw={'width_ratios':[1.45,1]},layout='constrained')
    blue='#116779'; dark='#303a40'; orange='#ad4a1f'
    wall_r=np.linspace(d.rc-d.rounding,d.rc,400)
    wall_z=d.rounding-np.sqrt(np.maximum(0.,d.rounding**2-(wall_r-d.rc)**2))
    ax.fill_between(np.r_[wall_r,rv[-1]],np.r_[wall_z,0],wv[-1],color='#e0e4e4',edgecolor='#879493',linewidth=.8)
    ax.streamplot(rv,wv,np.ma.masked_invalid(jr),np.ma.masked_invalid(jw),color=blue,density=.93,linewidth=1.0,arrowsize=1.1)
    ax.plot(r,cap,color=orange,lw=2.2,label='Conversion surface')
    ax.scatter([1],[d.lip],color=dark,s=18,zorder=5)
    ax.annotate('Lip',(1,d.lip),xytext=(1.12,.1),arrowprops={'arrowstyle':'-','color':dark},color=dark)
    ax.set(xlim=(0,rv[-1]),ylim=(wv[0],wv[-1]),xlabel='Radius along the brane  $r/R$',ylabel='Depth toward the bulk  $w/R$',title='Flow through the rounded mouth')
    ax.set_aspect('equal'); ax.legend(loc='lower left',fontsize=8,framealpha=.95)
    bx.plot(r,q,color=blue,lw=2.4);bx.axhline(0,color=dark,lw=.7)
    points=[0,.5,.9,.99];bx.scatter(points,np.interp(points,r,q),color=blue,s=18)
    bx.plot(1,0,marker='o',mfc='white',mec=orange,ms=6,clip_on=False)
    bx.annotate('Tangent at the lip',(1,0),xytext=(.57,.039),arrowprops={'arrowstyle':'-','color':orange},color=orange,fontsize=8)
    bx.set(xlim=(0,1.025),ylim=(-.006,.198),xlabel='Radius on the opening  $r/R$',ylabel='Inward number current per projected volume',title='Current crossing the surface')
    bx.grid(axis='y',alpha=.22);bx.spines[['right','top']].set_visible(False)
    fig.savefig(output,dpi=160);plt.close(fig)
    print('figure',output,'grid',state['n'],'far',state['far'],'rate',state['rate'],'samples',len(q))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('run',type=Path);parser.add_argument('output',type=Path)
    args=parser.parse_args();draw(args.run,args.output)

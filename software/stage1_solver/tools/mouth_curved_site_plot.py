"""Plot measured half-order displacement across the mouth from the guarded readout."""
import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def draw(run, output):
    outcome=json.loads((run/'outcome.json').read_text())
    assert outcome['exitCode']==outcome['childOutcome']['exitCode']==0 and outcome['limitsVerified']
    rows=[r for r in map(json.loads,(run/'stdout').read_text().splitlines())
          if r['test']=='front_shift' and r['site']=='half']
    fig,ax=plt.subplots(figsize=(6.7,3.8),layout='constrained')
    for branch,color,label,marker in [('converting','#185C70','Curved conversion layer','o'),
        ('uniform_restoring_curvature','#A75423','Same shape, uniform restoring stiffness','s')]:
        series=[r for r in rows if r['branch']==branch]
        assert len(series)==6
        ax.semilogy([r['radial_fraction'] for r in series],[r['mean_depth_shift_per_h'] for r in series],
                    color=color,marker=marker,label=label,linewidth=1.9,markersize=5)
    ax.set(xlabel='Fraction of local mouth radius',ylabel='Mean depth shift per unit site coupling',
           title='Response to one site at the wall',ylim=(8e-8,17))
    ax.set_xticks([0,.25,.5,.75,1],['0 (axis)','0.25','0.5','0.75','1 (wall)'])
    ax.grid(axis='y',alpha=.22)
    ax.legend(loc='lower right',frameon=False)
    fig.savefig(output,dpi=180)
    plt.close(fig)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    draw(args.run,args.output)

"""Format the measured CP10 force, work, depth and ranked-cell tables as Markdown."""
import argparse
import json
from pathlib import Path


def tables(record, unfolded=None):
    lines=['# Checkpoint 10 — signed spatial tables for radius 4','',
           'All differences are the spatial-label result minus the original material result. '
           'The full depth force includes the common field load. Stored energy and the order-shift component '
           'use the same frozen material quadrature points. Degree refers to the Gauss rule used in each direction.','']
    def heading(title):
        lines.extend(['## '+title,''])
    def table(names):
        lines.extend(['| '+' | '.join(names)+' |','|'+'|'.join(['---']+['---:']*(len(names)-1))+'|'])
    def row(*items):
        lines.append('| '+' | '.join(map(str,items))+' |')
    def signed(x):
        return f'{x:+.8f}'
    levels=record['levels']
    heading('Raw values')
    table(['Degree','Grid','Full force, material','Full force, spatial','Order work, material',
           'Order work, spatial','Crystal energy, material','Crystal energy, spatial'])
    for g in ('6','8'):
        d=record['degrees'][g]
        for i,n in enumerate(levels):
            s=d['totals'][i]; f=d['full_depth'][i]
            row(g,f'{n}×',signed(f['native_total']),signed(f['spatial_total']),signed(s[2]),signed(s[3]),
                f'{s[0]:.10f}',f'{s[1]:.10f}')
    lines.append('')
    heading('Spatial buckets')
    lines.extend(['A is the ordered interior outside both neighborhoods. B contains cells touching the '
        'physical weld and their immediate neighbors. C contains cells touching any folded prism or the '
        'boundary between active and fully melted cells, with one layer of neighbors. '
        'A neighbor shares a face, edge or vertex on the native spatial grid. '
        'Cells in both neighborhoods are assigned to C. The last table in this section reports the effect of assigning those cells to B.',''])
    for g in ('6','8'):
        d=record['degrees'][g]
        lines.extend([f'**Degree {g}**',''])
        table(['Grid','Signed quantity','A','B','C','Sum'])
        for i,n in enumerate(levels):
            for text,c in (('Stored energy',4),('Order work',5),('Full depth force',None)):
                vals=[b[i][c] for b in d['primary']] if c is not None else [b[0] for b in d['full_depth'][i]['primary']]
                row(f'{n}×',text,*map(signed,vals),signed(sum(vals)))
        lines.append('')
    lines.extend(['**Degree 8 with weld-first assignment in overlap cells**',''])
    table(['Grid','Quantity','A','B','C'])
    d=record['degrees']['8']
    for i,n in enumerate(levels):
        for text,c in (('Stored energy',4),('Order work',5),('Full depth force',None)):
            v=[b[i][c] for b in d['weld_priority']] if c is not None else [b[0] for b in d['full_depth'][i]['weld_priority']]
            row(f'{n}×',text,*map(signed,v))
    lines.append('')
    heading('Original material depth and weld side')
    lines.extend(['Each interval is a fraction of the saved material depth. Core and collar refer to the side '
                  'of the original bonded material point. Rows beginning at 2/24, 3/24 and 4/24 cover the requested band.',''])
    for g in ('6','8'):
        d=record['degrees'][g]
        for i,n in enumerate(levels):
            lines.extend([f'**Degree {g}, {n}× grid**',''])
            table(['Material interval','Core ΔE','Collar ΔE','Core ΔFχ','Collar ΔFχ'])
            values=[]
            for j,b in enumerate(d['depth_side']):
                v=[b[side][i][c] for c in (4,5) for side in (0,1)]
                values.append(v)
                row(f'{j}/24 – {j+1}/24',*map(signed,v))
            row('**2/24 – 5/24 subtotal**',*map(signed,map(sum,zip(*values[2:5]))))
            row('**All depths**',*map(signed,map(sum,zip(*values))))
            lines.append('')
    heading('Ten largest net cell contributions on the 4× grid')
    lines.extend(['The degree-8 cells are ranked by the magnitude of their net signed order work difference. '
                  'Local fields and physical position refer to the quadrature point with the largest absolute weighted '
                  'contribution in that cell. The two determinants are forward volume ratios. The polar angle is in radians.',''])
    table(['Rank','Fine cell (x,t,a)','Bucket','Cell ΔFχ','r','w','θ','χ','−∂χ/∂ln L'])
    for i,p in enumerate(d['top10'],1):
        row(i,'('+', '.join(map(str,p['fine_cell']))+')',p['bucket'],f"{p['delta_F']:+.10f}",
            *[f'{p[k]:.7g}' for k in ('r','w','theta','chi','minus_dchi')])
    lines.append('')
    table(['Rank','Material side','J material','J spatial','W material','W spatial'])
    for i,p in enumerate(d['top10'],1):
        row(i,p['material_side'],*[f'{p[k]:.8g}' for k in ('J_L','J_E','W_L','W_E')])
    lines.extend(['','Full-precision values, overlap controls, the degree-6 ranked cells and the historical '
                  'degree-3 energy check are in the [measurement record](../../../research/pde_ledger_v3/_measurements/mouth_r4_bucket_localization.json).',''])
    if unfolded is not None:
        heading('Unfolded radius-4 control')
        lines.extend(['The native material deformation was halved at the original nodes. Both evaluations '
                      'use the same held shape, order, carrier and stored core/collar metrics. Each entry '
                      'is a spatial-minus-material error on that control.',''])
        table(['Degree','Grid','Stored energy','Log depth','Scaled log exit','Mean width',
               'Depth asymmetry','Waist profile','Polar–equatorial'])
        for g in ('6','8'):
            for level,values in zip(levels,unfolded['results'][g]['spatial_minus_material']):
                row(g,f'{level}×',*map(signed,values))
        lines.extend(['','The energy column uses the frozen native material points. '
                      'The six forces use the same motion of each order knot in both representations. '
                      'The [unfolded measurement record](../../../research/pde_ledger_v3/_measurements/mouth_r4_bucket_unfolded.json) '
                      'also gives the original and spatial forces separately, the global orientation bounds, '
                      'and a near-static energy check from the signed depth quadrature.',''])
    return '\n'.join(lines)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('measurement',type=Path)
    p.add_argument('output',type=Path)
    p.add_argument('--unfolded',type=Path)
    args=p.parse_args()
    with args.output.open('x') as stream:
        stream.write(tables(json.loads(args.measurement.read_bytes()),
                            json.loads(args.unfolded.read_bytes()) if args.unfolded else None))

"""Render accepted states of the radius-four survival run from its pinned field maps."""
import argparse
import hashlib
import json
from pathlib import Path
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib import patheffects
import numpy as np
from PIL import Image
from scipy.interpolate import RegularGridInterpolator


INK = "#132b3c"
SUBTLE = "#415869"
BLUE = "#168c91"
RUST = "#b64b2b"
PAPER = "#faf9f5"
ORDER = LinearSegmentedColormap.from_list("shear_order", ["#0c2746", "#225276", "#398d9d", "#9ad5b6", "#f7e69a"])
DENSITY = plt.get_cmap("magma")


def read(path):
    return json.loads(path.read_bytes())


def sha(path):
    with path.open("rb") as f:
        return hashlib.file_digest(f, "sha256").hexdigest()


def states(run):
    pin = run / "pin/software/stage1_solver/tools"
    sys.path.insert(0, str(pin))
    from mouth_r4_fold_order import all_low_after, axes_model
    from mouth_sealed_shape_relax import Meridian
    from mouth_joined_angular_work import bump
    from mouth_wound_carrier import straight_vortex

    previous = Path(read(run/"coarse/invocation.json")["command"][-3]) / "coarse"
    original = Path(read(previous/"invocation.json")["command"][-5]) / "coarse"
    reference_depth = read(original/"endpoint.json")["depth"]
    f0 = straight_vortex(1, 28.)
    logs = [json.loads(s) for s in (run/"coarse/stdout").read_text().splitlines() if s.endswith("}")]
    accepted = {d["iteration"]:d for d in logs if d["kind"]=="accepted"}
    records = {d["iteration"]:d for d in logs if d["kind"]=="residuals"}
    baseline = read(previous/"latest.json")
    assert baseline["sha256"] == accepted[min(accepted)]["sha256"]
    records.setdefault(baseline["iteration"], baseline)

    series = []
    for i, saved in sorted(accepted.items()):
        path = run/"coarse"/saved["file"]
        if sha(path) != saved["sha256"]:
            raise ValueError(f"Changed checkpoint {path}")
        with np.load(path) as f:
            if int(f["radius"])!=4 or int(f["winding"])!=1:
                raise ValueError(f"Unexpected state {path}")
            shape = Meridian(4.,f["shape_params"],extra=float(f["shape_extra"]))
            shift = float(f["x"][5])
            fields = f["field"]
        if abs(saved["depth"]-shape.length)>1e-12 or abs(saved["exit_radius"]-shape.exit)>1e-12:
            raise ValueError("Recorded wall dimensions disagree with checkpoint")
        m = axes_model(shape,reference_depth)
        shaft = fields[2,m.ids[1,0]]
        front = all_low_after(m.z[1]/shape.length,shaft,0.)
        if front is not None:
            front *= shape.length
        current = records.get(i,{}).get("current",{})
        independently_logged = current.get("shaft_order",{}).get("terminal_from_fraction",{}).get("disk")
        if independently_logged and not np.isclose(front,independently_logged["0.0"]*shape.length,atol=1e-12,rtol=0):
            raise ValueError("Zero-order contour differs from logged full-shaft result")
        # z-axis of the ordinary brane is the winding axis. Its equatorial
        # slice contains x and y; w is perpendicular to both ordinary axes.
        x = np.linspace(-4.65,4.65,701)
        w = np.linspace(-.43,1.08,341)
        xx,ww = np.meshgrid(x,w)
        radius = np.abs(xx)
        b = shape.radius(ww)[0]
        delta = -.5*shift*bump(ww,shape.length)[0]
        inside = radius <= b+delta
        xi = radius / b
        for _ in range(5):
            xi -= (b*xi+delta*xi**3-radius)/(b+3*delta*xi**2)
        xi = np.where(inside,xi,(radius-b-delta)/(28.-b-delta))
        if np.max(abs(np.where(inside,b*xi+delta*xi**3,b+delta+(28.-b-delta)*xi)-radius))>1e-10:
            raise ValueError("Physical coordinate inversion failed")
        slab = np.where(ww<0,0,np.where(ww<=shape.length,1,2))
        section = np.empty((*xx.shape,2))
        for (layer,side), ids in m.ids.items():
            select = (slab==layer) & (inside==(side==0))
            ids = ids[:,-1,:]
            data = np.stack((fields[2,ids],fields[0,ids]),axis=-1)
            sampler = RegularGridInterpolator((m.rr[side],m.z[layer]),data)
            section[select] = sampler((np.clip(xi[select],0,1),ww[select]))
        order = section[:,:,0]
        density = (f0(radius)[0]*section[:,:,1])**2
        if np.min(order)<-1e-12 or np.max(order)>1+1e-12 or not np.isfinite(density).all():
            raise ValueError("Invalid section")
        if np.max(abs(density[:,len(x)//2]))>1e-12:
            raise ValueError("One-quantum density should vanish on the axis")
        wall_depth = np.linspace(0,shape.length,301)
        wall_radius = shape.radius(wall_depth)[0]-.5*shift*bump(wall_depth,shape.length)[0]
        if not np.allclose(np.interp(saved["depth_samples"],wall_depth,wall_radius),saved["equator_radius"],rtol=0,atol=2e-5):
            raise ValueError("Equatorial wall does not reproduce saved radii")
        series.append(dict(index=i,path=path,sha=saved["sha256"],shape=shape,shift=shift,wall=wall_radius,
                           w=wall_depth,x=x,z=w,order=order,density=density,front=front,current=current))
    return series


def geometry_view(ax,s):
    """Orthographic equatorial section. Projection preserves length without vertical exaggeration."""
    from matplotlib.collections import PolyCollection
    from matplotlib.colors import to_rgb
    az, elev = np.deg2rad([-66.,18.])
    def project(r,phi,w):
        across = r*np.sin(phi-az)
        along = r*np.cos(phi-az)
        return np.stack(np.broadcast_arrays(across, -np.sin(elev)*along-w*np.cos(elev),
                                            np.cos(elev)*along-w*np.sin(elev)),axis=-1)
    angles = np.deg2rad(np.linspace(-12,246,131))
    w = np.linspace(0,s['shape'].length,35)
    z,p = np.meshgrid(w,angles,indexing='ij')
    r = np.interp(z,s['w'],s['wall'])
    pos = project(r,p,z)
    quads = np.stack((pos[:-1,:-1],pos[:-1,1:],pos[1:,1:],pos[1:,:-1]),axis=2)
    # Shading distinguishes near and far walls. Only the lower two panels encode fields.
    light = .5+.5*np.cos((p[:-1,:-1]+p[:-1,1:])/2-np.deg2rad(120))
    dark, bright = np.array(to_rgb('#32626e')), np.array(to_rgb('#b9d8ce'))
    colors = dark+(bright-dark)*light[...,None]
    order = np.argsort(quads[:,:,:,2].mean(axis=2).ravel())
    polygons = quads.reshape(-1,4,3)[order,:,:2]
    rgb = colors.reshape(-1,3)[order]
    ax.add_collection(PolyCollection(polygons,facecolors=rgb,edgecolors='face',linewidths=.35))
    for level,rad in ((0.,s['shape'].neck),(s['shape'].length,s['shape'].exit)):
        rim = project(rad,angles,level)
        ax.plot(rim[:,0],rim[:,1],color=INK,lw=1.25)
    for angle in (angles[0],angles[-1]):
        edge=project(s['wall'],angle,s['w'])
        ax.plot(edge[:,0],edge[:,1],color=INK,lw=.95)
    ghost = project(s['shape'].neck,np.deg2rad(np.linspace(-114,-12,70)),0.)
    ax.plot(ghost[:,0],ghost[:,1],color='#789191',ls=(0,(3,3)),lw=.85)
    for label,rad,angle,level,at in (
            ('Entrance',s['shape'].neck,24,0.,(4.72,.62)),
            ('Exit',s['shape'].exit,215,s['shape'].length,(-4.72,-1.48))):
        peg=project(rad,np.deg2rad(angle),level)
        ax.annotate(label,peg[:2],at,color=SUBTLE,ha='center',va='center',fontsize=9,
                    arrowprops=dict(arrowstyle='-',color=SUBTLE,lw=.8,shrinkA=2,shrinkB=1))
    ax.set(xlim=(-5.65,5.65),ylim=(-2.07,1.59),aspect='equal')
    ax.set_axis_off()


def field_panel(ax,s,array,cmap,norm,kind):
    ax.pcolormesh(s["x"],s["z"],array,cmap=cmap,norm=norm,shading="nearest",rasterized=True)
    for sign in (-1,1):
        wall,=ax.plot(sign*s["wall"],s["w"],lw=1.35,color="#fafbf6")
        wall.set_path_effects([patheffects.Stroke(linewidth=2.7,foreground="#20313c"),patheffects.Normal()])
    if kind=="order" and s["front"] is not None:
        r=np.interp(s["front"],s["w"],s["wall"])
        ax.hlines(s["front"],-r,r,color="#e5c184",lw=1.4,ls=(0,(4,3)))
    ax.set(xlim=(s["x"][0],s["x"][-1]),ylim=(s["z"][-1],s["z"][0]),
           xticks=[-4,-2,0,2,4],yticks=[-.4,0,.4,.8,1.],xlabel="Distance across the mouth",ylabel="Depth into bulk,  w")
    ax.set_aspect(2.7, adjustable="box")
    ax.tick_params(length=3,color="#9faeb0",labelsize=9)
    for spine in ax.spines.values(): spine.set_color("#9faeb0")


def render(s,series,density_norm,output):
    from matplotlib.lines import Line2D
    fig=plt.figure(figsize=(12.4,8.3),facecolor=PAPER)
    fig.text(.06,.94,f'Radius-four mouth  /  checkpoint {s["index"]}',size=24,weight="bold",color=INK)
    fig.text(.06,.904,"This slice shows two brane directions and depth into the bulk (w). It passes through the center of the vortex.",size=11,color=SUBTLE)
    fig.text(.06,.839,"Mouth geometry",weight="bold",size=14,color=INK)
    fig.text(.06,.813,"Equatorial cutaway  ·  true proportions",size=10,color=SUBTLE)
    geom=fig.add_axes([.047,.570,.435,.231],facecolor=PAPER)
    geometry_view(geom,s)
    # The shallow shape is literal: place dimensional readouts beneath it.
    fig.text(.06,.543,f'Entrance radius  4.0000     Exit radius  {s["shape"].exit:.4f}     Depth  {s["shape"].length:.4f}',size=10,color=INK)
    fig.text(.555,.839,"Outer wall near the exit",weight="bold",size=14,color=INK)
    ax=fig.add_axes([.603,.595,.322,.186],facecolor=PAPER)
    minimum=min(t['shape'].exit for t in series)
    maximum=max(t['shape'].exit for t in series)
    for t in (series[0],s) if s is not series[0] else (s,):
        base=t is series[0] and s is not series[0]
        ax.plot(t['wall'],t['w'],lw=2.2,color=RUST if base else BLUE,ls='--' if base else '-',
                zorder=2 if base else 3,label=f'Checkpoint {t["index"]}')
        ax.plot(t['shape'].exit,t['shape'].length,marker='o',ms=5,mfc=PAPER if base else BLUE,
                mec=RUST if base else BLUE)
    ax.set(xlim=(minimum-.04,maximum+.24),ylim=(max(t['shape'].length for t in series)+.014,.67),xlabel='Radius',ylabel='Depth, w')
    ax.grid(True,lw=.6,alpha=.27);ax.spines[['top','right']].set_visible(False)
    ax.legend(frameon=False,fontsize=9,loc='upper left')
    ax.tick_params(labelsize=9)
    fig.text(.06,.474,"Crystalline order",weight="bold",size=14,color=INK)
    fig.text(.546,.474,"Coherent carrier density",weight="bold",size=14,color=INK)
    order_norm=Normalize(0,1)
    for left,data,cmap,norm,kind in ((.07,s['order'],ORDER,order_norm,'order'),(.558,s['density'],DENSITY,density_norm,'density')):
        panel=fig.add_axes([left,.180,.39,.265])
        field_panel(panel,s,data,cmap,norm,kind)
        cax=fig.add_axes([left+.055,.120,.27,.012])
        sm=matplotlib.cm.ScalarMappable(norm=norm,cmap=cmap)
        cb=fig.colorbar(sm,cax=cax,orientation="horizontal",ticks=[0, .25, .5, .75, 1.] if kind=='order' else None)
        cb.ax.tick_params(labelsize=8,length=2);cb.outline.set_visible(False)
        cb.set_label('Order χ   (0 = melted,  1 = ordered)' if kind=='order' else 'Density n   (undisturbed fluid = 1)',size=9)
    fig.add_artist(Line2D([.064,.088],[.046,.046],transform=fig.transFigure,color="#987441",lw=1.6,ls=(0,(4,3))))
    fig.text(.096,.040,f'Entire shaft melted below w = {s["front"]:.4f}' if s['front'] is not None else 'A fully melted cross-section does not persist through the exit.',size=9,color=SUBTLE)
    fig.text(.546,.040,'Sections stretch depth 2.7×. White lines trace the modeled wall.',size=9,color=SUBTLE)
    fig.savefig(output,dpi=156,facecolor=PAPER)
    plt.close(fig)


def main(run,destination):
    destination.mkdir(parents=True,exist_ok=True)
    series=states(run)
    max_n=max(float(s['density'].max()) for s in series)
    density_norm=Normalize(0,max(1.,np.ceil(max_n*10)/10))
    plt.rcParams.update({'font.family':'DejaVu Sans','text.color':INK,'axes.labelcolor':SUBTLE,
                         'xtick.color':SUBTLE,'ytick.color':SUBTLE,'axes.edgecolor':SUBTLE,
                         'font.size':10})
    frames=[]
    for s in series:
        path=destination/f'checkpoint_{s["index"]:02d}.png'
        render(s,series,density_norm,path)
        with Image.open(path) as im:
            frames.append(im.convert('RGB').resize((1116,747),Image.Resampling.LANCZOS))
        print(s['index'],s['sha'],s['shape'].length,s['shape'].exit,s['front'],s['density'].min(),s['density'].max(),flush=True)
    if len(frames)>1:
        # A shared palette keeps the same field colors from flickering.
        palette=Image.fromarray(np.concatenate([np.asarray(im.resize((372,249))) for im in frames],axis=0)).quantize(colors=256)
        small=[im.quantize(palette=palette,dither=Image.Dither.NONE) for im in frames]
        small[0].save(destination/'checkpoint_sequence.gif',save_all=True,append_images=small[1:],
                      duration=1800,loop=0,optimize=False,disposal=2)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('run',type=Path)
    p.add_argument('destination',type=Path)
    args=p.parse_args()
    main(args.run.resolve(),args.destination.resolve())

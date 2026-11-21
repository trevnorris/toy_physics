"""Plot the tested 4D throat sections and the axial current required by continuity."""
import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def read_run(folder):
    outcome = json.loads((folder/"outcome.json").read_text())
    assert outcome["limitsVerified"] and outcome["exitCode"] == outcome["childOutcome"]["exitCode"] == 0
    rows = {row["test"]: row for row in map(json.loads, (folder/"stdout").read_text().splitlines())}
    with np.load(folder/"field.npz") as fields:
        points, r, h = fields["points"], fields["r"], fields["h"]
    s = rows["carrier"]
    inside = (points[:,1]>s["neck_w"]) & (points[:,1]<s["exit_w"])
    z, index = np.unique(points[inside,1], return_inverse=True)
    a = np.zeros(len(z))
    np.maximum.at(a, index, points[inside,0])
    arc_z = np.linspace(h[-1], s["neck_w"], 100)
    rc = s["neck"] + s["neck_w"]
    arc_r = rc-np.sqrt(np.maximum(0, s["neck_w"]**2-(arc_z-s["neck_w"])**2))
    return s, rows["sections"]["values"], np.r_[arc_z,z,s["exit_w"]], np.r_[arc_r,a,s["exit_radius"]], r, h


def draw(first, second, dest):
    fig, (wall, flow) = plt.subplots(1, 2, figsize=(9.4,4.8), layout="constrained")
    colors = ["#006b69", "#a14d28"]
    for data, color, label, style in zip((first,second), colors, ("Earlier flare","Delayed flare"), ("-","--"), strict=True):
        s, cuts, depth, radius, r, h = data
        for sign in (-1,1):
            wall.plot(sign*radius,depth,color=color,ls=style,lw=2,label=label if sign==1 else None)
        area = 4*np.pi*radius**3/3
        flow.plot(depth, s["feed"]/area, color=color, ls=style, lw=2, label=label)
        cut_z = [cut["w"] for cut in cuts]
        cut_avg = [cut["integrated_jw"]/(4*np.pi*cut["radius"]**3/3) for cut in cuts]
        flow.scatter(cut_z,cut_avg,s=21,edgecolors=color,facecolors="white",zorder=3)
    s, _, _, _, r, h = first
    wall.plot(np.r_[-r[::-1],r[1:]],np.r_[h[::-1],h[1:]],color="#343a43",lw=1.6)
    wall.annotate("opening", (0,h[0]),xytext=(0,-15),textcoords="offset points",ha="center",fontsize=9)
    wall.text(0,s["neck_w"],"waist",ha="center",va="center",fontsize=9)
    wall.text(0,s["exit_w"]+.07,"into the bulk",ha="center",va="bottom",fontsize=9)
    wall.set(xlabel="Distance across a brane diameter (mouth radii)",ylabel="Depth along $w$ (mouth radii)",
             title="The two receiver walls",xlim=(-1.72,1.72),ylim=(-.29,3.32),aspect="equal")
    flow.set(xlabel="Depth along $w$ (mouth radii)",ylabel="Mean inward current per 3-volume",
             title="Current through an axial section")
    flow.legend(loc="upper right", frameon=False)
    flow.text(.98,.55,"Lines from total intake\nPoints from the solved field",
              ha="right",transform=flow.transAxes,fontsize=9,color="#4b535a")
    flow.set_ylim(bottom=0)
    for ax in (wall,flow):
        ax.spines[["top","right"]].set_visible(False)
        ax.set_axisbelow(True)
        ax.grid(alpha=.16)
    fig.savefig(dest,dpi=180)
    plt.close(fig)
    for s, _, _, _, _, _ in (first,second):
        print(json.dumps(dict(profile=s["profile"], neck_radius=s["neck"],
             mean_axial_at_neck=s["feed"]/(4*np.pi*s["neck"]**3/3),
             mean_axial_at_bore=s["feed"]/(4*np.pi*s["bore"]**3/3),
             mean_axial_at_exit=s["feed"]/(4*np.pi*s["exit_radius"]**3/3),
             waist_to_bore_current=(s["bore"]/s["neck"])**3)))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("cubic",type=Path)
    parser.add_argument("delayed",type=Path)
    parser.add_argument("output",type=Path)
    args = parser.parse_args()
    draw(read_run(args.cubic),read_run(args.delayed),args.output)

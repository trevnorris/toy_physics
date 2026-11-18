"""Plot signed displacement and the two momentum readings of the continued benchmark."""

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("fields",type=Path)
    parser.add_argument("figure",type=Path)
    args = parser.parse_args()
    fig,(motion,budget) = plt.subplots(1,2,figsize=(10,3.9),layout="constrained")
    colors = ("#245383","#b25422")
    for sign,style in (("positive","-"),("negative","--")):
        with np.load(args.fields/f"benchmark_continued_128_{sign}.npz") as fields:
            for i,field in enumerate(("Wall","Wave envelope")):
                label = f"{field}, {'+' if sign=='positive' else '−'} drive"
                motion.plot(fields["z"],fields["displacement"][i,:,-1],style,
                            color=colors[i],label=label)
                if sign=="positive":
                    budget.plot(fields["times"],fields["stored"][i],color=colors[i],label=f"{field}, stored")
                    budget.plot(fields["times"],fields["cut_impulse"][i],":",lw=2.5,color=colors[i],
                                label=f"{field}, crossed")
            cut = float(fields["cut_depth"])
    motion.axvline(cut,color="#4a5159",lw=1)
    motion.text(cut,0.04,"Cut ",transform=motion.get_xaxis_transform(),va="bottom",ha="right")
    motion.set(xlabel="Depth from rim",ylabel="Sideways displacement",title="Motion at time 1")
    budget.set(xlabel="Time",ylabel="Sideways momentum",title="Positive drive")
    for ax in (motion,budget):
        ax.ticklabel_format(axis="y",style="sci",scilimits=(-2,-2),useMathText=True)
        ax.spines[["top","right"]].set_visible(False)
        ax.grid(axis="y",alpha=0.22)
        ax.legend(fontsize=8,frameon=False)
    fig.suptitle("Wall and wave in the test segment · benchmark inputs, solver units")
    fig.savefig(args.figure,dpi=180)
    plt.close(fig)


if __name__ == "__main__":
    main()

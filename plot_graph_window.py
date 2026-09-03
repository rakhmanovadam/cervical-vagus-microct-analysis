"""
Draw a readable window of the fascicle tracking graph.

The graph is not a network blob: it is a few thousand long chains running down
the nerve in parallel, splitting and merging. So position each node by the data
rather than by a generic layout algorithm - slice number down the vertical axis,
real physical position across the horizontal - and show only a short window.

    python plot_graph_window.py                 # default window
    python plot_graph_window.py --start 27 --mm 3
"""
import argparse, glob, os
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
MM_PER_PX = 0.0114

ap = argparse.ArgumentParser()
ap.add_argument("--start", type=float, default=44.0, help="window start, mm along nerve")
ap.add_argument("--mm", type=float, default=4.0, help="window length in mm")
a = ap.parse_args()

nodes = pd.read_csv(os.path.join(HERE, "graph_nodes.csv"))
f0 = nodes["frame"].min()
nodes["mm"] = (nodes["frame"] - f0) * MM_PER_PX

# Rebuild the edges. Only edges whose BOTH ends fall in the window are drawn.
import xml.etree.ElementTree as ET
ns = {"g": "http://graphml.graphdrawing.org/xmlns"}
root = ET.parse(sorted(glob.glob(os.path.join(HERE, "*RawFascicleTracking*.graphml")))[0]).getroot()
E = [(e.get("source"), e.get("target")) for e in root.find("g:graph", ns).findall("g:edge", ns)]

win = nodes[(nodes["mm"] >= a.start) & (nodes["mm"] < a.start + a.mm)]
pos = win.set_index("node_id")[["ellipse_center-0", "mm"]]
keep = set(pos.index)

fig, ax = plt.subplots(figsize=(10, 6.5))
for s, t in E:
    if s in keep and t in keep:
        x0, y0 = pos.loc[s]; x1, y1 = pos.loc[t]
        ax.plot([x0 * MM_PER_PX, x1 * MM_PER_PX], [y0, y1],
                color="0.55", lw=0.5, zorder=1)

# Colour by role so the branch points are the thing you notice.
plain = win[~win["is_split"] & ~win["is_merge"]]
ax.scatter(plain["ellipse_center-0"] * MM_PER_PX, plain["mm"],
           s=3, color="#4878a8", zorder=2, label="fascicle, continuing")
sp = win[win["is_split"]]
ax.scatter(sp["ellipse_center-0"] * MM_PER_PX, sp["mm"],
           s=34, color="#d1495b", zorder=3, label="splits below this point")
mg = win[win["is_merge"] & ~win["is_split"]]
ax.scatter(mg["ellipse_center-0"] * MM_PER_PX, mg["mm"],
           s=34, color="#2a9d8f", marker="^", zorder=3, label="merges into this point")

ax.invert_yaxis()          # superior at the top, matching the nerve
ax.set_xlabel("position across the nerve (mm)")
ax.set_ylabel("distance along nerve, superior to inferior (mm)")
ax.set_title("Fascicle tracking graph, %.1f-%.1f mm\n%d nodes, %d branch points"
             % (a.start, a.start + a.mm, len(win),
                int(win["is_split"].sum() + win["is_merge"].sum())))
ax.legend(loc="upper right", fontsize=8, framealpha=0.9)
fig.tight_layout()
out = os.path.join(HERE, "fig4_graph_window.png")
fig.savefig(out, dpi=150)
print("wrote fig4_graph_window.png: %d nodes, %d splits, %d merges"
      % (len(win), int(win["is_split"].sum()), int(win["is_merge"].sum())))

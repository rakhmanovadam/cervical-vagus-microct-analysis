"""
Flythrough animation of the fascicular plexus in one human cervical vagus nerve.

DATA SOURCE
    Reconstructing Human Vagal Anatomy (REVA), Case Western Reserve University and
    Duke University, published on the NIH SPARC portal via Pennsieve.
    DOI 10.26275/rqkx-w7yx, CC-BY-4.0. Subject SR042, left cervical trunk, CL1.

WHAT IT DRAWS
    Travelling from the top of the neck downwards, one frame per 5 slices. Each
    fascicle is drawn as its own measured ellipse, at its own measured position
    and size, coloured by its track_id so a single fascicle keeps one colour for
    as long as it exists. Fascicles splitting or merging on that slice are
    outlined in black. The dashed outline is the segmented nerve boundary, which
    is missing on 62% of slices and is labelled as such when absent.

WHY
    The static plots reduce each cross-section to a single number, which discards
    where the fascicles actually sit. Spatial arrangement is the thing the REVA
    project is trying to establish, so it is worth keeping.

HOW TO RUN
    python make_animation.py            build the cache if needed, then render
    python make_animation.py --cache    rebuild the cache only

    Needs ffmpeg on the PATH (brew install ffmpeg / apt install ffmpeg).
    Writes fascicle_plexus.mp4, and graph_nodes.csv as a cache.
"""

import argparse
import glob
import os
import sys
import xml.etree.ElementTree as ET

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")          # render to file, never try to open a window
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse                  # the fascicle shapes
from matplotlib.animation import FuncAnimation, FFMpegWriter

HERE = os.path.dirname(os.path.abspath(__file__))
MM_PER_PX = 0.0114             # 11.4 micrometres per pixel, as millimetres
CACHE = os.path.join(HERE, "graph_nodes.csv")
OUT = os.path.join(HERE, "fascicle_plexus.mp4")

STEP = 5        # render every 5th slice: 5538 slices becomes 1108 frames
FPS = 30        # frames per second, so 1108 frames is about 37 seconds
DPI = 130       # output resolution

# The nerve's apparent caliber varies from about 1.7 mm to 6.4 mm along this
# segment. The wide end is the jugular foramen region, where the segmented
# boundary encloses the ganglion and surrounding structures rather than the trunk
# alone. A fixed field of view would either waste most of the frame or clip that
# region, so the camera follows the nerve and prints the current scale.
ZOOM_PAD = 1.35                        # field of view as a multiple of nerve width
ZOOM_MIN_MM, ZOOM_MAX_MM = 1.8, 9.0    # limits on the half-width, in mm
SMOOTH = 51                            # slices to average over, so the camera
                                       # glides instead of jittering frame to frame


def find_one(pattern):
    """Return the single file in HERE matching a wildcard pattern."""
    hits = sorted(glob.glob(os.path.join(HERE, pattern)))
    if not hits:
        sys.exit("no file matching %r in %s" % (pattern, HERE))
    return hits[0]


def build_cache():
    """Read the 45 MB GraphML once and flatten it to a small CSV.

    Parsing that much XML takes minutes, so doing it inside the render loop would
    be wasteful. Every node already carries its own geometry (position, axes,
    angle) plus a track_id, so the morphology CSVs are not needed for drawing.
    """
    path = find_one("*RawFascicleTracking*.graphml")
    print("parsing %s ..." % os.path.basename(path))

    # Every tag in this file lives in the GraphML XML namespace, so searches have
    # to be namespace-aware. This dict maps the short prefix "g" onto that URL.
    ns = {"g": "http://graphml.graphdrawing.org/xmlns"}
    root = ET.parse(path).getroot()

    # The file declares attribute names once, then refers to them by short ids
    # like "v_area". Build {id: readable name} so the values can be labelled.
    keys = {k.get("id"): k.get("attr.name") for k in root.findall("g:key", ns)}
    graph = root.find("g:graph", ns)

    # --- nodes: one per fascicle per slice -------------------------------
    rows = []
    for n in graph.findall("g:node", ns):
        # Dict comprehension: for each <data> child, map its readable attribute
        # name to its text value.
        d = {keys.get(x.get("key")): x.text for x in n.findall("g:data", ns)}
        d["node_id"] = n.get("id")
        rows.append(d)
    nodes = pd.DataFrame(rows)

    def truthy(s):
        """GraphML writes booleans as the text "true"/"false". In Python,
        bool("false") is True, so the raw strings cannot be tested directly."""
        return s.astype(str).str.strip().str.lower().isin(["true", "1", "yes"])

    # --- edges: one per slice-to-slice link ------------------------------
    e_src, e_tgt, e_split, e_merge = [], [], [], []
    for e in graph.findall("g:edge", ns):
        d = {keys.get(x.get("key")): x.text for x in e.findall("g:data", ns)}
        e_src.append(e.get("source"))     # the slice-above end of the arrow
        e_tgt.append(e.get("target"))     # the slice-below end
        e_split.append(d.get("split"))
        e_merge.append(d.get("merge"))
    edges = pd.DataFrame({"source": e_src, "target": e_tgt,
                          "split": e_split, "merge": e_merge})

    # A split belongs to the node the arrows LEAVE; a merge to the node they
    # ARRIVE at. Flagging the node rather than the edge is what stops a
    # three-way split from being counted as three separate events.
    split_nodes = set(edges.loc[truthy(edges["split"]), "source"])
    merge_nodes = set(edges.loc[truthy(edges["merge"]), "target"])
    # .isin(set) tests membership for every row at once.
    nodes["is_split"] = nodes["node_id"].isin(split_nodes)
    nodes["is_merge"] = nodes["node_id"].isin(merge_nodes)

    # XML gives everything as text, so convert the numeric columns.
    # errors="coerce" turns anything unparseable into NaN rather than raising.
    num = ["area", "equivalent_diameter", "ellipse_center-0", "ellipse_center-1",
           "ellipse_major_axis", "ellipse_minor_axis", "ellipse_angle",
           "frame", "track_id"]
    for c in num:
        nodes[c] = pd.to_numeric(nodes[c], errors="coerce")

    nodes = nodes[["node_id", "frame", "track_id", "area", "equivalent_diameter",
                   "ellipse_center-0", "ellipse_center-1", "ellipse_major_axis",
                   "ellipse_minor_axis", "ellipse_angle", "is_split", "is_merge"]]
    nodes.to_csv(CACHE, index=False)     # index=False omits pandas' row numbers
    print("wrote %s: %d nodes, %d frames, %d tracks, %d split, %d merge"
          % (os.path.basename(CACHE), len(nodes), nodes["frame"].nunique(),
             nodes["track_id"].nunique(), nodes["is_split"].sum(),
             nodes["is_merge"].sum()))
    return nodes


def load_landmarks(offset_mm):
    """The five anatomical landmarks, shifted into this segment's coordinates."""
    lev = pd.read_csv(find_one("*vagal_levels*.csv")).dropna(subset=["axis-0"]).copy()
    lev["dist_mm"] = lev["axis-0"] * MM_PER_PX - offset_mm
    lev["short"] = (lev["name"].str.replace("left level of ", "", regex=False)
                    .str.replace("the ", "", regex=False))
    return lev.sort_values("dist_mm")


def track_color(track_ids):
    """Give each fascicle lineage a stable colour.

    A fascicle has to keep the same colour for its whole life or the eye cannot
    follow it. Hashing the track_id into a fixed palette gives a colour that
    depends only on the id, so it does not matter which frame renders first.
    """
    # np.vstack stacks three 20-colour palettes into one 60-colour array.
    base = np.vstack([plt.get_cmap("tab20")(np.linspace(0, 1, 20)),
                      plt.get_cmap("tab20b")(np.linspace(0, 1, 20)),
                      plt.get_cmap("tab20c")(np.linspace(0, 1, 20))])
    # 2654435761 is Knuth's multiplicative hash constant. Multiplying by it and
    # taking the remainder scatters consecutive ids across the palette, so
    # neighbouring fascicles rarely end up the same colour.
    idx = (track_ids.astype(np.int64) * 2654435761) % len(base)
    return base[idx]                     # fancy indexing: one colour per id


def load_nerve_outline():
    """The nerve's outer boundary per slice, drawn behind the fascicles.

    The nerve file uses center_x / center_y while the graph uses
    ellipse_center-0 / -1, and nothing states how they correspond. Both mappings
    were tested by asking which one places the fascicles inside the nerve:
    center_x <-> ellipse_center-0 contains 99.2% of fascicle centres, the swap
    contains 6.7%. That test also confirms width=major, height=minor, angle=angle
    for matplotlib's Ellipse.
    """
    n = pd.read_csv(find_one("*NerveMorph*.csv"))
    n = n[n["segment"] == "CL1"]
    # Rows where the boundary could not be traced are dropped here, but their
    # absence is announced on screen rather than passed over silently.
    return n.dropna(subset=["center_x", "center_y", "major_axis"]).set_index("index")


def render(nodes):
    """Draw every frame and encode them into an MP4."""
    nodes = nodes.dropna(subset=["frame", "ellipse_center-0", "ellipse_center-1"])
    f0 = nodes["frame"].min()            # first slice number in this segment
    offset_mm = f0 * MM_PER_PX           # distance from the top of the whole scan
    nodes = nodes.assign(dist_mm=(nodes["frame"] - f0) * MM_PER_PX)

    lev = load_landmarks(offset_mm)
    nerve = load_nerve_outline()

    # Group the nodes by slice once, up front. Looking up a prepared dict per
    # frame is far faster than filtering the whole table 1108 times.
    by_frame = {f: g for f, g in nodes.groupby("frame")}
    frames = sorted(by_frame)[::STEP]    # [::STEP] takes every STEPth item
    counts = nodes.groupby("dist_mm").size().sort_index()

    # --- camera path -----------------------------------------------------
    # Follow the nerve's own centre and width. Where the nerve row is blank, fall
    # back to the average position of the fascicles on that slice.
    cam = pd.DataFrame(index=sorted(by_frame))
    cam["cx"] = nerve["center_x"].reindex(cam.index)
    cam["cy"] = nerve["center_y"].reindex(cam.index)
    # This file's "major_axis" is not reliably the larger of the two: on slice
    # 2097 major is 132 px and minor is 333 px. Take the max explicitly.
    cam["span"] = nerve[["major_axis", "minor_axis"]].max(axis=1).reindex(cam.index)
    fallback = nodes.groupby("frame")[["ellipse_center-0", "ellipse_center-1"]].mean()
    cam["cx"] = cam["cx"].fillna(fallback["ellipse_center-0"])
    cam["cy"] = cam["cy"].fillna(fallback["ellipse_center-1"])
    # ffill then bfill: carry the last known width forward, then fill any
    # remaining gap at the very start by reaching backwards.
    cam["span"] = cam["span"].ffill().bfill()
    # .rolling(N, center=True).mean() replaces each value with the average of the
    # N values centred on it, which smooths the camera motion.
    cam = cam.rolling(SMOOTH, center=True, min_periods=1).mean()
    # np.clip keeps the zoom inside sensible limits.
    cam["half"] = np.clip(cam["span"] * MM_PER_PX * ZOOM_PAD / 2,
                          ZOOM_MIN_MM, ZOOM_MAX_MM)

    # --- figure layout ---------------------------------------------------
    fig = plt.figure(figsize=(7.2, 8.4))
    # A 2-row grid: the cross-section view on top, three times as tall as the
    # position tracker underneath.
    gs = fig.add_gridspec(2, 1, height_ratios=[3.0, 1.0], hspace=0.30,
                          left=0.11, right=0.97, top=0.88, bottom=0.09)
    ax = fig.add_subplot(gs[0])          # the nerve cross-section
    axc = fig.add_subplot(gs[1])         # the count-vs-distance tracker

    ax.set_aspect("equal")               # circles look like circles
    ax.set_xlabel("mm")
    ax.set_ylabel("mm")

    axc.plot(counts.index, counts.values, lw=0.5, color="0.6")
    axc.set_xlim(counts.index.min(), counts.index.max())
    axc.set_ylim(0, counts.max() * 1.30)
    axc.set_xlabel("distance along nerve, superior to inferior (mm)")
    axc.set_ylabel("fascicles")
    for _, r in lev.iterrows():
        if counts.index.min() <= r["dist_mm"] <= counts.index.max():
            axc.axvline(r["dist_mm"], color="0.5", lw=0.8, ls=":")
            axc.text(r["dist_mm"], axc.get_ylim()[1], " " + r["short"],
                     rotation=90, va="top", ha="left", fontsize=6, color="0.35")
    marker = axc.axvline(0, color="#d62728", lw=1.4)   # the moving position line

    fig.suptitle("Human cervical vagus, fascicular plexus\n"
                 "REVA SR042 left cervical trunk, microCT 11.4 µm",
                 fontsize=11, y=0.975)
    banner = fig.text(0.54, 0.912, "", ha="center", va="center", fontsize=9.5)
    scale = fig.text(0.955, 0.862, "", ha="right", va="center", fontsize=8,
                     color="0.45")

    # Shapes drawn this frame, kept so they can be removed before the next one.
    patches = []

    def draw(f):
        """Draw one frame. FuncAnimation calls this once per slice number."""
        # Erase last frame. Without this, every frame's shapes accumulate and
        # the video ends as a smear. This is the single most common mistake in
        # matplotlib animations.
        for p in patches:
            p.remove()
        patches.clear()

        g = by_frame[f]
        cx, cy, half = cam.loc[f, "cx"], cam.loc[f, "cy"], cam.loc[f, "half"]
        ax.set_xlim(-half, half)
        ax.set_ylim(-half, half)

        # The nerve boundary, if it was traceable on this slice.
        if f in nerve.index:
            r = nerve.loc[f]
            outline = Ellipse(((r["center_x"] - cx) * MM_PER_PX,
                               (r["center_y"] - cy) * MM_PER_PX),
                              width=r["major_axis"] * MM_PER_PX,
                              height=r["minor_axis"] * MM_PER_PX,
                              angle=r["angle"], facecolor="none",
                              edgecolor="0.45", lw=1.0, ls="--", zorder=0)
            ax.add_patch(outline)
            patches.append(outline)

        # The fascicles. Positions are converted to millimetres relative to the
        # camera centre, so the drawing stays put while the nerve wanders.
        cols = track_color(g["track_id"].fillna(0).values)
        for (_, r), col in zip(g.iterrows(), cols):
            event = bool(r["is_split"]) or bool(r["is_merge"])
            e = Ellipse(((r["ellipse_center-0"] - cx) * MM_PER_PX,
                         (r["ellipse_center-1"] - cy) * MM_PER_PX),
                        width=r["ellipse_major_axis"] * MM_PER_PX,
                        height=r["ellipse_minor_axis"] * MM_PER_PX,
                        angle=r["ellipse_angle"],
                        facecolor=col, alpha=0.85,
                        edgecolor="k" if event else "none",   # outline branchings
                        lw=1.8 if event else 0, zorder=2)
            ax.add_patch(e)
            patches.append(e)

        # Update the readouts.
        d = (f - f0) * MM_PER_PX
        marker.set_xdata([d, d])
        below = lev[lev["dist_mm"] <= d]
        where = below["short"].iloc[-1] if len(below) else "above first landmark"
        n_ev = int(g["is_split"].sum() + g["is_merge"].sum())
        banner.set_text("%.2f mm along nerve   |   %d fascicles   |   below %s%s"
                        % (d, len(g), where, "   |   branching" if n_ev else ""))
        # 62% of slices have no traceable nerve boundary. Say so, rather than
        # letting the dashed outline silently vanish and look like an absence.
        scale.set_text("field %.1f mm%s"
                       % (2 * half,
                          "" if f in nerve.index else
                          "   ·   nerve boundary not traced"))
        return patches + [marker, banner, scale]

    # FuncAnimation calls draw() once per entry in `frames`. blit=False redraws
    # the whole figure each time, which is slower but avoids artefacts from the
    # changing axis limits.
    anim = FuncAnimation(fig, draw, frames=frames, interval=1000 / FPS, blit=False)
    writer = FFMpegWriter(fps=FPS, bitrate=2600,
                          metadata={"title": "REVA SR042 cervical vagus plexus"})
    print("rendering %d frames -> %s" % (len(frames), os.path.basename(OUT)))
    anim.save(OUT, writer=writer, dpi=DPI)
    print("done: %.1f MB, %.0f s at %d fps"
          % (os.path.getsize(OUT) / 1e6, len(frames) / FPS, FPS))


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--cache", action="store_true",
                    help="rebuild the GraphML cache and stop")
    args = ap.parse_args()

    # Rebuild the cache if asked, or if it does not exist yet.
    if args.cache or not os.path.exists(CACHE):
        nodes = build_cache()
        if args.cache:
            sys.exit(0)
    else:
        nodes = pd.read_csv(CACHE)
    render(nodes)

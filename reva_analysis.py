"""
Fascicle morphometry of one human cervical vagus nerve, from public REVA microCT data.

DATA SOURCE
    Reconstructing Human Vagal Anatomy (REVA), Case Western Reserve University and
    Duke University. Published on the NIH SPARC portal via Pennsieve.
    DOI 10.26275/rqkx-w7yx, licensed CC-BY-4.0.
    Subject SR042, left cervical trunk, segment CL1 only.

WHAT THIS SCRIPT DOES
    Reads the four released tables and reproduces basic fascicle morphometry from
    them. It does not touch the microCT images. REVA's team segmented and tracked
    the fascicles; this script consumes their measurements.

HOW TO RUN
    python reva_analysis.py --stage 1     (then 2, 3, 4, 5)

    Stage 1  inspect     print structure and missingness of every file
    Stage 2  count       fascicles per cross-section along the nerve -> fig1
    Stage 3  diameter    verify the diameter formula, then measure sizes -> fig2
    Stage 4  events      count split/merge events -> fig3
    Stage 5  levels      map anatomical landmarks onto the segment

    Stages 1, 3 and 5 are checks. Stages 2 and 4 are measurements. The checks come
    first on purpose: if the units or the coordinate system are misunderstood, every
    measurement downstream is wrong in a way that still looks plausible.
"""

# --- standard library -------------------------------------------------------
import argparse   # parses the --stage flag from the command line
import glob       # finds files by wildcard pattern, e.g. "*FasMorph*.csv"
import os         # builds file paths that work on any operating system
import sys        # used here only for sys.exit() when a file is missing
import xml.etree.ElementTree as ET   # reads the GraphML file, which is XML

# --- third-party libraries --------------------------------------------------
import numpy as np    # fast maths on whole arrays of numbers at once
import pandas as pd   # spreadsheet-like tables ("DataFrames") in Python
import matplotlib

# matplotlib normally tries to open a GUI window to draw into. "Agg" switches it
# to a non-interactive backend that writes straight to a PNG file. This line must
# come BEFORE importing pyplot, or the choice is ignored.
matplotlib.use("Agg")
import matplotlib.pyplot as plt


# ============================================================================
# CONFIGURATION
# Every number that could be wrong lives here, in one place, so a reader can
# check the assumptions without reading the whole script.
# ============================================================================

# __file__ is this script's own path. abspath makes it absolute, dirname strips
# the filename off. Result: the folder this script lives in. Using this instead
# of the current working directory means the script works no matter where you
# run it from.
HERE = os.path.dirname(os.path.abspath(__file__))

# microCT voxel size for this dataset: one pixel is 11.4 micrometres across.
# To convert a length in pixels to millimetres you MULTIPLY by MM_PER_PX.
# To convert an area in square pixels you multiply by MM_PER_PX squared.
UM_PER_PX = 11.4
MM_PER_PX = UM_PER_PX / 1000.0          # 0.0114 mm per pixel
MM2_PER_PX2 = MM_PER_PX ** 2            # 0.00012996 mm^2 per square pixel

# The released files contain two segments, CL1 and CL2. The tracking graph and
# the landmark file are both CL1-only, so this analysis is restricted to CL1.
# CL2 is 151 mm and is left unanalyzed.
GRAPH_SEGMENT = "CL1"

# Upadhye et al. 2022 (J Neural Eng 19(5), doi:10.1088/1741-2552/ac9643) report
# 17.8 +/- 6.1 fascicle split/merge events per centimetre of nerve, in eight
# mid-cervical vagus nerves. Their cross-sections were sampled every 100 um.
# Both of those matter for the comparison, so both are recorded here.
BENCHMARK_EVENTS_PER_CM = 17.8
BENCHMARK_SD = 6.1
BENCHMARK_SAMPLING_UM = 100.0


def find_one(pattern):
    """Locate exactly one file matching a wildcard pattern, inside HERE.

    Filenames in this release are not perfectly consistent (FasMorph vs
    FasMorphology, and the subject prefix varies), so matching by pattern is
    more robust than hardcoding names.
    """
    # os.path.join glues folder and pattern with the right separator for the OS.
    # glob.glob expands the "*" wildcards into a list of real paths.
    # sorted() makes the result deterministic rather than filesystem-order.
    hits = sorted(glob.glob(os.path.join(HERE, pattern)))
    if not hits:
        # An empty list is "falsy" in Python, so `not hits` means "nothing found".
        sys.exit("no file matching %r in %s" % (pattern, HERE))
    if len(hits) > 1:
        print("warning: %d matches for %r, using %s" % (len(hits), pattern, hits[0]))
    return hits[0]


# Resolve the four input files once, at import time, so a missing file fails
# immediately with a clear message rather than halfway through a stage.
FAS_CSV = find_one("*FasMorph*.csv")          # one row per fascicle per slice
NERVE_CSV = find_one("*NerveMorph*.csv")      # one row per slice, whole-nerve outline
LEVELS_CSV = find_one("*vagal_levels*.csv")   # five anatomical landmarks
GRAPHML = find_one("*RawFascicleTracking*.graphml")   # fascicle tracking graph


# ============================================================================
# LOADERS
# Unit conversion happens exactly once, here, in a column whose name states its
# units. Converting at each point of use is how unit bugs get introduced.
# ============================================================================

def load_fascicles():
    """One row per fascicle per slice, with millimetre columns added."""
    df = pd.read_csv(FAS_CSV)
    # df["new"] = df["old"] * k creates a new column by multiplying an existing
    # one element-wise. pandas does this across all 46,195 rows in one operation
    # (this is "vectorised") rather than looping in Python, which would be slow.
    df["area_mm2"] = df["area"] * MM2_PER_PX2
    df["diam_mm"] = df["equivalent_diameter"] * MM_PER_PX
    return df


def load_nerve():
    """One row per slice describing the nerve's outer boundary.

    Careful: 3,458 of 5,538 rows in this file are blank, because the boundary
    could not be traced where the nerve touches another structure. Anything
    computed from this file is computed on the 37.6% that survived.
    """
    df = pd.read_csv(NERVE_CSV)
    # Note the column is spelled eq_diameter here and equivalent_diameter in the
    # fascicle file. Same quantity, two spellings, in files from the same release.
    df["area_mm2"] = df["area"] * MM2_PER_PX2
    df["diam_mm"] = df["eq_diameter"] * MM_PER_PX
    return df


def segment_offset_mm(df):
    """Distance from the top of the whole scan to the start of this segment.

    The `index` column holds the global slice number, which for CL1 starts at
    2097 rather than 0. Multiplying that by the pixel size gives the offset
    between the whole-scan ruler and this segment's own ruler.
    """
    return df["index"].min() * MM_PER_PX


def load_landmarks(offset_mm):
    """The five anatomical landmarks, converted into this segment's coordinates.

    The landmark file gives positions in the whole-scan frame (23.93 to 86.16 mm).
    The morphology files restart at 0 for each segment. Subtracting the offset
    moves the landmarks into the morphology files' frame so the two can be
    plotted together. Stage 5 verifies the offset before this is trusted.
    """
    # .dropna(subset=[...]) removes rows where that column is empty.
    # .copy() makes an independent copy, so later edits don't emit pandas'
    # SettingWithCopyWarning about modifying a view of another DataFrame.
    lev = pd.read_csv(LEVELS_CSV).dropna(subset=["axis-0"]).copy()

    # The landmark file has three coordinate columns. axis-0 ranges 2099-7557,
    # matching the slice numbers, so it is the along-the-nerve axis. axis-1 and
    # axis-2 range 856-1238, matching positions inside a single image, so they
    # are across-the-slice and are not needed here.
    lev["z_mm"] = lev["axis-0"] * MM_PER_PX
    lev["dist_global_mm"] = lev["z_mm"] - offset_mm

    # Shorten the names for use as plot labels. .str.replace() applies a string
    # replacement to every row at once. regex=False treats the pattern as plain
    # text rather than a regular expression.
    lev["short"] = (lev["name"].str.replace("left level of ", "", regex=False)
                    .str.replace("the ", "", regex=False))
    return lev.sort_values("z_mm")


def draw_landmarks(ax, lev):
    """Draw dotted vertical lines at each landmark on an existing plot."""
    # .iterrows() yields (index, row) pairs. The underscore is the conventional
    # name for a value you must accept but don't intend to use.
    for _, r in lev.iterrows():
        # Skip landmarks that fall outside the visible x-range of the plot.
        if not (ax.get_xlim()[0] <= r["dist_global_mm"] <= ax.get_xlim()[1]):
            continue
        # zorder=0 puts these behind the data line rather than over it.
        ax.axvline(r["dist_global_mm"], color="0.55", lw=0.8, ls=":", zorder=0)
        ax.text(r["dist_global_mm"], ax.get_ylim()[1], " " + r["short"],
                rotation=90, va="top", ha="left", fontsize=7, color="0.35")


def load_graph():
    """Read the GraphML tracking file into two tables: nodes and edges.

    GraphML is XML. Each <node> is one fascicle on one slice. Each <edge> says
    "this fascicle here is the same fascicle as that one on the next slice down",
    and carries flags saying whether that link is a split or a merge.

    The file declares its attribute names once in a <key> block and then refers
    to them by short ids like "v_area", so the ids have to be translated back to
    readable names.
    """
    # XML namespaces: every tag in this file is really
    # {http://graphml.graphdrawing.org/xmlns}node, so searches must be
    # namespace-aware. This dict maps a short prefix "g" onto that URL.
    ns = {"g": "http://graphml.graphdrawing.org/xmlns"}
    tree = ET.parse(GRAPHML)
    root = tree.getroot()

    # Build {key id: (readable name, "node" or "edge")}.
    keys = {}
    for k in root.findall("g:key", ns):
        keys[k.get("id")] = (k.get("attr.name"), k.get("for"))

    graph = root.find("g:graph", ns)
    nodes, edges = [], []

    for n in graph.findall("g:node", ns):
        rec = {"node_id": n.get("id")}
        for d in n.findall("g:data", ns):
            # keys.get(x, default) returns the default if x is missing, which
            # keeps an unexpected attribute from crashing the parse.
            name = keys.get(d.get("key"), (d.get("key"), None))[0]
            rec[name] = d.text
        nodes.append(rec)

    for e in graph.findall("g:edge", ns):
        # "source" is the slice-above end of the arrow, "target" the slice-below.
        rec = {"source": e.get("source"), "target": e.get("target")}
        for d in e.findall("g:data", ns):
            name = keys.get(d.get("key"), (d.get("key"), None))[0]
            rec[name] = d.text
        edges.append(rec)

    # pd.DataFrame(list_of_dicts) turns a list of records into a table, using
    # the dict keys as column names.
    return keys, pd.DataFrame(nodes), pd.DataFrame(edges)


def truthy(s):
    """Convert a column of GraphML text into real True/False values.

    GraphML declares split and merge as attr.type="boolean" but still writes them
    as the strings "true" and "false". In Python, bool("false") is True, because
    any non-empty string is truthy. Testing the raw column would therefore flag
    every single edge as a split.
    """
    # .astype(str) guards against pandas having already parsed some values.
    # .str.strip() removes stray whitespace, .str.lower() normalises case,
    # .isin([...]) returns True where the value matches any listed spelling.
    return s.astype(str).str.strip().str.lower().isin(["true", "1", "yes"])


# ============================================================================
# STAGE 1 - INSPECT
# No arithmetic. Print what is actually in the files, so later stages are built
# on observed structure rather than assumed structure.
# ============================================================================

def stage1():
    print("=" * 70)
    print("STAGE 1 - inspect")
    print("=" * 70)

    for label, path in [("fascicles", FAS_CSV), ("nerve", NERVE_CSV),
                        ("levels", LEVELS_CSV)]:
        df = pd.read_csv(path)
        print("\n--- %s: %s" % (label, os.path.basename(path)))
        # df.shape is a (rows, columns) tuple; the % operator unpacks it.
        print("rows: %d  cols: %d" % df.shape)
        print("columns:", list(df.columns))
        print(df.dtypes.to_string())

        # .isna() gives a True/False table of missing cells; .sum() counts the
        # Trues per column (True counts as 1). Then keep only non-zero counts.
        na = df.isna().sum()
        na = na[na > 0]
        if len(na):
            print("MISSING VALUES - do NOT drop these silently. A blank row in the")
            print("nerve file means the outer boundary could not be traced, usually")
            print("where the nerve touches another structure. That is information.")
            print(na.to_string())
        else:
            print("no missing values")

        if "segment" in df.columns:
            print("segments present:",
                  sorted(df["segment"].dropna().unique().tolist()))

        if "dist_global" in df.columns:
            d = df["dist_global"].dropna()
            print("dist_global mm: min %.4f  max %.4f  unique slices %d"
                  % (d.min(), d.max(), d.nunique()))
            # A units trap-check. np.diff gives differences between consecutive
            # sorted values. If the median gap equals the pixel size, the column
            # is already in millimetres and MM_PER_PX is right. If it were 1.0,
            # the column would be slice indices and every conversion below would
            # be wrong by a factor of 87.
            step = np.diff(np.sort(d.unique()))
            print("median slice spacing: %.5f mm (expect %.5f)"
                  % (np.median(step), MM_PER_PX))
        print(df.head(3).to_string())

    print("\n--- graph: %s" % os.path.basename(GRAPHML))
    keys, nodes, edges = load_graph()
    print("declared keys (id -> name, scope):")
    for kid, (name, scope) in keys.items():
        print("  %-8s %-20s %s" % (kid, name, scope))
    print("nodes: %d  columns: %s" % (len(nodes), list(nodes.columns)))
    print("edges: %d  columns: %s" % (len(edges), list(edges.columns)))
    print(nodes.head(3).to_string())
    print(edges.head(3).to_string())

    # Print the distinct values of each edge attribute. This is how the split /
    # merge / identity flags were confirmed rather than assumed.
    for col in edges.columns:
        if col in ("source", "target"):
            continue
        print("edge attr %r values: %s"
              % (col, edges[col].value_counts(dropna=False).head(6).to_dict()))


# ============================================================================
# STAGE 2 - FASCICLE COUNT ALONG THE NERVE
# ============================================================================

def stage2():
    print("=" * 70)
    print("STAGE 2 - fascicle count along the nerve")
    print("=" * 70)

    fas = load_fascicles()
    # Boolean masking: the expression inside the brackets produces a True/False
    # value per row, and df[mask] keeps only the True rows.
    fas = fas[fas["segment"] == GRAPH_SEGMENT]

    # The whole measurement, in one line. Each row is one fascicle on one slice,
    # so grouping rows by their slice and counting group members gives the number
    # of fascicles visible at that depth. .size() counts rows per group.
    counts = fas.groupby("dist_global").size().sort_index()

    print("slices analyzed: %d" % len(counts))
    print("length of segment: %.2f mm" % (counts.index.max() - counts.index.min()))
    print("fascicle count  min %d  median %.1f  mean %.2f  max %d"
          % (counts.min(), counts.median(), counts.mean(), counts.max()))
    print("total fascicle cross-sections measured: %d" % len(fas))

    # Break the count down by anatomical band, using landmarks as bin edges.
    lev = load_landmarks(segment_offset_mm(fas))
    print("\nmean count between consecutive landmarks:")
    edges_mm = [counts.index.min()] + lev["dist_global_mm"].tolist() + \
        [counts.index.max()]
    labels = ["segment start"] + lev["short"].tolist()
    # zip(list[:-1], list[1:]) is the standard idiom for consecutive pairs:
    # it pairs each element with the one after it.
    for lo, hi, lab in zip(edges_mm[:-1], edges_mm[1:], labels):
        # Half-open interval (>= lo, < hi) so no slice is counted in two bands.
        seg = counts[(counts.index >= lo) & (counts.index < hi)]
        if len(seg):
            print("  %6.2f - %6.2f mm  mean %5.1f   (below %s)"
                  % (lo, hi, seg.mean(), lab))

    print("\nNOTE: the high counts at the top of this segment are most likely the")
    print("jugular and nodose ganglia, not fascicles. A ganglion is a mass of")
    print("nerve cell bodies; segmenting it yields many apparent sub-units.")

    # --- plot -----------------------------------------------------------
    # plt.subplots() returns a figure (the whole image) and an axes (one plot
    # area inside it). figsize is in inches; dpi at save time sets the pixels.
    fig, ax = plt.subplots(figsize=(9, 4))
    ax.plot(counts.index, counts.values, lw=0.6, color="#1f77b4")
    ax.axhline(counts.median(), ls="--", lw=1, color="k",
               label="median = %.0f" % counts.median())
    ax.set_xlabel("distance along nerve, superior to inferior (mm)")
    ax.set_ylabel("fascicles per cross-section")
    ax.set_title("SR042 left cervical vagus (%s): fascicle count vs distance"
                 % GRAPH_SEGMENT)
    ax.legend(loc="upper right")
    draw_landmarks(ax, lev)
    fig.tight_layout()   # shrinks margins so labels are not clipped
    fig.savefig(os.path.join(HERE, "fig1_fascicle_count.png"), dpi=150)
    print("wrote fig1_fascicle_count.png")


# ============================================================================
# STAGE 3 - DIAMETER DEFINITION AND SIZE DISTRIBUTIONS
# ============================================================================

def stage3():
    print("=" * 70)
    print("STAGE 3 - equivalent diameter and size distributions")
    print("=" * 70)

    fas = load_fascicles()
    fas = fas[fas["segment"] == GRAPH_SEGMENT]
    nerve = load_nerve()
    nerve = nerve[nerve["segment"] == GRAPH_SEGMENT]

    # --- the check ------------------------------------------------------
    # A fascicle is an irregular blob. "Equivalent diameter" gives it one size
    # number by asking: if this blob were a circle of the same area, how wide
    # would that circle be? Area = pi*r^2, so diameter = 2*sqrt(area/pi).
    #
    # Recomputing it from the area column and comparing to the supplied column
    # proves both are in the same units (pixels). If they disagreed, the units
    # would be misunderstood and every millimetre figure below would be wrong.
    expected = 2.0 * np.sqrt(fas["area"] / np.pi)
    resid = np.abs(expected - fas["equivalent_diameter"])
    print("max |2*sqrt(area/pi) - equivalent_diameter| = %.3e px" % resid.max())
    if resid.max() < 1e-6:
        print("definition confirmed: equivalent_diameter = 2*sqrt(area/pi), in pixels")
    else:
        # Stop rather than print plausible-looking but wrong numbers.
        print("DEFINITION MISMATCH - stop and check units before continuing")
        return

    print("\nfascicle diameter (mm): median %.4f  IQR %.4f-%.4f  max %.4f"
          % (fas["diam_mm"].median(),
             fas["diam_mm"].quantile(0.25), fas["diam_mm"].quantile(0.75),
             fas["diam_mm"].max()))
    print("fascicle area (mm^2):   median %.5f  max %.5f"
          % (fas["area_mm2"].median(), fas["area_mm2"].max()))

    # Everything below this line comes from the nerve file, which is 62% blank.
    n_have = nerve["area_mm2"].notna().sum()
    print("\n--- whole-nerve figures, computed on %d of %d slices (%.1f%%) ---"
          % (n_have, len(nerve), 100.0 * n_have / len(nerve)))
    print("These are BIASED. The blank rows are not randomly distributed: they")
    print("are 92-100%% blank through the jugular foramen region and 26%% blank")
    print("distally, so this is mostly a sample of the lower half of the segment.")
    print("nerve diameter (mm):    median %.4f  min %.4f  max %.4f"
          % (nerve["diam_mm"].median(), nerve["diam_mm"].min(),
             nerve["diam_mm"].max()))
    print("nerve area (mm^2):      median %.4f" % nerve["area_mm2"].median())

    # Fascicular area fraction: how much of the nerve's cross-section is
    # fascicle. Dividing an area by an area cancels the units, so if the pixel
    # conversion were wrong this number would be absurd rather than plausible.
    # That makes it a useful independent sanity check on the unit handling.
    frac = fas.groupby("dist_global")["area_mm2"].sum() / \
        nerve.set_index("dist_global")["area_mm2"]
    frac = frac.dropna()
    print("fascicular area fraction: median %.3f  (range %.3f-%.3f)"
          % (frac.median(), frac.min(), frac.max()))

    fig, axes = plt.subplots(1, 2, figsize=(11, 4))
    axes[0].hist(fas["diam_mm"], bins=80, color="#1f77b4")
    axes[0].set_xlabel("fascicle equivalent diameter (mm)")
    axes[0].set_ylabel("cross-sections")
    axes[0].set_title("fascicle size distribution")
    axes[1].plot(nerve["dist_global"], nerve["diam_mm"], lw=0.6, color="#d62728")
    axes[1].set_xlabel("distance along nerve (mm)")
    axes[1].set_ylabel("nerve equivalent diameter (mm)")
    axes[1].set_title("whole-nerve caliber (62% of slices missing)")
    fig.tight_layout()
    fig.savefig(os.path.join(HERE, "fig2_diameters.png"), dpi=150)
    print("wrote fig2_diameters.png")


# ============================================================================
# STAGE 4 - SPLIT AND MERGE EVENTS
# ============================================================================

def stage4():
    print("=" * 70)
    print("STAGE 4 - split and merge events")
    print("=" * 70)

    keys, nodes, edges = load_graph()
    print("edge attribute columns:",
          [c for c in edges.columns if c not in ("source", "target")])

    n_split_edges = int(truthy(edges["split"]).sum())
    n_merge_edges = int(truthy(edges["merge"]).sum())
    print("edges total: %d" % len(edges))
    print("edges flagged split: %d" % n_split_edges)
    print("edges flagged merge: %d" % n_merge_edges)

    # --- counting rule --------------------------------------------------
    # An event is a NODE, not an edge. A fascicle splitting three ways emits
    # three split-flagged edges but is still one event, so counting edges
    # overcounts by however many multi-way branchings there are.
    #
    # A split belongs to the node the arrows LEAVE, so use "source".
    # A merge belongs to the node the arrows ARRIVE at, so use "target".
    # Swapping these would count the daughters instead of the branch point.
    nodes = nodes.set_index("node_id")
    frame = pd.to_numeric(nodes["frame"], errors="coerce")

    # .unique() collapses repeated node ids, which is what turns three arrows
    # from one three-way split into a single event.
    split_nodes = set(edges.loc[truthy(edges["split"]), "source"].unique())
    merge_nodes = set(edges.loc[truthy(edges["merge"]), "target"].unique())

    # 47 nodes are flagged as both a split and a merge: several fascicles arrive
    # and more than one leaves, at the same point. Whether that is one event or
    # two is a judgement call. The union (each node counted once) is used here;
    # the sum is printed alongside so a reader can see both.
    both = split_nodes & merge_nodes
    n_events = len(split_nodes | merge_nodes)

    print("splitting nodes: %d  (from %d split-flagged edges)"
          % (len(split_nodes), n_split_edges))
    print("merging nodes:   %d  (from %d merge-flagged edges)"
          % (len(merge_nodes), n_merge_edges))
    print("nodes flagged as both: %d" % len(both))
    print("branching events: %d  (counting both-nodes once; %d if counted twice)"
          % (n_events, len(split_nodes) + len(merge_nodes)))

    span_mm = (frame.max() - frame.min()) * MM_PER_PX
    print("segment span: %.2f mm" % span_mm)

    # --- the comparison -------------------------------------------------
    # Upadhye 2022 report events per CENTIMETRE OF NERVE: "Over the middle 1 cm
    # of all eight nerves, there were 17.8 +/- 6.1 merging and splitting events."
    # So nerve length is the denominator. Dividing by total fascicle path length
    # instead answers a different question and is not comparable to their figure.
    rate = n_events / (span_mm / 10.0)
    print("\nevents per cm of nerve: %.1f      (Upadhye 2022: %.1f +/- %.1f)"
          % (rate, BENCHMARK_EVENTS_PER_CM, BENCHMARK_SD))
    print("ratio: %.1fx their rate" % (rate / BENCHMARK_EVENTS_PER_CM))

    # Position of each event along the nerve, for the per-band table and figure.
    ev_nodes = sorted(split_nodes | merge_nodes)
    ev_mm = (frame.loc[ev_nodes] - frame.min()) * MM_PER_PX
    ev_mm = ev_mm.dropna()

    # --- why the rates differ, part 1: sampling density ------------------
    # Upadhye analysed a cross-section every 100 um. This dataset has one every
    # 11.4 um. A fascicle that splits and rejoins within 100 um is invisible at
    # their sampling and visible here. Binning events into windows of a given
    # width and counting occupied windows imitates coarser sampling.
    print("\neffect of sampling density:")
    for step_um in [UM_PER_PX, 50.0, BENCHMARK_SAMPLING_UM, 200.0]:
        binw = step_um / 1000.0
        # np.floor(x / w) assigns each position to a bin; .nunique() counts how
        # many distinct bins contain at least one event.
        occupied = np.floor(ev_mm / binw).astype(int).nunique()
        print("  sampled every %6.1f um: %5.1f events/cm"
              % (step_um, occupied / (span_mm / 10.0)))

    # --- why the rates differ, part 2: anatomical position ---------------
    # Their nerves were mid-cervical with 6.6 +/- 2.8 fascicles. This segment is
    # upper cervical and its proximal half is dominated by ganglion.
    fas = load_fascicles()
    fas = fas[fas["segment"] == GRAPH_SEGMENT]
    lev = load_landmarks(segment_offset_mm(fas))
    cnt = fas.groupby("dist_global").size()

    print("\nrate by anatomical band:")
    print("  %-34s %9s %11s %12s" % ("band", "fascicles", "events/cm", "at 100um/cm"))
    bounds = lev["dist_global_mm"].tolist()
    names = lev["short"].tolist()
    for lo, hi, lab in zip(bounds[:-1], bounds[1:], names[:-1]):
        length_cm = (hi - lo) / 10.0
        in_band = ev_mm[(ev_mm >= lo) & (ev_mm < hi)]
        coarse = np.floor(in_band / (BENCHMARK_SAMPLING_UM / 1000.0)).astype(int).nunique()
        c = cnt[(cnt.index >= lo) & (cnt.index < hi)]
        print("  %-34s %9.1f %11.1f %12.1f"
              % ("below " + lab, c.mean(), len(in_band) / length_cm,
                 coarse / length_cm))

    print("\nA third difference cannot be tested here: Upadhye traced fascicles by")
    print("hand and counted a split only once the daughters had their own")
    print("perineurium. These events come from automated tracking, which has no")
    print("such judgement. Some unknown share of the excess may be the algorithm.")

    fig, ax = plt.subplots(figsize=(9, 3.5))
    ax.hist(ev_mm, bins=60, color="#2ca02c")
    draw_landmarks(ax, lev)
    ax.set_xlabel("distance along nerve (mm)")
    ax.set_ylabel("branching events per bin")
    ax.set_title("split/merge events, SR042 left cervical vagus (%s)" % GRAPH_SEGMENT)
    fig.tight_layout()
    fig.savefig(os.path.join(HERE, "fig3_events.png"), dpi=150)
    print("wrote fig3_events.png")


# ============================================================================
# STAGE 5 - ANATOMICAL LANDMARKS
# ============================================================================

def stage5():
    print("=" * 70)
    print("STAGE 5 - anatomical landmarks")
    print("=" * 70)

    nerve = load_nerve()
    nerve = nerve[nerve["segment"] == GRAPH_SEGMENT]

    # --- recovering an undocumented coordinate offset --------------------
    # Problem: the landmark file measures from the top of the whole scan, while
    # dist_global restarts at 0 for each segment. The release documents no way
    # to relate the two, so landmarks cannot be placed on the morphology data.
    #
    # Observation: the morphology files carry `index`, which starts at 2097 and
    # not 0. That looks like a global slice number.
    #
    # Hypothesis: dist_global = (index - index.min()) * MM_PER_PX
    #
    # This block tests that hypothesis on every row before using it. Guessing a
    # relationship like this and not checking it would silently place every
    # landmark in the wrong place, and the figure would still look fine.
    off = segment_offset_mm(nerve)
    resid = np.abs((nerve["index"] - nerve["index"].min()) * MM_PER_PX
                   - nerve["dist_global"]).max()
    print("index range: %d - %d" % (nerve["index"].min(), nerve["index"].max()))
    print("max |(index - index_min)*%.4f - dist_global| = %.3e mm"
          % (MM_PER_PX, resid))
    if resid > 1e-9:
        print("OFFSET UNVERIFIED - landmark mapping below is not trustworthy")
    else:
        print("verified: dist_global = (index - %d) * %.4f"
              % (nerve["index"].min(), MM_PER_PX))
    print("segment offset = %.4f mm; segment spans dist_global 0.00 - %.2f mm\n"
          % (off, nerve["dist_global"].max()))

    lev = load_landmarks(off)
    for _, r in lev.iterrows():
        inside = 0 <= r["dist_global_mm"] <= nerve["dist_global"].max()
        print("%-55s  z = %7.2f mm   dist_global = %7.2f mm  %s"
              % (r["name"], r["z_mm"], r["dist_global_mm"],
                 "" if inside else "(outside segment)"))

    print("\nAll five landmarks fall inside the segment, which therefore runs from")
    print("the superior border of the jugular foramen to the greater horn of the")
    print("hyoid. The dissection endpoints were chosen anatomically.")


# A dictionary mapping the --stage number onto the function that runs it. Storing
# functions in a dict like this avoids a long if/elif chain and makes the valid
# choices available to argparse automatically.
STAGES = {1: stage1, 2: stage2, 3: stage3, 4: stage4, 5: stage5}

# This guard means the code below runs only when the file is executed directly,
# not when it is imported by another script.
if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--stage", type=int, required=True, choices=sorted(STAGES),
                    help="which stage to run (1-5)")
    args = ap.parse_args()
    STAGES[args.stage]()      # look up the function, then call it

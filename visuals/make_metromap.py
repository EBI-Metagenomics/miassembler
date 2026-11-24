#!/usr/bin/env python3
"""
Generate the miassembler metro map (nf-core workflow-schematic style).

Usage:
    python3 visuals/make_metromap.py [--font MavenPro.ttf] [--out visuals/miassembler_metromap.svg]

If --font is given, the font is embedded in the SVG (base64) so the image renders
identically everywhere (GitHub, browsers). Maven Pro is the font recommended by the
nf-core brand guidelines: https://fonts.google.com/specimen/Maven+Pro

Render the PNG with a headless browser, e.g.:
    "/Applications/Google Chrome.app/Contents/MacOS/Google Chrome" --headless \
        --screenshot=visuals/miassembler_metromap.png --window-size=3260,1130 \
        --hide-scrollbars visuals/miassembler_metromap.svg
"""

import argparse
import base64
from pathlib import Path

W, H = 3260, 1130
SHIFT_X = 30      # trims the empty left margin before the input icons

# ---------------------------------------------------------------------------
# Lines (one per read type), ordered top -> bottom inside the trunk
# ---------------------------------------------------------------------------
LINES = {
    "short":  {"color": "#1C7ED6", "label": "Short reads (Illumina)"},
    "ont_hq": {"color": "#FAB005", "label": "ONT, high quality"},
    "ont_lq": {"color": "#F76707", "label": "ONT, low quality"},
    "pb_lq":  {"color": "#AE3EC9", "label": "PacBio, low quality (CLR)"},
    "hifi":   {"color": "#0CA678", "label": "PacBio, high quality (HiFi)"},
}
LONG = ["ont_hq", "ont_lq", "pb_lq", "hifi"]
ALL = ["short"] + LONG

GAP = 14          # distance between parallel tracks
LW = 9            # line width
R = 13            # station radius
TRUNK_TOP = 402   # y of first track in the shared trunk
SHORT_Y = 250     # y of the short-reads row
LONG_TOP = 590    # y of the first long-reads track

TY = {k: TRUNK_TOP + i * GAP for i, k in enumerate(ALL)}
LY = {k: LONG_TOP + i * GAP for i, k in enumerate(LONG)}

X_START, X_SPLIT0, X_SPLIT1 = 230, 560, 650
X_MERGE0, X_MERGE1, X_END = 2620, 2710, 3122
FORK = (1520, 1580, 1920, 1980)   # short-read assembler fork: out0, out1, in0, in1
FORK_X = 1750
FORK_DY = 45

FONT = "'Maven Pro', 'Helvetica Neue', Arial, sans-serif"

# ---------------------------------------------------------------------------
# Stages (grey boxes)
# ---------------------------------------------------------------------------
STAGES = [
    (250, 580, ["1", "Input"]),
    (600, 1480, ["2", "Reads QC & decontamination"]),
    (1500, 2075, ["3", "Assembly & polishing"]),
    (2095, 2440, ["4", "Contig QC &", "decontamination"]),
    (2460, 2920, ["5", "Coverage & statistics"]),
    (2940, 3110, ["6", "Reporting"]),
]
BOX_TOP, BOX_BOTTOM = 60, 750

# ---------------------------------------------------------------------------
# Stations: (x, lines, tool name, description lines, optional, row)
#   row: "trunk" | "short" | "long" | "fork_top" | "fork_bottom"
# ---------------------------------------------------------------------------
STATIONS = [
    # 1. Input
    (352, ALL, "fetch-tool", ["Download reads", "from ENA"], True, "trunk"),
    (508, ALL, "FastQC", ["Raw read", "quality"], False, "trunk"),
    # 2. Reads QC - short
    (740, ["short"], "fastp", ["Trimming &", "filtering"], False, "short"),
    (1050, ["short"], "bwa-mem2", ["Human, PhiX &", "host removal"], True, "short"),
    (1395, ["short"], "FastQC", ["Processed", "read quality"], False, "short"),
    # 2. Reads QC - long
    (710, LONG, "fastplong", ["Filtering &", "quality level"], False, "long"),
    (885, ["ont_hq", "ont_lq"], "Porechop_ABI", ["Adapter", "trimming"], False, "long"),
    (1050, ["ont_hq", "ont_lq", "pb_lq"], "Chopper", ["Filtering &", "lambda removal"], False, "long"),
    (1235, LONG, "minimap2", ["Human &", "host removal"], True, "long"),
    (1395, LONG, "FastQC", ["Processed", "read quality"], False, "long"),
    # 3. Assembly - short (fork: alternative assemblers)
    (FORK_X, ["short"], "metaSPAdes", [], False, "fork_top"),
    (FORK_X, ["short"], "MEGAHIT", [], False, "fork_bottom"),
    # 3. Assembly - long
    (1600, ["hifi"], "HiFiAdapterFilt", ["HiFi adapter", "removal"], False, "long"),
    (1750, LONG, "Flye", ["Assembly"], False, "long"),
    (1880, ["ont_lq", "pb_lq"], "Racon", ["Polishing"], False, "long"),
    (2008, ["ont_lq"], "Medaka", ["Polishing"], False, "long"),
    # 4. Contig QC
    (2178, ["short"], "SeqKit", ["Length", "filtering"], False, "short"),
    (2350, ["short"], "minimap2", ["Human, PhiX &", "host removal"], True, "short"),
    (2178, LONG, "minimap2", ["Human &", "host removal"], True, "long"),
    (2350, ["ont_lq", "pb_lq"], "proovframe", ["Frameshift", "correction"], False, "long"),
    # 5. Coverage & stats
    (2540, ["short"], "bwa-mem2", ["Coverage"], False, "short"),
    (2540, LONG, "minimap2", ["Coverage"], False, "long"),
    (2820, ALL, "metaQUAST", ["Assembly", "statistics"], False, "trunk"),
    # 6. Reporting
    (3030, ALL, "MultiQC", ["Run & study", "reports"], False, "trunk"),
]


def esc(s):
    return s.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")


def path_through(points):
    """Polyline through waypoints; y changes are drawn as smooth S-curves."""
    d = f"M{points[0][0]},{points[0][1]}"
    for (x0, y0), (x1, y1) in zip(points, points[1:]):
        if y0 == y1:
            d += f" L{x1},{y1}"
        else:
            xm = (x0 + x1) / 2
            d += f" C{xm},{y0} {xm},{y1} {x1},{y1}"
    return d


def line_paths():
    out = []
    for k in LONG:
        pts = [(X_START, TY[k]), (X_SPLIT0, TY[k]), (X_SPLIT1, LY[k]),
               (X_MERGE0, LY[k]), (X_MERGE1, TY[k]), (X_END, TY[k])]
        out.append((k, path_through(pts)))
    f0, f1, f2, f3 = FORK
    common = [(X_START, TY["short"]), (X_SPLIT0, TY["short"]), (X_SPLIT1, SHORT_Y), (f0, SHORT_Y)]
    tail = [(f3, SHORT_Y), (X_MERGE0, SHORT_Y), (X_MERGE1, TY["short"]), (X_END, TY["short"])]
    up = [(f0, SHORT_Y), (f1, SHORT_Y - FORK_DY), (f2, SHORT_Y - FORK_DY), (f3, SHORT_Y)]
    down = [(f0, SHORT_Y), (f1, SHORT_Y + FORK_DY), (f2, SHORT_Y + FORK_DY), (f3, SHORT_Y)]
    out.append(("short", path_through(common + up[1:] + tail[1:])))
    out.append(("short", path_through(down)))
    return out


def text(x, y, s, size=22, weight=400, anchor="middle", fill="#000"):
    return (f'<text x="{x}" y="{y}" font-family="{FONT}" font-size="{size}" '
            f'font-weight="{weight}" text-anchor="{anchor}" fill="{fill}">{esc(s)}</text>')


def station_y(lines, row):
    if row == "trunk":
        ys = [TY[k] for k in lines]
    elif row == "long":
        ys = [LY[k] for k in lines]
    elif row == "fork_top":
        ys = [SHORT_Y - FORK_DY]
    elif row == "fork_bottom":
        ys = [SHORT_Y + FORK_DY]
    else:
        ys = [SHORT_Y]
    return min(ys), max(ys)


def label_anchors(row):
    """y for the tool-name baseline (above) and first description line (below)."""
    if row == "trunk":
        return TY[ALL[0]] - 30, TY[ALL[-1]] + 44
    if row == "long":
        return LY[LONG[0]] - 30, LY[LONG[-1]] + 44
    if row == "fork_top":
        return SHORT_Y - FORK_DY - 30, None   # name above the upper branch
    if row == "fork_bottom":
        return None, SHORT_Y + FORK_DY + 46   # name below the lower branch
    return SHORT_Y - 30, SHORT_Y + 44


def station(x, lines, name, desc, optional, row):
    y0, y1 = station_y(lines, row)
    s = []
    s.append(f'<rect x="{x - R}" y="{y0 - R}" width="{2 * R}" height="{y1 - y0 + 2 * R}" '
             f'rx="{R}" fill="#fff" stroke="#000" stroke-width="4"/>')
    if not optional:   # mandatory steps carry a black dot
        r = 5.5
        s.append(f'<rect x="{x - r}" y="{y0 - r}" width="{2 * r}" height="{y1 - y0 + 2 * r}" '
                 f'rx="{r}" fill="#000"/>')
    top, bottom = label_anchors(row)
    if row == "fork_top":
        s.append(text(x, top, name, 26, 700))
    elif row == "fork_bottom":
        s.append(text(x, bottom, name, 26, 700))
    else:
        s.append(text(x, top, name, 26, 700))
        for i, d in enumerate(desc):
            s.append(text(x, bottom + i * 28, d, 24))
    return "\n".join(s)


def file_icon(x, y, label, caption, scale=1.0):
    """nf-core style document icon; (x, y) = top-left corner."""
    w, h, fold = 62 * scale, 80 * scale, 18 * scale
    p = (f"M{x},{y} H{x + w - fold} L{x + w},{y + fold} V{y + h} H{x} Z")
    s = [f'<path d="{p}" fill="#fff" stroke="#000" stroke-width="5" stroke-linejoin="round"/>',
         f'<path d="M{x + w - fold},{y} V{y + fold} H{x + w}" fill="none" stroke="#000" '
         f'stroke-width="4" stroke-linejoin="round"/>',
         f'<rect x="{x - 8}" y="{y + h * 0.42}" width="{w + 16}" height="{26 * scale}" rx="4" fill="#000"/>',
         text(x + w / 2, y + h * 0.42 + 19 * scale, label, 18 * scale, 700, fill="#fff")]
    for i, c in enumerate(caption):
        s.append(text(x + w / 2, y + h + 30 + i * 28, c, 24))
    return "\n".join(s)


def build(font_b64=None):
    o = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" viewBox="0 0 {W} {H}">']
    if font_b64:
        o.append("<defs><style>@font-face{font-family:'Maven Pro';font-weight:400 900;"
                 f"src:url(data:font/ttf;base64,{font_b64}) format('truetype');}}</style></defs>")
    o.append(f'<rect width="{W}" height="{H}" fill="#fff"/>')
    o.append(f'<g transform="translate({-SHIFT_X},0)">')

    # stage boxes
    for x0, x1, (_num, *title) in STAGES:
        o.append(f'<rect x="{x0}" y="{BOX_TOP}" width="{x1 - x0}" height="{BOX_BOTTOM - BOX_TOP}" '
                 f'rx="24" fill="#E9E9E9"/>')
        for i, t in enumerate(title):
            o.append(text(x0 + 24, BOX_TOP + 44 + i * 30, t, 26, 700, "start", "#8A8A8A"))

    # lines (short drawn first so long lines sit on top in the trunk)
    for k, d in line_paths():
        o.append(f'<path d="{d}" fill="none" stroke="{LINES[k]["color"]}" stroke-width="{LW}" '
                 f'stroke-linejoin="round"/>')
    # description between the alternative short-read assemblers
    o.append(text(FORK_X, SHORT_Y + 8, "Assembly", 24))


    # input terminus (triangle) and output terminus (bar)
    ty0, ty1 = TY[ALL[0]] - 12, TY[ALL[-1]] + 12
    o.append(f'<path d="M{X_START},{ty0} L{X_START + 26},{(ty0 + ty1) / 2} L{X_START},{ty1} Z" fill="#000"/>')
    o.append(f'<rect x="{X_END}" y="{ty0}" width="11" height="{ty1 - ty0}" fill="#000"/>')

    # stations
    for st in STATIONS:
        o.append(station(*st))

    # inputs / outputs
    o.append(file_icon(110, 280, "CSV", ["Samplesheet"]))
    o.append(file_icon(110, 470, "ENA", ["Study / run", "accession"]))
    o.append(file_icon(3165, 200, "FASTA", ["Cleaned", "contigs"], 1.0))
    o.append(file_icon(3165, 390, "TSV", ["Coverage"], 1.0))
    o.append(file_icon(3165, 550, "HTML", ["QC reports"], 1.0))

    # ---------------- bottom: title + legend ----------------
    lt = 790
    o.append(text(STAGES[0][0], lt + 70, "ebi-metagenomics/", 44, 700, "start", "#555"))
    o.append(text(STAGES[0][0], lt + 140, "miassembler", 84, 800, "start", "#000"))
    o.append(text(STAGES[0][0], lt + 190, "Metagenomic assembly of short and long reads", 26, 400, "start", "#444"))

    # line legend
    bw = 1230
    bx = STAGES[-1][1] - bw   # right edge aligned with the last stage box
    o.append(f'<rect x="{bx}" y="{lt}" width="{bw}" height="300" rx="24" fill="#E9E9E9"/>')
    o.append(text(bx + 30, lt + 50, "Legend", 32, 700, "start"))
    for i, k in enumerate(ALL):
        y = lt + 100 + i * 40
        o.append(f'<line x1="{bx + 30}" y1="{y}" x2="{bx + 100}" y2="{y}" stroke="{LINES[k]["color"]}" '
                 f'stroke-width="{LW + 3}"/>')
        o.append(text(bx + 120, y + 8, LINES[k]["label"], 26, 400, "start"))

    # station legend
    sx = bx + 500
    items = [
        ("circle", "Mandatory"),
        ("optional", "Optional (depends on parameters / samplesheet)"),
        ("in", "Input"),
        ("out", "Output"),
    ]
    for i, (kind, lab) in enumerate(items):
        y = lt + 100 + i * 40
        cx = sx + 20
        if kind == "circle":
            o.append(f'<circle cx="{cx}" cy="{y}" r="{R}" fill="#fff" stroke="#000" stroke-width="4"/>')
            o.append(f'<circle cx="{cx}" cy="{y}" r="5.5" fill="#000"/>')
        elif kind == "optional":
            o.append(f'<circle cx="{cx}" cy="{y}" r="{R}" fill="#fff" stroke="#000" stroke-width="4"/>')
        elif kind == "pill":
            o.append(f'<rect x="{cx - R}" y="{y - R - 4}" width="{2 * R}" height="{2 * R + 8}" rx="{R}" '
                     f'fill="#fff" stroke="#000" stroke-width="4"/>')
            o.append(f'<rect x="{cx - 5.5}" y="{y - 9.5}" width="11" height="19" rx="5.5" fill="#000"/>')
        elif kind == "in":
            o.append(f'<path d="M{cx - 10},{y - 14} L{cx + 12},{y} L{cx - 10},{y + 14} Z" fill="#000"/>')
        else:
            o.append(f'<rect x="{cx - 5}" y="{y - 15}" width="10" height="30" fill="#000"/>')
        o.append(text(sx + 52, y + 8, lab, 26, 400, "start"))

    o.append("</g>")
    o.append("</svg>")
    return "\n".join(o)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--font", type=Path, help="TTF to embed (Maven Pro)")
    ap.add_argument("--out", type=Path, default=Path(__file__).parent / "miassembler_metromap.svg")
    a = ap.parse_args()
    font = base64.b64encode(a.font.read_bytes()).decode() if a.font else None
    a.out.write_text(build(font))
    print(f"wrote {a.out}")


if __name__ == "__main__":
    main()

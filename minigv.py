#!/usr/bin/env python3
"""
minIGV - minimal, figure-ready IGV snapshots.

IGV desktop (File > Save PNG/SVG Image) writes SVG through Batik. The file
carries a lot of screen chrome that you usually delete by hand in
Illustrator: grey panel backgrounds, silver panel borders, track separator
lines, the chromosome ideogram with its red "current view" box, the
"249 kb" span bar, Java-only fonts (Dialog), ".bw" suffixes on track names,
and no page size. This script removes or fixes those with rules based on
IGV's fixed drawing style, so it works on any locus / any set of tracks.

Only the Python standard library is needed.

Usage
-----
    minigv snapshot.svg                           # writes snapshot.clean.svg
    minigv snapshot.svg -o figure.svg --drop-ruler
    minigv *.svg --report                         # batch, print what was removed

(installed with `pip install git+https://github.com/yu-kuai/minIGV`;
`python minigv.py ...` works the same from a checkout)

In Python / Jupyter:

    from minigv import clean_svg
    fig = clean_svg("snapshot.svg", drop_ruler=True)      # same options as the flags
    fig                    # inline preview (last line of a cell); also fig.plot(width="70%")
    fig.compare()          # IGV original above the cleaned version
    print(fig)             # what was removed / changed
    fig.save("figure.svg")

    from minigv import interactive, list_tracks
    interactive("snapshot.svg")                       # checkboxes per step and per track,
                                                      # live preview, Save button (ipywidgets),
                                                      # prints the matching clean_svg(...) call
    list_tracks("snapshot.svg")                       # names for drop_tracks= / --drop-track

By default everything below is applied; turn pieces off with the --keep-*
and --no-* flags (see --help).

  backgrounds  near-white panel / name-row background rectangles
  borders      silver panel frames, the name-panel divider, track separators
  ideogram     cytoband ideogram, band labels and the red current-view box
  span         the "<N> kb" span bar (line + arrowheads + label)
  font         Java logical fonts (Dialog, SansSerif...) -> --font (Arial)
  compact      move the locus label down next to the ruler once the ideogram
               is gone
  gaps         close blank bands taller than --min-gap (IGV pads each panel to a
               fixed height) down to --gap
  crop         set width/height/viewBox to the remaining content (+ --margin)

Optional: --drop-ruler (coordinate ticks + labels), --drop-locus,
--drop-track NAME (whole track; repeatable; --list-tracks shows the names),
--tidy-names (join IGV's wrapped track names, drop .bw/.bam/...; an extension is
kept when dropping it would give two tracks the same name), --rename FILE
(CSV/TSV old,new pairs, substring replace in track names).

Notes
-----
* Written for IGV desktop (Batik) SVGs. igv.js / IGV-Web SVGs are drawn
  differently and are not handled.
* Track data, gene models, data-range labels ("[0 - 60]") and gene names are
  never touched: removal rules only match full-panel near-white rectangles,
  long light-grey lines, IGV's frame paths, cytoband-labelled groups and the
  span bar.
* Some Python builds (e.g. some Homebrew 3.12+/3.14) ship a broken pyexpat;
  if you see "No module named expat", run with another interpreter such as
  /usr/bin/python3 or conda's python.
"""

__version__ = "0.1.0"

import argparse
import csv
import math
import os
import re
import sys
import xml.etree.ElementTree as ET
from bisect import bisect_right
from functools import lru_cache

SVG_NS = "http://www.w3.org/2000/svg"
XLINK_NS = "http://www.w3.org/1999/xlink"
ET.register_namespace("", SVG_NS)
ET.register_namespace("xlink", XLINK_NS)
NS = "{%s}" % SVG_NS

JAVA_FONTS = {"dialog", "dialoginput", "sansserif", "serif", "monospaced"}
TRACK_EXT = re.compile(r"\.(bw|bigwig|bigWig|bam|cram|bed|bedgraph|bedGraph|bb|bigbed|tdf|wig|gz|narrowPeak|broadPeak|seg)$", re.I)
EXT_FRAGMENT = re.compile(r"^(bw|bigwig|bam|cram|bed|bedgraph|bb|bigbed|tdf|wig|gz)$", re.I)
CYTOBAND = re.compile(r"^[pq](\d+|[A-H]\d*)(\.\d+)?$")      # human p15.4 / q13.2, mouse qA1 / qB5.1
COORD_LABEL = re.compile(r"^[\d,.]+\s?(bp|kb|Mb|Gb)$", re.I)
LOCUS = re.compile(r"^[\w.\-]+:[\d,]+-[\d,]+$")
LATEX_GREEK = {
    r"$\alpha$": "α", r"$\beta$": "β", r"$\gamma$": "γ", r"$\delta$": "δ",
    r"$\epsilon$": "ε", r"$\kappa$": "κ", r"$\mu$": "μ", r"$\sigma$": "σ",
}
NAMED_COLORS = {
    "white": (255, 255, 255), "black": (0, 0, 0), "silver": (192, 192, 192),
    "gray": (128, 128, 128), "grey": (128, 128, 128), "lightgray": (211, 211, 211),
    "lightgrey": (211, 211, 211), "darkgray": (169, 169, 169), "red": (255, 0, 0),
}


# --------------------------------------------------------------------------- geometry

@lru_cache(maxsize=None)
def parse_transform(s):
    """Return the affine matrix (a, b, c, d, e, f) of an SVG transform list."""
    m = (1.0, 0.0, 0.0, 1.0, 0.0, 0.0)
    for name, args in re.findall(r"(\w+)\s*\(([^)]*)\)", s or ""):
        v = [float(x) for x in re.split(r"[\s,]+", args.strip()) if x]
        if name == "translate":
            t = (1, 0, 0, 1, v[0], v[1] if len(v) > 1 else 0)
        elif name == "scale":
            t = (v[0], 0, 0, v[1] if len(v) > 1 else v[0], 0, 0)
        elif name == "matrix" and len(v) == 6:
            t = tuple(v)
        elif name == "rotate":
            r = math.radians(v[0])
            cos, sin = math.cos(r), math.sin(r)
            t = (cos, sin, -sin, cos, 0, 0)
            if len(v) == 3:
                t = mul(mul((1, 0, 0, 1, v[1], v[2]), t), (1, 0, 0, 1, -v[1], -v[2]))
        else:
            continue
        m = mul(m, t)
    return m


def mul(m, n):
    a, b, c, d, e, f = m
    A, B, C, D, E, F = n
    return (a * A + c * B, b * A + d * B, a * C + c * D, b * C + d * D,
            a * E + c * F + e, b * E + d * F + f)


def apply(m, x, y):
    a, b, c, d, e, f = m
    return a * x + c * y + e, b * x + d * y + f


def bbox_of(points):
    xs = [p[0] for p in points]
    ys = [p[1] for p in points]
    return min(xs), min(ys), max(xs), max(ys)


def num(el, key, default=0.0):
    try:
        return float(el.get(key, default))
    except (TypeError, ValueError):
        return default


# --------------------------------------------------------------------------- style

@lru_cache(maxsize=None)
def parse_style(s):          # cached: callers must not mutate the returned dict
    out = {}
    for part in (s or "").split(";"):
        if ":" in part:
            k, v = part.split(":", 1)
            out[k.strip()] = v.strip()
    return out


@lru_cache(maxsize=None)
def color(v):
    """'rgb(1,2,3)' / '#abc' / named -> (r, g, b); None for none/unknown."""
    if not v or v in ("none", "transparent"):
        return None
    v = v.strip().lower()
    if v in NAMED_COLORS:
        return NAMED_COLORS[v]
    m = re.match(r"rgb\(\s*(\d+)\s*,\s*(\d+)\s*,\s*(\d+)\s*\)", v)
    if m:
        return tuple(int(x) for x in m.groups())
    m = re.match(r"#([0-9a-f]{6}|[0-9a-f]{3})$", v)
    if m:
        h = m.group(1)
        if len(h) == 3:
            h = "".join(ch * 2 for ch in h)
        return tuple(int(h[i:i + 2], 16) for i in (0, 2, 4))
    return None


def is_light(c, thr):
    return c is not None and min(c) >= thr


def is_grey(c, lo, hi):
    return c is not None and max(c) - min(c) <= 8 and lo <= c[0] <= hi


def is_red(c):
    return c is not None and c[0] >= 150 and c[1] <= 60 and c[2] <= 60


# --------------------------------------------------------------------------- model

class Node:
    __slots__ = ("el", "parent", "tag", "m", "style", "box", "text")

    def __init__(self, el, parent, tag, m, style):
        self.el, self.parent, self.tag, self.m, self.style = el, parent, tag, m, style
        self.box = None
        self.text = None

    @property
    def fill(self):
        return color(self.style.get("fill"))

    @property
    def stroke(self):
        return color(self.style.get("stroke"))

    @property
    def font_size(self):
        m = re.match(r"([\d.]+)", self.style.get("font-size", "12"))
        return float(m.group(1)) * abs(self.m[3] or 1) if m else 12.0


def clip_boxes(root):
    """id -> local bbox of each <clipPath> (Batik writes userSpaceOnUse rectangles)."""
    out = {}
    for cp in root.iter(NS + "clipPath"):
        pts = []
        for sh in cp:
            t = sh.tag.replace(NS, "")
            if t == "rect":
                x, y = num(sh, "x"), num(sh, "y")
                pts += [(x, y), (x + num(sh, "width"), y + num(sh, "height"))]
            elif t in ("path", "polygon"):
                v = [float(x) for x in re.findall(r"-?\d*\.?\d+(?:e-?\d+)?", sh.get("d") or sh.get("points") or "")]
                pts += [(v[i], v[i + 1]) for i in range(0, len(v) - 1, 2)]
        if pts and cp.get("id"):
            out[cp.get("id")] = bbox_of(pts)
    return out


def build(root):
    nodes, parent_of = [], {}
    clips = clip_boxes(root)

    def walk(el, parent, m, style):
        tag = el.tag.replace(NS, "")
        if tag in ("defs", "clipPath", "mask", "pattern", "linearGradient", "radialGradient"):
            return
        m = mul(m, parse_transform(el.get("transform")))
        st = dict(style)
        for k in ("fill", "stroke", "font-size", "font-family", "font-weight"):
            if el.get(k):
                st[k] = el.get(k)
        st.update(parse_style(el.get("style")))
        n = Node(el, parent, tag, m, st)
        n.box = element_box(n)
        # clip to the element's own clip-path; fully clipped-away art counts as absent
        cm = re.match(r"url\(#([^)]+)\)", st.get("clip-path", ""))
        if n.box and cm and cm.group(1) in clips:
            cb = clips[cm.group(1)]
            cb = bbox_of([apply(m, cb[0], cb[1]), apply(m, cb[2], cb[3])])
            b = (max(n.box[0], cb[0]), max(n.box[1], cb[1]), min(n.box[2], cb[2]), min(n.box[3], cb[3]))
            # lines are zero-thickness, so "touching" counts as visible
            n.box = b if (b[0] <= b[2] + 1e-6 and b[1] <= b[3] + 1e-6) else None
        if tag == "text":
            n.text = "".join(el.itertext()).strip()
        nodes.append(n)
        for c in list(el):
            parent_of[c] = el
            walk(c, n, m, st)

    walk(root, None, (1, 0, 0, 1, 0, 0), {})
    return nodes, parent_of


def element_box(n):
    el, m, t = n.el, n.m, n.tag
    try:
        if t == "rect":
            x, y, w, h = num(el, "x"), num(el, "y"), num(el, "width"), num(el, "height")
            return bbox_of([apply(m, x, y), apply(m, x + w, y + h), apply(m, x + w, y), apply(m, x, y + h)])
        if t == "line":
            return bbox_of([apply(m, num(el, "x1"), num(el, "y1")), apply(m, num(el, "x2"), num(el, "y2"))])
        if t in ("polygon", "polyline"):
            v = [float(x) for x in re.split(r"[\s,]+", (el.get("points") or "").strip()) if x]
            pts = [apply(m, v[i], v[i + 1]) for i in range(0, len(v) - 1, 2)]
            return bbox_of(pts) if pts else None
        if t == "path":
            v = [float(x) for x in re.findall(r"-?\d*\.?\d+(?:e-?\d+)?", el.get("d") or "")]
            pts = [apply(m, v[i], v[i + 1]) for i in range(0, len(v) - 1, 2)]
            return bbox_of(pts) if pts else None
        if t == "text":
            x, y = apply(m, num(el, "x"), num(el, "y"))
            txt = "".join(el.itertext())
            fs = float(re.match(r"([\d.]+)", n.style.get("font-size", "12") + "x").group(1) or 12)
            fs *= abs(m[3] or 1)
            return (x, y - fs * 0.8, x + 0.56 * fs * len(txt.strip()), y + fs * 0.25)
    except (ValueError, AttributeError, IndexError):
        return None
    return None


def area(b):
    return (b[2] - b[0]) * (b[3] - b[1]) if b else 0


def visible(n):
    if n.tag == "text":
        return bool(n.text) and n.style.get("fill") != "none"
    if n.tag in ("rect", "line", "polygon", "polyline", "path"):
        if n.tag == "line":
            return n.style.get("stroke", "black") != "none"
        return n.fill is not None or n.stroke is not None
    return False


# --------------------------------------------------------------------------- cleaner

class Cleaner:
    def __init__(self, tree, args):
        self.tree, self.root, self.args = tree, tree.getroot(), args
        self.nodes, self.parent_of = build(self.root)
        self.by_el = {n.el: n for n in self.nodes}
        self.removed = set()
        self.dead = set()          # removed elements and everything inside them
        self.pending = {}          # parent -> children to detach (see flush)
        self.report = {}
        self.texts = [n for n in self.nodes if n.tag == "text" and n.text and n.box]   # box None = clipped away
        self._layout()

    # -- layout of the IGV window -------------------------------------------
    def _layout(self):
        frames = [n for n in self.nodes if n.tag == "rect" and n.fill is None and is_grey(n.stroke, 100, 245)
                  and n.box and area(n.box) > 5000]
        canvas = [n.box for n in self.nodes if n.tag == "rect" and n.box and area(n.box) > 5000]
        canvas += [f.box for f in frames]
        if canvas:
            self.canvas = (min(b[0] for b in canvas), min(b[1] for b in canvas),
                           max(b[2] for b in canvas), max(b[3] for b in canvas))
        else:
            all_boxes = [n.box for n in self.nodes if n.box]
            self.canvas = bbox_of([(b[0], b[1]) for b in all_boxes] + [(b[2], b[3]) for b in all_boxes])
        coords = [t for t in self.texts if COORD_LABEL.match(t.text)]
        loci = [t for t in self.texts if LOCUS.match(t.text)]
        left = [t.box[0] for t in coords + loci]
        self.data_left = (min(left) - 6) if left else 150.0
        # ruler tick labels share one baseline; the span label ("249 kb") stands alone
        rows = {}
        for t in coords:
            rows.setdefault(round(t.box[3]), []).append(t)
        self.tick_labels = max(rows.values(), key=len) if rows else []
        if len(self.tick_labels) < 2:
            self.tick_labels = []
        self.span_labels = [t for t in coords if t not in self.tick_labels]
        self.locus = loci[0] if loci else None
        tops = sorted(f.box[1] for f in frames)
        hdr = [f for f in frames if f.box[1] <= (tops[0] + 2 if tops else 0)]
        if self.tick_labels:
            ty = self.tick_labels[0].box[3]
            below = [f.box[3] for f in frames if f.box[1] <= ty <= f.box[3]]
            self.header_bottom = min(below) if below else ty + 20
        elif hdr:
            self.header_bottom = hdr[0].box[3]
        else:
            self.header_bottom = None
        self.tracks = self._find_tracks()

    def _find_tracks(self):
        """Track rows as [{'name', 'y0', 'y1'}], top to bottom.

        IGV paints one near-white background box per track in the name column;
        those give exact row bounds. Without them, rows are split halfway between
        consecutive track names."""
        top = self.header_bottom or 0
        names = [t for t in self.texts if t.box[0] < self.data_left and t.box[1] >= top - 1]
        if names:                                   # track names are left-aligned in the name column;
            x_min = min(t.box[0] for t in names)    # this drops gene labels clipped at the data edge
            names = [t for t in names if t.box[0] <= x_min + 3]
        cands = [n for n in self.nodes if n.tag == "rect" and n.box and is_light(n.fill, 235)
                 and n.box[0] <= 2 and n.box[2] <= self.data_left + 2 and n.box[3] - n.box[1] >= 8
                 and n.box[1] >= top - 1]
        inside = lambda a, b: a.box[1] >= b.box[1] - 0.5 and a.box[3] <= b.box[3] + 0.5 and a is not b
        rows = [r for r in cands if sum(inside(o, r) for o in cands) < 2]      # drop whole-panel boxes
        rows = sorted({(round(r.box[1], 2), round(r.box[3], 2)) for r in rows})
        # each name belongs to the tightest row that holds it (rows of stacked panels can overlap)
        owner = {}
        for t in names:
            cy = (t.box[1] + t.box[3]) / 2
            fits = [r for r in rows if r[0] <= cy < r[1]]
            if fits:
                owner.setdefault(min(fits, key=lambda r: r[1] - r[0]), []).append(t)
        out = []
        for y0, y1 in rows:
            lab = sorted(owner.get((y0, y1), []), key=lambda t: t.box[3])
            if lab:
                out.append({"name": "".join("".join(t.el.itertext()) for t in lab).strip(), "y0": y0, "y1": y1})
        if out:
            return out
        # fallback: group name lines that sit within ~1.8 line heights, split rows at midpoints
        groups = []
        for t in sorted(names, key=lambda t: t.box[3]):
            if groups and t.box[3] - groups[-1][-1].box[3] <= t.font_size * 1.8:
                groups[-1].append(t)
            else:
                groups.append([t])
        centers = [((g[0].box[1] + g[-1].box[3]) / 2) for g in groups]
        for i, g in enumerate(groups):
            y0 = (centers[i - 1] + centers[i]) / 2 if i else top
            y1 = (centers[i] + centers[i + 1]) / 2 if i + 1 < len(groups) else self.canvas[3]
            out.append({"name": "".join("".join(t.el.itertext()) for t in g).strip(), "y0": y0, "y1": y1})
        return out

    def drop_tracks(self, names):
        """Remove whole tracks (name, data, scale label) by their IGV name."""
        want = {n.strip().lower() for n in names}
        bands = sorted((t["y0"], t["y1"]) for t in self.tracks if t["name"].strip().lower() in want)
        if not bands:
            return
        starts = [b[0] for b in bands]
        for n in self.nodes:
            if n.box is None or n.tag not in ("rect", "line", "polygon", "polyline", "path", "text"):
                continue
            cy = (n.box[1] + n.box[3]) / 2
            i = bisect_right(starts, cy) - 1
            if i >= 0 and cy < bands[i][1] and n.el not in self.dead:
                self.remove(n, "elements of dropped tracks")
        self.report["tracks dropped"] = len(bands)

    def remove(self, node, why):
        if node.el in self.removed:
            return
        self.removed.add(node.el)
        self.dead.update(node.el.iter())
        p = self.parent_of.get(node.el)
        if p is not None:                 # detach later in one pass per parent (list.remove is O(n))
            self.pending.setdefault(p, set()).add(node.el)
        self.report[why] = self.report.get(why, 0) + 1

    def alive(self, n):
        return n.el not in self.dead

    # -- rules ------------------------------------------------------------------
    def backgrounds(self):
        for n in self.nodes:
            if n.tag == "rect" and is_light(n.fill, 235) and n.box:
                w, h = n.box[2] - n.box[0], n.box[3] - n.box[1]
                if area(n.box) >= 3000 or max(w, h) >= 50:      # panels, name rows, divider column
                    self.remove(n, "background rectangles")

    def borders(self):
        for n in self.nodes:
            if not n.box:
                continue
            w, h = n.box[2] - n.box[0], n.box[3] - n.box[1]
            if n.tag == "rect" and n.fill is None and is_grey(n.stroke, 100, 245) and max(w, h) >= 50:
                self.remove(n, "panel borders")
            elif n.tag == "path" and is_grey(n.fill, 100, 245) and (n.el.get("d") or "").count("M") >= 2 and area(n.box) >= 3000:
                self.remove(n, "panel borders")      # IGV draws frames as an outer+inner rectangle path
            elif n.tag == "line" and is_light(n.stroke, 180) and math.hypot(w, h) >= 100:
                self.remove(n, "divider / separator lines")

    def ideogram(self):
        bands = [t for t in self.texts if CYTOBAND.match(t.text)]
        if len(bands) >= 3:
            grp = self._common_ancestor([b.el for b in bands])
            protect = {t.el for t in self.tick_labels + ([self.locus] if self.locus else [])}
            g = self.by_el.get(grp)
            if g is not None and g.tag == "g" and not (self._descendants(grp) & protect):
                ys = [b.box for b in bands]
                y0, y1 = min(b[1] for b in ys) - 25, max(b[3] for b in ys) + 5
                self.remove(g, "ideogram")
            else:   # fall back: drop band labels and every shape in the band row
                y0 = min(b.box[1] for b in bands) - 25
                y1 = max(b.box[3] for b in bands) + 5
                for n in self.nodes:
                    if n.box and n.tag in ("rect", "polygon", "path", "line", "text") and y0 <= n.box[1] and n.box[3] <= y1 \
                            and n.box[0] >= self.data_left and n.el not in protect:
                        self.remove(n, "ideogram")
        else:
            y0 = y1 = None
        # red current-view box (and red centromere marks if they sat outside the group)
        for n in self.nodes:
            if n.box and (is_red(n.fill) or is_red(n.stroke)) and self.header_bottom and n.box[3] <= self.header_bottom:
                if y0 is None or (y0 - 5 <= n.box[1] and n.box[3] <= y1 + 5):
                    self.remove(n, "current-view box")

    def span(self):
        for lab in self.span_labels:
            self.remove(lab, "span bar")
            cy = (lab.box[1] + lab.box[3]) / 2
            for n in self.nodes:
                if not n.box or n.box[0] < self.data_left - 2:
                    continue
                ncy = (n.box[1] + n.box[3]) / 2
                flat = n.box[3] - n.box[1] <= 8
                if n.tag in ("line", "polygon") and flat and abs(ncy - cy) <= 6:
                    self.remove(n, "span bar")

    def ruler(self):
        if not self.tick_labels or self.header_bottom is None:
            return
        y0 = min(t.box[1] for t in self.tick_labels) - 2
        for n in self.nodes:
            if n.box and n.box[0] >= self.data_left - 2 and n.box[1] >= y0 and n.box[3] <= self.header_bottom + 1 \
                    and n.tag in ("line", "text", "rect", "polygon"):
                if n.tag == "text" and not COORD_LABEL.match(n.text or ""):
                    continue
                self.remove(n, "ruler")

    def drop_locus(self):
        if self.locus:
            self.remove(self.locus, "locus label")

    def compact(self):
        """After the ideogram is gone, sit the locus label just above the ruler."""
        if not self.locus or not self.alive(self.locus):
            return
        fs = self.locus.font_size
        ticks = [t for t in self.tick_labels if self.alive(t)]
        if ticks:
            target = min(t.box[1] for t in ticks) - 0.6 * fs   # baseline just above the ruler labels
        elif self.header_bottom is not None:
            target = self.header_bottom - 0.4 * fs             # ruler dropped: sit at the header's bottom
        else:
            return
        cur = self.locus.box[3] - fs * 0.25
        dy = target - cur
        if dy > fs:                                    # only ever move it down into the freed space
            self._shift(self.locus, dy)
            self.report["locus label moved down"] = 1

    def fonts(self, family):
        def fix(el):
            st = el.get("style")
            if st and "font-family" in st:
                def rep(m):
                    fam = m.group(1).strip().strip("'\"")
                    return "font-family:'%s'" % family if fam.lower() in JAVA_FONTS else m.group(0)
                el.set("style", re.sub(r"font-family:\s*([^;]+)", rep, st))
            if el.get("font-family") and el.get("font-family").strip("'\"").lower() in JAVA_FONTS:
                el.set("font-family", family)
        for el in self.root.iter():
            fix(el)
        self.report["fonts -> %s" % family] = 1

    def names(self, strip_ext, rename, merge=True):
        name_texts = [t for t in self.texts if self.alive(t) and t.box[0] < self.data_left and t.box[3] > (self.header_bottom or 0)]
        name_texts.sort(key=lambda t: (round(t.box[0]), t.box[3]))
        def tidy(s):
            if strip_ext:
                s = TRACK_EXT.sub("", s)
            for old, rep in rename:                      # case-insensitive (IGV/ENCODE names vary)
                s = re.sub(re.escape(old), lambda _m, r=rep: r, s, flags=re.I)
            return s

        # merge IGV's wrapped names (same x, next baseline within ~1.8 line heights),
        # but only when the tidied result fits the name column - otherwise keep the wrap
        i = 0
        while merge and i < len(name_texts) - 1:
            a, b = name_texts[i], name_texts[i + 1]
            room = self.data_left - a.box[0] - 4
            if abs(a.box[0] - b.box[0]) < 1 and 0 < b.box[3] - a.box[3] <= a.font_size * 1.8 \
                    and (EXT_FRAGMENT.match(b.text) or 0.56 * a.font_size * len(tidy(a.text + b.text)) <= room):
                merged = "".join(a.el.itertext()).lstrip() + b.text   # IGV keeps the space it wrapped at
                shift = (b.box[3] - a.box[3]) / 2
                self._set_text(a, merged)
                self._shift(a, shift)
                self.remove(b, "wrapped track names merged")
                name_texts.pop(i + 1)
                continue
            i += 1
        # stripping extensions must not make two different tracks share a name
        # (e.g. H3K4me3.bw signal + H3K4me3.bed peaks): keep the extension for those
        cand = {t.el: tidy(t.text) for t in name_texts}
        seen = {}
        for t in name_texts:
            seen.setdefault(cand[t.el], set()).add(t.text)
        for t in name_texts:
            new = cand[t.el]
            if strip_ext and len(seen[new]) > 1:
                no_strip = t.text
                for old, rep in rename:
                    no_strip = re.sub(re.escape(old), lambda _m, r=rep: r, no_strip, flags=re.I)
                new = no_strip
                self.report["extension kept (name clash)"] = self.report.get("extension kept (name clash)", 0) + 1
            if new != t.text:
                self._set_text(t, new)
                self.report["track names changed"] = self.report.get("track names changed", 0) + 1

    def close_gaps(self, keep, min_gap):
        """Shrink empty horizontal bands (IGV pads panels to a fixed height) to `keep` pt.

        Every drawn element below a gap is moved up. Batik clips each element in
        its own user space, so the clip region travels with it."""
        leaves = [n for n in self.nodes if n.box and n.tag in ("rect", "line", "polygon", "polyline", "path", "text")
                  and self.alive(n)]
        spans = sorted((n.box[1], n.box[3]) for n in leaves if visible(n))
        if not spans:
            return
        gaps, cur_end = [], spans[0][1]
        for y0, y1 in spans[1:]:
            if y0 - cur_end > min_gap:
                gaps.append((cur_end, y0))
            cur_end = max(cur_end, y1)
        if not gaps:
            return
        total = 0.0
        for n in leaves:
            top = n.box[1]
            dy = -sum((g1 - g0) - keep for g0, g1 in gaps if g1 <= top + 1e-6)
            if dy:
                self._shift(n, dy)
        total = sum((g1 - g0) - keep for g0, g1 in gaps)
        self.report["blank gaps closed"] = "%d (%.0f pt)" % (len(gaps), total)

    def simplify(self):
        self.flush()
        """Fewer objects, same picture: IGV draws one <rect> per bar (tens of thousands)
        and many zero-size lines with stroke:none. Delete the invisible lines and merge
        each run of sibling rects that share style/transform into a single <path>."""
        merged = dropped = 0
        for parent in list(self.root.iter()):
            kids = list(parent)
            if len(kids) < 2:
                continue
            new, run, key = [], [], None

            def flush():
                nonlocal merged
                if len(run) >= 2:
                    d = []
                    for r in run:
                        x, y, w, h = num(r, "x"), num(r, "y"), num(r, "width"), num(r, "height")
                        d.append("M%g %gh%gv%gh%gZ" % (x, y, w, h, -w))
                    path = ET.Element(NS + "path", {"d": "".join(d)})
                    for k, v in run[0].attrib.items():
                        if k not in ("x", "y", "width", "height", "rx", "ry"):
                            path.set(k, v)
                    new.append(path)
                    merged += len(run)
                else:
                    new.extend(run)
                run.clear()

            for k in kids:
                tag = k.tag.replace(NS, "")
                st = k.get("style", "")
                if tag == "line" and re.search(r"(^|;)\s*stroke:\s*none", st):
                    dropped += 1               # a line without a stroke draws nothing
                    continue
                kk = (st, k.get("transform"), k.get("fill"), k.get("stroke")) if tag == "rect" and not k.get("rx") \
                    and len(k) == 0 and "stroke:none" in st.replace(" ", "") else None
                if kk is not None and kk == key:
                    run.append(k)
                    continue
                flush()
                key = kk
                if kk is not None:
                    run.append(k)
                else:
                    new.append(k)
            flush()
            if merged or dropped:
                parent[:] = new
        if dropped:
            self.report["invisible lines removed"] = dropped
        if merged:
            self.report["bar rectangles merged"] = merged

    def crop(self, margin):
        boxes = [n.box for n in self.nodes if self.alive(n) and n.box and visible(n)]
        if not boxes:
            return
        cx0, cy0, cx1, cy1 = self.canvas
        x0 = max(min(b[0] for b in boxes), cx0) - margin
        y0 = max(min(b[1] for b in boxes), cy0) - margin
        x1 = min(max(b[2] for b in boxes), cx1) + margin
        y1 = min(max(b[3] for b in boxes), cy1) + margin
        w, h = x1 - x0, y1 - y0
        self.root.set("viewBox", "%g %g %g %g" % (x0, y0, w, h))
        self.root.set("width", "%gpt" % w)
        self.root.set("height", "%gpt" % h)
        self.report["page set to"] = "%g x %g pt" % (w, h)

    # -- helpers ----------------------------------------------------------------
    def _shift(self, node, dy):
        """Move a text down by dy (page units). Batik clips each text to its band
        (clip-path in the element's own user space), so prepend a translate to the
        element itself: the clip region moves with it."""
        el = node.el
        sy = node.m[3] or 1
        el.set("transform", ("translate(0,%g) " % (dy / sy) + el.get("transform", "")).strip())
        node.m = mul(node.m, (1, 0, 0, 1, 0, dy / sy))
        if node.box:
            node.box = (node.box[0], node.box[1] + dy, node.box[2], node.box[3] + dy)

    def _set_text(self, node, new):
        el = node.el
        for c in list(el):
            el.remove(c)
        el.text = new
        node.text = new

    def _descendants(self, el):
        return set(el.iter())

    def _common_ancestor(self, els):
        def chain(e):
            out = []
            while e is not None:
                out.append(e)
                e = self.parent_of.get(e)
            return out[::-1]
        chains = [chain(e) for e in els]
        lca = None
        for level in zip(*chains):
            if all(x is level[0] for x in level):
                lca = level[0]
            else:
                break
        return lca

    def flush(self):
        for p, gone in self.pending.items():
            p[:] = [c for c in p if c not in gone]
        self.pending.clear()

    def prune_empty_groups(self):
        self.flush()
        changed = True
        while changed:
            changed = False
            for el in list(self.root.iter()):
                for c in list(el):
                    if c.tag == NS + "g" and len(c) == 0 and not (c.text or "").strip():
                        el.remove(c)
                        changed = True


# --------------------------------------------------------------------------- CLI

def read_rename(path):
    pairs = []
    with open(path, newline="", encoding="utf-8-sig") as fh:
        sample = fh.read(4096)
        fh.seek(0)
        if "\t" in sample:
            rows = csv.reader(fh, delimiter="\t")
        elif "," in sample:
            rows = csv.reader(fh)
        else:
            rows = (line.split() for line in fh)
        for i, row in enumerate(rows):
            row = [c.strip() for c in row if c.strip()]
            if len(row) < 2:
                continue
            old, new = row[0], row[1]
            if i == 0 and old.lower() in ("old", "from", "target", "targets", "name"):
                continue
            for k, v in LATEX_GREEK.items():
                new = new.replace(k, v)
            if old != new:
                pairs.append((old, new))
    pairs.sort(key=lambda p: -len(p[0]))   # longest first, so substrings don't clobber
    return pairs


def clean(path, out, args, rename):
    tree, c = _run(path, args, rename)
    tree.write(out, encoding="utf-8", xml_declaration=True)
    return c.report


def _run(path, args, rename):
    tree = ET.parse(path)
    c = Cleaner(tree, args)
    if not args.keep_backgrounds:
        c.backgrounds()
    if not args.keep_borders:
        c.borders()
    if not args.keep_ideogram:
        c.ideogram()
    if not args.keep_span:
        c.span()
    if args.drop_ruler:
        c.ruler()
    if getattr(args, "drop_tracks", None):
        c.drop_tracks(args.drop_tracks)
    if args.drop_locus:
        c.drop_locus()
    elif not args.keep_ideogram and not args.no_compact:
        c.compact()
    if args.font.lower() != "keep":
        c.fonts(args.font)
    if getattr(args, "tidy_names", False):
        c.names(not args.keep_ext, rename)
    elif rename:                                   # an explicit rename table applies on its own
        c.names(False, rename, merge=False)
    if not args.keep_gaps:
        c.close_gaps(args.gap, args.min_gap)
    c.prune_empty_groups()
    if not args.no_crop:
        c.crop(args.margin)
    if not getattr(args, "no_simplify", False):
        c.simplify()
    c.flush()
    return tree, c


# --------------------------------------------------------------------------- Python / Jupyter API

class CleanedSVG:
    """Result of clean_svg(). In Jupyter it displays itself; also .plot(), .compare(), .save()."""

    def __init__(self, source, svg, report, canvas):
        self.source, self.svg, self.report, self._canvas = source, svg, report, canvas

    # -- display -----------------------------------------------------------------
    def plot(self, width="100%"):
        """Show the cleaned figure inline (Jupyter)."""
        from IPython.display import HTML, display
        display(HTML(self._html(self.svg, width)))

    def compare(self, width="100%"):
        """Show IGV's original above the cleaned version (Jupyter)."""
        from IPython.display import HTML, display
        with open(self.source, encoding="utf-8") as fh:
            orig = fh.read()
        x0, y0, x1, y1 = self._canvas
        orig = re.sub(r"<svg\b", '<svg viewBox="%g %g %g %g"' % (x0, y0, x1 - x0, y1 - y0), orig, count=1)
        cap = 'style="font:bold 12px sans-serif;margin:8px 0 2px"'
        display(HTML("<div %s>IGV original</div>%s<div %s>cleaned</div>%s"
                     % (cap, self._html(orig, width, border=True), cap, self._html(self.svg, width, border=True))))

    def _repr_html_(self):
        return self._html(self.svg, "100%")

    def __repr__(self):
        lines = ["CleanedSVG(%s)" % os.path.basename(self.source)]
        lines += ["    %-30s %s" % (k, v) for k, v in self.report.items()]
        return "\n".join(lines)

    # -- output ------------------------------------------------------------------
    def save(self, path):
        if os.path.abspath(path) == os.path.abspath(self.source):
            raise ValueError("refusing to overwrite the IGV original: %s" % path)
        with open(path, "w", encoding="utf-8") as fh:
            fh.write(self.svg)
        return path

    @staticmethod
    def _html(svg, width, border=False):
        # Every IGV SVG names its clip paths clipPath1, clipPath2...; several inline
        # SVGs in one notebook page would steal each other's clips, so make ids unique.
        CleanedSVG._n = getattr(CleanedSVG, "_n", 0) + 1
        pre = "igv%d_" % CleanedSVG._n
        svg = re.sub(r'<\?xml[^>]*\?>|<!DOCTYPE[^>]*>', "", svg)
        svg = re.sub(r'\bid="([^"]+)"', lambda m: 'id="%s%s"' % (pre, m.group(1)), svg)
        svg = re.sub(r'url\(#([^)]+)\)', lambda m: "url(#%s%s)" % (pre, m.group(1)), svg)
        svg = re.sub(r'href="#([^"]+)"', lambda m: 'href="#%s%s"' % (pre, m.group(1)), svg)
        # scale to the cell width, keep the aspect ratio
        def fix(m):
            tag = re.sub(r'\s(width|height)="[^"]*"', "", m.group(0))
            return tag[:-1] + ' width="100%" preserveAspectRatio="xMinYMin meet">'
        svg = re.sub(r"<svg\b[^>]*>", fix, svg, count=1)
        frame = "border:1px solid #ddd;" if border else ""
        return '<div style="width:%s;%sbackground:white">%s</div>' % (width, frame, svg)


def clean_svg(path, rename=None, font="Arial", margin=4, gap=10, min_gap=20,
              keep_gaps=False, keep_backgrounds=False, keep_borders=False, keep_ideogram=False,
              keep_span=False, tidy_names=False, keep_ext=False, compact=True, crop=True,
              drop_ruler=False, drop_locus=False, drop_tracks=None, simplify=True):
    """Clean an IGV desktop SVG and return a CleanedSVG (displays inline in Jupyter).

    rename: optional CSV/TSV path, dict {old: new}, or list of (old, new) pairs
            (substring replace in track names).
    tidy_names: opt-in - join IGV's wrapped track names and drop .bw/.bam/...
    The other arguments mirror the command-line flags.

        from minigv import clean_svg
        fig = clean_svg("snapshot.svg", drop_ruler=True)
        fig                      # preview
        fig.compare()            # IGV original vs cleaned
        fig.save("figure.svg")
    """
    if isinstance(rename, str):
        pairs = read_rename(rename)
    elif isinstance(rename, dict):
        pairs = sorted(rename.items(), key=lambda p: -len(p[0]))
    else:
        pairs = sorted(rename or [], key=lambda p: -len(p[0]))
    args = argparse.Namespace(
        font=font, margin=margin, gap=gap, min_gap=min_gap, keep_gaps=keep_gaps,
        keep_backgrounds=keep_backgrounds, keep_borders=keep_borders, keep_ideogram=keep_ideogram,
        keep_span=keep_span, keep_ext=keep_ext, tidy_names=tidy_names, no_compact=not compact,
        no_crop=not crop, drop_ruler=drop_ruler, drop_locus=drop_locus, drop_tracks=list(drop_tracks or []),
        no_simplify=not simplify)
    tree, c = _run(path, args, pairs)
    svg = ET.tostring(tree.getroot(), encoding="unicode")
    return CleanedSVG(path, '<?xml version="1.0" encoding="utf-8"?>\n' + svg, c.report, c.canvas)


def list_tracks(path):
    """Track names in IGV's order (use them with drop_tracks= / --drop-track)."""
    c = Cleaner(ET.parse(path), argparse.Namespace())
    return [t["name"] for t in c.tracks]


def call_line(path, out=None, **kw):
    """One-line Python call (only non-default options) for these settings."""
    import inspect
    defaults = {k: v.default for k, v in inspect.signature(clean_svg).parameters.items() if k != "path"}
    args = [repr(path)] + ["%s=%r" % (k, v) for k, v in kw.items()
                           if k in defaults and v != defaults[k] and not (k == "drop_tracks" and not v)]
    call = "clean_svg(%s)" % ", ".join(args)
    return call + (".save(%r)" % out if out else "")


def cli_line(path, out=None, **kw):
    """Equivalent command line for these settings."""
    import shlex
    flags = {"keep_backgrounds": "--keep-backgrounds", "keep_borders": "--keep-borders",
             "keep_ideogram": "--keep-ideogram", "keep_span": "--keep-span", "keep_gaps": "--keep-gaps",
             "drop_ruler": "--drop-ruler", "drop_locus": "--drop-locus", "tidy_names": "--tidy-names",
             "keep_ext": "--keep-ext"}
    parts = ["minigv", shlex.quote(path)]
    if out:
        parts += ["-o", shlex.quote(out)]
    for k, f in flags.items():
        if kw.get(k):
            parts.append(f)
    if kw.get("compact") is False:
        parts.append("--no-compact")
    if kw.get("crop") is False:
        parts.append("--no-crop")
    if kw.get("simplify") is False:
        parts.append("--no-simplify")
    if kw.get("font", "Arial") != "Arial":
        parts += ["--font", shlex.quote(kw["font"])]
    for k, f, d in (("gap", "--gap", 10), ("min_gap", "--min-gap", 20), ("margin", "--margin", 4)):
        if k in kw and kw[k] != d:
            parts += [f, "%g" % kw[k]]
    for t in kw.get("drop_tracks") or []:
        parts += ["--drop-track", shlex.quote(t)]
    return " ".join(parts)


def interactive(path, width="100%"):
    """Jupyter UI: tick boxes for each cleaning step and each track, live preview, save.
    Below the preview it prints the clean_svg(...) call (and command line) that
    reproduces the current settings, for scripting.

        from minigv import interactive
        interactive("snapshot.svg")
    """
    import html as _html
    import ipywidgets as W
    from IPython.display import HTML, display

    tracks = list_tracks(path)
    step = lambda d, v=True: W.Checkbox(value=v, description=d, indent=False, layout=W.Layout(width="auto"))
    show = {
        "backgrounds": step("grey backgrounds", False),
        "borders": step("panel borders / separators", False),
        "ideogram": step("chromosome ideogram", False),
        "span": step('span bar ("249 kb")', False),
        "ruler": step("coordinate ruler"),
        "locus": step("locus label (chr:start-end)"),
    }
    layout_opts = {
        "compact": step("move locus label down"),
        "keep_gaps": step("keep IGV's empty space", False),
        "crop": step("crop page to content"),
    }
    font = W.Dropdown(options=["Arial", "Helvetica", "keep"], value="Arial", description="font", layout=W.Layout(width="200px"))
    gap = W.FloatSlider(value=10, min=0, max=40, step=1, description="gap (pt)", continuous_update=False)
    track_boxes = [W.Checkbox(value=True, description=t, indent=False, layout=W.Layout(width="auto")) for t in tracks]
    preview_w = W.IntSlider(value=100, min=30, max=100, step=5, description="preview %", continuous_update=False)
    out_path = W.Text(value=os.path.splitext(path)[0] + ".clean.svg", description="save to", layout=W.Layout(width="70%"))
    save_btn = W.Button(description="Save SVG", button_style="success", icon="save")
    status = W.HTML()
    code = W.HTML()
    preview = W.Output(layout=W.Layout(border="1px solid #ddd"))
    state = {}

    def label(text):
        return W.HTML("<b>%s</b>" % text)

    def settings():
        return dict(
            keep_backgrounds=show["backgrounds"].value, keep_borders=show["borders"].value,
            keep_ideogram=show["ideogram"].value, keep_span=show["span"].value,
            drop_ruler=not show["ruler"].value, drop_locus=not show["locus"].value,
            compact=layout_opts["compact"].value, keep_gaps=layout_opts["keep_gaps"].value,
            crop=layout_opts["crop"].value, font=font.value, gap=gap.value,
            drop_tracks=[b.description for b in track_boxes if not b.value])

    def show_code(*_):
        kw = settings()
        py = "from minigv import clean_svg; " + call_line(path, out_path.value, **kw)
        sh = cli_line(path, out_path.value, **kw)
        box = ("<div style='font:12px monospace;background:#f6f8fa;border:1px solid #ddd;padding:6px 8px;"
               "white-space:pre-wrap;word-break:break-all;user-select:all'>%s</div>")
        code.value = ("<b>Python</b> (same settings, for scripts)" + box % _html.escape(py) +
                      "<b>Command line</b>" + box % _html.escape(sh))

    def render(*_):
        status.value = "<i>cleaning…</i>"
        show_code()
        try:
            fig = clean_svg(path, **settings())
        except Exception as e:                       # keep the UI alive and show what went wrong
            status.value = "<span style='color:#b00'>error: %s</span>" % e
            return
        state["fig"] = fig
        redraw()
        status.value = "page %s" % fig.report.get("page set to", "unchanged")

    def redraw(*_):
        if "fig" in state:
            preview.clear_output(wait=True)
            with preview:
                display(HTML(state["fig"]._html(state["fig"].svg, "%d%%" % preview_w.value)))

    def save(_):
        if "fig" in state:
            try:
                state["fig"].save(out_path.value)
                status.value = "saved → %s" % out_path.value
            except Exception as e:
                status.value = "<span style='color:#b00'>%s</span>" % e

    for w in list(show.values()) + list(layout_opts.values()) + track_boxes + [font, gap]:
        w.observe(render, names="value")
    preview_w.observe(redraw, names="value")
    out_path.observe(show_code, names="value")
    save_btn.on_click(save)

    controls = W.HBox([W.VBox([label("show")] + list(show.values())),
                       W.VBox([label("layout")] + list(layout_opts.values()) + [font, gap]),
                       W.VBox([label("tracks")] + track_boxes)],
                      layout=W.Layout(gap="24px", align_items="flex-start"))
    ui = W.VBox([controls, W.HBox([preview_w, status]), preview, W.HBox([out_path, save_btn]), code])
    render()
    display(ui)


def main(argv=None):
    p = argparse.ArgumentParser(description="minIGV: remove IGV screen chrome from IGV desktop SVG snapshots for figures.",
                                formatter_class=argparse.RawDescriptionHelpFormatter, epilog=__doc__.split("Usage")[1])
    p.add_argument("--version", action="version", version="minIGV " + __version__)
    p.add_argument("svg", nargs="+", help="IGV desktop SVG file(s)")
    p.add_argument("-o", "--output", help="output file (single input) or directory (several inputs); "
                                          "default: <name>.clean.svg next to the input")
    p.add_argument("--rename", metavar="FILE", help="optional CSV/TSV with old,new track-name pairs (substring replace)")
    p.add_argument("--font", default="Arial", help="font replacing Java's Dialog/SansSerif (default Arial; 'keep' to leave)")
    p.add_argument("--margin", type=float, default=4, help="margin around content when cropping, in pt (default 4)")
    p.add_argument("--gap", type=float, default=10, help="blank space left where an empty band is closed, in pt (default 10)")
    p.add_argument("--min-gap", type=float, default=20, help="only close blank bands taller than this, in pt (default 20)")
    p.add_argument("--keep-gaps", action="store_true", help="don't close empty space inside IGV panels")
    p.add_argument("--keep-backgrounds", action="store_true")
    p.add_argument("--keep-borders", action="store_true")
    p.add_argument("--keep-ideogram", action="store_true")
    p.add_argument("--keep-span", action="store_true", help='keep the "<N> kb" span bar')
    p.add_argument("--tidy-names", action="store_true",
                   help="join IGV's wrapped track names and drop .bw/.bam/... (off by default)")
    p.add_argument("--keep-ext", action="store_true", help="with --tidy-names: keep .bw/.bam/...")
    p.add_argument("--no-compact", action="store_true", help="don't move the locus label down")
    p.add_argument("--no-crop", action="store_true")
    p.add_argument("--no-simplify", action="store_true",
                   help="keep IGV's one-rect-per-bar structure (default merges bars into paths: same look, far fewer objects)")
    p.add_argument("--drop-ruler", action="store_true", help="also remove the coordinate ruler")
    p.add_argument("--drop-locus", action="store_true", help="also remove the chr:start-end label")
    p.add_argument("--drop-track", action="append", dest="drop_tracks", metavar="NAME",
                   help="remove a whole track by its IGV name (repeatable); see --list-tracks")
    p.add_argument("--list-tracks", action="store_true", help="print the track names IGV wrote and exit")
    p.add_argument("--report", action="store_true", help="print what was changed")
    args = p.parse_args(argv)

    if args.list_tracks:
        for path in args.svg:
            print(path)
            for name in list_tracks(path):
                print("   ", name)
        return 0
    rename = read_rename(args.rename) if args.rename else []
    many = len(args.svg) > 1
    if many and args.output and not os.path.isdir(args.output):
        p.error("with several inputs, -o must be an existing directory")
    for path in args.svg:
        stem = os.path.splitext(os.path.basename(path))[0]
        if args.output and not many:
            out = args.output
        elif args.output:
            out = os.path.join(args.output, stem + ".clean.svg")
        else:
            out = os.path.join(os.path.dirname(path) or ".", stem + ".clean.svg")
        if os.path.abspath(out) == os.path.abspath(path):
            p.error("refusing to overwrite the input: %s" % path)
        rep = clean(path, out, args, rename)
        print("%s -> %s" % (path, out))
        if args.report:
            for k, v in rep.items():
                print("    %-30s %s" % (k, v))


if __name__ == "__main__":
    sys.exit(main())

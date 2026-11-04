#!/usr/bin/env python3
"""Render a FOSSIL figure scene as a static SVG (issue #47).

Usage: mesh2svg.py SCENE.fig OUT.svg [--theme dark|paper] [--edges none|faint|strong]

The scene is written by a docs example program (marker ``!fig NAME`` in docs/examples/src/*.f90): FOSSIL computes
every coordinate, this script only projects, sorts, shades and draws. Text format, one record per line, ``#`` starts
a comment:

  size W H                      canvas in px (default 640 400)
  camera AZ EL                  view direction in degrees: azimuth about +z from +x, elevation above the xy plane
                                (default -60 25); the scene is fitted to the canvas, orthographic
  zoom X Y Z R                  frame the sphere of radius R about (X, Y, Z) instead of the whole scene; what falls
                                outside is clipped to the card
  caption TEXT                  a caption band under the drawing
  mesh STYLE N                  then N lines "x1 y1 z1 x2 y2 z2 x3 y3 z3 [v]": triangles, flat-shaded (two-sided
                                Lambert); with the styles field (v in [-1, 1], diverging: teal, cream, ember) and frac
                                (v in [0, 1], sequential: cream to sepia) the colour is v's, and the triangles of one
                                cell show no seam (cells are written as two triangles)
  segs STYLE N                  then N lines "x1 y1 z1 x2 y2 z2": segments
  pts STYLE N                   then N lines "x y z": points
  label STYLE X Y Z TEXT        a label anchored at a point, drawn on top

STYLE is a role of the palette below (body, ghost, patch, hole, probe, foot, inside, outside, box, accent, cell), or
field / frac for a mesh coloured by a value.
Meshes, segments and points are depth-sorted together (painter's algorithm, by the centroid of each primitive), so the
overlays are hidden by the body where they pass behind it; a "ghost" mesh is translucent. The same scene gives the same
bytes (fixed rounding, stable sort).
"""

from __future__ import annotations

import math
import sys

THEMES = {
    # background, ink (text), caption band, edge colour
    "dark": {"bg": "#1b1714", "ink": "#efe6da", "band": "#231e1a", "edge": "#120f0d"},
    "paper": {"bg": "#fdfaf5", "ink": "#3b2414", "band": "#f6efe4", "edge": "#6e4626"},
}
# role: (base colour, kind-specific options)
ROLES = {
    "body": "#c9a27a",      # sepia stone
    "ghost": "#9c8b7a",     # the body before an operation
    "patch": "#e8590c",     # what an operation added (ember)
    "hole": "#ff6b3d",      # boundary edges
    "probe": "#4fb3bf",     # query points
    "foot": "#f2d16b",      # closest points
    "inside": "#4fb3bf",
    "outside": "#e8590c",
    "box": "#8c7b6b",       # bounding boxes, grids
    "accent": "#f2d16b",
    "cell": "#8a5a33",
}
ON_SURFACE = ("hole", "foot")  # segments lying on a mesh, and every point: a depth bias keeps them over its facets
RANK = {"mesh": 0, "segs": 1, "pts": 2}
LIGHT = (-0.45, 0.55, 0.70)  # in view space: from the upper left, towards the viewer


def colormap(style: str, v: float):
    """field: v in [-1, 1], teal - cream - ember; frac: v in [0, 1], cream - sepia."""
    def mix(a, b, t):
        a, b = hexrgb(a), hexrgb(b)
        return tuple(a[i] + (b[i] - a[i]) * t for i in range(3))
    if style == "frac":
        return mix("#f3e3c9", "#8a5a33", max(0.0, min(1.0, v)))
    v = max(-1.0, min(1.0, v))
    return mix("#efe0c4", "#2f8f9d", -v) if v < 0 else mix("#efe0c4", "#d9480f", v)


def hexrgb(c: str):
    return tuple(int(c[i:i + 2], 16) for i in (1, 3, 5))


def rgbhex(r, g, b) -> str:
    return "#%02x%02x%02x" % tuple(max(0, min(255, round(v))) for v in (r, g, b))


def read(path: str):
    scene = {"size": (640, 400), "camera": (-60.0, 25.0), "caption": "", "items": [], "labels": [], "zoom": None}
    lines = [ln.split("#", 1)[0].strip() if not ln.lstrip().startswith("label") and not ln.lstrip().startswith(
        "caption") else ln.strip() for ln in open(path)]
    lines = [ln for ln in lines if ln]
    i = 0
    while i < len(lines):
        tok = lines[i].split()
        key = tok[0]
        if key == "size":
            scene["size"] = (int(tok[1]), int(tok[2]))
        elif key == "camera":
            scene["camera"] = (float(tok[1]), float(tok[2]))
        elif key == "zoom":
            scene["zoom"] = tuple(float(x) for x in tok[1:5])
        elif key == "caption":
            scene["caption"] = lines[i].split(None, 1)[1] if len(tok) > 1 else ""
        elif key in ("mesh", "segs", "pts"):
            style, n = tok[1], int(tok[2])
            width = {"mesh": 9, "segs": 6, "pts": 3}[key]
            for row in lines[i + 1:i + 1 + n]:
                cols = row.split()
                v = [float(x) for x in cols[:width]]
                pts = [tuple(v[k:k + 3]) for k in range(0, width, 3)]
                val = float(cols[width]) if len(cols) > width else 0.0
                scene["items"].append((key, style, pts, val))
            i += n
        elif key == "label":
            style = tok[1]
            xyz = tuple(float(x) for x in tok[2:5])
            text = lines[i].split(None, 5)[5] if len(tok) > 5 else ""
            scene["labels"].append((style, xyz, text))
        else:
            raise SystemExit(f"mesh2svg: unknown record '{key}' in {path}")
        i += 1
    return scene


def basis(az: float, el: float):
    a, e = math.radians(az), math.radians(el)
    d = (math.cos(e) * math.cos(a), math.cos(e) * math.sin(a), math.sin(e))  # from the scene towards the viewer
    right = (-math.sin(a), math.cos(a), 0.0)
    up = (d[1] * right[2] - d[2] * right[1], d[2] * right[0] - d[0] * right[2], d[0] * right[1] - d[1] * right[0])
    return right, up, d


def dot(u, v):
    return u[0] * v[0] + u[1] * v[1] + u[2] * v[2]


def sub(u, v):
    return (u[0] - v[0], u[1] - v[1], u[2] - v[2])


def cross(u, v):
    return (u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2], u[0] * v[1] - u[1] * v[0])


def f1(v: float) -> str:
    s = f"{v:.1f}"
    return "0" if s in ("0.0", "-0.0") else s.rstrip("0").rstrip(".")


def render(scene, theme: str, edges: str) -> str:
    th = THEMES[theme]
    w, h = scene["size"]
    band = 26 if scene["caption"] else 0
    pad = 18
    right, up, view = basis(*scene["camera"])
    light = LIGHT
    ln = math.sqrt(dot(light, light))
    light = tuple(c / ln for c in light)

    def proj(p):
        return dot(p, right), dot(p, up), dot(p, view)

    projected = [(kind, style, [proj(p) for p in pts], pts, val) for kind, style, pts, val in scene["items"]]
    xs = [q[0] for _, _, qs, _, _ in projected for q in qs] + [proj(l[1])[0] for l in scene["labels"]]
    ys = [q[1] for _, _, qs, _, _ in projected for q in qs] + [proj(l[1])[1] for l in scene["labels"]]
    if not xs:
        raise SystemExit("mesh2svg: empty scene")
    x0, x1, y0, y1 = min(xs), max(xs), min(ys), max(ys)
    if scene["zoom"]:
        zx, zy, _ = proj(scene["zoom"][:3])
        r = scene["zoom"][3]
        x0, x1, y0, y1 = zx - r, zx + r, zy - r, zy + r
        projected = [it for it in projected  # the primitives whose box overlaps the view (a long segment across it too)
                     if min(q[0] for q in it[2]) <= x1 + r and max(q[0] for q in it[2]) >= x0 - r
                     and min(q[1] for q in it[2]) <= y1 + r and max(q[1] for q in it[2]) >= y0 - r]
    scale = min((w - 2 * pad) / max(x1 - x0, 1e-30), (h - band - 2 * pad) / max(y1 - y0, 1e-30))
    ox = (w - scale * (x1 - x0)) / 2 - scale * x0
    oy = (h - band + scale * (y1 - y0)) / 2 + scale * y0

    def px(q):
        return ox + scale * q[0], oy - scale * q[1]

    zs = [q[2] for _, _, qs, _, _ in projected for q in qs]
    bias = 0.01 * (max(zs) - min(zs))  # primitives lying on the surface (ON_SURFACE) are drawn over their facets

    def depth(k):
        kind, style, qs, _, _ = projected[k]
        z = sum(q[2] for q in qs) / len(qs)
        z += bias if kind == "pts" or (kind == "segs" and style in ON_SURFACE) else 0.0
        return (z, RANK[kind], k)  # ties (a view from above: one depth) put segments, then points, over the meshes

    order = sorted(range(len(projected)), key=depth)
    out = []
    for k in order:
        kind, style, qs, pts, val = projected[k]
        base = colormap(style, val) if style in ("field", "frac") else hexrgb(ROLES.get(style, ROLES["body"]))
        if kind == "mesh":
            n = cross(sub(pts[1], pts[0]), sub(pts[2], pts[0]))
            nn = math.sqrt(dot(n, n))
            if nn == 0.0:
                continue
            nv = (dot(n, right) / nn, dot(n, up) / nn, dot(n, view) / nn)
            lam = abs(dot(nv, light))
            shade = 0.30 + 0.70 * lam
            col = rgbhex(*(c * shade for c in base))
            p = [px(q) for q in qs]
            d = "M" + "L".join(f"{f1(a)} {f1(b)}" for a, b in p) + "Z"
            if style == "ghost":
                out.append(("ghost", f'<path d="{d}" fill="{col}"/>'))
            elif style in ("field", "frac"):
                out.append(("cell", f'<path d="{d}" fill="{col}" stroke="{col}"/>'))
            else:
                out.append(("solid", f'<path d="{d}" fill="{col}"' + (f' stroke="{col}"' if edges == "none" else "")
                            + "/>"))
        elif kind == "segs":
            (a, b), (c, e) = px(qs[0]), px(qs[1])
            wdt = "2.2" if style == "hole" else "1.4"
            out.append(("free", f'<path d="M{f1(a)} {f1(b)}L{f1(c)} {f1(e)}" stroke="{rgbhex(*base)}" '
                                f'stroke-width="{wdt}" stroke-linecap="round"/>'))
        else:
            a, b = px(qs[0])
            out.append(("free", f'<circle cx="{f1(a)}" cy="{f1(b)}" r="3" fill="{rgbhex(*base)}" '
                                f'stroke="{th["bg"]}" stroke-width="1"/>'))
    labels = []
    for style, xyz, text in scene["labels"]:
        a, b = px(proj(xyz))
        labels.append(f'<text x="{f1(a + 6)}" y="{f1(b - 6)}" fill="{ROLES.get(style, th["ink"])}" '
                   f'font-family="ui-monospace,SFMono-Regular,Menlo,monospace" font-size="12" paint-order="stroke" '
                   f'stroke="{th["bg"]}" stroke-width="3">{esc(text)}</text>')
    cap = ""
    if band:
        cap = (f'<text x="{pad}" y="{h - 9}" fill="{th["ink"]}" font-family="ui-monospace,SFMono-Regular,Menlo,'
               f'monospace" font-size="13">{esc(scene["caption"])}</text>')
    return (f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {w} {h}" width="{w}" height="{h}" role="img" '
            f'aria-label="{esc(scene["caption"])}"><rect width="{w}" height="{h}" rx="10" fill="{th["bg"]}"/>'
            f'<clipPath id="card"><rect width="{w}" height="{h - band}" rx="10"/></clipPath>'
            f'<g clip-path="url(#card)">{groups(out, th, edges)}{"".join(labels)}</g>{cap}</svg>\n')


def groups(out, th, edges) -> str:
    """Runs of consecutive primitives of one kind in one <g> carrying their shared attributes (the depth order stays)."""
    attrs = {
        "solid": {"none": 'stroke-width="0.6" stroke-linejoin="round"',
                  "faint": f'stroke="{th["edge"]}" stroke-opacity="0.35" stroke-width="0.4" stroke-linejoin="round"',
                  "strong": f'stroke="{th["edge"]}" stroke-opacity="0.8" stroke-width="0.6" stroke-linejoin="round"'}[edges],
        "ghost": 'fill-opacity="0.18" stroke="#9c8b7a" stroke-opacity="0.25" stroke-width="0.4"',
        "free": 'fill="none"',
        "cell": 'stroke-width="0.9" stroke-linejoin="round"',  # covers the seam of a cell's two triangles
    }
    res, run, kind = [], [], None
    for k, svg in out + [(None, "")]:
        if k != kind:
            if run:
                res.append(f'<g {attrs[kind]}>{"".join(run)}</g>')
            run, kind = [], k
        run.append(svg)
    return "".join(res)


def esc(s: str) -> str:
    return s.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;").replace('"', "&quot;")


def main(argv):
    args = [a for a in argv[1:] if not a.startswith("--")]
    opts = dict(a[2:].split("=", 1) for a in argv[1:] if a.startswith("--") and "=" in a)
    if len(args) != 2:
        raise SystemExit(__doc__)
    svg = render(read(args[0]), opts.get("theme", "dark"), opts.get("edges", "faint"))
    with open(args[1], "w") as fh:
        fh.write(svg)


if __name__ == "__main__":
    main(sys.argv)

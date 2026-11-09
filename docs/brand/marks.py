#!/usr/bin/env python3
"""FOSSIL brand marks: the ammonite in the Sumi hand, in FOSSIL's sepia.

The mark is drawn the way ADAM, FLUME, PRISM and llm-wiki are drawn (the "Sumi family"):

* the brush: every stroke is a filled outline around a centreline, tapered by a width profile and
  wavering slightly (seeded, so the output is byte-for-byte reproducible);
* the shell: one logarithmic spiral growing 2.15 times per turn, drawn as a single stroke that starts
  thin at the core and swells to the aperture; nine septa (chamber walls) cross the last whorl;
* the seal: a 10 x 10 square in the upper right, carved with a refined quadrant;
* the ground: a balanced quadtree (16 -> 8 -> 4) drawn in the mesh colour, refined to the finest
  level exactly under the ink. It is never decoration: it is regenerated from the strokes.

Writes, in the output directory (default docs/public):
  logo.svg          the full mark for light pages (Limestone): ink, mesh and seal on paper
  logo-dark.svg     the full mark for dark pages (Shale)
  favicon.svg       below 48 px: no mesh, no seal, thicker strokes, on a rounded Limestone tile
  favicon-dark.svg  the same on a Shale tile (the VitePress nav logo in the dark theme)
  icons/<slug>.svg, icons/<slug>-dark.svg
                    the 18 feature-card marks of the landing page (docs/index.md), Limestone and Shale ink on a
                    transparent ground, drawn for 48 px with the same brush: 2-6 card-size strokes, no mesh, no seal,
                    the seal colour only as an accent that carries the meaning (the kintsugi crack, the query point,
                    the contact, the die). Slugs: stl, repair, processing, analysis, distance, level-sets, curves,
                    ib, amr, cut-cells, cut-solver, moving, collision, blocks, overset, body-fitted, gpu, portable.

Use logo*.svg at 48 px and above, favicon*.svg below. Never hand-edit the SVGs: change this file
and regenerate. PNG favicons (16, 32, 180 px) are rasterised from favicon.svg, e.g.
  rsvg-convert -w 32 -h 32 docs/public/favicon.svg -o docs/public/favicon-32.png

Usage: python3 docs/brand/marks.py [output directory, default docs/public]
"""
import math
import os
import random
import sys

# Limestone (light) and Shale (dark): FOSSIL's sepia in the Sumi roles.
PALETTE = {
    'light': dict(ground='#fdfaf5', ink='#3b2414', mesh='#806a52', seal='#8a5a33'),
    'dark': dict(ground='#1b1714', ink='#f2ede6', mesh='#8f7d69', seal='#e2c29b'),
}
LABEL = 'FOSSIL logo: an ammonite in one brush stroke, with its septa, on an AMR ground'

SEAL = (50, 4, 10)            # x, y, side of the seal (viewBox 0 0 64 64)
FAVICON_THICKEN = 1.7         # stroke width factor below 48 px


# ------------------------------------------------------------------------------------------ brush
def resample(pts, step=0.45):
    """The polyline resampled at a uniform arc-length step."""
    out = [pts[0]]
    acc = 0.0
    for (x0, y0), (x1, y1) in zip(pts, pts[1:]):
        d = math.hypot(x1 - x0, y1 - y0)
        if d == 0:
            continue
        t = step - acc
        while t <= d:
            out.append((x0 + (x1 - x0) * t / d, y0 + (y1 - y0) * t / d))
            t += step
        acc = d - (t - step)
    if math.hypot(out[-1][0] - pts[-1][0], out[-1][1] - pts[-1][1]) > 1e-6:
        out.append(pts[-1])
    return out


def profile(kind, t):
    """Width multiplier along the stroke, t in [0, 1]."""
    if kind == 'load':        # loaded press, long release
        a = min(1, t / 0.12)
        return math.sin(a * math.pi / 2) * (0.25 + 0.75 * (1 - t) ** 0.85) + 0.03
    if kind == 'grow':        # thin start, swelling, quick lift
        b = min(1, (1 - t) / 0.10)
        return (0.12 + 0.88 * t ** 0.9) * math.sin(b * math.pi / 2) + 0.03
    if kind == 'press':       # even line: set down, a slow release, lifted (lines and loops of the card marks)
        a = min(1, t / 0.10)
        b = min(1, (1 - t) / 0.18)
        return (0.42 + 0.58 * math.sin(a * math.pi / 2) * math.sin(b * math.pi / 2)) * (1 - 0.18 * t)
    raise ValueError(kind)


class Stroke:
    def __init__(self, pts, width, kind, waver, seed, step=0.45, digits=2):
        self.c = resample(pts, step)
        self.digits = digits
        self.w, self.kind, self.waver, self.seed = width, kind, waver, seed

    def outline(self, scale=1.0):
        """SVG path data of the filled stroke."""
        c, n = self.c, len(self.c)
        rnd = random.Random(self.seed)
        ph = [rnd.uniform(0, 2 * math.pi) for _ in range(4)]
        left, right = [], []
        for i, (x, y) in enumerate(c):
            j0, j1 = max(0, i - 1), min(n - 1, i + 1)
            tx, ty = c[j1][0] - c[j0][0], c[j1][1] - c[j0][1]
            d = math.hypot(tx, ty) or 1
            nx, ny = -ty / d, tx / d
            t = i / (n - 1)
            hw = 0.5 * self.w * scale * profile(self.kind, t)
            hl = hw * (1 + self.waver * (0.6 * math.sin(11 * t + ph[0]) + 0.4 * math.sin(29 * t + ph[1])))
            hr = hw * (1 + self.waver * 0.5 * math.sin(17 * t + ph[3]))
            off = self.waver * 0.35 * self.w * math.sin(7 * t + ph[2]) * 0.3
            left.append((x + nx * (hl + off), y + ny * (hl + off)))
            right.append((x - nx * (hr - off), y - ny * (hr - off)))
        return path(left + right[::-1], self.digits)


def path(pts, digits=2):
    """SVG path data of a closed outline."""
    return 'M' + ' L'.join(f'{x:.{digits}f} {y:.{digits}f}' for x, y in pts) + ' Z'


# --------------------------------------------------------------------------------------- the shell
def ammonite():
    """The spiral and its septa, as strokes."""
    cx, cy = 30.5, 35.5
    b = math.log(2.15) / (2 * math.pi)        # whorl expansion per turn
    rmax = 22.5
    th_end = 3.1 * 2 * math.pi
    a = rmax / math.exp(b * th_end)
    th0 = 0.55 * 2 * math.pi

    def r(th):
        return a * math.exp(b * th)

    def ang(th):                               # counter-clockwise on screen, aperture lower left
        return -th + math.radians(200)

    spiral = [(cx + r(t) * math.cos(ang(t)), cy + r(t) * math.sin(ang(t)))
              for t in [th0 + (th_end - th0) * i / 400 for i in range(401)]]
    strokes = [Stroke(spiral, 4.6, 'grow', waver=0.10, seed=3)]
    for k in range(9):                         # septa across the last whorl
        t = th_end - 0.32 - k * 0.52
        if t - 2 * math.pi < th0:
            break
        r_out, r_in = r(t) * 0.92, r(t - 2 * math.pi) * 1.18
        if r_out - r_in < 1.5:
            continue
        pts = []
        for i in range(21):
            u = i / 20
            rr = r_out + (r_in - r_out) * u
            aa = ang(t) + 0.10 * math.sin(math.pi * u)
            pts.append((cx + rr * math.cos(aa), cy + rr * math.sin(aa)))
        strokes.append(Stroke(pts, 1.25 + 0.10 * (8 - k), 'load', waver=0.12, seed=20 + k))
    return strokes


# -------------------------------------------------------------------------------------- the ground
def quadtree(strokes):
    """Leaves (x, y, h) of a balanced quadtree refined 16 -> 8 -> 4 under the ink."""
    pts = [p for s in strokes for p in s.c]

    def touched(x, y, h, pad):
        return any(x - pad <= px <= x + h + pad and y - pad <= py <= y + h + pad for px, py in pts)

    def split(x, y, h):
        hh = h // 2
        return {(x, y, hh), (x + hh, y, hh), (x, y + hh, hh), (x + hh, y + hh, hh)}

    leaves = {(x, y, 16) for x in range(0, 64, 16) for y in range(0, 64, 16)}
    for target, pad in ((8, 0.4), (4, 0.0)):
        new = set()
        for (x, y, h) in leaves:
            new |= split(x, y, h) if h > target and touched(x, y, h, pad) else {(x, y, h)}
        leaves = new

    def leaf_at(px, py):
        for (x, y, h) in leaves:
            if x <= px < x + h and y <= py < y + h:
                return (x, y, h)
        return None

    changed = True                             # 2:1 balance across faces
    while changed:
        changed = False
        new = set(leaves)
        for (x, y, h) in list(leaves):
            for px, py in ((x - 0.1, y + h / 2), (x + h + 0.1, y + h / 2), (x + h / 2, y - 0.1), (x + h / 2, y + h + 0.1)):
                if 0 <= px < 64 and 0 <= py < 64:
                    lf = leaf_at(px, py)
                    if lf and lf[2] > 2 * h and lf in new:
                        new.discard(lf)
                        new |= split(*lf)
                        changed = True
        leaves = new
    return leaves


def ground(strokes):
    """The quadtree as lines; an edge shared by two levels keeps the coarser (thicker) width."""
    width = {16: 0.38, 8: 0.22, 4: 0.14}
    segs = {}
    for (x, y, h) in quadtree(strokes):
        for seg in (((x, y), (x + h, y)), ((x, y + h), (x + h, y + h)), ((x, y), (x, y + h)), ((x + h, y), (x + h, y + h))):
            segs[seg] = max(segs.get(seg, 0), width[h])
    return ''.join(f'<line x1="{x1}" y1="{y1}" x2="{x2}" y2="{y2}" stroke-width="{w}"/>'
                   for ((x1, y1), (x2, y2)), w in sorted(segs.items()))


# ---------------------------------------------------------------------------------------- the seal
def seal(fill, carve):
    x, y, s = SEAL
    out = [f'<rect x="{x}" y="{y}" width="{s}" height="{s}" rx="1.00" fill="{fill}"/>']
    for (a, b, w) in ((55.30, 9.30, 1.25), (57.15, 9.30, 1.25), (55.30, 11.15, 1.25),
                      (57.15, 11.15, .33), (58.07, 11.15, .33), (57.15, 12.08, .33), (58.07, 12.08, .33)):
        out.append(f'<rect x="{a:.2f}" y="{b:.2f}" width="{w}" height="{w}" fill="{carve}"/>')
    out.append(f'<rect x="51.60" y="5.60" width="3.10" height="3.10" fill="none" stroke="{carve}" stroke-width="0.70"/>')
    return ''.join(out)


# ------------------------------------------------------------------------------- the feature cards
# One mark per landing-page card (docs/index.md), drawn with the same brush for ~48 px: no mesh, no seal, card-size
# strokes, an accent in the seal colour where it carries the meaning. viewBox 0 0 64 64, transparent ground.
W_MAIN, W_MID, W_FINE = 5.6, 4.0, 2.9      # stroke widths of the card marks (viewBox units; 48 px = 0.75 px/unit)


class Dot:
    """A brush press: a filled, slightly wavering round blot of radius r (or an ellipse rx, ry)."""

    def __init__(self, cx, cy, r, seed, ry=None, waver=0.06):
        self.cx, self.cy, self.rx, self.ry, self.seed, self.waver = cx, cy, r, ry or r, seed, waver

    def outline(self, scale=1.0):
        rnd = random.Random(self.seed)
        ph = [rnd.uniform(0, 2 * math.pi) for _ in range(3)]
        pts = []
        for i in range(48):
            a = 2 * math.pi * i / 48
            f = 1 + self.waver * (0.6 * math.sin(3 * a + ph[0]) + 0.4 * math.sin(5 * a + ph[1]))
            pts.append((self.cx + self.rx * f * math.cos(a), self.cy + self.ry * f * math.sin(a)))
        return path(pts, 1)


class Blot:
    """A filled closed outline given by its points (the solid bodies of the card marks)."""

    def __init__(self, pts):
        self.pts = pts

    def outline(self, scale=1.0):
        return path(self.pts, 1)


def line(p, q, bow=0.0, n=24):
    """From p to q, bowed sideways by bow (a fraction of the length) like a hand-drawn line."""
    (x0, y0), (x1, y1) = p, q
    nx, ny = -(y1 - y0), x1 - x0
    return [(x0 + (x1 - x0) * u + nx * bow * 4 * u * (1 - u), y0 + (y1 - y0) * u + ny * bow * 4 * u * (1 - u))
            for u in (i / n for i in range(n + 1))]


def arc(cx, cy, rx, a0, a1, ry=None, rot=0.0, n=96):
    """An elliptic arc from angle a0 to a1 (degrees, screen: y down, so positive angles turn clockwise)."""
    ry = ry or rx
    c, s = math.cos(math.radians(rot)), math.sin(math.radians(rot))
    out = []
    for i in range(n + 1):
        a = math.radians(a0 + (a1 - a0) * i / n)
        x, y = rx * math.cos(a), ry * math.sin(a)
        out.append((cx + c * x - s * y, cy + s * x + c * y))
    return out


def bezier(p0, p1, p2, p3, n=64):
    out = []
    for i in range(n + 1):
        u = i / n
        a, b, c, d = (1 - u) ** 3, 3 * u * (1 - u) ** 2, 3 * u * u * (1 - u), u ** 3
        out.append((a * p0[0] + b * p1[0] + c * p2[0] + d * p3[0], a * p0[1] + b * p1[1] + c * p2[1] + d * p3[1]))
    return out


def rounded(pts, r):
    """A polyline with its interior corners filleted (radius r, at least the half width: the outline never folds)."""
    out = [pts[0]]
    for a, b, c in zip(pts, pts[1:], pts[2:]):
        u1 = (a[0] - b[0], a[1] - b[1])
        u2 = (c[0] - b[0], c[1] - b[1])
        l1, l2 = math.hypot(*u1), math.hypot(*u2)
        u1, u2 = (u1[0] / l1, u1[1] / l1), (u2[0] / l2, u2[1] / l2)
        th = math.acos(max(-1.0, min(1.0, u1[0] * u2[0] + u1[1] * u2[1])))
        t = min(r / math.tan(th / 2), 0.45 * l1, 0.45 * l2)
        t1, t2 = (b[0] + u1[0] * t, b[1] + u1[1] * t), (b[0] + u2[0] * t, b[1] + u2[1] * t)
        out += [((1 - u) ** 2 * t1[0] + 2 * u * (1 - u) * b[0] + u * u * t2[0],
                 (1 - u) ** 2 * t1[1] + 2 * u * (1 - u) * b[1] + u * u * t2[1]) for u in (i / 12 for i in range(13))]
    out.append(pts[-1])
    return out


def square(x0, y0, x1, y1, r=2.6, rot=0.0):
    """A square drawn in one stroke from the middle of its left side, round and back past the start."""
    cx, cy = (x0 + x1) / 2, (y0 + y1) / 2
    c, s = math.cos(math.radians(rot)), math.sin(math.radians(rot))
    corners = [(x0, cy), (x0, y0), (x1, y0), (x1, y1), (x0, y1), (x0, cy + 0.10 * (y0 - y1))]
    pts = [(cx + c * (x - cx) - s * (y - cy), cy + s * (x - cx) + c * (y - cy)) for x, y in corners]
    return rounded(pts, r)


def gap(p, q, g0, g1):
    """The straight stroke from p to q, stopped g0 short of p and g1 short of q."""
    d = math.hypot(q[0] - p[0], q[1] - p[1])
    u = ((q[0] - p[0]) / d, (q[1] - p[1]) / d)
    return line((p[0] + u[0] * g0, p[1] + u[1] * g0), (q[0] - u[0] * g1, q[1] - u[1] * g1))


def foot(curve, p):
    """The point of the sampled curve closest to p."""
    return min(curve, key=lambda q: (q[0] - p[0]) ** 2 + (q[1] - p[1]) ** 2)


def card_marks():
    """slug -> (alt text, ink strokes, accent strokes)."""
    m = {}
    seeds = iter(range(100, 10000))

    def st(pts, w, kind='press', waver=0.10, src=None):  # coarser outlines than the logo's: drawn for 48 px
        return Stroke(pts, w, kind, waver, next(src or seeds), step=0.8, digits=1)

    def own(base, used):
        """Seeds of a redrawn mark, from its own range; the shared ones its first drawing used are skipped, so the
        marks after it keep their bytes."""
        for _ in range(used):
            next(seeds)
        return iter(range(base, base + 100))

    # STL in, STL out: two facets sharing an edge
    b, d = (27, 10), (40, 55)                                    # the shared edge
    m['stl'] = ('two triangular facets sharing an edge', [
        st(line((6.6, 46.4), (28, 8.4), 0.02), W_MID, 'load'),
        st(line((25.5, 9.4), (57.5, 22.6), -0.03), W_MID, 'load'),
        st(line((57, 20.4), (39.4, 56.6), -0.02), W_MID, 'load'),
        st(line((41.6, 55.4), (6.4, 43.6), -0.03), W_MID, 'load'),
        st(line(b, d, 0.0), W_MAIN, 'load')], [])

    # Repair dirty meshes: kintsugi, a bowl mended along its crack in the seal colour
    m['repair'] = ('a bowl mended along its crack in the seal colour (kintsugi)', [
        st(line((6, 23), (58, 23), 0.01), W_MID, 'load'),
        st(arc(32, 23, 23, 180, 0, ry=27), W_MAIN, 'grow'),
        st(line((24, 55), (40, 55)), W_FINE, 'load')], [
        st(rounded([(37, 24), (31, 31), (37, 37), (30.5, 44), (34, 50.5)], 1.8), W_FINE + 0.3, 'press', 0.06)])

    # Geometry processing: a square and a circle overlapping (a boolean)
    m['processing'] = ('a square and a circle overlapping, as in a boolean', [
        st(square(8, 24, 40, 56), W_MID),
        st(arc(40, 24, 16, 205, 205 - 352), W_MAIN, 'load')], [])

    # Geometry analysis: a pebble and the dimension line that measures it
    m['analysis'] = ('a pebble measured by a dimension line', [
        st(arc(32, 25, 20, 165, 165 - 352, ry=13.5, rot=-8), W_MAIN, 'load'),
        st(line((14, 51), (50, 51), 0.01), W_MID, 'press'),
        st(line((11.5, 43), (11.5, 59)), W_FINE, 'load'),
        st(line((52.5, 43), (52.5, 59)), W_FINE, 'load')], [])

    # Signed distance: a query point (a ring), the perpendicular to its closest point on a wall, the right angle there
    sd = own(20000, 5)
    oc, orad = (6, 74), 42.0                                      # the wall: a convex arc
    wall = arc(*oc, orad, -101, -15)
    p = (53, 9)
    h = math.hypot(p[0] - oc[0], p[1] - oc[1])
    nx, ny = (p[0] - oc[0]) / h, (p[1] - oc[1]) / h               # the normal through p
    tx, ty = ny, -nx                                              # the tangent, toward the square's side
    b0 = (oc[0] + (orad + 2.7) * nx, oc[1] + (orad + 2.7) * ny)   # the foot, on the wall's ink
    q = 5.8
    m['distance'] = ('a query point as a ring, the perpendicular stroke to its closest point on a curved wall, and '
                     'the right-angle mark at the foot', [
        st(wall, W_MAIN, 'load', src=sd),
        st(gap(p, b0, 6.6, 0.0), W_FINE + 0.3, 'press', src=sd),
        st(rounded([(b0[0] + q * tx, b0[1] + q * ty), (b0[0] + q * (tx + nx), b0[1] + q * (ty + ny)),
                    (b0[0] + q * nx, b0[1] + q * ny)], 0.9), W_FINE - 0.5, 'press', 0.05, src=sd)], [
        st(arc(p[0], p[1], 4.5, 200, 200 - 356, n=48), W_FINE + 0.1, 'press', 0.05, src=sd)])

    # Level sets: karesansui, the raked gravel around a stone
    def rock(d, k=1.0):
        r0 = 7.5
        out = []
        for i in range(181):
            a = 2 * math.pi * i / 180 + 0.3
            amp = r0 / (r0 + d)
            r = r0 * (1 + amp * (0.16 * math.cos(2 * a + 0.4) + 0.09 * math.cos(3 * a + 1.3))) + d
            out.append((33 + r * math.cos(a), 33 + 0.86 * r * math.sin(a)))
        return out
    m['level-sets'] = ('a stone with the raked gravel round it in level lines (karesansui)', [
        Blot(rock(0)),
        st(rock(7.5)[8:190 - 14], W_FINE),
        st(rock(15)[70:] + rock(15)[1:58], W_FINE),
        st(rock(22.5)[130:] + rock(22.5)[1:118], W_FINE)], [])

    # 2-D curves: one closed loop, an aerofoil section
    def foil():
        out = []
        for i in range(121):
            u = i / 120                                          # trailing edge, upper side, nose, lower side, back
            x = 0.5 * (1 + math.cos(2 * math.pi * u))
            yt = 5 * 0.22 * (0.2969 * math.sqrt(x) - 0.1260 * x - 0.3516 * x ** 2 + 0.2843 * x ** 3 - 0.1036 * x ** 4)
            y = yt if u < 0.5 else -yt
            y += 0.05 * math.sin(math.pi * x)                    # camber
            X, Y = 52 * (1 - x) - 26, -52 * y                     # nose on the left
            c, s = math.cos(math.radians(-12)), math.sin(math.radians(-12))
            out.append((32 + c * X - s * Y, 33 + s * X + c * Y))
        return out
    m['curves'] = ('one closed brush loop, an aerofoil section', [st(foil(), W_MAIN, 'press')], [])

    # Immersed-boundary queries: strata cut by a body's boundary, one point on it
    body = arc(66, 33, 25, 108, 252)
    strata = []
    for k, (y, w) in enumerate(((17, 0.6), (27, -0.5), (37, 0.5), (47, -0.4))):
        x_end = 66 - math.sqrt(max(0.0, 25 ** 2 - (y - 33) ** 2)) - 3.2
        strata.append(st([(x, y + w * math.sin(0.25 * x + k)) for x in [6 + (x_end - 6) * i / 40 for i in range(41)]],
                         W_FINE + 0.3, 'load'))
    m['ib'] = ('layered strata cut by the boundary of a body, one point marked on it', [
        st(body, W_MAIN, 'press')] + strata, [Dot(66 - math.sqrt(25 ** 2 - 64), 25, 3.6, next(seeds))])

    # AMR tagging: a square refined toward a stroke
    m['amr'] = ('a square refined twice toward a brush stroke', [
        st(square(8, 8, 56, 56), W_MID),
        st(line((32, 8), (32, 56)), W_FINE),
        st(line((8, 32), (56, 32)), W_FINE),
        st(line((20, 32), (20, 56)), W_FINE - 0.4),
        st(line((8, 44), (32, 44)), W_FINE - 0.4),
        st(bezier((4, 34), (18, 38), (26, 48), (30, 60)), W_MAIN, 'load')], [])

    # Exact cut cells: a cell cut by a curve, the solid side hatched
    cut = bezier((5, 25), (22, 28), (36, 40), (40, 59))
    hatch = []
    for c in (50, 59, 68):                                       # lines x + y = c, from the cell edge to the curve
        x0 = max(10.0, c - 54.0)
        p0 = (x0 + 3.2, c - x0 - 3.2)
        ts = [i / 200 * 40 for i in range(201)]
        t_hit = min(ts, key=lambda t: min((q[0] - p0[0] - t) ** 2 + (q[1] - p0[1] + t) ** 2 for q in cut))
        t_end = t_hit - 3.4
        hatch.append(st(line(p0, (p0[0] + t_end, p0[1] - t_end)), W_FINE + 0.3, 'load'))
    m['cut-cells'] = ('a grid cell cut by a curve, its solid side hatched', [
        st(square(10, 10, 54, 54), W_MID),
        st(cut, W_MAIN, 'load')] + hatch, [])

    # Cut-cell solver support: a 2 x 2 patch of cells, a wall leaving a sliver in the corner cell, the sliver linked
    # across the shared face into its larger neighbour (cell merging)
    cs = own(20100, 8)
    g0, g1, g2 = 10, 32, 54                                       # grid lines
    wc, wr = (-13, 77), 53.5                                      # the wall: a flat arc about a point off the corner
    cin, win = 1.3, 2.6                                           # the sliver stays clear of the grid's and wall's ink
    ro = wr + win
    a0 = -math.degrees(math.asin((wc[1] - g1 - cin) / ro))        # offset wall meets y = g1 + cin
    a1 = -math.degrees(math.acos((g1 - cin - wc[0]) / ro))        # offset wall meets x = g1 - cin
    sliver = [(wc[0] + ro * math.cos(math.radians(a)), wc[1] + ro * math.sin(math.radians(a)))
              for a in (a0 + (a1 - a0) * i / 16 for i in range(17))] + [(g1 - cin, g1 + cin)]
    s0 = (sum(x for x, _ in sliver) / len(sliver), sum(y for _, y in sliver) / len(sliver))
    tip = (47, 41)
    d = math.hypot(tip[0] - s0[0], tip[1] - s0[1])
    ux, uy = (tip[0] - s0[0]) / d, (tip[1] - s0[1]) / d
    barbs = []
    for sgn in (1, -1):                                           # the arrowhead: two barbs 35 degrees off the shaft
        c, s_ = math.cos(math.radians(35)), sgn * math.sin(math.radians(35))
        bx, by = -(c * ux - s_ * uy), -(s_ * ux + c * uy)
        barbs.append(st(line((tip[0] + 7 * bx, tip[1] + 7 * by), tip), W_FINE, 'load', src=cs))
    m['cut-solver'] = ('a 2 by 2 patch of cells cut by a wall, the small sliver left in the corner cell linked by an '
                       'arrow across the shared face into its larger neighbour', [
        st(line((g0, 6), (g0, 58)), W_FINE, src=cs), st(line((g1, 6), (g1, 58)), W_FINE, src=cs),
        st(line((g2, 6), (g2, 58)), W_FINE, src=cs),
        st(line((6, g0), (58, g0)), W_FINE, src=cs), st(line((6, g1), (58, g1)), W_FINE, src=cs),
        st(line((6, g2), (58, g2)), W_FINE, src=cs),
        st(arc(*wc, wr, -72, -18), W_MAIN, 'load', src=cs),
        st(line(s0, (tip[0] - 1.0 * ux, tip[1] - 1.0 * uy), -0.06), W_MID, 'press', src=cs)] + barbs, [
        Blot(sliver)])

    # Moving bodies: a body and its trail
    trail = []
    for y, x0 in ((24, 14), (32, 5), (40, 14)):
        x1 = 44 - math.sqrt(max(0.0, 13 ** 2 - (y - 32) ** 2)) - 3.6
        trail.append(st(line((x0, y), (x1, y), 0.0), W_MID, 'grow'))
    m['moving'] = ('a round body with the trail of its motion', [
        st(arc(44, 32, 13, 200, 200 - 352), W_MAIN, 'load')] + trail, [])

    # Collision detection: two bodies just touching
    ca, ra, rb = (22, 37), 14.5, 11.0
    u = (math.cos(math.radians(-32)), math.sin(math.radians(-32)))
    cb = (ca[0] + (ra + rb + 1.0) * u[0], ca[1] + (ra + rb + 1.0) * u[1])
    cp = (ca[0] + (ra + 0.5) * u[0], ca[1] + (ra + 0.5) * u[1])
    sparks = []
    for ang in (-60, 0, 60):
        v = (math.cos(math.radians(-32 - 90 + ang)), math.sin(math.radians(-32 - 90 + ang)))
        sparks.append(st(line((cp[0] + 3.2 * v[0], cp[1] + 3.2 * v[1]), (cp[0] + 12.0 * v[0], cp[1] + 12.0 * v[1])),
                         W_MID, 'load'))
    m['collision'] = ('two round bodies just touching, the contact marked', [
        st(arc(*ca, ra, 150, 150 - 352), W_MAIN, 'load'),
        st(arc(*cb, rb, 330, 330 - 352), W_MID, 'load')], sparks)

    # Structured blocks: a bent grid
    def bent(uu, vv):
        return (10 + 44 * uu + 7 * math.sin(math.pi * vv), 54 - 44 * vv + 6 * math.sin(math.pi * uu))
    m['blocks'] = ('a bent grid of a few curvilinear lines', [
        st([bent(0, i / 40) for i in range(41)], W_MID, 'load'),
        st([bent(1, i / 40) for i in range(41)], W_MID, 'load'),
        st([bent(i / 40, 0) for i in range(41)], W_MID, 'load'),
        st([bent(i / 40, 1) for i in range(41)], W_MID, 'load'),
        st([bent(0.5, i / 40) for i in range(41)], W_FINE, 'press'),
        st([bent(i / 40, 0.5) for i in range(41)], W_FINE, 'press')], [])

    # Overset assembly: two grids overlapping
    def grid_lines(cx, cy, h, rot, w):
        c, s = math.cos(math.radians(rot)), math.sin(math.radians(rot))

        def tr(x, y):
            return (cx + c * x - s * y, cy + s * x + c * y)
        return [st(line(tr(0, -h), tr(0, h)), w), st(line(tr(-h, 0), tr(h, 0)), w)]
    m['overset'] = ('two square grids overlapping, one turned', [
        st(square(7, 23, 41, 57), W_MID)] + grid_lines(24, 40, 17, 0, W_FINE) + [
        st(square(29, 7, 57, 35, rot=24), W_MID)] + grid_lines(43, 21, 14, 24, W_FINE), [])

    # Body-fitted grids: an O-grid fitted to an elongated body: offset rings and lines along the body normals
    ea, eb, rot = 12.5, 5.0, math.radians(-22)

    def onorm(t, d):
        x, y = ea * math.cos(t), eb * math.sin(t)
        nx, ny = eb * math.cos(t), ea * math.sin(t)
        h = math.hypot(nx, ny)
        x, y = x + d * nx / h, y + d * ny / h
        return (32 + math.cos(rot) * x - math.sin(rot) * y, 32 + math.sin(rot) * x + math.cos(rot) * y)

    def oring(d, t0):
        return [onorm(t0 + 2 * math.pi * 0.978 * i / 120, d) for i in range(121)]
    m['body-fitted'] = ('an O-grid fitted to an elongated solid body: offset rings and lines along its normals', [
        Blot([onorm(2 * math.pi * i / 96, 0) for i in range(96)]),
        st(oring(8.5, 2.2), W_FINE, 'press'),
        st(oring(18.0, 5.3), W_MID, 'load')] + [
        st([onorm(t, d) for d in (2.6 + 18.5 * i / 20 for i in range(21))], W_FINE, 'press')
        for t in (0.0, 0.95, 2.2, math.pi, 4.1, 5.35)], [])

    # On the GPU, bit for bit: a chip with its pins, the die in the seal colour
    pins = []
    for v in (27, 37):
        pins += [st(line((v, 15), (v, 6)), W_FINE, 'load'), st(line((v, 49), (v, 58)), W_FINE, 'load'),
                 st(line((15, v), (6, v)), W_FINE, 'load'), st(line((49, v), (58, v)), W_FINE, 'load')]
    m['gpu'] = ('a chip with its pins, the die in the seal colour', [
        st(square(15, 15, 49, 49, r=3.0), W_MAIN)] + pins, [
        Blot(rounded([(32, 24.5), (24.5, 24.5), (24.5, 39.5), (39.5, 39.5), (39.5, 24.5), (32, 24.5)], 1.6))])

    # Free and portable: a ginkgo leaf, its blade fanned brush presses with the notch between them
    base = (32, 47)
    blades = []
    for ang, rr in ((-160, 22), (-135, 25), (-111, 25.5), (-69, 25.5), (-45, 25), (-20, 22)):
        a_ = math.radians(ang)
        tip = (base[0] + rr * math.cos(a_), base[1] + rr * math.sin(a_))
        blades.append(st(line(tip, (base[0] + 2 * math.cos(a_), base[1] + 2 * math.sin(a_)), 0.02), 9.0, 'load', 0.12))
    m['portable'] = ('a ginkgo leaf, its blade fanned brush presses with a notch', blades + [
        st(bezier(base, (32, 52), (34, 57), (38, 61)), W_MID, 'grow')], [])
    return m


def card_svg(label, ink, accent, theme):
    p = PALETTE[theme]
    out = (f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 64 64" width="48" height="48" role="img" '
           f'aria-label="{label}"><g fill="{p["ink"]}">' + ''.join(f'<path d="{s.outline()}"/>' for s in ink) + '</g>')
    if accent:
        out += f'<g fill="{p["seal"]}">' + ''.join(f'<path d="{s.outline()}"/>' for s in accent) + '</g>'
    return out + '</svg>\n'


# ------------------------------------------------------------------------------------------- files
def logo(strokes, theme):
    p = PALETTE[theme]
    ink = ''.join(f'<path d="{s.outline()}"/>' for s in strokes)
    return (f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 64 64" width="256" height="256" role="img" aria-label="{LABEL}">'
            f'<g stroke="{p["mesh"]}" stroke-linecap="round" opacity="0.55">{ground(strokes)}</g>'
            f'<g fill="{p["ink"]}">{ink}</g>{seal(p["seal"], p["ground"])}</svg>\n')


def favicon(strokes, theme):
    p = PALETTE[theme]
    xs = [x for s in strokes for x, _ in s.c]
    ys = [y for s in strokes for _, y in s.c]
    pad = max(s.w for s in strokes) * FAVICON_THICKEN * 0.6
    x0, y0, x1, y1 = min(xs) - pad, min(ys) - pad, max(xs) + pad, max(ys) + pad
    side = max(x1 - x0, y1 - y0) * 1.16
    vx, vy = (x0 + x1 - side) / 2, (y0 + y1 - side) / 2
    ink = ''.join(f'<path d="{s.outline(FAVICON_THICKEN)}"/>' for s in strokes)
    return (f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="{vx:.2f} {vy:.2f} {side:.2f} {side:.2f}" width="64" height="64" role="img" aria-label="FOSSIL">'
            f'<rect x="{vx:.2f}" y="{vy:.2f}" width="{side:.2f}" height="{side:.2f}" rx="{side * 0.1875:.2f}" fill="{p["ground"]}"/>'
            f'<g fill="{p["ink"]}">{ink}</g></svg>\n')


def main():
    out = sys.argv[1] if len(sys.argv) > 1 else os.path.join('docs', 'public')
    os.makedirs(out, exist_ok=True)
    strokes = ammonite()
    files = {
        'logo.svg': logo(strokes, 'light'),
        'logo-dark.svg': logo(strokes, 'dark'),
        'favicon.svg': favicon(strokes, 'light'),
        'favicon-dark.svg': favicon(strokes, 'dark'),
    }
    for slug, (label, ink, accent) in card_marks().items():
        files[os.path.join('icons', slug + '.svg')] = card_svg(label, ink, accent, 'light')
        files[os.path.join('icons', slug + '-dark.svg')] = card_svg(label, ink, accent, 'dark')
    os.makedirs(os.path.join(out, 'icons'), exist_ok=True)
    for name, text in files.items():
        with open(os.path.join(out, name), 'w', encoding='utf-8') as fh:
            fh.write(text)
        print(f'{os.path.join(out, name)}  {len(text)} bytes')


if __name__ == '__main__':
    main()

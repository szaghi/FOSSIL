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
    raise ValueError(kind)


class Stroke:
    def __init__(self, pts, width, kind, waver, seed):
        self.c = resample(pts)
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
        pts = left + right[::-1]
        return 'M' + ' L'.join(f'{x:.2f} {y:.2f}' for x, y in pts) + ' Z'


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
    for name, text in files.items():
        with open(os.path.join(out, name), 'w', encoding='utf-8') as fh:
            fh.write(text)
        print(f'{os.path.join(out, name)}  {len(text)} bytes')


if __name__ == '__main__':
    main()

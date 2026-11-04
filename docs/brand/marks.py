#!/usr/bin/env python3
"""FOSSIL brand marks: an ammonite meshed in triangles, carved by one light (issue #43).

The shell is a logarithmic spiral r = R exp(b (theta - theta_max)) growing 1.8 times per turn. Each
whorl is joined to the one inside it by a strip of quads, every quad split into two triangles, the
way a mesher would triangulate it; the innermost turn is a triangle fan. Every facet is flat (no
gradients) and takes a colour between a shadow and a ridge tone by how it faces a light from the
upper left: the flat mesh reads as a carved relief.

Writes, in docs/public:
  logo.svg       the mark for light pages (chocolate #4a2c16 to cream #e2c29b, core #3b2414)
  logo-dark.svg  the mark for dark pages (caramel #8a5a33 to ivory #fbecd2, core #5e3a1f)
  favicon.svg    the light mark (it reads on light and dark browser tabs)

Usage: python3 docs/brand/marks.py [output directory, default docs/public]
"""
import math
import os
import sys

LIGHT = dict(shadow='#4a2c16', ridge='#e2c29b', core='#3b2414')
DARK = dict(shadow='#8a5a33', ridge='#fbecd2', core='#5e3a1f')
LIGHT_DIR = (-0.6, -0.8)   # upper left; SVG y points down
EXPANSION = 1.8            # whorl growth per turn
TURNS = 3.5
STEP_DEG = 15              # segment angle of the mesh
RADIUS = 27.0
ROTATION = -1.2            # puts the aperture at the lower left
CENTRE = (32.0, 32.0)


def fmt(x):
    """0.1 units: below a pixel up to 640 px."""
    return f'{x:.1f}'.rstrip('0').rstrip('.')


def mix(a, b, t):
    """Colour between a and b (t in [0, 1])."""
    ca = [int(a[i:i + 2], 16) for i in (1, 3, 5)]
    cb = [int(b[i:i + 2], 16) for i in (1, 3, 5)]
    t = max(0.0, min(1.0, t))
    return '#' + ''.join(f'{round(ca[i] + (cb[i] - ca[i]) * t):02x}' for i in range(3))


def facets():
    """Triangles of the shell: (points, side) with side 0 outer face, 1 inner face, 2 core."""
    cx, cy = CENTRE
    b = math.log(EXPANSION) / (2 * math.pi)
    tmax = 2 * math.pi * TURNS

    def point(t):
        r = RADIUS * math.exp(b * (t - tmax))
        return (cx + r * math.cos(t + ROTATION), cy + r * math.sin(t + ROTATION))

    tris = []
    n = int(round((tmax - 2 * math.pi) / math.radians(STEP_DEG)))
    ts = [2 * math.pi + k * (tmax - 2 * math.pi) / n for k in range(n + 1)]
    for t0, t1 in zip(ts[:-1], ts[1:]):
        o0, o1, i0, i1 = point(t0), point(t1), point(t0 - 2 * math.pi), point(t1 - 2 * math.pi)
        tris.append(((o0, o1, i1), 0))
        tris.append(((o0, i1, i0), 1))
    for k in range(6):   # the innermost turn: a fan from the centre
        t0, t1 = 2 * math.pi * k / 6, 2 * math.pi * (k + 1) / 6
        tris.append((((cx, cy), point(t0), point(t1)), 2))
    return tris


def shade(points, side):
    """Lightness in [0, 1]: how the facet faces the light (outer faces point away from the centre)."""
    cx, cy = CENTRE
    mx = sum(p[0] for p in points) / 3
    my = sum(p[1] for p in points) / 3
    a = math.atan2(my - cy, mx - cx)
    nx, ny = (math.cos(a), math.sin(a)) if side == 0 else (-math.cos(a), -math.sin(a))
    lx, ly = LIGHT_DIR
    ln = math.hypot(lx, ly)
    return 0.5 + 0.5 * (nx * lx + ny * ly) / ln


def mark(palette, label):
    """The SVG of the mark (viewBox 64)."""
    body = []
    for points, side in facets():
        colour = palette['core'] if side == 2 else mix(palette['shadow'], palette['ridge'], shade(points, side))
        pts = ' '.join(f'{fmt(x)},{fmt(y)}' for x, y in points)
        body.append(f'<polygon points="{pts}" fill="{colour}" stroke="{colour}"/>')
    # a hairline stroke in the fill colour seals the anti-aliasing seams between facets
    return (f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 64 64" width="256" height="256" role="img" '
            f'aria-label="{label}">'
            '<g stroke-width="0.35" stroke-linejoin="round">' + ''.join(body) + '</g></svg>\n')


def main():
    out = sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.path.dirname(__file__), '..', 'public')
    os.makedirs(out, exist_ok=True)
    label = 'FOSSIL logo: an ammonite meshed in triangles, carved by light'
    files = {'logo.svg': mark(LIGHT, label), 'logo-dark.svg': mark(DARK, label), 'favicon.svg': mark(LIGHT, label)}
    for name, text in files.items():
        with open(os.path.join(out, name), 'w') as fh:
            fh.write(text)
        print(f'{name}: {len(text)} bytes')


if __name__ == '__main__':
    main()

#!/usr/bin/env python3
"""Render the animation data of a docs example as animated SVGs (issue #43).

Usage: anim2svg.py DATA OUTDIR

DATA is written by docs/examples/src/anim_airfoil.f90: a header (cells NX NY, frames NF, section points NP; the grid box
X0 X1 Y0 Y1), then per frame the fluid fraction of every cell (row by row), the signed distance at every grid node and
the section in world coordinates. Two images are written:

  airfoil-cut-cells.svg   the cut cells of each frame, shaded by their solid fraction, under the section
  airfoil-level-sets.svg  iso-distance contours of each frame (marching squares on the node distances)

The animation is SMIL with discrete frames (no script), so it plays wherever the image is shown, a GitHub README
included; the same data give the same bytes. Called by scripts/docs_examples.sh for the "!anim NAME" marker.
"""

from __future__ import annotations

import sys

CELL = 14  # cell size (px)
PAD = 16  # padding around the grid (px)
TITLE = 26  # height of the caption band (px)
FRAME = 0.12  # seconds per frame
BG = "#1b1714"  # the FOSSIL dark theme
GRID = "#3a312a"
INK = "#efe6da"
CREAM = "#e2c29b"
SOLID = "#8a5a33"
EMBER = "#e8590c"
LEVELS = [0.05, 0.1, 0.2, 0.3, 0.45, 0.6]


def read(path: str):
    tokens = open(path).read().split()
    pos = 0

    def take(n, conv=float):
        nonlocal pos
        out = [conv(t) for t in tokens[pos:pos + n]]
        pos += n
        return out

    nx, ny, nf, npt = take(4, int)
    x0, x1, y0, y1 = take(4)
    frames = []
    for _ in range(nf):
        assert tokens[pos] == "frame"
        pos += 2
        alpha = take(nx * ny)
        dist = take((nx + 1) * (ny + 1))
        sec = take(2 * npt)
        frames.append((alpha, dist, list(zip(sec[0::2], sec[1::2]))))
    return nx, ny, (x0, x1, y0, y1), frames


def fmt(v: float) -> str:
    return f"{v:.1f}".rstrip("0").rstrip(".")


def frame_group(k: int, n: int, body: str) -> str:
    """A group visible during frame k of n only (discrete SMIL)."""
    a, b = k / n, (k + 1) / n
    if k == 0:
        vals, times = "visible;hidden", f"0;{b:.4f}"
    else:
        vals, times = "hidden;visible;hidden", f"0;{a:.4f};{b:.4f}"
    anim = (f'<animate attributeName="visibility" values="{vals}" keyTimes="{times}" calcMode="discrete" '
            f'dur="{n * FRAME:.2f}s" repeatCount="indefinite"/>')
    return f'<g visibility="{"visible" if k == 0 else "hidden"}">{anim}{body}</g>'


def canvas(nx, ny, caption, frames_svg):
    w, h = nx * CELL + 2 * PAD, ny * CELL + 2 * PAD + TITLE
    grid = [f'<rect x="{PAD}" y="{PAD}" width="{nx * CELL}" height="{ny * CELL}" fill="none" stroke="{GRID}"/>']
    for i in range(1, nx):
        grid.append(f'M{PAD + i * CELL} {PAD}v{ny * CELL}')
    for j in range(1, ny):
        grid.append(f'M{PAD} {PAD + j * CELL}h{nx * CELL}')
    lines = f'<path d="{"".join(grid[1:])}" stroke="{GRID}" stroke-width="0.6"/>'
    text = (f'<text x="{PAD}" y="{h - 9}" fill="{INK}" font-family="ui-monospace,SFMono-Regular,Menlo,monospace" '
            f'font-size="13">{caption}</text>')
    return (f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {w} {h}" width="{w}" height="{h}" role="img" '
            f'aria-label="{caption}"><rect width="{w}" height="{h}" rx="10" fill="{BG}"/>{grid[0]}{lines}'
            f'{"".join(frames_svg)}{text}</svg>\n')


def main():
    data, outdir = sys.argv[1], sys.argv[2]
    nx, ny, (x0, x1, y0, y1), frames = read(data)
    sx = nx * CELL / (x1 - x0)
    sy = ny * CELL / (y1 - y0)

    def px(x, y):
        return PAD + (x - x0) * sx, PAD + (y1 - y) * sy

    def section(sec, fill):
        pts = " ".join(f"{fmt(px(x, y)[0])},{fmt(px(x, y)[1])}" for x, y in sec)
        return f'<polygon points="{pts}" fill="{fill}" stroke="{CREAM}" stroke-width="1.6" stroke-linejoin="round"/>'

    # cut cells: solid cells in sepia, cut cells in ember, more opaque the more solid they are
    out = []
    for k, (alpha, _, sec) in enumerate(frames):
        cells = []
        for j in range(ny):
            for i in range(nx):
                a = alpha[j * nx + i]
                if a >= 1.0:
                    continue
                x, y = PAD + i * CELL, PAD + (ny - 1 - j) * CELL
                if a <= 0.0:
                    cells.append(f'<rect x="{x}" y="{y}" width="{CELL}" height="{CELL}" fill="{SOLID}"/>')
                else:
                    cells.append(f'<rect x="{x}" y="{y}" width="{CELL}" height="{CELL}" fill="{EMBER}" '
                                 f'fill-opacity="{0.25 + 0.65 * (1 - a):.2f}"/>')
        out.append(frame_group(k, len(frames), "".join(cells) + section(sec, "none")))
    with open(f"{outdir}/airfoil-cut-cells.svg", "w") as fh:
        fh.write(canvas(nx, ny, "cut cells of a pitching NACA 0012: solid fraction per cell", out))

    # level sets: marching squares on the node distances
    out = []
    for k, (_, dist, sec) in enumerate(frames):
        paths = []
        for lev_i, lev in enumerate(LEVELS):
            seg = []
            for j in range(ny):
                for i in range(nx):
                    c = [(i, j), (i + 1, j), (i + 1, j + 1), (i, j + 1)]
                    d = [dist[jj * (nx + 1) + ii] for ii, jj in c]
                    pts = []
                    for e in range(4):
                        da, db = d[e], d[(e + 1) % 4]
                        if (da - lev) * (db - lev) < 0:
                            t = (lev - da) / (db - da)
                            (ia, ja), (ib, jb) = c[e], c[(e + 1) % 4]
                            gx = x0 + (ia + t * (ib - ia)) * (x1 - x0) / nx
                            gy = y0 + (ja + t * (jb - ja)) * (y1 - y0) / ny
                            pts.append(px(gx, gy))
                    for p, q in zip(pts[0::2], pts[1::2]):
                        seg.append(f"M{fmt(p[0])} {fmt(p[1])}L{fmt(q[0])} {fmt(q[1])}")
            opacity = 0.95 - 0.12 * lev_i
            paths.append(f'<path d="{"".join(seg)}" fill="none" stroke="{CREAM}" stroke-opacity="{opacity:.2f}" '
                         f'stroke-width="1.2" stroke-linecap="round"/>')
        out.append(frame_group(k, len(frames), "".join(paths) + section(sec, SOLID)))
    with open(f"{outdir}/airfoil-level-sets.svg", "w") as fh:
        fh.write(canvas(nx, ny, "signed distance of a pitching NACA 0012: level sets 0.05 to 0.6", out))


if __name__ == "__main__":
    main()

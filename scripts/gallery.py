#!/usr/bin/env python3
"""Write the figure gallery of the landing page and of the README from one list (issue #47).

Usage: gallery.py            rewrite the gallery regions of docs/index.md and README.md
       gallery.py --check    exit 1 if either region differs from what this list writes (CI: check_doc_snippets.sh)

Both pages carry a region between the lines ``<!-- gallery:start -->`` and ``<!-- gallery:end -->``: the two airfoil
animations, then one figure per feature card of the landing page, each computed by a program of docs/examples/src and
linked to the page that explains it. The regions are generated, so the README cannot show other figures than the site:
edit this list, never the regions.
"""

from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
SITE = "https://szaghi.github.io/FOSSIL"
BASE = "/FOSSIL"  # the site's base (docs/.vitepress/config.mts): raw HTML links are not rewritten by VitePress
START, END = "<!-- gallery:start -->", "<!-- gallery:end -->"

ANIMATIONS = [
    ("airfoil-cut-cells", "animation: the exact cut cells of a NACA 0012 section pitching and translating through a "
                          "48 by 24 grid, each cell shaded by its solid fraction"),
    ("airfoil-level-sets", "animation: iso-distance contours from 0.05 to 0.6 around the same pitching section"),
]
ANIMATION_CAPTION = ("A NACA 0012 section pitching and translating through a Cartesian grid, every frame computed by "
                     "FOSSIL: its exact cut cells and the level sets of its signed distance, with the pose passed to "
                     "the query instead of rebuilding the tree.")
# (feature, image, link, caption): the order of the feature cards
CARDS = [
    ("STL in, STL out", "tut01-dragon", "/manual/tutorial/01-first-surface",
     "the dragon loaded, 6588 facets, and the box analyze computes"),
    ("Repair dirty meshes", "anim-holes", "/manual/tutorial/02-repair",
     "holes closed one by one by fill_holes, refined and faired"),
    ("Signed distance", "tut04-probes", "/manual/tutorial/04-distance",
     "probes and their closest points on the dragon"),
    ("Level sets", "tut06-field", "/manual/tutorial/06-distance-field",
     "the signed-distance field on a stretched grid"),
    ("2-D curves", "tut07-field", "/manual/tutorial/07-curves",
     "a section of the dragon and its 2-D signed distance"),
    ("Immersed-boundary queries", "tut08-features", "/manual/tutorial/08-ib-queries",
     "closest points on a face, an edge or a vertex"),
    ("AMR tagging", "tut09-tags", "/manual/tutorial/09-amr-tagging",
     "cut, inside and outside cells of a block, in signed bands"),
    ("Exact cut cells", "tut10-cells2d", "/manual/tutorial/10-cut-cells",
     "the fluid fraction of the cut cells of a section"),
    ("Moving bodies", "tut11-swept", "/manual/tutorial/11-moving-bodies",
     "the area a moving section sweeps in a time step"),
    ("Mass properties", "anim-mass-properties", "/guide/advanced/mass-properties",
     "a T-handle flipping about its intermediate axis, from its inertia tensor"),
    ("Surface deviation", "anim-deviation", "/guide/advanced/surface-deviation",
     "a decimation sweep coloured by its certified distance from the original"),
    ("Feature edges", "anim-feature-tagging", "/guide/advanced/feature-edges",
     "cells touching the sharp edges of a part, as the grid refines"),
    ("Thickness and gaps", "anim-thin-regions", "/guide/advanced/thickness",
     "a wing's walls thinner than the cell, shrinking to the trailing edge"),
    ("Geometry processing", "anim-boolean", "/manual/cookbook#unite-two-bodies",
     "two cubes, then their union, intersection and differences"),
    ("On the GPU, bit for bit", "tut12-slice", "/manual/tutorial/12-gpu",
     "the device code path on a slice: 0 differences from the host"),
    ("Free and portable", "tour", "/manual/tutorial/01-first-surface",
     "the tutorial programs, run one after the other"),
]


def site_region() -> str:
    img = "./examples/images"
    out = [START, "", '<div class="fossil-casts">']
    out += [f'<img src="{img}/{name}.svg" alt="{alt}">' for name, alt in ANIMATIONS]
    out += ["</div>", "", f'<p class="fossil-gallery-caption">{ANIMATION_CAPTION}</p>', "",
            '<div class="fossil-gallery">']
    for title, name, link, caption in CARDS:
        out.append(f'<a class="fossil-card" href="{BASE}{link}">'
                   f'<img src="{img}/{name}.svg" alt="{title}: {caption}" loading="lazy">'
                   f"<strong>{title}</strong><span>{caption}</span></a>")
    out += ["</div>", "", END]
    return "\n".join(out)


def readme_region() -> str:
    img = "docs/examples/images"
    out = [START, "", '<p align="center">']
    out += [f'<img src="{img}/{name}.svg" alt="{alt}" width="420">' for name, alt in ANIMATIONS]
    out += [f"<br><sub>{ANIMATION_CAPTION}</sub></p>", "", "<table>"]
    for i in range(0, len(CARDS), 3):
        out.append("<tr>")
        for title, name, link, caption in CARDS[i:i + 3]:
            out.append(f'<td align="center" width="33%"><a href="{SITE}{link}"><img src="{img}/{name}.svg" '
                       f'alt="{title}: {caption}" width="260"></a><br><b>{title}</b><br><sub>{caption}</sub></td>')
        out.append("</tr>")
    out += ["</table>", "", END]
    return "\n".join(out)


def splice(path: Path, region: str, check: bool) -> bool:
    text = path.read_text()
    a, b = text.find(START), text.find(END)
    if a < 0 or b < a:
        raise SystemExit(f"gallery: no {START} ... {END} region in {path}")
    new = text[:a] + region + text[b + len(END):]
    if check:
        return new == text
    path.write_text(new)
    return True


def main(argv) -> int:
    check = "--check" in argv
    missing = [n for n, _ in ANIMATIONS] + [c[1] for c in CARDS]
    missing = [n for n in missing if not (ROOT / "docs/examples/images" / f"{n}.svg").exists()]
    if missing:
        print(f"gallery: missing images: {', '.join(missing)}")
        return 1
    ok = splice(ROOT / "docs/index.md", site_region(), check)
    ok = splice(ROOT / "README.md", readme_region(), check) and ok
    if check and not ok:
        print("gallery: docs/index.md or README.md is out of date: run scripts/gallery.py")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))

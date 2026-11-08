#!/usr/bin/env python3
"""Write the figure gallery of the landing page and of the README from one list (issues #47, #77).

Usage: gallery.py            rewrite the gallery regions of docs/index.md and README.md
       gallery.py --check    exit 1 if either region differs from what this list writes, if an image is missing, or if
                             the cards below do not match the feature cards of docs/index.md (CI: check_doc_snippets.sh)

Both pages carry a region between the lines ``<!-- gallery:start -->`` and ``<!-- gallery:end -->``: the hero
animations, two pairs, then one figure per feature card of the landing page, each computed by a program of
docs/examples/src and linked to the page that explains it. CARDS matches the ``features:`` cards of the frontmatter of
docs/index.md one to one, in the same order, with the same title and link: every card has a figure and every figure a
card. The regions are generated, so the README cannot show other figures than the site: edit this list, never the
regions.
"""

from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
SITE = "https://szaghi.github.io/FOSSIL"
BASE = "/FOSSIL"  # the site's base (docs/.vitepress/config.mts): raw HTML links are not rewritten by VitePress
START, END = "<!-- gallery:start -->", "<!-- gallery:end -->"

# the hero animations: pairs of (image, alt text), each pair with its caption
ANIMATIONS = [
    ([("airfoil-cut-cells", "animation: the exact cut cells of a NACA 0012 section pitching and translating through a "
                            "48 by 24 grid, each cell shaded by its solid fraction"),
      ("airfoil-level-sets", "animation: iso-distance contours from 0.05 to 0.6 around the same pitching section")],
     "A NACA 0012 section pitching and translating through a Cartesian grid, every frame computed by FOSSIL: its exact "
     "cut cells and the level sets of its signed distance, with the pose passed to the query instead of rebuilding the "
     "tree."),
    ([("anim-overset", "animation: an O-grid around a cylinder crossing a Cartesian background, its holes, receptors "
                       "and orphans recomputed in every frame"),
      ("anim-hyperbolic", "animation: an O-grid marched out of a NACA 0012 section, three layers per frame, up to 24 "
                          "layers")],
     "Grids around bodies: an overset assembly recomputed as the cylinder's O-grid moves through the background, with "
     "the orphans that appear where it leaves it, and a body-fitted O-grid marched out of the section by hyperbolic "
     "extrusion."),
]
# (feature, image, link, caption): the feature cards of docs/index.md, in their order
CARDS = [
    ("STL in, STL out", "tut01-dragon", "/manual/tutorial/01-first-surface",
     "the dragon loaded, 6588 facets, and the box analyze computes"),
    ("Repair dirty meshes", "anim-holes", "/manual/tutorial/02-repair",
     "holes closed one by one by fill_holes, refined and faired"),
    ("Geometry processing", "anim-boolean-csg", "/manual/tutorial/12-geometry",
     "a block drilled nine times, each step checked by an independent boolean"),
    ("Geometry analysis", "anim-deviation", "/manual/tutorial/12-geometry",
     "a decimation sweep coloured by its certified distance from the original"),
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
    ("Cut-cell solver support", "anim-cut-solver", "/guide/advanced/cut-solver",
     "the small pieces of a moving ellipse and their merge targets, every frame"),
    ("Moving bodies", "tut11-swept", "/manual/tutorial/11-moving-bodies",
     "the area a moving section sweeps in a time step"),
    ("Collision detection", "anim-ccd", "/manual/tutorial/14-contact",
     "conservative advancement to the time of impact: plain bound, 25 steps; directional, 5"),
    ("Structured blocks", "anim-block-walk", "/manual/tutorial/13-blocks-overset",
     "points tracked through an O-grid by stencil walks, across its seam"),
    ("Overset assembly", "ck-overset", "/manual/tutorial/13-blocks-overset",
     "an O-grid over a background: the holes, the receptors around them, the O-grid's receptors"),
    ("Body-fitted grids", "ck-wall-spacing", "/manual/tutorial/15-body-fitted",
     "an O-grid marched out of a NACA 0012 section by hyperbolic extrusion"),
    ("On the GPU, bit for bit", "tut16-slice", "/manual/tutorial/16-gpu",
     "the device code path on a slice: 0 differences from the host"),
    ("Free and portable", "tour", "/guide/installation",
     "the tutorial programs, run one after the other"),
]


def site_region() -> str:
    img = "./examples/images"
    out = [START, ""]
    for pair, caption in ANIMATIONS:
        out += ['<div class="fossil-casts">']
        out += [f'<img src="{img}/{name}.svg" alt="{alt}">' for name, alt in pair]
        out += ["</div>", "", f'<p class="fossil-gallery-caption">{caption}</p>', ""]
    out += ['<div class="fossil-gallery">']
    for title, name, link, caption in CARDS:
        out.append(f'<a class="fossil-card" href="{BASE}{link}">'
                   f'<img src="{img}/{name}.svg" alt="{title}: {caption}" loading="lazy">'
                   f"<strong>{title}</strong><span>{caption}</span></a>")
    out += ["</div>", "", END]
    return "\n".join(out)


def readme_region() -> str:
    img = "docs/examples/images"
    out = [START, ""]
    for pair, caption in ANIMATIONS:
        out += ['<p align="center">']
        out += [f'<img src="{img}/{name}.svg" alt="{alt}" width="420">' for name, alt in pair]
        out += [f"<br><sub>{caption}</sub></p>", ""]
    out += ["<table>"]
    for i in range(0, len(CARDS), 3):
        out.append("<tr>")
        for title, name, link, caption in CARDS[i:i + 3]:
            out.append(f'<td align="center" width="33%"><a href="{SITE}{link}"><img src="{img}/{name}.svg" '
                       f'alt="{title}: {caption}" width="260"></a><br><b>{title}</b><br><sub>{caption}</sub></td>')
        out.append("</tr>")
    out += ["</table>", "", END]
    return "\n".join(out)


def frontmatter_cards() -> list[tuple[str, str]]:
    """(title, link) of the ``features:`` cards of the frontmatter of docs/index.md, in their order."""
    text = (ROOT / "docs/index.md").read_text().split("\n---", 1)[0]
    cards: list[list[str]] = []
    in_features = False
    for line in text.splitlines():
        if not line.startswith(" "):
            in_features = line.startswith("features:")
            continue
        if not in_features:
            continue
        item = line.strip()
        if item.startswith("- "):
            cards.append(["", ""])
            item = item[2:]
        key, _, value = item.partition(":")
        if cards and key in ("title", "link"):
            cards[-1][0 if key == "title" else 1] = value.strip().strip('"')
    return [(t, ln) for t, ln in cards]


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
    status = 0
    names = [n for pair, _ in ANIMATIONS for n, _ in pair] + [c[1] for c in CARDS]
    missing = [n for n in names if not (ROOT / "docs/examples/images" / f"{n}.svg").exists()]
    if missing:
        print(f"gallery: missing images (run scripts/docs_examples.sh): {', '.join(missing)}")
        status = 1 if check else 0
    cards = [(c[0], c[2]) for c in CARDS]
    if frontmatter_cards() != cards:
        print("gallery: CARDS do not match the feature cards of docs/index.md one to one (title, link, order)")
        status = 1
    ok = splice(ROOT / "docs/index.md", site_region(), check)
    ok = splice(ROOT / "README.md", readme_region(), check) and ok
    if check and not ok:
        print("gallery: docs/index.md or README.md is out of date: run scripts/gallery.py")
        status = 1
    return status


if __name__ == "__main__":
    sys.exit(main(sys.argv))

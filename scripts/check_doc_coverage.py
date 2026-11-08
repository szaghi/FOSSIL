#!/usr/bin/env python3
"""Check that the documentation covers every feature with figures and run examples (issue #77).

Usage: check_doc_coverage.py            report the gaps, exit 0 unless an included file is missing
       check_doc_coverage.py --strict   exit 1 on any gap too (scripts/check_doc_snippets.sh runs it so)

Rules:
  - (issue #85) every level-2 section of a reference page or tutorial chapter that explains something (not API,
    Example, Known limitations, ...) shows a generated figure; every tutorial chapter shows an animation and at least
    two runs; the reference pages of iterative or time-dependent features (ANIMATED) show an animation; every entry
    point of the feature map appears in the snippet of a cookbook recipe;
  - every reference page (docs/guide/usage.md, docs/guide/advanced/*.md) and every tutorial chapter shows at least one
    generated figure (docs/examples/images/*.svg) and at least one run example (docs/examples/output/*.ansi);
  - every row of the feature map (docs/guide/features.md, tables with a Tutorial column) links a tutorial chapter or a
    cookbook recipe (a Cookbook column, when the table has one);
  - the feature cards of the landing page (docs/index.md frontmatter) and the gallery tiles (scripts/gallery.py CARDS)
    are the same list, in the same order;
  - every snippet, output and image a page includes exists (always an error: a missing file is a broken page);
  - every example program in docs/examples/src is shown by some page (its snippet, one of its outputs or figures).
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DOCS = ROOT / "docs"
EX = DOCS / "examples"
INCLUDE = re.compile(r"(?:<<<\s*@/|\]\()(?:\.\./)*(?:\./)?(examples/(?:snippets|output|images)/[A-Za-z0-9_.\-]+?)(?:\{[^}]*\})?(?:\)|\s|$)")
INCLUDE_HTML = re.compile(r'src="\./(examples/images/[A-Za-z0-9_.\-]+)"')


def pages() -> list[Path]:
    out = [p for p in DOCS.rglob("*.md") if "node_modules" not in p.parts and ".vitepress" not in p.parts]
    return sorted(out)


def includes(text: str) -> list[str]:
    return INCLUDE.findall(text) + INCLUDE_HTML.findall(text)


# sections that need no figure (issue #85): reference material, not explanation
SECTION_EXEMPT_REF = re.compile(r"^(API|Example|See also|References?|Left out|Known limitations|Limitations|"
                                r"Design|Sources?|Status|Notation|Proved and measured|Measured|Results|"
                                r"Compilers|Further reading|Related|Benchmarks?|Performance|Changelog|Command line|Validation|Limits|"
                                r"Where the library|Placement rules|Output for inspection|Recommendation)\b", re.I)
SECTION_EXEMPT_TUT = re.compile(r"^(Running it|Run it|Next|What next|Summary|Recap|Build|Building)\b", re.I)
# reference pages of iterative or time-dependent features: they show an animation (issue #85)
ANIMATED = {"isotropic-remesh", "decimate", "smoothing", "hole-filling", "alpha-wrap", "booleans", "distance-field",
            "cut-cells", "moving-bodies", "cut-solver", "curvilinear-blocks", "overset", "collision", "body-fitted",
            "surface-deviation", "thickness", "feature-edges", "mass-properties", "body-distance", "amr-tagging",
            "ib-queries"}
# feature-map entry points that are not routines a recipe calls
ENTRY_SKIP = {"_all", "_any", "_device", "fossil_cut_cells_3d", "pose", "geometry", "quality", "validity", "compute"}


def sections(text: str) -> list[tuple[str, str]]:
    """The level-2 sections of a page: (title, body), code fences ignored."""
    out, title, body, fence = [], None, [], False
    for line in text.splitlines():
        if line.startswith("```"):
            fence = not fence
        if not fence and line.startswith("## "):
            if title is not None:
                out.append((title, "\n".join(body)))
            title, body = line[3:].strip(), []
        elif title is not None:
            body.append(line)
    if title is not None:
        out.append((title, "\n".join(body)))
    return out


def is_animation(inc: str) -> bool:
    """A generated image with SMIL frames."""
    path = DOCS / inc
    return path.suffix == ".svg" and path.exists() and "<animate" in path.read_text()


def entry_points() -> list[str]:
    """The routine names of the feature map's Entry point column."""
    names: set[str] = set()
    for line in (DOCS / "guide/features.md").read_text().splitlines():
        cells = line.split("|")
        if not line.startswith("| ") or len(cells) < 4 or "`" not in cells[2]:
            continue
        for code in re.findall(r"`([^`]+)`", cells[2]):
            names.update(re.findall(r"[a-z_][a-z0-9_]*", code.split("(")[0].split("%")[-1]))
    return sorted(n for n in names if n not in ENTRY_SKIP)


def main() -> int:
    strict = "--strict" in sys.argv[1:]
    gaps: list[str] = []
    errors: list[str] = []
    shown: set[str] = set()

    for page in pages():
        text = page.read_text()
        for inc in includes(text):
            shown.add(Path(inc).name)
            if not (DOCS / inc).exists():
                errors.append(f"{page.relative_to(ROOT)}: missing {inc}")

    ref = [DOCS / "guide/usage.md"] + sorted((DOCS / "guide/advanced").glob("*.md"))
    tut = sorted((DOCS / "manual/tutorial").glob("*.md"))
    for page in ref + tut:
        text = page.read_text()
        inc = includes(text)
        figs = [i for i in inc if i.startswith("examples/images/")]
        runs = [i for i in inc if i.startswith("examples/output/")]
        name = page.relative_to(ROOT)
        missing = [what for what, n in (("figure", figs), ("run example", runs)) if not n]
        if missing:
            gaps.append(f"{name}: no {' and no '.join(missing)}")
        # issue #85: a figure in every section that explains something
        exempt = SECTION_EXEMPT_TUT if page in tut else SECTION_EXEMPT_REF
        for title, body in sections(text):
            if exempt.search(title):
                continue
            if not any(i.startswith("examples/images/") for i in includes(body)):
                gaps.append(f"{name}: section '{title}' has no figure")
        anims = [f for f in figs if is_animation(f)]
        if page in tut:
            if not anims:
                gaps.append(f"{name}: no animation")
            if len(set(runs)) < 2:
                gaps.append(f"{name}: fewer than two runs")
        elif page.stem in ANIMATED and not anims:
            gaps.append(f"{name}: no animation (an iterative or time-dependent feature)")

    # issue #85: every entry point of the feature map shown by a recipe
    cook = (DOCS / "manual/cookbook.md").read_text()
    recipe_code = "\n".join((DOCS / i).read_text() for i in includes(cook)
                            if i.startswith("examples/snippets/") and (DOCS / i).exists())
    for entry in entry_points():
        if not re.search(rf"\b{re.escape(entry)}\b", recipe_code):
            gaps.append(f"docs/manual/cookbook.md: no recipe shows `{entry}`")

    rows = 0
    for table in re.split(r"\n\s*\n", (DOCS / "guide/features.md").read_text()):
        lines = [l for l in table.splitlines() if l.startswith("|")]
        if len(lines) < 3:
            continue
        head = [c.strip() for c in lines[0].strip("|").split("|")]
        if "Tutorial" not in head:
            continue
        cols = [head.index(c) for c in ("Tutorial", "Cookbook") if c in head]
        for line in lines[2:]:
            cells = [c.strip() for c in line.strip("|").split("|")]
            rows += 1
            if not any(cells[c] for c in cols if c < len(cells)):
                gaps.append(f"docs/guide/features.md: no tutorial or recipe for '{cells[0][:70]}'")

    front = (DOCS / "index.md").read_text().split("---")[1]
    cards = re.findall(r"^\s+title:\s*(.+)$", front.split("features:")[1], re.M)
    cards = [c.strip().strip('"') for c in cards]
    sys.dont_write_bytecode = True  # no __pycache__ in scripts/
    sys.path.insert(0, str(ROOT / "scripts"))
    import gallery  # noqa: E402
    tiles = [c[0] for c in gallery.CARDS]
    if cards != tiles:
        gaps.append(f"landing: cards {cards} differ from the gallery tiles {tiles}")

    for src in sorted((EX / "src").glob("*.f90")):
        text = src.read_text()
        names = {f"{src.stem}.f90"} | {f"{src.stem}-{r}.f90" for r in re.findall(r"^\s*!region\s+(\S+)", text, re.M)}
        names |= {f"{i}.ansi" for i in re.findall(r"^\s*!run\s+(?:-s\s+)?(\S+)", text, re.M)}
        names |= {f"{n}.svg" for n in re.findall(r"^\s*!(?:fig|image)\s+(\S+)", text, re.M)}
        names |= {f"{n.split()[0]}.svg" for n in re.findall(r"^\s*!cast\s+(.+)$", text, re.M)}
        if src.stem == "anim_airfoil":
            names |= {"airfoil-cut-cells.svg", "airfoil-level-sets.svg"}
        if not names & shown:
            gaps.append(f"{src.relative_to(ROOT)}: shown by no page")

    for e in errors:
        print(f"ERROR {e}")
    for g in gaps:
        print(f"GAP   {g}")
    print(f"doc coverage: {len(ref)} reference pages, {len(tut)} tutorial chapters, {rows} feature rows, "
          f"{len(gaps)} gaps, {len(errors)} missing files{' (strict)' if strict else ''}")
    return 1 if errors or (strict and gaps) else 0


if __name__ == "__main__":
    sys.exit(main())

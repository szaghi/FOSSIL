#!/usr/bin/env bash
# Build and run the documentation examples, regenerating everything the pages include from them (issue #43).
#
#   docs/examples/src/*.f90    the example programs (hand-written), with marker comments:
#                                !run ID COMMAND        a run shown in the pages (!run -s: with its exit status)
#                                !region NAME ... !endregion NAME   a part of the program included on its own
#                                !image ID              the output of the run ID also as an image, images/ID.svg
#                                !cast NAME ID ID ...   the runs ID as one animated image, images/NAME.svg: a
#                                                       terminal session typing each command, then showing its output
#                                !anim NAME             the program wrote NAME.dat: scripts/anim2svg.py draws its
#                                                       geometry animations into images/ (NAME-*.svg)
#                                !fig NAME              the program wrote NAME.fig, a scene of FOSSIL's geometry:
#                                                       scripts/mesh2svg.py draws it into images/NAME.svg (issue #47)
#   docs/examples/lib/*.f90    modules shared by the programs (fossil_fig: the scene writer), compiled first
#   docs/examples/files/       input files, copied into the directory where the runs happen, with the test bodies
#                              src/tests/{dragon,bunny,cube-01,naca0012-ascii}.stl
#   docs/examples/snippets/    generated: <program>.f90 without the markers, and <program>-<region>.f90
#   docs/examples/images/      generated: <ID>.svg and the casts (scripts/ansi2svg.py)
#   docs/examples/output/      generated: <ID>.ansi, "$ COMMAND" then its output, the run directory shown as /home/user
#
# Ported from FoXy. The library is rebuilt from scratch by FoBiS (mode static-gnu, into static/) and the examples are
# built by the same compiler, gfortran or $FC: the outputs must not depend on the compiler (use explicit formats).
# The runs happen in one scratch directory, in the order of the files and of the lines (a run may read what an
# earlier one wrote), with HOME pointing there and a minimal environment, so nothing outside it is read or written.
#
# Usage: bash scripts/docs_examples.sh            (FC=gfortran-14 bash scripts/docs_examples.sh: another compiler)
set -euo pipefail

root=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
ex=$root/docs/examples
build=$root/build/docs-examples # git-ignored
run_dir=$build/run
home=/home/user                 # how the run directory is shown

fc=${FC:-gfortran}
mkdir -p "$root/build"
(cd "$root" && fobis clean --mode static-gnu && fobis build --mode static-gnu --fc "$fc") \
  > "$root/build/docs-examples.log" 2>&1 || {
  cat "$root/build/docs-examples.log"; echo "docs_examples: library build failed" >&2; exit 1; }
lib=$root/static/libfossil.a
mods=$root/static/mod
rm -rf -- "$build"
mkdir -p "$build/bin" "$build/mod" "$run_dir" "$ex/snippets" "$ex/output" "$ex/images"
rm -f -- "$ex"/snippets/*.f90 "$ex"/output/*.ansi "$ex"/images/*.svg

markers='^ *!(run|region|endregion|image|cast|anim|fig) '
# snippets: the whole program and each region, without the markers, dedented
dedent() { awk '{l[NR]=$0; if ($0 ~ /[^ ]/) {match($0, /^ */); if (m == "" || RLENGTH < m) m = RLENGTH}}
                END {for (i = 1; i <= NR; i++) print substr(l[i], m + 1)}' "$1"; }
for src in "$ex"/src/*.f90; do
  name=$(basename "$src" .f90)
  grep -Ev "$markers" "$src" > "$ex/snippets/$name.f90" || true
  for region in $(sed -n 's/^ *!region \([A-Za-z0-9_-]*\).*/\1/p' "$src"); do
    awk -v r="$region" '$1 == "!endregion" && $2 == r {on = 0}
                        on && $0 !~ /^ *!(run|region|endregion|image|cast|anim|fig) / {print}
                        $1 == "!region" && $2 == r {on = 1}' "$src" > "$build/region.f90"
    dedent "$build/region.f90" > "$ex/snippets/$name-$region.f90"
  done
done

# shared modules, then the programs
objs=()
for src in "$ex"/lib/*.f90; do
  [ -e "$src" ] || continue
  "$fc" -c -I"$mods" -J"$build/mod" -o "$build/$(basename "$src" .f90).o" "$src"
  objs+=("$build/$(basename "$src" .f90).o")
done
for src in "$ex"/src/*.f90; do
  "$fc" -I"$mods" -I"$build/mod" -J"$build/mod" -o "$build/bin/$(basename "$src" .f90)" "$src" "${objs[@]}" "$lib"
done

# runs, in the order of the files and of the lines
if [ -d "$ex/files" ]; then cp -R "$ex/files/." "$run_dir/"; fi
for body in dragon.stl bunny.stl cube-01.stl naca0012-ascii.stl; do cp "$root/src/tests/$body" "$run_dir/"; done # test bodies
run() { # run [-s] ID COMMAND
  local show=0 status=0
  if [ "$1" = -s ]; then show=1; shift; fi
  local id=$1; shift
  local cmd="$*"
  {
    printf '$ %s\n' "$cmd"
    (cd "$run_dir" && env -i HOME="$run_dir" PATH="$build/bin:/usr/bin:/bin" SHELL=/bin/bash LC_ALL=C \
                     LD_LIBRARY_PATH="${LD_LIBRARY_PATH:-}" GFORTRAN_UNBUFFERED_PRECONNECTED=y bash -c "$cmd" 2>&1) \
      || status=$?
    if [ $show = 1 ]; then printf '[exit status %d]\n' "$status"; fi
  } | sed -e "s|$run_dir|$home|g" > "$ex/output/$id.ansi"
}
for src in "$ex"/src/*.f90; do
  while IFS= read -r line; do
    line=${line#*!run }
    if [ "${line%% *}" = -s ]; then line=${line#-s }; run -s "${line%% *}" "${line#* }"
    else run "${line%% *}" "${line#* }"; fi
  done < <(grep -E '^ *!run ' "$src" || true)
done
# images
for id in $(sed -n 's/^ *!image \([A-Za-z0-9_-]*\).*/\1/p' "$ex"/src/*.f90); do
  python3 "$root/scripts/ansi2svg.py" "$ex/output/$id.ansi" "$ex/images/$id.svg"
done
while read -r name ids; do
  captures=()
  for id in $ids; do captures+=("$ex/output/$id.ansi"); done
  python3 "$root/scripts/ansi2svg.py" --cast "$ex/images/$name.svg" "${captures[@]}"
done < <(sed -n 's/^ *!cast \(.*\)/\1/p' "$ex"/src/*.f90)
for name in $(sed -n 's/^ *!anim \([A-Za-z0-9_-]*\).*/\1/p' "$ex"/src/*.f90); do
  python3 "$root/scripts/anim2svg.py" "$run_dir/$name.dat" "$ex/images"
done
for name in $(sed -n 's/^ *!fig \([A-Za-z0-9_-]*\).*/\1/p' "$ex"/src/*.f90); do
  python3 "$root/scripts/mesh2svg.py" "$run_dir/$name.fig" "$ex/images/$name.svg"
done
echo "docs_examples: $(ls "$ex"/src/*.f90 | wc -l) programs, $(ls "$ex"/output/*.ansi | wc -l) runs," \
     "$(grep -h '^ *!fig ' "$ex"/src/*.f90 | wc -l) figures"

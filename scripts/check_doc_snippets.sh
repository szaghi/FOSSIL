#!/usr/bin/env bash
# Extract every ```fortran code block from docs/guide/api/*.md and try to
# compile each one against the installed FOSSIL library. Fails on the first
# block that does not compile.
#
# Assumes:
#   - The static library has been built via `fobis build --mode static-gnu`
#     so that the library archive and its module files exist.
#   - PENF / VecFor modules are reachable via the build's module path.
#
# Usage: scripts/check_doc_snippets.sh
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
DOCS_DIR="$ROOT/docs/guide/api"
WORK="$(mktemp -d)"
FAILLOG="$WORK/.failures"
: > "$FAILLOG"
trap 'rm -rf "$WORK"' EXIT

# Locate the most recently built static library and its module directory.
LIB="$(find "$ROOT" -maxdepth 4 -type f -name 'libfossil*.a' \
        2>/dev/null | head -n1 || true)"
if [[ -z "${LIB:-}" ]]; then
  echo "error: libfossil.a not found; run 'fobis build --mode static-gnu' first" >&2
  exit 2
fi
LIB_DIR="$(dirname "$LIB")"
MOD_DIRS=()
while IFS= read -r d; do MOD_DIRS+=("-I" "$d"); done < <(
  find "$ROOT" -maxdepth 4 -type d -name 'mod*' 2>/dev/null
)

for md in "$DOCS_DIR"/*.md; do
  [[ -e "$md" ]] || continue
  # Split the file into snippets; awk emits one .f90 path per program block.
  snippets=$(
    awk '
      /^```fortran$/ { in_block = 1; next }
      /^```$/        { if (in_block) print "---SNIPPET-END---"; in_block = 0; next }
      in_block       { print }
    ' "$md" | awk -v base="$(basename "$md" .md)" -v work="$WORK" '
      BEGIN { idx = 0; out = "" }
      /^---SNIPPET-END---$/ {
        if (out != "") {
          idx++
          fn = sprintf("%s/%s_%02d.f90", work, base, idx)
          print out > fn
          close(fn)
          print fn
        }
        out = ""
        next
      }
      { out = out $0 "\n" }
    '
  )
  for f in $snippets; do
    # Skip non-program snippets (e.g., ASCII tree previews).
    if ! grep -qE '^program ' "$f"; then
      continue
    fi
    if ! gfortran "${MOD_DIRS[@]}" -c "$f" -o "$f.o" -J "$WORK" 2> "$f.err"; then
      echo "------------------------------------------------------------"
      echo "FAIL: $(basename "$f")"
      echo "------------------------------------------------------------"
      sed -n '1,40p' "$f"
      echo "--- gfortran stderr ---"
      cat "$f.err"
      echo "x" >> "$FAILLOG"
    fi
  done
done

failed=$(wc -l < "$FAILLOG" | tr -d ' ')
if [[ "$failed" -gt 0 ]]; then
  echo "snippet check: $failed failure(s)"
  exit 1
fi
echo "snippet check: all blocks compiled."

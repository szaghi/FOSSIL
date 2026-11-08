call body%find_self_intersections(pairs=crossings, status=status)
print '(A,I0,A,I0,A,I0)', 'pairs:    ', size(crossings), ', crossing ', count(crossings%kind == SELFINT_CROSSING), &
      ', sheets ', count(crossings%kind == SELFINT_SHEETS)

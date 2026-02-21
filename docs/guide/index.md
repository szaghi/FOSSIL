---
title: About FOSSIL
---

# About FOSSIL

**FOSSIL** (FOrtran Stereo Litography parser) is a pure Fortran 2003+ library for reading, writing, and manipulating [STL (stereolithography)](https://en.wikipedia.org/wiki/STL_(file_format)) mesh files. It provides three OOP-designed types — `file_stl_object`, `surface_stl_object`, and `facet_object` — that together handle the full lifecycle from file I/O to advanced geometry operations.

Modern computational science and engineering tools often require STL mesh handling. FOSSIL fills that gap with a clean, user-friendly API — no C or C++ dependencies, no external mesh libraries.

The library is designed with a test-driven approach; every surface operation is exercised by an automated test program in `src/tests/`.

## Authors

- Stefano Zaghi — [@szaghi](https://github.com/szaghi)

Contributions are welcome — see the [Contributing](contributing) page.

## Copyrights

FOSSIL is distributed under a multi-licensing system:

| Use case | License |
|----------|---------|
| FOSS projects | [GPL v3](http://www.gnu.org/licenses/gpl-3.0.html) |
| Closed source / commercial | [BSD 2-Clause](http://opensource.org/licenses/BSD-2-Clause) |
| Closed source / commercial | [BSD 3-Clause](http://opensource.org/licenses/BSD-3-Clause) |
| Closed source / commercial | [MIT](http://opensource.org/licenses/MIT) |

> Anyone interested in using, developing, or contributing to FOSSIL is welcome — pick the license that best fits your needs.
